#!/usr/local/bin/perl -w
#
# make_agp_file.pl
# v0.3
#
# dl
#
# Usage : make_agp_file.pl


# v0.3
# 010815 : dl  : PP version

#################################################################################
# Initialise variables                                                          #
#################################################################################


$|=1;
use Socket;
#use strict;
use vars qw ($debug $seq_len $sv_acc $sv_ver);
use lib '/wormsrv2/scripts';
use Wormbase;

 ##############################
 # Script variables (run)     #
 ##############################

my $maintainers = "dl1\@sanger.ac.uk krb\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate     = `date +%y%m%d`;   chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $version     = &get_cvs_version("$0");

my @gff_files = ('I','II','III','IV','V','X');
my $outdir    = "/wormsrv2/autoace/CHROMOSOMES";
my $datadir   = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";

#$debug = 1;

#################################################################################
# Main Loop                                                                     #
#################################################################################

foreach my $chromosome (@gff_files) {

    my $i = 1;
    my ($start,$stop,$clone,$acc,$gap_span,$offset,$span,$gpseq,$gspan,$last_stop,$last_start);
#    @reports = "";

    $file = "$outdir/CHROMOSOME_$chromosome.agp";
    $last_stop = 2;

    &usage(1,$chromosome) unless ("-e $datadir/CHROMOSOME_${chromosome}.clone_acc.gff");

    # read data from gff file
    open (GFF, "<$datadir/CHROMOSOME_$chromosome.clone_acc.gff");
    while (<GFF>) {
	
	$seq_len = "";
	$sv_acc  = "";
	$sv_ver  = "";
	
	# parse from each line
	($start, $stop) = (/(\d+)\s+(\d+)/); 
	($clone, $acc) = (/Sequence \"(\S+)\" acc\=(\S+)/);

        # catch gaps
	if ($last_stop < $start) {
	    $gap_span = $start - $stop{$i-1};
	    
	    if ($gap_span == 1) {
		push (@report,"Putative butt-end join\n");
	    }
	    else {
#		print "[ " . ($stop{$i-1}+1) . " : " . ($start-1) ."] so insert padding characters over the gap\n" if ($debug);
		push (@report, "$chromosome\t$stop{$i-1}\t" . ($start - 1) . "\t$i\tN\t$gap_span\n");
		$start{$i}  = $stop{$i-1} + 1;
		$stop{$i}   = $start - 1;
		$span{$i}   = $gap_span -1;
		$acc{$i}    = "gap";
		$clone{$i}  = "gap";
		$last_start = $start;
		$last_stop  = $stop;
		$i++;
	    }
	}
	
        # subsumed under previous clone
	if (($start > $last_start) && ($stop < $last_stop)) {
	    push (@report,"$clone is redundant\n");
	    next;
	}

	$clone{$i} = $clone;
	$start{$i} = $start;
	$stop{$i}  = $stop;
	$acc{$i}   = $acc;
	$span{$i}  = $stop - $start + 1;

#	printf "%8s [%8d => %8d : %6d] $acc{$i}\n", $clone{$i}, $start{$i}, $stop{$i}, $span{$i};
	&getseq($acc,$chromosome);	
	$ver{$i}    = $sv_ver;
#	$len{$i}    = $seq_len;
	$last_stop  = $stop;
	$last_start = $start;
	$i++;
	
    }
    
    close GFF;
    
    $start{$i} = $stop{$i-1} + 1;
    $limit = $i;
    
    print "\n";

    open (LOG, ">${file}.log");

    # write report lines (redundant clones, butt-ends, gaps) to the log file
    foreach (@report) {
	print LOG $_;
    }
    

    open (OUT, ">$file");
    for ($i=1; $i<$limit;$i++) {
	$span2get = $start{$i+1} - $start{$i};
	
	if ($clone{$i} eq "gap") {
	    printf LOG "%3d %8s\t[%8d => %8d] [%8d] : Pad with %6d bp of '-'s}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
	    $unique_stop = $start{$i} + $span2get;
	    print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tN\t$span2get\n";
	} 
	else {
	    printf LOG "%3d %8s\t[%8d => %8d] [%8d] : Extract %6d bp from accession $acc{$i}, version $ver{$i}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
	    $unique_stop = $start{$i} + $span2get ;
	    print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tF\t$acc{$i}.$ver{$i}\t1\t$span2get\t+\n";
	}
    }
    close LOG;
    close OUT;
}

 ##############################
 # hasta luego                #
 ##############################

exit(0);


#################################################################################
### Subroutines                                                               ###
#################################################################################


sub getseq {

    my ($acc,$chromosome) = @_;
    my $absent;
    
    my $querycontent1 = "-e+[embl-acc:'$acc']";
    my $querycontent2 = "-e+[emblnew-acc:'$acc']";
    
    my $request1 = "/srs6bin/cgi-bin/wgetz?$querycontent1";
    my $request2 = "/srs6bin/cgi-bin/wgetz?$querycontent2";
    
    my $server = "srs.ebi.ac.uk";

    if (!defined(open_TCP(*F,$server,80))) {
        print "Error connecting to server \n";
        exit(-1);
    }

    # get from emblnew
    print F "GET $request2 HTTP/1.0\n\n";
    print F "Accept: */*\n";
    print F "User-Agent: socketsrs/1.0\n\n";
    
    undef ($absent);
    while (my $return_line=<F>) {
	if ($return_line =~ /Your query retrieved no entries/) {
#	    print "Query returned no entry : '$request2'\n" if ($debug);
	    $absent = 1;
	    last;
	}
    }
    close F;
    
    # else get from embl
    if ($absent) {
	if (!defined(open_TCP(*G,$server,80))) {
	    print "Error connecting to server \n";
	    exit(-1);
	} 
	
	print G "GET $request1 HTTP/1.0\n\n";
	print G "Accept: */*\n";
	print G "User-Agent: socketsrs/1.0\n\n";
    }
    else {
	if (!defined(open_TCP(*G,$server,80))) {
	    print "Error connecting to server \n";
	    exit(-1);
	} 
	
	print G "GET $request2 HTTP/1.0\n\n";
	print G "Accept: */*\n";
	print G "User-Agent: socketsrs/1.0\n\n";
    }

# Parsing annotation
    
    while (my $return_line=<G>) {

    # Sequence version
    # SV   AL032679.2
	if ($return_line =~ /SV\s+(\S+)/) {
	    ($sv_acc,$sv_ver) = split (/\./, $1);
	}
    
    # Sequence
    # SQ   Sequence 2679 BP; 882 A; 510 C; 415 G; 872 T; 0 other;

	if ($return_line =~ /^SQ\s+Sequence\s(\d+)\sBP/) {
	    $seq_len = $1;
	    return ($sv_acc,$sv_ver,$seq_len);
	    last;
	}
    }
    close G;
}


 ########################################################
 # Output: successful network connection in file handle #
 ########################################################

sub open_TCP {
    my ($FS,$dest,$port) = @_;
    my $proto = getprotobyname ('tcp');
    socket ($FS,PF_INET,SOCK_STREAM,$proto);
    my $sin = sockaddr_in($port,inet_aton($dest));
    connect($FS,$sin) || return undef;
    my $old_fh = select($FS);
    $| = 1;
    select($old_fh);
}

sub usage {
    my $error = shift;
    my $chromosome = shift;
    if ($error == 1){ 
	# No gff file to work from
	print "The gff file '$datadir/CHROMOSOME_${chromosome}.clone_acc.gff' doesn't exist.\n";
	exit(0);
    }
    elsif ($error == 0) {
        # Normal help menu
	exec ('perldoc',$0);
    }
}


__END__

=pod

=head2   NAME - make_agp_file.pl

=head1 USAGE

=over 4

=item make_agp_file.pl [-options]

=back

make_agp_file produces agp format lists for each chromosome. These files
delineate the yellow brick road (aka golden path) through a tiling set
of EMBL/GenBank entries (with sequence versions).

make_agp_file requires gff files with accessions (this is made after running
GFFsplitter late in the build procedure).

make_agp_file.pl mandatory arguments:

=over 4

=item none,

=back

make_agp_file.pl  OPTIONAL arguments:

=over 4

=back

=cut

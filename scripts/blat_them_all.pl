#!/usr/local/bin/perl5.6.0
#
# blat_them_all.pl
# kj2
#
# wrapper to generate blat data by
# - getting sequence out of autoace
# - blating it against all ESTs/mRNAs
# - mapping it back to autoace
# - producing confirmed introns
# - producing virtual objects to put the data into
#
# -e : run everything for ESTs
# -m : run everything for mRNAs
# -x : run everything for parasitic nematode ESTs
# -o : run everything for miscellaneous peptides (worm non RNA division coding sequence, not HTG)
#
# -n : dump dna/chromosome data from autoace
# -b : blating (autoace.fa, chromosome.ace already present)
# -s : sorting/mapping (blat2ace.pl)
# -v : produce the virtual objects
#
# 16.10.01 Kerstin Jekosch
# 17.10.01 kj: modified to get everything onto wormsrv2 and to include an mRNA and parasitic nematode blatx option
# 26.11.01 kj: runs everything for autoace AND camace now
# 13.11.01 kj: added some file removers to tidy up at the end
# 14.11.01 kj: added option to just produce virtual objects
# 01.02.02 dl: added option to search miscPep file
# 01.02.02 dl: uncommented report logging & mail routines
# 01.02.02 dl: routine to convert '-' -> 'N' needs to be within the same BLOCK as the tace command
#            : else you get a zero length fasta file each time and the confirm intron routines fail
# 02.02.21 dl: typos in the naming of the confirmed_intron virtual objects
# 02.04.08 dl: old style logging for autoace.fa check, prevented complete run of subs
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2002-11-15 18:21:35 $

use strict;
use Ace;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Std;
use IO::Handle;
use vars qw($opt_e $opt_m $opt_x $opt_b $opt_s $opt_o $opt_v $opt_n $opt_h $opt_d $opt_c);
$|=1;


###############
# directories #
###############

our $dbdir = "/wormsrv2/autoace";
our $seq   = "/wormsrv2/autoace/BLAT/autoace.fa";               
our $chrom = "/wormsrv2/autoace/BLAT/chromosome.ace";
our $blat  = "/wormsrv2/autoace/BLAT";



my $bin     = "/wormsrv2/scripts";
my $db;

our %homedb;
our $blatex  = '/nfs/disk100/wormpub/blat/blat';
#our $giface  = '/nfs/disk100/wormpub/ACEDB/bin.ALPHA/giface';
our $giface  = '/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_5/giface';
our $data;
our %word = (
	     EST      => 'BLAT_EST',
	     mRNA     => 'BLAT_mRNA',
	     EMBL     => 'BLAT_EMBL',
	     NEMATODE => 'BLATX_NEMATODE',
	     );

########################################
# command-line options & ramifications #
########################################

getopts('emxbsvonhdc');

if ($opt_c) {
    $dbdir = "/wormsrv1/camace";
    $seq   = "/wormsrv1/camace/BLAT/autoace.fa";               
    $chrom = "/wormsrv1/camace/BLAT/chromosome.ace";
    $blat  = "/wormsrv1/camace/BLAT";
}


# Help pod documentation
&usage(0) if ($opt_h);

(our $debug = 1) if ($opt_d);

# Exit if no option choosen [n|b|s|v]
&usage(1) unless ($opt_n || $opt_b || $opt_s || $opt_v); 

# Exit if no data type choosen [EST|mRNA|EMBL|NEMATODE]
&usage(2) unless ($opt_e || $opt_m || $opt_o || $opt_x || $opt_n); 

# Exit if multiple data types choosen [EST|mRNA|EMBL|NEMATODE]
&usage(3) if (($opt_e + $opt_m + $opt_o + $opt_x) > 1);

# assign type variable
($data = 'EST')      if ($opt_e);
($data = 'mRNA')     if ($opt_m);
($data = 'EMBL')     if ($opt_o);
($data = 'NEMATODE') if ($opt_x);

my $query;
$query      = '/nfs/disk100/wormpub/analysis/ESTs/C.elegans_nematode_ESTs'     if ($opt_e); # EST data set
$query      = '/nfs/disk100/wormpub/analysis/ESTs/C.elegans_nematode_mRNAs'    if ($opt_m); # mRNA data set
$query      = '/nfs/disk100/wormpub/analysis/ESTs/non_C.elegans_nematode_ESTs' if ($opt_x); # ParaNem EST data set
$query      = '/nfs/disk100/wormpub/analysis/ESTs/C.elegans_nematode_miscPep'  if ($opt_o); # Other CDS data set

#############
# LOG stuff #
#############

my $maintainer = "All";
$maintainer = "dl1\@sanger.ac.uk" if ($debug);

my $rundate    = `date +%y%m%d`;   chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;
my $WS_version = &get_wormbase_version_name;

my $logfile = "/wormsrv2/logs/blat_them_all.${WS_version}.$rundate.$$";
system ("/bin/touch $logfile");
open (LOG,">>$logfile") or die ("Could not create logfile\n");
LOG->autoflush();
open (STDOUT,">>$logfile");
STDOUT->autoflush();
open (STDERR,">>$logfile"); 
STDERR->autoflush();

print LOG "# blat_them_all\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";
print LOG "WormBase version : ${WS_version}\n";
print LOG "\n";
print LOG "======================================================================\n";
print LOG " -e : perform blat for ESTs\n"                     if ($opt_e);
print LOG " -m : perform blat for mRNAs\n"                    if ($opt_m);
print LOG " -o : perform blat for miscPep\n"                  if ($opt_o);
print LOG " -x : perform blatx for parasitic nematode ESTs\n" if ($opt_x);
print LOG "======================================================================\n";
print LOG "\n";
print LOG "Starting blat process .. \n\n";

###########################
# Main loop               #
###########################

# Write sequence data from autoace
if ($opt_n) {
    
    # This should check if it needs to (re)dump the DNA.
    # HOW??? - 

    &dump_dna;
}  

# CHECK: 
&usage(11) unless (-e "/wormsrv2/autoace/BLAT/superlinks.ace");
    
# assign contigs to laboratory
%homedb = &which_db;

# BLAT the query data type 
if ($opt_b) {

    # CHECK: does the autoace.fa exist
    # exit if autoace.fa file is absent
    &usage(5) unless (-e "$seq");
    
    # CHECK: how old is the autoace.fa file ?
    # exit if autoace.fa file created prior to start of (re)build 
    &usage(6) if ( (-M "/wormsrv2/autoace/logs/A1:Build_in_progress" < -M "/wormsrv2/autoace/BLAT/autoace.fa") && (!$opt_c) );
   
    print LOG "running blat and putting the result in $blat/${data}_out.psl\n";

    # BLAT system call
    if ($opt_x) {
	system("$blatex $seq $query -t=dnax -q=dnax $blat/${data}_out.psl") && die "Blat failed\n";
    }
#    elsif ($opt_o) {
#	system("$blatex $seq $query -q=prot -t=dnax $blat/${data}_out.psl") && die "Blat failed\n";
#    }
    else {
	system("$blatex $seq $query $blat/${data}_out.psl") && die "Blat failed\n";
    }
}

# map to autoace #
if ($opt_s) {

    if ($opt_e) {

	# map BEST hits for whole genome
	print "Mapping blat data to autoace\n";
	
	unless ($opt_c) {
	    system("$bin/blat2ace.pl -ei") && die "Mapping failed\n"; 
	}
	else {
	    system("$bin/blat2ace.pl -eiz") && die "Mapping failed\n"; 
	}
    }

    if ($opt_m) {
	print "Mapping blat data to autoace\n";
	unless ($opt_c) {
	    system("$bin/blat2ace.pl -mi") && die "Mapping failed\n"; 
	}
	else {
	    system("$bin/blat2ace.pl -miz") && die "Mapping failed\n"; 
	}
    }

    if ($opt_o) {
	print "Mapping blat data to autoace\n";
	unless ($opt_c) {
	    system("$bin/blat2ace.pl -oi") && die "Mapping failed\n"; 
	}
	else {
	    system("$bin/blat2ace.pl -oiz") && die "Mapping failed\n"; 
	}
    }

    if ($opt_x) {
	print "Mapping blat data to autoace\n";
	system("$bin/blat2ace.pl -x") && die "Mapping failed\n"; 
    }

    # produce confirmed introns #
    if ($opt_e) {
	print "Producing EST confirmed introns in databases\n\n";

	print "Confirm intron data in autoace\n";
        &confirm_introns('autoace','EST');

	print "Confirm intron data in camace\n";
	&confirm_introns('camace', 'EST');

	print "Confirm intron data in stlace\n";	
	&confirm_introns('stlace', 'EST');
    }
    if ($opt_m) {
	print "Producing mRNA confirmed introns in databases\n";

	print "Confirm intron data in autoace\n";
        &confirm_introns('autoace','mRNA');

	print "Confirm intron data in camace\n";
	&confirm_introns('camace', 'mRNA');

	print "Confirm intron data in stlace\n";	
	&confirm_introns('stlace', 'mRNA');
    }
    if ($opt_o) {
	print "Producing EMBL confirmed introns in databases\n";

	print "Confirm intron data in autoace\n";
        &confirm_introns('autoace','EMBL');

	print "Confirm intron data in camace\n";
	&confirm_introns('camace', 'EMBL');

	print "Confirm intron data in stlace\n";	
	&confirm_introns('stlace', 'EMBL');
    }
}

# produce files for the virtual objects #

if ($opt_v) {

    print LOG "Producing $data files for the virtual objects\n\n";
    &virtual_objects_blat($data);
}


close LOG

# mail logfile to maintainer
&mail_maintainer("WormBase Report: blat_them_all ",$maintainer,$logfile);

##############################
# hasta luego                #
##############################

exit(0);

#################################################################################
### Subroutines                                                               ###
#################################################################################


###################################################
# dump_dna : get data out of autoace              #
###################################################

# tace query for chromosome DNA files and chromosome link files
# not run for blatting (-b), mapping (-s) and virtual_obj generation (-v)

sub dump_dna {

    local (*CHANGE,*NEW);

my $command;

    unless ($opt_c) {
	$command  = "query find Sequence \"CHROMOSOME*\"\nshow -a -f /wormsrv2/autoace/BLAT/chromosome.ace\n";
	$command .= "follow subsequence\nshow -a -f /wormsrv2/autoace/BLAT/superlinks.ace\n";
	$command .= "dna -f /wormsrv2/autoace/BLAT/autoace.first\nquit\n";
    }
    else {
	$command  = "query find Sequence \"SUPERLINK*\"\nshow -a -f /wormsrv1/camace/BLAT/superlinks.ace\n";
	$command .= "dna -f /wormsrv1/camace/BLAT/autoace.first\nquit\n";
    }

    my ($sequence,$name);
    
    # tace dump chromosome consensus sequences
    system("echo '$command' | $giface $dbdir") && &usage(5);
    

    # move '-'s into 'n's => blat excludes '-'
    # [021115 dl] Do we need this anymore? no more '-' in sequence
    unless ($opt_c) {
	open (CHANGE, "</wormsrv2/autoace/BLAT/autoace.first");
    }
    else {
	open (CHANGE, "</wormsrv1/camace/BLAT/autoace.first");
    }
    open (NEW, ">$seq");
    while (<CHANGE>) {
	chomp;
	$sequence = $_;
	$sequence =~ tr/-/n/;
	print NEW "$sequence\n";
    }
    close CHANGE;

    # remove intermediary sequence file
    unless ($opt_c) {
	unlink ('/wormsrv2/autoace/BLAT/autoace.first') if (-e '/wormsrv2/autoace/BLAT/autoace.first');
    }
    else {
	unlink ('/wormsrv1/camace/BLAT/autoace.first') if (-e '/wormsrv1/camace/BLAT/autoace.first');
    }
}


sub which_db {
    
    local (*LINK);
    my (%homedb,$name);
    
    print LOG "Assign LINK* objects to laboratory\n\n";
    # deal with the superlink objects
    open (LINK, "</wormsrv2/autoace/BLAT/superlinks.ace") || die "whoops $!";
    while (<LINK>) {
	if (/^Sequence\s+\:\s+\"(\S+)\"/) {
	    $name = $1;
	    next;
	}
	if (/^From_Laboratory\s+\"(\S+)\"/) {
	    $homedb{$name} = $1;
	    print LOG "assigning $1 to $name\n";
	    undef ($name);
 	    next;
	}
    }
    close LINK;
    
    print LOG "\n";
   
    # return hash for (super)link objects
    return (%homedb);

}

#############################
# virtual object generation #
#############################

sub virtual_objects_blat {
    
    my ($data) = shift;
    local (*OUT_autoace_homol,*OUT_camace_homol,*OUT_stlace_homol);
    local (*OUT_autoace_feat,*OUT_camace_feat,*OUT_stlace_feat);
    my ($name,$length,$total,$first,$second,$m,$n);

    # autoace
    open (OUT_autoace_homol, ">$blat/virtual_objects.autoace.$word{$data}.ace") or die "$!";
    open (OUT_autoace_feat,  ">$blat/virtual_objects.autoace.ci.$data.ace")     or die "$!";
    # camace
    open (OUT_camace_homol,  ">$blat/virtual_objects.camace.$word{$data}.ace")  or die "$!";
    open (OUT_camace_feat,   ">$blat/virtual_objects.camace.ci.$data.ace")      or die "$!";
    # stlace
    open (OUT_stlace_homol,  ">$blat/virtual_objects.stlace.$word{$data}.ace")  or die "$!";
    open (OUT_stlace_feat,   ">$blat/virtual_objects.stlace.ci.$data.ace")      or die "$!";
    
    open (ACE, "<$chrom") || die &usage(11);
    while (<ACE>) {
	if (/Subsequence\s+\"(\S+)\" (\d+) (\d+)/) {
	    $name   = $1;
	    $length = $3 - $2 + 1;
	    $total = int($length/100000) +1;
	    
	    # autoace
	    print OUT_autoace_homol "Sequence : \"$name\"\n";
	    print OUT_autoace_feat  "Sequence : \"$name\"\n";
	    # camace
	    print OUT_camace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
	    print OUT_camace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
	    # stlace
	    print OUT_stlace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
	    print OUT_stlace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
	    
	    for ($n = 0; $n <= $total; $n++) {
		$m      = $n + 1;
		$first  = ($n*100000) + 1;
		$second = $first + 149999;
		if (($length - $first) < 100000) {
		    $second = $length;
		    # autoace
		    print OUT_autoace_homol "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n";
		    print OUT_autoace_feat  "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
		    # camace
		    print OUT_camace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
		    print OUT_camace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
		    # stlace
		    print OUT_stlace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
		    print OUT_stlace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
		    last;
		}					
		else {
		    ($second = $length) if ($second >  $length);
		    # autoace
		    print OUT_autoace_homol "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n";
		    print OUT_autoace_feat  "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
		    # camace
		    print OUT_camace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
		    print OUT_camace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
		    # stlace
		    print OUT_stlace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
		    print OUT_stlace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
		}
	    }
	    print OUT_autoace_homol "\n";
	    print OUT_autoace_feat  "\n";
	    print OUT_camace_homol  "\n" if ($homedb{$name} eq "HX");
	    print OUT_camace_feat   "\n" if ($homedb{$name} eq "HX");
	    print OUT_stlace_homol  "\n" if ($homedb{$name} eq "RW");
	    print OUT_stlace_feat   "\n" if ($homedb{$name} eq "RW");
	}
    }
    close ACE;
    close OUT_autoace_homol;
    close OUT_autoace_feat;
    close OUT_camace_homol;
    close OUT_camace_feat;
    close OUT_stlace_homol;
    close OUT_stlace_feat;

    # clean up if you are dealing with parasitic nematode conensus data
    if ($data eq "NEMATODE") {
	unlink ("$blat/virtual_objects.autoace.ci.$data.ace");
	unlink ("$blat/virtual_objects.camace.ci.$data.ace");
	unlink ("$blat/virtual_objects.stlace.ci.$data.ace");
    }

}


sub confirm_introns {

    my ($db,$data) = @_;
    local (*GOOD,*BAD,*SEQ);

    # open the output files
    open (GOOD, ">$blat/$db.good_introns.$data.ace") or die "$!";
    open (BAD,  ">$blat/$db.bad_introns.$data.ace")  or die "$!";

    my ($link,@introns,$dna,$switch,$tag);

    ($tag = "cDNA") if ($opt_m || $opt_o);
    ($tag = "EST")  if ($opt_e);
    
    $/ = "";
    open(CI,  "<$blat/${db}.ci.${data}.ace")      or die "Cannot open $blat/$db.ci.$data.ace $!\n";
    while (<CI>) {
	next unless /^\S/;
	if (/Sequence : \"(\S+)\"/) {
	    $link = $1;
	    print "Sequence : $link\n";
	    @introns = split /\n/, $_;
	    
	    # get the link sequence #
	    print "Extracting DNA sequence for $link\n";
	    undef ($dna);

	    open(SEQ, "<$blat/autoace.fa") || &usage(5);
	    $switch = 0;
	    $/ = "\n";
	    
	    # added shortcuts next & last to speed this section

	    while (<SEQ>) {
		if (/^\>$link$/) {
		    $switch = 1;
		    next;
		}
		elsif (/^(\w+)$/) {
		    if ($switch == 1) {
			chomp;
			$dna .= $1;
		    }			
		}
		elsif ($switch == 1) {
		    $switch = 0;
		    last;
		}
		else { 
		    $switch = 0;
		}
	    }
	    close SEQ;
	    
	    print "DNA sequence is " . length($dna) . " bp long.\n";

	    # evaluate introns #
	    $/ = "";
	    foreach my $test (@introns) {
		if ($test =~ /Confirmed_intron/) {
		    my @f = split / /, $test;
		    
		    #######################################
		    # get the donor and acceptor sequence #
		    #######################################
		    
		    my ($first,$last,$start,$end,$pastfirst,$prelast);
		    if ($f[1] < $f[2]) {
			($first,$last,$pastfirst,$prelast) = ($f[1]-1,$f[2]-1,$f[1],$f[2]-2);
		    }
		    else {
			($first,$last,$pastfirst,$prelast) = ($f[2]-1,$f[1]-1,$f[2],$f[1]-2);
		    }	

		    $start = substr($dna,$first,2);
		    $end   = substr($dna,$prelast,2);

#		    print "Coords start $f[1] => $start, end $f[2] => $end\n";
		    
		    ##################
		    # map to S_Child #
		    ##################
		    
		    my $lastvirt = int((length $dna) /100000) + 1;
		    my ($startvirt,$endvirt,$virtual);
		    if ((int($first/100000) + 1 ) > $lastvirt) {
			$startvirt = $lastvirt;
		    }
		    else {
			$startvirt = int($first/100000) + 1;
		    }
		    if ((int($last/100000) + 1 ) > $lastvirt) {
			$endvirt = $lastvirt;
		    }
		    else {
			$endvirt = int($first/100000) + 1;
		    }
		    
		    if ($startvirt == $endvirt) { 
			$virtual = "Confirmed_intron_EST:" .$link."_".$startvirt unless ($opt_m);
			$virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($opt_m);
		    }
		    elsif (($startvirt == ($endvirt - 1)) && (($last%100000) <= 50000)) {
			$virtual = "Confirmed_intron_EST:" .$link."_".$startvirt unless ($opt_m);
			$virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($opt_m);
		    }
		
		    #################
		    # check introns #
		    #################
		    
		    my $firstcalc = int($f[1]/100000);
		    my $seccalc   = int($f[2]/100000);
		    print STDERR "Problem with $test\n" unless (defined $firstcalc && defined $seccalc); 
		    my ($one,$two);
		    if ($firstcalc == $seccalc) {
			$one = $f[1]%100000;
			$two = $f[2]%100000;
		    }
		    elsif ($firstcalc == ($seccalc-1)) {
			$one = $f[1]%100000;
			$two = $f[2]%100000 + 100000;
			print STDERR "$virtual: $one $two\n";
		    }
		    elsif (($firstcalc-1) == $seccalc) {
			$one = $f[1]%100000 + 100000;
			$two = $f[2]%100000;
			print STDERR "$virtual: $one $two\n";
		    } 
		    print STDERR "Problem with $test\n" unless (defined $one && defined $two); 
		    
		    if ( ( (($start eq 'gt') || ($start eq 'gc')) && ($end eq 'ag')) ||
			 (  ($start eq 'ct') && (($end eq 'ac') || ($end eq 'gc')) ) ) {	 
			print GOOD "Feature_data : \"$virtual\"\n";
			print GOOD "Confirmed_intron $one $two $tag\n\n";
		    }  	
		    else {
			print BAD "Feature_data : \"$virtual\"\n";
			print BAD "Confirmed_intron $one $two $tag\n\n";		
		    }
		}
	    }
	}
    }
    close CI;

    close GOOD;
    close BAD;

}




sub usage {
    my $error = shift;

    if ($error == 1) {
	# no option supplied
	print "\nNo process option choosen [-n|b|s|v]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    elsif ($error == 2) {
	# No data-type choosen
	print "\nNo data option choosen [-e|m|o|x]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    elsif ($error == 3) {
	# 'Multiple data-types choosen
	print "\nMultiple data option choosen [-e|m|o|x]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    elsif ($error == 4) {
	# 'autoace.fa' file exists
	print "\nThe WormBase 'autoace.fa' file exists already\n";
	print "Remove this file before starting the process again\n\n";
	exit(0);
    }
    elsif ($error == 5) {
	# 'autoace.fa' file is not there or unreadable
	print "\nThe WormBase 'autoace.fa' file you does not exist or is non-readable.\n";
	print "Check File: '/wormsrv2/autoace/BLAT/autoace.fa'\n\n";
	exit(0);
    }
    elsif ($error == 6) {
	# BLAT failure
	print "\BLAT failure.\n";
	print "Whoops! you're going to have to start again.\n\n";
	exit(0);
    }
    elsif ($error == 11) {
	# 'superlinks.ace' file is not there or unreadable
	print "\nThe WormBase 'superlinks.ace' file you does not exist or is non-readable.\n";
	print "Check File: '/wormsrv2/autoace/BLAT/superlinks.ace'\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }
}

__END__

=pod

=head2   NAME - blat_them_all.pl

=head1 USAGE

=over 4

=item  blat_them_all.pl -options

=back

wrapper to generate blat data by getting sequence out of autoace, blating it against all ESTs/mRNAs, 
mapping it back to autoace/camace, producing confirmed introns, producing virtual objects to put the data into

blat_them_all mandatory arguments:

=over 4

=item -e   run everything for ESTs

=back

or

=item -m   run everything for mRNAs

=back

or

=item -o   run everything for miscPep

=back

or

=item -x   run everything for nematode EST consensus sequences

=back

blat_them_all optional arguments

=item -b   start with blating (autoace.fa, chromosome.ace already present)

=item -s   start later with sorting/mapping (use blat2ace.pl on whatever_out.psl)

=back

=cut




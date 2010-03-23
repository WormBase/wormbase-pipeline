#!/software/worm/bin/perl -w
#
# remap_wig.pl
#
# by ar2
#
# This is to remap wig plot files between genome release where coordinates have changed.
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2010-03-23 10:00:43 $

use strict;
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/bioperl-live';
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Modules::Remap_Sequence_Change;
use Bio::Graphics::Wiggle;
use Bio::Graphics::Wiggle::Loader;

my ($help, $debug, $test, $store, $species, $wormbase);
my $wig;
my ($old, $new, $chrom, $out, $split, $dirty);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "store:s"    => \$store,
	    "species:s"  => \$species,
	    "wig:s"      => \$wig,
	    "old:i"      => \$old,
	    "new:i"      => \$new,
	    "chrom:s"    => \$chrom,
	    "split"      => \$split, # remapping needs to be done per chrom so file may need splitting first.
	    "out:s"      => \$out,   # output file
	    "dirty"      => \$dirty
            );

die "missing data.  Valid file, old and new releases are required\n" unless ((-e $wig) and $old and $new);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species
                             );
}

#&convert_to_wdb($wig);


my $log = Log_files->make_build_log($wormbase);
my @mapping_data = Remap_Sequence_Change::read_mapping_data($old, $new, $wormbase->species);

if( $split ) { 
    my %chroms; # file handle for each sequence
    my $track_desc;
    open(FH,"<$wig") or $log->log_and_die("cant open $wig: $!\n");
    my $printer;
    my $linecount =0;
    while(<FH>) {
	$linecount++;
	if(/track/) {
	    $track_desc = $_;
	    next;
	}
	if(/chrom=(\w+)\s/) {
	    if(defined $printer) {
		close $printer;
		print STDERR $linecount,"\n";
		
	    }
	    open CH,">>${wig}_$1" or $log->log_and_die("cant open ${wig}_$1: $!\n");
	    $printer = *CH;
	    $chroms{$1} = "${wig}_$1";
	}
	print $printer $_;
    }    
    close $printer;
    print STDERR $linecount,"\n";

    #stickem back together and put track desc at top
    $out = "$wig.".$old."_".$new unless $out;
    open (OUT,">$out") or $log->log_and_die("cant open $out : $!\n");
    print OUT $track_desc;
    close OUT or $log->log_and_die("error closing $out : $!\n");
    
    foreach my $seq (keys %chroms){
	&remap_file($chroms{$seq},$seq,$chroms{$seq}."_remap");
	$wormbase->run_command("cat ".$chroms{$seq}."_remap >> $out",$log);
	$wormbase->run_command("rm -f ".$chroms{$seq}."_remap",$log) unless $dirty;    #cleanup
	$wormbase->run_command("rm -f ".$chroms{$seq},$log) unless $dirty;             #cleanup
    }    
}
else {
    &remap_file($wig,$chrom,$out);
}

sub remap_file {
    my $wig = shift; #input
    my $chrom = shift;
    my $out = shift; #output

    open(WIG,"<$wig") or $log->log_and_die("cant open wig file, $wig $!\n");
    $out = "$wig.".$old."_".$new unless $out;
    open(NEW,">$out") or $log->log_and_die("cant write output to $out $!\n");

    while(<WIG>) {
	unless(/^[a-z]/) {
	    my @data = split;
	    if( scalar @data > 2) {
		#BED file
		my ($chrom, $coord1, $coord2, $score) = @data;
		print NEW join("\t",("$chrom",Remap_Sequence_Change::remap_BED($chrom,$coord1, $coord2 ,$old, $new, \@mapping_data),$score)),"\n"
		}
	    else {
		# wig file
		my ($coord, $score) = @data;
		print NEW Remap_Sequence_Change::remap_wig($chrom,$coord,$old, $new, \@mapping_data)."\t$score\n";
	    }
	}
	else {
	    print NEW;
	}
    }
}

sub convert_to_wdb {
    my $file = shift;
    my $fh = IO::File->new($file);
    my $wigger = Bio::Graphics::Wiggle->new("$file.wib",1);
    while(<$fh>){ 
	my @data = split(/\s+/,$_);
	if($data[0] =~ /^\d+$/){
	    $wigger->set_value($data[0]=>$data[1]);
	}
    }
}

$log->mail;
exit;

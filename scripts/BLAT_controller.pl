#!/usr/local/bin/perl5.8.0 -w
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2006-04-25 12:40:23 $

use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Log_files;
use Storable;

my ($test, $database, $debug);
my ($mask, $dump, $run, $postprocess, $ace, $load, $process, $virtual);
my @types;
my $all;
my $store;

GetOptions (
	    'debug:s'     => \$debug,
	    'test'        => \$test,
	    'database:s'  => \$database,
	    'store:s'     => \$store,

	    'mask'        => \$mask,
	    'dump'        => \$dump,
	    'process'     => \$process,
	    'virtual'     => \$virtual,
	    'run'         => \$run,
	    'postprocess' => \$postprocess,
	    'ace'         => \$ace,
	    'load'        => \$load,

	    'types:s'     => \@types,
	    'all'         => \$all
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}


my $log = Log_files->make_build_log($wormbase);

my $wormpub = $wormbase->wormpub;
my $EST_dir = "$wormpub/analysis/ESTs";

$database = $wormbase->autoace unless $database;
my $blat_dir = $wormbase->blat;

#make sure passed type is valid.
@types = split(/,/,join(',',@types));
my @alltypes = qw( est mrna ncrna ost nematode embl tc1 washu nembase );
foreach my $t(@types) {
  $log->log_and_die("invalid type passed to $0 : $t\n") unless (grep {$t eq $_} @alltypes);
}
@types = @alltypes if ($all or !(@types));



if( $mask ) {
  my $opts = '-'.join(" -",@types);
  $log->write_to("running transcriptmasker.pl -$opts\n");
  $wormbase->run_script("transcriptmasker.pl -$opts", $log);

  # copy non-masked sequence type to database BLAT dir
  $log->write_to("copying nematode_ESTs, TC1s and embl_cds to $EST_dir\n");

  $wormbase->run_command("scp $EST_dir/other_nematode_ESTs      $blat_dir/", $log);
  $wormbase->run_command("scp $EST_dir/nembase_nematode_contigs $blat_dir/", $log);
  $wormbase->run_command("scp $EST_dir/washu_nematode_contigs   $blat_dir/", $log);
  $wormbase->run_command("scp $EST_dir/elegans_TC1s             $blat_dir/", $log);
  $wormbase->run_command("scp $EST_dir/elegans_embl_cds         $blat_dir/", $log);
}


# Now make blat target database using autoace (will be needed for all possible blat jobs)
# This also makes a backup copy of the old psl files (in case you need them to refer to)

&dump_dna if $dump;

# run all blat jobs on distributed cluster ( eg cbi1 )
if ( $run ) {
  foreach my $type ( @types ){
    $log->write_to("blatting $type\n");
    $wormbase->run_script("batch_BLAT.pl -$type -blat", $log);
  }
}

if( $postprocess ) {
  # merge psl files and convert to ace format
  $log->write_to("merging PSL files \n");
  system("cat $blat_dir/nematodeEST_* | perl -ne 'print if (/^[0-9]/)' > $blat_dir/nematode_out.psl"); # /d causes compiler warning (?)
  system("cat $blat_dir/elegansEST_*  | perl -ne 'print if (/^[0-9]/)' > $blat_dir/est_out.psl");
}

if ( $process or $virtual ) {
  foreach my $option (@types ) {
    #create virtual objects
    $wormbase->run_script("blat_them_all.pl -virtual -$option", $log) if $virtual;
    $wormbase->run_script("blat_them_all.pl -process -$option", $log) if $process;
  }
}

#if( $ace ) {
#  foreach my $type ( @types ) {
#    my $cmd = "blat2ace.pl -$type";
#    $cmd .= " -intron" if (( $type eq "mrna") or ($type eq "est") or ($type eq "ost") );
#    $wormbase->run_script("$cmd", $log);
#  }
#}

if( $load ) {
  foreach my $type (@types){
    $log->write_to("loading BLAT data - $type\n");

    # virtual objs
    my $file =  "$blat_dir/virtual_objects.autoace.blat.$type.ace";
    $wormbase->load_to_database( $database, $file,"virtual_objects_$type");

    # Don't need to add confirmed introns from nematode data (because there are none!)
    unless ( ($type eq "nematode") || ($type eq "washu") || ($type eq "nembase") || ($type eq "tc1") || ($type eq "embl")|| ($type eq "ncrna") ) {
      $file = "$blat_dir/virtual_objects.autoace.ci.$type.ace"; 
      $wormbase->load_to_database($database, $file, "blat_confirmed_introns_$type");

      $file = "$blat_dir/autoace.good_introns.$type.ace";
      $wormbase->load_to_database($database, $file, "blat_good_introns_$type");
    }

    # BLAT results
    $file = "$blat_dir/autoace.blat.$type.ace";
    $wormbase->load_to_database($database, $file, "blat_${type}_data");
  }
}

$log->mail;
exit(0);


###############################################################################################################


#############################################################################
# dump_dna                                                                  #
# gets data out of autoace/camace, runs tace query for chromosome DNA files #
# and chromosome link files.                                                #
#############################################################################

sub dump_dna {

  $log->write_to("dumping DNA from $database\n");
  local (*CHANGE,*NEW);

  my $command;
  my $giface = $wormbase->giface;

  $command  = "query find Sequence \"CHROMOSOME*\"\n";
  $command .= "show -a -f $blat_dir/chromosome.ace\n";
  $command .= "follow Subsequence\n";
  $command .= "show -a -f $blat_dir/superlinks.ace\n";
  $command .= "dna -f $blat_dir/autoace.first\nquit\n";

  # tace dump chromosomal DNA and superlinks file
  $wormbase->run_command("echo '$command' | $giface $database", $log);

  # Change '-'s in chromosome sequences into 'n's because blat excludes '-'
  # Not strictly needed anymore but left in for safety
  my $sequence;

  open (CHANGE, "<$blat_dir/autoace.first");
  open (NEW, ">$blat_dir/autoace.fa");
  while (<CHANGE>) {
    chomp;
    next unless (/\w+/); # remove blank lines acedb puts in at start of seq dumps.
    $sequence = $_;
    $sequence =~ tr/-/n/;
    print NEW "$sequence\n";
  }
  close(CHANGE);
  close(NEW);

  # remove intermediary sequence file
  unlink ("${blat_dir}/autoace.first") if (-e "${blat_dir}/autoace.first");

}

__END__

=pod

=head1 NAME - BLAT_controller.pl

=head2 USAGE

BLAT_controller.pl is a wrapper for the several scripts involved in running the various BLAT jobs in the WormBase build

BLAT_controller.pl  arguments:

=over 4

=item 

-debug user  =  debug mode to specify who gets log mail

=item 

-test        =  use the test build

=item 

-store string  =  load a previously serialised Wormbase object ( from Wormbase.pm ) to maintain test and debug setting when called from autoace_builder.pl

=item 

-database string = run on database other than build


 * script action options

=item 

-mask        runs transcriptmasker.pl  to mask ESTs polyA TSLs etc

=item 

-dump        dumps the target DNA sequences to

=item 

-process     runs blat_them_all.pl -process which runs blat2ace.pl to convert psl -> ace files

=item 

-virtual     runs blat_them_all.pl -virtual to create the virtual object the BLAT results hang off

=item 

-run         runs batch_BLAT.pl which shatters ESTs and nematode and submits bsub jobs to cbi1

=item 

-postprocess merges the psl files from the shattered EST and nematode files

=item 

-load       loads the BLAT acefiles in to the specified database

=item 

-types string  allows the specification of certain types to run BLAT. Valid types : est mrna ncrna ost nematode embl tc1 washu nembase

=item 

-all        runs all of the above tpyes.

=cut

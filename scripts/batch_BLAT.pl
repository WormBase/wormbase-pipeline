#!/usr/local/bin/perl5.8.0 -w
#
# batch_BLAT.pl
#
# Anthony Rogers
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2007-05-23 13:33:11 $
 

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $verbose, $est, $mrna, $ncrna, $ost, $nematode, $washu, $nembase, $embl, $tc1, $all, $export, $no_bsub);
my ($blat, $process, $virtual);
my ($store, $test);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "est"        => \$est,
	    "mrna"       => \$mrna,
	    "ncrna"      => \$ncrna,
	    "ost"        => \$ost,
	    "nematode"   => \$nematode,
	    "nembase"    => \$nembase,
	    "washu"      => \$washu,
	    "embl"       => \$embl,
	    "tc1"        => \$tc1,
	    "all"        => \$all,
	    "dump"       => \$export,
	    "no_bsub"    => \$no_bsub,
	    "blat"       => \$blat,
	    "process"    => \$process,
	    "virtual"    => \$virtual,
	    "store:s"    => \$store,
	    "test"       => \$test
	    );

if( $help ) { system ('perldoc',$0); exit(0);}

# turn everything on for -all option
if( $all ) {   
    $est      = 1;
    $mrna     = 1;
    $ncrna    = 1;
    $ost      = 1;
    $nematode = 1;
    $nembase  = 1;
    $washu    = 1;
    $embl     = 1;
    $tc1      = 1;
}
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

$log->log_and_die("failed do run blat and process / virtual at same time\n") if ($blat and ( $process or $virtual ));

if ( $blat ) {
    
    # ESTs (~2 hours) 
    &split_run( "est" ) if $est;
    
    # OSTs (5min) 
    &run_bsub("elegans_OSTs.masked", "ost_out.psl") if $ost;
    
    # TC1s (10min)
    &run_bsub("elegans_TC1s", "tc1_out.psl") if $tc1;
    
    # mRNA (2 hours)
    &run_bsub("elegans_mRNAs.masked", "mrna_out.psl", " -fine") if $mrna;
    
    # EMBL (5 mins) 
    &run_bsub("elegans_embl_cds", "embl_out.psl") if $embl;
    
    # ncRNAs (?) 
    &run_bsub("elegans_ncRNAs.masked", "ncrna_out.psl") if $ncrna;
    
    # Nembase contigs
    &split_run( "nembase_nematode_contigs", "nembase_out.psl", "-t=dnax -q=dnax" ) if $nembase;
    
    # WashU contigs
    &split_run( "washu_nematode_contigs", "washu_out.psl", "-t=dnax -q=dnax" ) if $washu;
    
    # splitting Nematode_ESTs
    &split_run( "nematode" ) if ( $nematode );
    
} # end $blat loop

if ( $process or $virtual ) {

  my @blat_jobs;
  push(@blat_jobs,"est")      if ( ($est)      || ($all) );
  push(@blat_jobs,"ost")      if ( ($ost)      || ($all) );
  push(@blat_jobs,"mrna")     if ( ($mrna)     || ($all) );
  push(@blat_jobs,"ncrna")    if ( ($ncrna)    || ($all) );
  push(@blat_jobs,"embl")     if ( ($embl)     || ($all) );
  push(@blat_jobs,"tc1")      if ( ($tc1)      || ($all) );
  push(@blat_jobs,"nembase")  if ( ($nembase)  || ($all) );
  push(@blat_jobs,"washu")    if ( ($washu)    || ($all) );
  push(@blat_jobs,"nematode") if ( ($nematode) || ($all) );

  # once the jobs have finished and been processed
  # run the main blat job

  foreach my $option (@blat_jobs ) {
      $wormbase->run_script("blat_them_all.pl -process -$option", $log);
      # also create virtual objects
      $wormbase->run_script("blat_them_all.pl -virtual -$option", $log);
  }
} # end loop for $process or $virtual 

$log->mail("All","BUILD REPORT: batch_BLAT.pl");

exit (0);


###########################
####### Subroutines #######
###########################


sub run_bsub
  {
    my $source = shift;
    my $output = shift;
    my $parameters = shift;

    my ($job_name) = $source =~ /_(\w+)/;
    $job_name = "BLAT_"."$job_name";
    my $blat       = "/software/worm/bin/blat/blat";
    my $autoace_fa = $wormbase->blat."/autoace.fa";
    #input and output in same place now
    my $EST_dir    = $wormbase->blat;
    my $output_dir = $wormbase->blat;

    my $error_dir = "$wormpub/BSUB_ERRORS";
    my $error     = "$error_dir/$output.err";

    my $bsub_cmd = "$blat $autoace_fa $EST_dir/$source ";
    $bsub_cmd .= $parameters if $parameters;
    $bsub_cmd .= " $output_dir/$output";

    $log->write_to("bsub job $bsub_cmd . . . ");

    #submit bsub command
    my $bsub_id = system("bsub -e $error -J $job_name $bsub_cmd") unless $no_bsub;

    $log->write_to("$bsub_id\n");
  }


sub split_run 
  {
    my $type = shift;
    my $shatter_tree;
    my $name_stem;
    my $source_file;
    my $opts;
    my $ESTdir = $wormbase->blat;
    my $shattered_dir;

    if( $type eq "est" ) {
      $shattered_dir = "shatteredESTs";
      $shatter_tree = "$ESTdir/$shattered_dir";
      $name_stem   = "$shatter_tree/elegansEST";
      $source_file = $wormbase->blat."/elegans_ESTs.masked";
      $opts = "";
    }
    elsif ($type eq "nematode" ) {
      $shattered_dir = "shattered_nematode";
      $shatter_tree = "$ESTdir/$shattered_dir";
      $name_stem   = "$shatter_tree/nematodeEST";
      $source_file = $wormbase->blat."/other_nematode_ESTs";
      $opts        = "-t=dnax -q=dnax"
    }
    elsif ($type eq "washu" ) {
      $shattered_dir = "shattered_washu";
      $shatter_tree = "$ESTdir/$shattered_dir";
      $name_stem   = "$shatter_tree/washu";
      $source_file = $wormbase->blat."/washu_nematode_contigs";
      $opts        = "-t=dnax -q=dnax"
    }
    elsif ($type eq "nembase" ) {
      $shattered_dir = "shattered_nembase";
      $shatter_tree = "$ESTdir/$shattered_dir";
      $name_stem   = "$shatter_tree/nembase";
      $source_file = $wormbase->blat."/nembase_nematode_contigs";
      $opts        = "-t=dnax -q=dnax"
    }
    else {
      $log->log_and_die("no type passed to split_run\n");
    }
    
	unless ( -e $shatter_tree ) {
	   mkdir( $shatter_tree ) ;
	   $wormbase->run_command("perl $ENV{'CVS_DIR'}/shatter $source_file 25000 $name_stem") and $log->log_and_die("cant shatter $source_file : $!\n");
	}
	else {
		$log->write_to("cDNA file already shattered - using existing files\n";
	}
    opendir(DIR,"$shatter_tree");

    while ( my $file = readdir( DIR ) ) {
      next unless $file =~ /\d+/;
      &run_bsub("$shattered_dir/${file}","$file.psl","$opts");
    }
  }

__END__

=pod

=head2   NAME - batch_BLAT.pl

=head1 USAGE

=over 4

=item  batch_BLAT.pl -options

=back

A wrapper script to hadle the submission of blat jobs to cbi1 cluster (or whatever cluster you are on)

If -nematode or -est is run it will run shatter to split the EST file into files with 25000 sequences each and submit one job per file.

All output is written to ~wormpub/BLAT

Error files are written to ~wormpub/BSUB_ERRORS

=over 4

=item options

  -help show this
  -debug send log only to specifed user

molecule specific runs
  -est
  -mrna
  -ncrna
  -ost
  -nematode
  -nembase
  -washu
  -embl
  -tc1
  -all all of the above

  -no_bsub dont actually submit the jobs ( for testing )

=back

END

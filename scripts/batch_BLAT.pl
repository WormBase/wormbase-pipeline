#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;

my ($help, $debug, $verbose, $est, $mrna, $ncrna, $ost, $nematode, $embl, $tc1, $all, $dump, $no_bsub);
my ($blat, $process, $virtual);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "est"        => \$est,
	    "mrna"       => \$mrna,
	    "ncrna"      => \$ncrna,
	    "ost"        => \$ost,
	    "nematode"   => \$nematode,
	    "embl"       => \$embl,
	    "tc1"        => \$tc1,
	    "all"        => \$all,
	    "dump"       => \$dump,
	    "no_bsub"    => \$no_bsub,
	    "blat"       => \$blat,
	    "process"    => \$process,
	    "virtual"    => \$virtual
	    );

if( $help ) { system ('perldoc',$0); exit(0);}

my $log = Log_files->make_build_log($debug);
my $wormpub = glob("~wormpub");

$log->log_and_die("cant do blat and process / virtual at same time\n") if ($blat and ( $process or $virtual ));

if ( $blat ) {

  # copy autoace.fa to cbi1
  $log->write_to("copying autoace.fa\n");

  if ( -e "$wormpub/BLAT/autoace.fa" and &new_file("$wormpub/BLAT/autoace.fa") ) {
    $log->write_to("using existing autoace.fa file\n\n");
  } else {
    system("scp wormsrv2:/wormsrv2/autoace/BLAT/autoace.fa $wormpub/BLAT/autoace.fa") and $log->log_and_die("cant copy autoace.fa\n");
  }

  #ESTs (~2 hours) 
  &run_bsub("elegans_ESTs.masked", "est_out.psl") if $est;

  #OSTs (5min) 
  &run_bsub("elegans_OSTs.masked", "ost_out.psl") if $ost;

  #TC1s (10min)
  &run_bsub("elegans_TC1s", "tc1_out.psl") if $tc1;

  #mRNA (2 hours)
  &run_bsub("elegans_mRNAs.masked", "mrna_out.psl", " -fine") if $mrna;

  #EMBL (5 mins) 
  &run_bsub("elegans_embl_cds", "embl_out.psl") if $embl;

  #ncRNAs (?) 
  &run_bsub("elegans_ncRNAs.masked", "ncrna_out.psl") if $ncrna;

  if ( $nematode ) {

    # splitting Nematode_ESTs
    my $shatter_dir = "${wormpub}/analysis/ESTs/shattered";
    system("shatter ${wormpub}/analysis/ESTs/other_nematode_ESTs 25000 $shatter_dir/nematodeEST") and $log->log_and_die("cant shatter $shatter_dir/nematodeEST\n");

    opendir(DIR,"$shatter_dir");

    while ( my $file = readdir( DIR ) ) {
      next unless $file =~ /\d+/;
      &run_bsub("shattered/${file}","$file.psl","-t=dnax -q=dnax");
    }
  }
}

if ( $process or $virtual ) {

  my @blat_jobs;
  push(@blat_jobs,"est")      if ( ($est)      || ($all) );
  push(@blat_jobs,"ost")      if ( ($ost)      || ($all) );
  push(@blat_jobs,"mrna")     if ( ($mrna)     || ($all) );
  push(@blat_jobs,"ncrna")    if ( ($ncrna)    || ($all) );
  push(@blat_jobs,"embl")     if ( ($embl)     || ($all) );
  push(@blat_jobs,"tc1")      if ( ($tc1)      || ($all) );
  push(@blat_jobs,"nematode") if ( ($nematode) || ($all) );

  # once the jobs have finished and been processed
  # run the main blat job
 
  foreach my $option (@blat_jobs ) {
    system("$wormpub/scripts/blat_them_all.pl -process -$option");
    # also create virtual objects
    system("$wormpub/scripts/blat_them_all.pl -virtual -$option");
  }
}


$log->mail;

exit (0);

sub run_bsub
  {
    my $source = shift;
    my $output = shift;
    my $parameters = shift;

    my ($job_name) = $source =~ /_(\w+)/;
    $job_name = "BLAT_"."$job_name";
    my $blat       = "$wormpub/bin.ALPHA/blat";
    my $autoace_fa = "$wormpub/BLAT/autoace.fa";
    my $EST_dir    = "$wormpub/analysis/ESTs";
    my $output_dir = "$wormpub/BLAT";

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

sub new_file
  {
    my $file = shift;
    my (@date) ="ls -ltr $file";

    # get the last block file in @date and split to find date
    my @line = split(/\s+/,$date[$#date]);
    
    my $file_age = sprintf("%.1f", -M $line[8]);

    return 0 if $file_age > 1;
    return 1;
}


__END__

=pod

=head2   NAME - batch_BLAT.pl

=head1 USAGE

=over 4

=item  batch_BLAT.pl -options

=back

A wrapper script to hadle the submission of blat jobs to cbi1 cluster (or whatever cluster you are on)

If -nematode is run it will run shatter to split the NematodeEST into files with 25000 sequences eachq and submit one job per file.

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
  -embl
  -tc1
  -all all of the above

  -no_bsub dont actually submit the jobs ( for testing )

=back

END

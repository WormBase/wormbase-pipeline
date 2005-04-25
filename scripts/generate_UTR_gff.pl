#!/usr/local/bin/perl5.8.0 -w
#
# prepare_UTR_GFF.pl
#
# moves GFF files to /nfs/disk100/wormpub so you can run make_UTR_GFF.pl on cbi1
use lib $ENV{'CVS_DIR'};
use strict;
use Getopt::Long;
use Wormbase;
use Log_files;

my $debug;
my ($prepare, $run);
my $scripts = shift | $ENV{'CVS_DIR'};
my $chromosome;

GetOptions (
	    'prepare' => \$prepare,
	    'run'     => \$run,
	    'debug:s' => \$debug,
	    'chromosome:s' => \$chromosome,
	    'chromosomes:s' => \$chromosome
	   );

my $log = Log_files->make_build_log($debug);

my $wormpub = glob("~wormpub");
my $datdir = "$wormpub/analysis/UTR";

my $GFFdir = "wormsrv2:/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";

my @chromosomes = $chromosome ? split(/,/,join(',',$chromosome)) : qw( I II III IV V X );

my $errors = 0;
if ( $prepare ) {
  $log->write_to("Copying GFF files from $GFFdir to $datdir \n");
  foreach my $chrom ( @chromosomes ) {

    # copy Coding_transcripts
    $errors++ if( system ("scp $GFFdir/CHROMOSOME_$chrom.Coding_transcript.gff $datdir/") );

    # copy Coding_exons
    $errors++ if( system ("scp $GFFdir/CHROMOSOME_$chrom.coding_exon.gff $datdir/") );

    # copy CDS
    $errors++ if( system ("scp $GFFdir/CHROMOSOME_$chrom.CDS.gff $datdir/") );
  }
}

$log->log_and_die("There were errors in the GFF copying so I stopping\n") unless ($errors == 0);

if( $run ) {
  $log->write_to("Submitting bsub jobs\n");
  foreach my $chrom ( @chromosomes ) {
    my $err  = "$datdir/$chrom.err.$$";
    my $bsub = "bsub -e $err \"$scripts/make_UTR_GFF.pl $chrom\"";
    print "$bsub\n";
    system("$bsub");
  }
}

$log->mail;

exit(0);

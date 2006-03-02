#!/usr/local/bin/perl5.8.0 -w
#
# GFF_method_dump.pl
#
# by Anthony Rogers
#
# Selectively dump GFF for certain acedb methods
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-03-02 16:08:47 $


use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Storable;

my ($help, $debug, $test, $quicktest, $database, @methods, @chromosomes, $dump_dir, @clones, $list );
my @sequences;
my $store;
GetOptions (
	    "help"          => \$help,
	    "debug=s"       => \$debug,
	    "test"          => \$test,
	    "store:s"       => \$store,
	    "quicktest"     => \$quicktest,
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,

	    # ive added method and methods for convenience
	    "method:s"      => \@methods,
	    "methods:s"     => \@methods,
	    "chromosomes:s" => \@chromosomes,
	    "chromosome:s"  => \@chromosomes,
	    "clone:s"       => \@clones,
	    "clones:s"      => \@clones,
	    "list:s"        => \$list
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



@methods     = split(/,/,join(',',@methods));
@chromosomes = split(/,/,join(',',@chromosomes));
@sequences = split(/,/,join(',',@clones)) if @clones;

&check_options;

my $giface = $wormbase->giface;

$database = $wormbase->autoace unless $database;
$dump_dir = "/tmp/GFF_CLASS" unless $dump_dir;

mkdir $dump_dir unless -e $dump_dir;

#make surdump_dir is writable
system("touch $dump_dir/dump") and die "cant write to $dump_dir\n";


# open database connection once
open (WRITEDB,"| $giface $database") or die "failed to open giface connection to $database\n";

foreach my $sequence ( @sequences ) {
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $command = "gif seqget $sequence +method $method; seqactions -hide_header; seqfeatures -version 2 -file $dump_dir/${sequence}_${method}.gff\n";

      print "$command";
      print WRITEDB $command;
    }
  }
  else { 
    my $command = "gif seqget $sequence; seqactions -hide_header; seqfeatures -version 2 -file $dump_dir/$sequence.gff\n";
    print "$command";
    print WRITEDB $command;
  }
}


close WRITEDB;

# remove write test
system("rm -f $dump_dir/dump");
exit(0);

sub check_options {


  unless($list or @clones ) {

    # -chromosomes
    my %chroms = qw(I 1 II 1 III 1 IV 1 V 1 X 1 MtDNA 1);
    my @chrom_keys = keys %chroms;
    unless (@chromosomes ) {
      @sequences= map("CHROMOSOME_"."$_",@chrom_keys);
      print "Dumping for all chromosomes\n";
    } 
    else {
      foreach (@chromosomes) {
	if ( $chroms{$_} ) {
	  push( @sequences,"CHROMOSOME_"."$_");
	}
	else {
	  die "$_ is not a valid chromosome\n";
	}
      }
    }
  }
  
  &process_list if $list;

  # -database
  if ( $database ){
    if( -e "$database" ) {
      if( -e "$database/wspec/models.wrm" ) {
	print "$database OK\nDumping @methods for chromosomes @chromosomes\n";
	return;
      }
    }
  }
  else {
    die "You must enter a valid database\n";
  }
  die "$database is not a valid acedb database\n";
}



sub process_list
  {
    open(LIST,"<$list") or die "bad list $list\n";
    while(<LIST>) {
      chomp;
      push(@sequences,$_);
    }
  }

=pod 

=head1 NAME - GFF_method_dump.pl

=head2 USAGE

=over 4

This script will GFF dump specified methods from a database

It is use by dump_gff_batch.pl so if you change it make sure it is still compatible !

=back

=item MANDATORY ARGS:

=over 4

-methods     Methods to dump eg curated,history (Comma separated list)

=back

=item OPTIONAL ARGS:

=over 4

-database    Database to dump from ( default /wormsrv2/autoace )

-chromosomes Chromosomes to dump as comma separated list eg I,II,X ( defaults to all )

=back

=item OUTPUT:

=over 4

A separate file is written for each method for each chromosome and is named 

CHROMOSOME_($chrom)_($method).gff

=back

=item EXAMPLES:

=over 4

GFF_method_dump.pl -database wormsrv2/autoace -chromosomes II,V -method curated,TRANSCRIPT

will GFF dump separate curated and TRANSCRIPT files for both chromosomes II and V ie

  CHROMOSOME_II_curated.gff
  CHROMOSOME_II_TRANSCRIPT.gff
  CHROMOSOME_V_curated.gff
  CHROMOSOME_V_TRANSCRIPT.gff

=back

=item REQUIREMENTS:

=over 4

Access to ~acedb for giface binary

=back

=item WARNING:

=over 4

At time of writing the version of acedb (4.9y) adds extra data to the GFF output. Eg if you ask for method curated you also get a load of ALLELES that have no method, so the GFF files produced should be treated with caution and may need post-processing.

=back

=cut

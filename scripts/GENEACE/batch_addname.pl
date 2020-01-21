#!/usr/local/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use NameDB_handler;
use Wormbase;
use Getopt::Long;
use Log_files;

=pod

=head batch_addname.pl

=item Options:

  -file      file containing list of GeneIDs and CGC name eg WBGene00008040 ttr-5
  -species   what species these are for - default = elegans
  -force     bypass CGC name validation check eg to add Cbr-cyp-33E1; use with care!
  -type      allows you to specify Sequence names to add (defaults to CGC if not speciefied)

e.g. perl batch_addname.pl -file genenames.txt -species briggsae

=cut

my ($help, $debug, $verbose, $store, $wormbase);
my ($file, $species, $force, $type);

GetOptions(
	   "help"       => \$help,
	   "debug=s"    => \$debug,
	   "verbose"    => \$verbose,
	   "store:s"    => \$store,
	   'file:s'     => \$file,
	   'species:s'  => \$species,
	   'force'      => \$force,
	   'type:s'     => \$type,
	  ) or die;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( 
			    -debug    => $debug,
			    -organism => $species,
			   );
}

$species = 'elegans' unless $species;
$type = 'CGC' unless $type;

if (($type ne "CGC") && ($type ne "Sequence")) {
  die "can't open $file : $!\n$type is not a valid type (CGC/Sequence)";
}

my $log = Log_files->make_log("NAMEDB:$file", $ENV{USER});

$log->write_to("loading $file to Nameserver\n\n");
$log->write_to("FORCE mode is ON!\n\n") if $force;

my $db = NameDB_handler->new($wormbase);

#open file and read

open (FILE,"<$file") or die "can't open $file : $!\n";
my $method = defined($force) ? 'force_name' : 'add_name';
my $count=0;
while(<FILE>) {
    my($id,$name) = split;
    my $success = $db->$method($id,$name,$type);
    if (defined $success) {
      $log->write_to("$id\t$name\tok\n");
      print "$id\t$name\tok\n";
    }
    else {
      $log->write_to("$id\t$name\tFAILED NAMEDB assignment\n");
      print "$id\t$name\tFAILED NAMEDB assignment\n";
    }
  $count++;
}

$log->write_to("=======================\nprocessed $count genes\n");
$log->mail;

















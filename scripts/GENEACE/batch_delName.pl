#!/usr/bin/env perl
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use Wormbase;
use NameDB_handler;
use Getopt::Long;
use Log_files;

use warnings;
use strict;

=pod

=head batch_delName.pl

=item Options:

  -file      file containing list of GeneIDs and CGC names eg WBGene00008040 ttr-5

e.g. perl batch_delName.pl -file genenames.txt

=cut

my ($help, $debug, $verbose, $store, $wormbase, $species);
my ($file);
GetOptions(
           "help"       => \$help,
           "debug=s"    => \$debug,
           "verbose"    => \$verbose,
           "store:s"    => \$store,
           'species:s'  => \$species,
	   'file:s'     => \$file,
	  ) or die;


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( 
                            -debug    => $debug,
                            -organism => $species,
                           );
}

my $log = Log_files->make_log("NAMEDB:$0:$file", $ENV{'USER'});


$log->write_to("loading $file to Nameserver\n\n");

my $db = NameDB_handler->new($wormbase);

# open file and read
# file needs to be in such a format:
#
# WBGene007 flunk-1
# WBGene123 nope-2

open (FILE,"<$file") or die "can't open $file : $!\n";
my $count=0;
while(<FILE>) {
    my($id,$name) = split;
    my $success = $db->delName($name);
    my $msg = defined $success ? 'deleted' : 'FAILED to delete';
    $log->write_to("$msg CGC name $name for $id\n");

    $count++;
}

$log->write_to("=======================\nprocessed $count genes\n");
$log->mail;

#!/usr/bin/perl
#
# dumps gene descriptions
#
# Options:
#   -format  <record || tab> (defaults to record)
#   -species <name>        WormBase species name
#   -store <storable file> pass a stored config
#   -debug <user>          send log mails to user
#   -test                  use the test database
#   -database <path>       use a custom database

use strict;

use Getopt::Long;
use IO::File;
use Storable;

use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/preFTPdumps/perSpecies";

use Wormbase;
use Log_files;
use Dumper;


my ($species,$format,$store,$debug,$test,$database,$outfile);
GetOptions(
     'species=s'  => \$species,
     'format=s'   => \$format,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'outfile=s'  => \$outfile,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);
$database = $wormbase->autoace if not defined $database;

$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $database) or $log->log_and_die("Could not connect to $database\n");

$outfile = $wormbase->reports . '/functional_descriptions.txt'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

my $separator = "=\n";
my $no_entry = 'none available';
$format ||= 'record';

$log->write_to("dumping ${\$wormbase->long_name} functional descriptions\n");
print $of "# ${\$wormbase->long_name} functional descriptions\n";
print $of "# WormBase version: " . $dbh->version . "\n";
print $of '# Generated: '.&get_date."\n";

# header
unless($format eq 'record') {
    print $of join("\t",
     qw/gene_id public_name molecular_name concise_description provisional_description detailed_description automated_description gene_class_description/
                  ),"\n";
}

my $i = $dbh->fetch_many(-query=>qq{find Gene Species="${\$wormbase->long_name}"});
while (my $gene = $i->next) {
    my $name = $gene->Public_name || 'not known';
    my $molecular_name = $gene->Molecular_name || 'not known';
    
    # Fetch each and any of the possible brief identifications
    my $concise  = $gene->Concise_description     || $no_entry;
    my $prov     = $gene->Provisional_description || $no_entry;
    my $detailed = $gene->Detailed_description    || $no_entry;
    my $automated = $gene->Automated_description  || $no_entry;
    my $gene_class = $gene->Gene_class;
    my $gene_class_description = $gene_class ? $gene_class->Description : 'not known';

    if ($format eq 'record') {
	      print $of join("\t",$gene,$name,$molecular_name) . "\n";
	      print $of rewrap("Concise description: $concise\n");
	      print $of rewrap("Provisional description: $prov\n");
	      print $of rewrap("Detailed description: $detailed\n");
	      print $of rewrap("Automated description: $automated\n");
	      print $of rewrap("Gene class description: $gene_class_description\n");
	      print $of $separator;
    } else {
	     print $of join("\t",$gene,$name,$molecular_name,$concise,$prov,$detailed,$automated,$gene_class_description),"\n";
    }
}

$dbh->close();
$of->close;
$log->mail;

exit(0);

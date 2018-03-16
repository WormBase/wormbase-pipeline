#!/usr/bin/env perl
#
# perl $CVS_DIR/GENEACE/fix_geneace_orthologs.pl -database ~/DATABASES/geneace/ -outfile Live_orthocheck.out
#
# Script to fix the analysis/evidence meta data that sometimes gets lost when manually moving orthologs/when genes are merged split/manually generated and loaded.
#
# Can use the -species option to only check the genes in that core species, but defaults to all Live genes with Orthologs
#


use strict;
use lib $ENV{"CVS_DIR"};
use Wormbase;
use Ace;
use Ace::Object;
use Carp;
use Getopt::Long;
use GENEACE::Geneace;
use File::Path;

my ($species, $database, $out, $debug, $wormbase);
GetOptions ("species:s"        => \$species,
	    "database:s" => \$database,
	    "outfile:s" => \$out,
	    "debug:s" => \$debug
    );

  $wormbase = Wormbase->new( -debug    => $debug,
                             -organism => $species,
      );

my $full_species_name;
if ($species){
    $full_species_name = $wormbase->full_name;
}

my $db = Ace->connect(-path => $database,) || die("Ace::Error");
my $query;
if ($species) {
    $query = "find Gene Species=\"${full_species_name}\"; Ortholog";
}
else {
    $query = "find Gene; Live; Ortholog";
}
my $genes = $db->fetch_many(-query => $query);
open (ACE,">$out") || die("Out::Error");
my $count;

while (my $gene = $genes->next){
  foreach my $o ($gene->Ortholog){
    # map {print "$gene ",$gene->Species," => $o ",$o->right," ($_)\n"}$o->right(3)
    unless ($o->right(3)){
     reverseO($o,$gene);
    }
  }
}
print "$count genes could not be fixed.\n";
close(ACE);
exit(0);

###############################################################################################

sub reverseO {
   my($g,$id) = @_;
   foreach my $or ($g->Ortholog){
     if ("$or" eq "$id"){
        my @e = $or->right(3);
        print ACE "Gene : $id\n";
	print "//ERROR: $g::$id has missing ortholog data\n" unless @e;
	$count++;
        map{ print ACE "Ortholog $g \"${\$g->Species}\" From_analysis $_\n"} @e;
        print ACE "\n";
     }
   }
}

#!/software/bin/perl -w
#
# promote_omim_data.pl                           
# 
# by Paul Davis                         
#
# This script promoted the OMIM disease data to the level of the gene.
# Script extended to also populate the Potential_model_for Human disease data utilising the DO ontology.
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2015-02-10 14:52:03 $      

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $noload, $database, $acefile, $flatfile, $obofile, $species);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
            "species=s"  => \$species,
	    "database:s" => \$database,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noload"     => \$noload,
	    "acefile:s"  => \$acefile,
	    "obofile:s"  => \$obofile,
            "flatfile:s" => \$flatfile,
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);


my $log = Log_files->make_build_log($wormbase);

my $tace            = $wormbase->tace;        # TACE PATH
my $WS_name         = $wormbase->get_wormbase_version_name();
my $species_name    = $wormbase->full_name;

$acefile = $wormbase->acefiles . "/omim_db_data.ace" if not defined $acefile;
$flatfile = $wormbase->ontology . "/disease_association.by_orthology." . $wormbase->get_wormbase_version_name . ".tsv.txt" 
    if not defined $flatfile;
$obofile = $wormbase->primaries . "/citace/temp_unpack_dir/home/citace/Data_for_${WS_name}/Data_for_Ontology/disease_ontology.${WS_name}.obo" 
    if not defined $obofile;
$database = $wormbase->autoace if not defined $database;

open (my $outfh, ">$acefile") or $log->log_and_die("Cannot write ace file $acefile");
open (my $outflat, ">$flatfile") or $log->log_and_die("Cannot write flat file $flatfile");

# Additional information is now required to add the Disease Ontology terms to genes where an OMIM disease ID has 
# been identified by human protein orthology.
my $omim2do = &gatherDOdata($obofile);

my (%gene2omim, %gene2do);

my $tmdef = &get_table_maker_def();
my $command = "Table-maker -p $tmdef\nquit\n";
$log->write_to("\nRetrieving OMIM, using Table-maker and query ${tmdef}...\n");
  
open (TACE, "echo '$command' | $tace $database | ") || die "Cannot query acedb. $command  $tace\n";
while (<TACE>) {
  chomp;
  s/\"//g;

  if (/^WBGene\d+/) {
    my ($gene, $ortholog, $species,  $analysis, $database, $database_field, $database_acc) = split(/\t/, $_);

    next if $ortholog !~ /^HGNC:/;

    $gene2omim{$gene}->{$database_field}->{$database_acc}++; 

    if (exists $omim2do->{$database_acc}) {
      foreach my $doid (keys %{$omim2do->{$database_acc}}) {
        $gene2do{$gene}->{$doid}->{$database_acc}->{$ortholog}->{$analysis} = 1;
      }
    }
  }
}

my %stats;

#
# Database lines; for these, will conservatively only add them if they have support from >= 2 methods
#
foreach my $gene (sort keys %gene2omim) {
      
  print $outfh "\nGene : \"$gene\"\n";
  my $gene_had_data = 0;

  foreach my $tp (keys %{$gene2omim{$gene}}) {
    foreach my $oid (sort keys %{$gene2omim{$gene}->{$tp}}) {
      if ($gene2omim{$gene}->{$tp}->{$oid} > 1) { 
        print $outfh "Database OMIM $tp $oid\n";
        $stats{database}->{$tp}++;
        $gene_had_data = 1;
      }
    }
  }

  $stats{database}->{genecount}++ if $gene_had_data;
}

#
# Potential_model_for lines
#
foreach my $gene (sort keys %gene2do) {
  $stats{doterm}->{genecount}++;

  print $outfh "\nGene : \"$gene\"\n";

  foreach my $do_id (keys %{$gene2do{$gene}}) {
    $stats{doterm}->{doterms}++;
    
    my (@list, %orths);
    foreach my $omim_id (sort keys %{$gene2do{$gene}->{$do_id}}) {
      push @list, "OMIM:$omim_id";
      foreach my $orth (keys %{$gene2do{$gene}->{$do_id}->{$omim_id}}) {
        $orths{$orth} = 1;
        my @methods = sort keys %{$gene2do{$gene}->{$do_id}->{$omim_id}->{$orth}};

        printf $outflat "%s\t%s\t%s\t%s\t%s\n", $gene, $do_id, "OMIM:$omim_id", $orth, join("|", @methods);
      }
      push @list, sort keys %orths;
    }
    my $list_str = join(",", @list);
    
    print $outfh "Potential_model \"$do_id\" \"Homo sapiens\" Inferred_automatically \"Inferred by orthology to human genes with OMIM annotation ($list_str)\"\n";
  }
}
close($outfh) or $log->log_and_die("Could not close acefile after writing\n");
close($outflat) or $log->log_and_die("Could not close flat file after writing\n");


if ($noload){
  $log->write_to("Output NOT loaded into ".$wormbase->autoace."\n");
} else {
  $log->write_to("loading $acefile to ".$wormbase->autoace."\n");
  $wormbase->load_to_database($wormbase->autoace,$acefile,'promote_omim_data.pl', $log);
}

# Close log files and exit
$log->write_to(sprintf("Wrote Database lines for %d genes, with %d gene entries and %d disease entries\n", 
                       $stats{database}->{genecount},
                       $stats{database}->{gene},
                       $stats{database}->{disease}));
$log->write_to(sprintf("Wrote Potential_model for %d genes, with %d DO term references\n", 
                       $stats{doterm}->{genecount},
                       $stats{doterm}->{doterms}));

$log->mail();
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub gatherDOdata {
  my ($obo_file) = @_;

  open (my $obo_fh, "<$obo_file") or $log->log_and_die("Can't open OBO file: $obo_file\n");
  my ($doid, $omimid, %omim2do);

  $log->write_to("Retrieving DO_term data, using the citace obo file...\n");
  while (<$obo_fh>) {
    # id: DOID:0050631
    if (/^id:\s+(DOID:\d+)/) {
      $doid = $1;
    }
    #xref: OMIM:203100
    if (/^xref:\s+OMIM:(\d+)/) {
      $omimid = $1;
      $omim2do{$omimid}->{$doid} = 1;
    }
  }
  my $countomim2do = (keys %omim2do);
  $log->write_to("Collected $countomim2do  OMIM::DO_terms\nFinished getting DO data...");
  if ($countomim2do < 4000) {
    $log->log_and_die("ERROR: $countomim2do OMIM::DO_term mappings collected which is less than the minimum 4000, this is bad!!\n");
  }
  return \%omim2do;
}


sub get_table_maker_def {

  my $def = '/tmp/omim.def';
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;

Sortcolumn 1

Colonne 1 
Width 16 
Optional 
Visible 
Class 
Class Gene
From 1
Condition Ortholog AND Species = "$species_name"


Colonne 2 
Width 40 
Mandatory 
Visible 
Class 
Class Gene 
From 1 
Tag Ortholog  

Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Species 
Right_of 2 
Tag  HERE 
Condition "Homo sapiens"

Colonne 4 
Width 12 
Mandatory 
Visible 
Class 
Class Analysis 
Right_of 3 
Tag  HERE  # From_analysis 

Colonne 5 
Width 12 
Mandatory 
Visible 
Class 
Class Database 
From 2 
Tag Database  
Condition "OMIM"
 
Colonne 6 
Width 12 
Mandatory 
Visible 
Class 
Class Database_field 
Right_of 5 
Tag HERE   
 
Colonne 7 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 6 
Tag HERE   

// End of these definitions
END

  print TMP $txt;
  close TMP;

  return $def;
}


sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}




__END__

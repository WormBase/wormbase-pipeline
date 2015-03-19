#!/usr/bin/env perl
#
# transfer_interpro_GO_terms
#
# Transfers GO_term annotations from protein domains up to the gene
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2015-03-19 12:17:46 $      

use strict;
use warnings;

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Ace;

#
# Still need to have blacklists, because the taxonomy constraints are
# no-where near complete
#
my @BAD_PHRASES =  (
  'sporulation',
  'forespore',
  'photosynthesis',
  'photosynthetic',
  'chlorophyll',
    );

my %BAD_ACCS =  (
  'GO:0009772' => 1,  # photosynthetic electron transport in photosystem II
  'GO:0045282' => 1,  # plasma membrane succinate dehydrogenase complex (only_in_taxon Bacteria)
  'GO:0009288' => 1,  # bacterial-type flagellum
  'GO:0007391' => 1,  # dorsal closure (only_in_taxon Insecta)
  'GO:0009103' => 1,  # lipopolysaccharide biosynthetic process (only in prokaryotes)
    );


##############################
# Script variables (run)     #
##############################

my ($debug, $store, $test, $species, $acefile, $noload, $database,
    $motif, $tmhmm );

##############################
# command-line options       #
##############################

GetOptions (
  "debug=s"           => \$debug,
  "test"              => \$test,
  "store:s"           => \$store,
  "database:s"        => \$database,
  "acefile:s"         => \$acefile,
  "noload"            => \$noload,
    );


# recreate configuration 
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
} else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test );
}

my $log = Log_files->make_build_log($wormbase);
my $tace      = $wormbase->tace;      # tace executable path

$database = $wormbase->autoace if not defined $database;
$acefile = $wormbase->acefiles."/propagated_GO_terms.ace" if not defined $acefile;

my %go_accs_to_terms = %{&get_GO_acc_to_term()};

my @annots = &transfer_from_motif();

my $annot_id = &get_next_GO_annotation_id();

$log->log_and_die("Did not find any annotations to transfer. Something probably went wrong\n")
    if not @annots;

open (my $acefh,">$acefile") or $log->log_and_die("cant open $acefile :$!\n");
foreach my $annot (@annots) {
  my ($gene, $go_id, $motif) = @$annot;
  printf $acefh "\nGO_annotation : \"%d\"\n", $annot_id++;
  print $acefh "Gene \"$gene\"\n";
  print $acefh "GO_term \"$go_id\"\n";
  print $acefh "GO_code \"IEA\"\n";
  print $acefh "Motif \"$motif\"\n";
  print $acefh "GO_reference \"Gene Ontology Consortium\" \"GO_REF\" \"0000002\"\n";
  
}
close($acefh) or $log->log_and_die("Could not cleanly close $acefile\n");

if (not $noload) {
  $wormbase->load_to_database($database,$acefile,'transfer_GO_terms', $log);
}

$log->mail();
exit(0);

my (%blacklist_rejections, %missing_go_def_rejections, %bad_phrase_rejections);

########################################################################################
sub transfer_from_motif {
  my (%res);
  my $def = &write_motif_def();
  
  my $query = $wormbase->table_maker_query($database, $def);
  while(<$query>) {
    chomp;
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;

    my($motif, $GO, $protein, $cds, $gene) = split(/\t/, $_);

    next if not defined $gene or not defined $cds or not defined $GO or not  defined $motif;
  
    if (not exists $go_accs_to_terms{$GO}) {
      $missing_go_def_rejections{$GO}++;
      next;
    }
    if (exists $BAD_ACCS{$GO}) {
      $blacklist_rejections{$GO}++;
      next;
    }  
    if (grep { $go_accs_to_terms{$GO} =~ /\b$_\$b/ } @BAD_PHRASES) {
      $bad_phrase_rejections{$GO}++;
      next;
    }
    
    $res{$gene}->{$GO}->{$motif} = 1;
  }
  unlink $def;
  
  foreach my $k (sort keys %missing_go_def_rejections) {
    $log->write_to(sprintf("Rejected %d annotations for %s because it has no definition in database\n", 
                           $missing_go_def_rejections{$k}, 
                           $k));
  }
  foreach my $k (sort keys %blacklist_rejections) {
    $log->write_to(sprintf("Rejected %d annotations for %s because it is on the blacklist\n", 
                           $blacklist_rejections{$k}, 
                           $k));
  }
  foreach my $k (sort keys %bad_phrase_rejections) {
    $log->write_to(sprintf("Rejected %d annotations for %d because definition contains a bad phrase\n",
                           $bad_phrase_rejections{$k}, 
                           $k));
  }
  
  my @annotations;
  
  foreach my $gene (sort keys %res) {
    foreach my $go_acc (sort keys %{$res{$gene}}) {
      foreach my $motif (sort keys %{$res{$gene}->{$go_acc}}) {
        push @annotations, [$gene, $go_acc, $motif];
      }
    }
  }
  $log->write_to(sprintf("Propagated %d annotations\n", scalar(@annotations)));
  
  return @annotations;
}

#################################################
sub transfer_from_tmhmm {

  my $mydb = Ace->connect(-path => $database ) 
      or $log->log_and_die("Connection failure " . Ace->error . "\n");

  my @annotations;

  my $query = "FIND Protein WHERE Species = \"".$wormbase->full_name."\"" 
      . " WHERE Feature AND NEXT = \"Tmhmm\";" 
      . " FOLLOW Corresponding_CDS;"
      . " FOLLOW Gene";
  my $genes = $mydb->fetch_many(-query => $query);
  while(my $gene = $genes->next){
    push @annotations, [$gene, "GO:0016021", "CBS", "TMHMM"];
  }

  $mydb->close or $log->log_and_die("Could not close Ace connection\n");

  return @annotations;
}

########################################################################################
sub get_next_GO_annotation_id {
  my $def = &write_GO_annot_def();
  my $query = $wormbase->table_maker_query($database, $def);

  my $max_id;

  while(<$query>) {
    chomp;
    s/\"//g;
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;

    my ($annot_id, $go_term) = split(/\t/, $_);
    next if $go_term !~ /^GO:/;

    $max_id = $annot_id if not defined $max_id or $max_id < $annot_id;
  }

  return $max_id + 100;
}


######################################################################
sub get_GO_acc_to_term {
  my %terms;

  # get the GO terms
  my $term_def = &write_GO_def();
  my $term_query = $wormbase->table_maker_query($database,$term_def);

  while(<$term_query>) {
    chomp;
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;

    my @data = split("\t",$_);
    my ( $GO, $term) = @data;
    $terms{$GO} = $term;
  }

  return \%terms;
}

############################################################
sub write_GO_def {
  my $txt = <<"ENDE";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class GO_term 
From 1 
  
Colonne 2 
Width 120
Mandatory
Visible 
Class 
Class Text 
From 1 
Tag Term

ENDE

  return &write_tm_def("GO_acc_to_term", $txt);
}

############################################################
sub write_motif_def {
  my $txt = <<"ENDE";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Motif 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class GO_term 
From 1 
Tag GO_term 
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Protein 
From 1 
Tag Pep_homol 
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Class 
Class CDS 
From 3 
Tag Corresponding_CDS 
 
Colonne 5 
Width 12 
Mandatory 
Visible 
Class 
Class Gene 
From 4 
Tag Gene 

ENDE

  return &write_tm_def("motif_to_go", $txt);
}

sub write_GO_annot_def {
  my $txt = <<"ENDE";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class GO_annotation 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class GO_term 
From 1 
Tag GO_term 
 
ENDE

  return &write_tm_def("current_go_annots", $txt);
}

###################################
sub write_tm_def {
  my ($fname, $string) = @_;
  
  my $file = "/tmp/$fname.def";

  open(my $fh, ">$file") or $log->log_and_die("Could not open $fname for writing\n");
  print $fh $string;
  close($fh) or $log->log_and_die("Could not close $fname after writing\n");

  return $file;
}


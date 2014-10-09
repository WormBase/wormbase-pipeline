#!/usr/bin/env perl
#
# inherit_GO_terms.pl
#
# map GO_terms to ?Sequence objects from ?Motif and ?Phenotype
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-10-09 16:04:08 $      

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
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
    $motif, $phenotype, $tmhmm, $trusted_papers_only);

##############################
# command-line options       #
##############################

GetOptions (
  "debug=s"           => \$debug,
  "test"              => \$test,
  "species:s"         => \$species,
  "store:s"           => \$store,
  "database:s"        => \$database,
  "acefile:s"         => \$acefile,
  "noload"            => \$noload,
  "phenotype"         => \$phenotype,
  "tmhmm"             => \$tmhmm,
  "motif"             => \$motif,
  "trustedpapersonly" => \$trusted_papers_only, 
    );


# recreate configuration 
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
} else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test, 
                             -organism => $species );
}

my $log = Log_files->make_build_log($wormbase);
my $tace      = $wormbase->tace;      # tace executable path

$database = $wormbase->autoace if not defined $database;
$acefile = $wormbase->acefiles."/inherited_GO_terms.ace" if not defined $acefile;

open (my $acefh,">$acefile") or $log->log_and_die("cant open $acefile :$!\n");

my %go_accs_to_terms = %{&get_GO_acc_to_term()};

&transfer_from_motif()          if ($motif);
&transfer_from_rnai_phenotype() if ($phenotype);
&transfer_from_tmhmm()          if ($tmhmm);

close($acefh) or $log->log_and_die("Could not cleanly close $acefile\n");

if (not $noload) {
  $wormbase->load_to_database($database,$acefile,'inherit_GO_terms', $log);
}

$log->mail();
exit(0);

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
  
    next if &is_bad_term($GO, $cds, $motif);
    
    $res{Gene}->{$gene}->{$GO}->{$motif} = 1;
    # $res{CDS}->{$cds}->{$GO}->{$motif} = 1;
  }

  foreach my $obj_type (sort keys %res) {
    foreach my $obj (sort keys %{$res{$obj_type}}) {
      print $acefh "\n$obj_type : \"$obj\"\n";
      foreach my $go_acc (sort keys %{$res{$obj_type}->{$obj}}) {
        foreach my $motif (sort keys %{$res{$obj_type}->{$obj}->{$go_acc}}) {
          print $acefh "GO_term \"$go_acc\" IEA inferred_automatically \"$motif\"\n";
        }
      }
    }
  }

  unlink $def;
}

#################################################
sub transfer_from_tmhmm {

  my $mydb = Ace->connect(-path => $database ) 
      or $log->log_and_die("Connection failure " . Ace->error . "\n");

  my $query = "FIND Protein WHERE Species = \"".$wormbase->full_name."\"" 
      . " WHERE Feature AND NEXT = \"Tmhmm\";" 
      . " FOLLOW Corresponding_CDS;"
      . " FOLLOW Gene";
  my $genes = $mydb->fetch_many(-query => $query);
  while(my $gene = $genes->next){
    print $acefh "\nGene : ".$gene->name."\n";
    print $acefh "GO_term \"GO:0016021\"\tIEA\tInferred_automatically\t\"CBS:TMHMM\"\n";
  }

  $mydb->close or $log->log_and_die("Could not close Ace connection\n");
}

##################################################
sub transfer_from_rnai_phenotype {
   
  my (%phen_names, %include_list);

  my $mydb = Ace->connect( -path=>$database ) or 
      $log->log_and_die("Connection error: " . Ace->error . "\n");

  my $phenos = $mydb->fetch_many(-class => 'Phenotype');
  while(my $pheno = $phenos->next){
    if ($pheno->Primary_name) {
      $phen_names{$pheno->name} = $pheno->Primary_name;
    }
  }
  $mydb->close or $log->log_and_die("Could not close Ace connection\n");

  READARRAY: while (<DATA>) {
    chomp;
    $include_list{$_} =1;
  }

  my $def = &write_phenotype_def();
  my $tm_query = $wormbase->table_maker_query($database,$def);
  while(<$tm_query>) {
    chomp;
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;

    my ( $rnai, $paper, $gene, $primary, $phen_id, $go ) = split("\t",$_);

    if (not defined $phen_id or
        not defined $go or
        not defined $gene or
        not defined $paper or
        $primary ne 'RNAi_primary' or
        $rnai !~ /WBRNAi/ or
        $phen_id !~ /^WBPhen/ or
        $paper !~ /^WBPaper/) {
      next;
    }
    next if $trusted_papers_only and not exists $include_list{$paper};

    my $phenotype = (exists $phen_names{$phen_id}) ? $phen_names{$phen_id} : $phen_id; 

    print $acefh "\nGene : \"$gene\"\n";
    print $acefh "GO_term \"$go\" IEA Inferred_automatically \"$phenotype ($phen_id|$rnai)\"\n";
  }

  unlink $def;
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

#######################################################################
sub is_bad_term {
  my ($go_acc, $gene, $entity) = @_;
  
  if (not exists $go_accs_to_terms{$go_acc}) {
    $log->write_to("Could not find term for GO-acc $go_acc, so will not be transferred to $gene via $entity\n");
    return 1;
  }
  
  if (exists $BAD_ACCS{$go_acc}) {
    $log->write_to("GO-acc $go_acc is on blacklist will not be transferred to $gene via $entity\n");
    return 1;
  }
  
  foreach my $phrase (@BAD_PHRASES) {
    my $go_term = $go_accs_to_terms{$go_acc};
    
    if ($go_term =~ /\b$phrase\b/) {
      $log->write_to("The invalid term '$phrase' was found in the term '$go_term' of GO-acc $go_acc and will not be transferred to $gene via $entity\n");
      return 1;
    }
  }
  
  return 0;
}


############################################# 
sub write_phenotype_def {
  my $txt = <<"ENDE";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class RNAi 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Paper 
From 1 
Tag Reference 
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Text 
Right_of 3 
Tag  HERE  # Inferred_automatically 
Condition RNAi_primary
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Phenotype 
From 1 
Tag Phenotype 
 
Colonne 6 
Width 12 
Mandatory 
Visible 
Class 
Class GO_term 
From 5 
Tag GO_term 

ENDE

  return &write_tm_def("Phenotype_to_GO", $txt);
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

###################################
sub write_tm_def {
  my ($fname, $string) = @_;
  
  my $file = "/tmp/$fname.def";

  open(my $fh, ">$file") or $log->log_and_die("Could not open $fname for writing\n");
  print $fh $string;
  close($fh) or $log->log_and_die("Could not close $fname after writing\n");

  return $file;
}


__DATA__
WBPaper00004402
WBPaper00004403
WBPaper00004540
WBPaper00004651
WBPaper00004769
WBPaper00005599
WBPaper00005654
WBPaper00006395
WBPaper00024497
WBPaper00024925
WBPaper00025054
WBPaper00026763
WBPaper00005736
WBPaper00026593
WBPaper00028783
WBPaper00029254
WBPaper00030951

__END__

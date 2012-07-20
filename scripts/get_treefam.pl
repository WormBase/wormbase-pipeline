#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Log_files;
use Getopt::Long;
use strict;
use Storable;
use DBI;

# 
# Perl script treefam_worm2.pl
# Written by Avril Coghlan (alc@sanger.ac.uk).
# 18-JAN-06. 
#
# For the TreeFam project.
#
# This perl script connects to the MYSQL database of
# TreeFam families and prints out a list of the C. elegans and
# C. briggsae genes in TreeFam families.
# 
# The output has the format:
# WORM_GENE NUMBER_OF_FAMILIES FAMILIES
# where WORM_GENE is the gene name, eg., R13F6.4 for a C. elegans
#                 gene or CBG100063 for a C. briggsae gene,
#       NUMBER_OF_FAMILIES is the number of TreeFam families that 
#                 WORM_GENE appears in,
#       FAMILIES is a list of the families that WORM_GENE is in. 
#
# The command-line format is:
# % perl <treefam_worm2.pl> 
#
#------------------------------------------------------------------#



my ($debug, $test,$store, $no_load, $species, $acefile, %worm, %family); 

my $treefam_host = "db.treefam.org";
my $treefam_port = 3308;
my $treefam_user = "anonymous";
my $treefam_pass = "";
my $treefam_dbnameprefix = "treefam_";
my $treefam_version;

GetOptions (
  "debug:s"        => \$debug,
  "test"           => \$test,
  "store:s"        => \$store,
  "species:s"      => \$species,
  "noload"         => \$no_load,
  "acefile=s"      => \$acefile,
  "dbhost=s"       => \$treefam_host,
  "dbuser=s"       => \$treefam_user, 
  "dbpass=s"       => \$treefam_pass, 
  "dbport=s"       => \$treefam_port,
  "dbnameprefix=s" => \$treefam_dbnameprefix,
  "dbversion=s"    => \$treefam_version,
    ) or die "Bad options\n";



my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
                             -organism => $species
      );
}
my $log = Log_files->make_build_log($wormbase);

my $species_tax_id = $wormbase->ncbi_tax_id;

$acefile = $wormbase->acefiles."/treefam.ace" if not defined $acefile;

if (not defined $treefam_version) {
  $treefam_version = "7";
} elsif (uc($treefam_version) eq 'LATEST') {
  $treefam_version = &get_latest_treefam_version();
} 


#------------------------------------------------------------------#

# GET THE LONG NAMES OF THE TREEFAM GENES FROM THE MYSQL DATABASE:

my $database  = $treefam_dbnameprefix . $treefam_version;

$log->write_to("connecting to treefam database : \tmysql:$database:$treefam_host:$treefam_port\n");
my $dbh  = DBI->connect_cached("dbi:mysql:$database:$treefam_host:$treefam_port", 
                               $treefam_user, 
                               $treefam_pass) or $log->log_and_die("Could not connect to database\n");

my $table_w  = 'genes';
# THIS TABLE HAS THE ID AND DISPLAY ID. SOMETIMES THE DISPLAY ID IS
# THE UNIPROT NAME, SOMETIMES THE ID IS:31234
my $st  = "SELECT ID, DISP_ID  from $table_w WHERE TAX_ID = $species_tax_id";
my $sth = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
$sth->execute or die "Cannot execute the query: $sth->errstr";
while ((my @array) = $sth->fetchrow_array) {
  my ($id, $disp_id) = @array;
  if ($treefam_version <= 7) {
    $worm{$id} = $id;
  } else {
    $worm{$disp_id} = $id;
  }

}
$dbh->disconnect();

#------------------------------------------------------------------#

if (keys %worm) {
  $dbh = DBI->connect_cached("dbi:mysql:$database:$treefam_host:$treefam_port", 
                             $treefam_user, 
                             $treefam_pass)  or $log->log_and_die("Could not connect to database\n");
  
  # FIRST READ IN TREEFAM-A AND THEN TREEFAM-B:  
  for (my $i = 1; $i <= 2; $i++) {
    if ($i == 1) { # LOOK AT TREEFAM-A:
      $table_w = 'fam_genes where FAM_TYPE="A"';
    } elsif ($i == 2) { # LOOK AT TREEFAM-B:
      $table_w = 'fam_genes where FAM_TYPE="B"';
    }
    
    # THE FIRST THREE COLUMNS IN THE TABLE famB_gene/famA_gene ARE THE TRANSCRIPT NAME, FAMILY NAME AND WHETHER THE
    # TRANSCRIPT IS IN THE SEED/FULL TREE:
    # eg., ENSMUST00000049178.2 TF105085 FULL
    
    my $st  = "SELECT ID, AC, FLAG FROM $table_w AND FLAG=\"FULL\""; 
    my $sth = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
    $sth->execute or die "Cannot execute the query: $sth->errstr";
    
    while ((my @array) = $sth->fetchrow_array) {
      my ($disp_id, $ac, $flag) = @array;
  
      if ($worm{$disp_id}) {
        my $id = $worm{$disp_id};
        # REMEMBER THE FAMILY THAT THIS WORM GENE IS IN:
        if (!$family{$id}){
          $family{$id} = $ac;
        }
        # as of treefam_4 they include an experimental clustering method 
        # used to populate treefam-c families.  These all begin TF5 and we only
        # want these if there is no A or B family
        else {
          $family{$id} = "$family{$id},$ac" unless($ac =~ /^TF5/);
        }
      }          
    }
  }

  $dbh->disconnect();
}

#------------------------------------------------------------------# 

# PRINT OUT A LIST OF THE WORM GENES THAT APPEAR IN TREEFAM, AND THE
# FAMILIES THAT THEY APPEAR IN:

# ar2 - convert to WORMPEP by going Gene_name->Gene->CDS->Protein
# TreeFam doesnt deal in isoforms

# gw3 - the ACE database used has been changed from WS155 to current_DB
# we can't use autoace because CDS->Protein is not yet in autoace at this stage of the Build

# klh 20120712 - Yes we can use autoace - we have already finished wormpep, so CDS->Protein connection
# are in place

my $db = Ace->connect( -path => $wormbase->autoace ) or $log->log_and_die(Ace->error,"\n");
open (OUT,">$acefile") or $log->log_and_die("Cannot write to $acefile\n");
foreach my $id (keys %family) {
  my $family = $family{$id};
  my @family = split(/\,/,$family); # THIS IS A LIST OF THE FAMILIES THAT A WORM GENE APPEARS IN.
  my $no_families = scalar(@family);

  my $gene = $id;
  
  my (@wormpep_obj);
  
  $gene = $gene =~ /(\S+\.\w+)\.\d+$/ ? $1 : $gene;
  if( $gene =~ /WBGene/ ) { # if a WBGeneID was used
    my $gene_obj = $db->fetch(Gene => "$gene");

    if ($gene_obj) {
      my @cds = $gene_obj->Corresponding_CDS;
      foreach my $cds (@cds) {
        if ($cds->Corresponding_protein) {
          push @wormpep_obj, $cds->Corresponding_protein;
        }
      }
    }
  } else { 
    # try to search for a CDS, else for a Gene_name
    my $cds_obj = $db->fetch(CDS=>"$gene");
    if ($cds_obj) {
      push @wormpep_obj, $cds_obj->Corresponding_protein;
    } else {
      my $gene_name = $db->fetch(Gene_name => "$gene");
      if ($gene_name) {      
        my $gene_obj = $gene_name->Molecular_name_for;
        $gene_obj ||= $gene_name->Other_name_for;
        $gene_obj ||= $gene_name->Sequence_name_for;
        
        if ($gene_obj) {
          my @cds = $gene_obj->Corresponding_CDS;
          foreach my $cds (@cds) {
            if ($cds->Corresponding_protein) {
              push @wormpep_obj, $cds->Corresponding_protein;
            }
          }
        }
      } 
    }
  }

  foreach my $wormpep (@wormpep_obj) {
    print OUT "\nProtein : \"",$wormpep->name,"\"\n";
    foreach (@family) { 
      print OUT "Database TREEFAM TREEFAM_ID $_\n";
    }
  }
}
close(OUT);

#------------------------------------------------------------------# 

$wormbase->load_to_database($wormbase->autoace, $acefile, 'treefam', $log) unless ($no_load or $test);

$log->mail;

exit(0);

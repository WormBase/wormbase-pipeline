#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Log_files;
use Getopt::Long;
use strict;
use Storable;

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

my ($debug, $test,$store, $no_load);
GetOptions (
	    "debug:s"       => \$debug,
	    "test"          => \$test,
	    "store:s"       => \$store,
	    "no_load"       => \$no_load
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
my $log = Log_files->make_build_log($wormbase);

# DECLARE MYSQL USERNAME AND HOST:

use DBI;

#------------------------------------------------------------------#

# GET THE LONG NAMES OF THE TREEFAM GENES FROM THE MYSQL DATABASE:

my %WORM                      = (); # HASH TABLE TO KEEP A LIST OF WORM GENES IN TREEFAM.
my $database                  = 'treefam_3';

$log->write_to("connecting to treefam database : \tmysql:treefam_3:db.treefam.org:3308\n");
my $dbh                       = DBI->connect("dbi:mysql:treefam_3:db.treefam.org:3308", 'anonymous', '') || return;
my $table_w                   = 'genes';
# THIS TABLE HAS THE ID AND DISPLAY ID. SOMETIMES THE DISPLAY ID IS
# THE UNIPROT NAME, SOMETIMES THE ID IS:
my $st                        = "SELECT ID, TAX_ID from $table_w"; 
my $sth                       = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
my$rv                        = $sth->execute or die "Cannot execute the query: $sth->errstr";
if ($rv >= 1)
{
   while ((my @array) = $sth->fetchrow_array)    {
      my $ID                  = $array[0];  # eg., AH3.1 for a C. elegans gene or
                                         #      WBGene00024691 for a C. briggsae gene
      my $TAX_ID              = $array[1];  # eg., 6239 for a C. elegans gene or
                                         #      6238 for a C. briggsae gene 
      if ($TAX_ID == 6239 || $TAX_ID == 6238) # IT IS A C. ELEGANS OR C. BRIGGSAE GENE
      {
         # REMEMBER THAT THIS WORM GENE IS IN TREEFAM:
         $WORM{$ID}        = $TAX_ID;  
      } 
   }
}
my $rc                        = $dbh->disconnect();
$rc                        = "";

#------------------------------------------------------------------#

# GET THE NAMES OF ALL THE TREEFAM FAMILIES AND THE GENES THAT ARE IN
# THEM FROM THE MYSQL DATABASE:

$database                  = 'treefam';
$dbh                       = DBI->connect("dbi:mysql:treefam_3:db.treefam.org:3308", 'anonymous', '') || return;
# FIRST READ IN TREEFAM-A AND THEN TREEFAM-B:

my %FAMILY                    = (); # HASH TABLE TO KEEP A RECORD OF THE TREEFAM FAMILIES THAT A WORM GENE IS IN.
for (my $i = 1; $i <= 2; $i++)
{
   # SPECIFY THE TABLE:

   if    ($i == 1) # LOOK AT TREEFAM-A:
   {
      $table_w             = 'famA_gene';
   }
   elsif ($i == 2) # LOOK AT TREEFAM-B:
   {
      $table_w             = 'famB_gene';
   }

   # THE FIRST THREE COLUMNS IN THE TABLE famB_gene/famA_gene ARE THE TRANSCRIPT NAME, FAMILY NAME AND WHETHER THE
   # TRANSCRIPT IS IN THE SEED/FULL TREE:
   # eg., ENSMUST00000049178.2 TF105085 FULL

   my $st                     = "SELECT ID, AC, FLAG from $table_w'"; 
   my $sth                    = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
   my $rv                     = $sth->execute or die "Cannot execute the query: $sth->errstr";
   if ($rv >= 1)
   {
      while ((my @array) = $sth->fetchrow_array)    {
         my $ID               = $array[0];  # eg., F40G9.2.1 for a C. elegans gene OR WBGene00027163 for a C. briggsae gene.
         my $AC               = $array[1];  # eg., TF105085, NAME OF THE TREEFAM FAMILY.
         my $FLAG             = $array[2];  # eg., FULL OR BOTH

 
         if (($FLAG eq 'FULL' || $FLAG eq 'BOTH')) # WE ARE ONLY INTERESTED IN THE FULL TREES.
         {
            # CHECK IF IT IS A C. BRIGGSAE OR C. ELEGANS GENE:
            if ($WORM{$ID})
            {
               # REMEMBER THE FAMILIES THAT THIS WORM GENE IS IN:
               if (!($FAMILY{$ID})) { $FAMILY{$ID} = $AC;                 }
               else                 { $FAMILY{$ID} = $FAMILY{$ID}.",".$AC;}
            }
         } 
      }
   }
}
$rc                        = $dbh->disconnect();
$rc                        = "";

#------------------------------------------------------------------# 

# PRINT OUT A LIST OF THE WORM GENES THAT APPEAR IN TREEFAM, AND THE
# FAMILIES THAT THEY APPEAR IN:

# ar2 - convert to WORMPEP by going Gene_name->Gene->CDS->Protein
# TreeFam doesnt deal in isoforms

# gw3 - the ACE database used has been changed from WS155 to current_DB
# we can't use autoace because CDS->Protein is not yet in autoace at this stage of the Build

my $db = Ace->connect( -path => $wormbase->database('current') ) or $log->log_and_die(Ace->error,"\n");
open (OUT,">".$wormbase->acefiles."/treefam.ace") or $log->log_and_die("cant write to ".$wormbase->acefiles." dir\n");
foreach my $ID (keys %FAMILY)
{
  my $family                 = $FAMILY{$ID};
  my @family                 = split(/\,/,$family); # THIS IS A LIST OF THE FAMILIES THAT A WORM GENE APPEARS IN.
  my $no_families            = $#family + 1; # THIS IS THE NUMBER OF FAMILIES THAT A WORM GENE APPEARS IN.
  #print "$ID $no_families $family\n"; 
  my $gene = $ID;
  my $gene_obj;
  $gene = $gene =~ /(\w+\.\w+)\.\d+/ ? $1 : $gene;
   if( $gene =~ /WBGene/ ) {
     $gene_obj = $db->fetch("query find Gene $gene");
   }
   else {
     $gene_obj = $db->fetch(Gene_name => "$gene");
   }
   if ($gene_obj) {
     $gene_obj = $gene_obj->Molecular_name_for;
   }
   if ( $gene_obj ) {
     if ( my $cds = $gene_obj->Corresponding_CDS ) {
       if ( my $wormpep = $cds->Corresponding_protein ) {
	 print OUT "\nProtein : \"",$wormpep->name,"\"\n";
	 foreach (@family) {
	   print OUT "Database TREEFAM TREEFAM_ID $_\n";
	 }
       }
     }
   }
}

#------------------------------------------------------------------# 

$wormbase->load_to_database($wormbase->autoace,$wormbase->acefiles."/treefam.ace",'treefam') unless $no_load;

$log->mail;
print STDERR "FINISHED.\n";
exit(0);

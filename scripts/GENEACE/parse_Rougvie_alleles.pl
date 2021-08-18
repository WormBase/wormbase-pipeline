#!/usr/local/bin/perl5.8.0 -w

use strict;
use Getopt::Long;
use lib -e $ENV{'CVS_DIR'};
use Wormbase;

my ($filein, $fileout, $wormbase, $USER);

GetOptions (
    'filein=s'         => \$filein,
    'fileout=s' => \$fileout,
    'user:s'    => \$USER,
    );

$wormbase = Wormbase->new("-organism" =>$species, -debug => $debug, -test => $test);
my $db;
my $tace            = $wormbase->tace;
my $database = "/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/geneace";
my $ace = Ace->connect (-path => $database,
			-program => $tace) || die "cannot connect to database at $database\n";

unless ((defined $filein) && ($fileout)) {
    print "You need to specif both -filein and -fileout\n";
}
print "$filein\n" if defined($filein);
print "$fileout\n" if defined($fileout);


my $curator;
if (defined $USER){ 
    if ($USER eq 'pad') {
        $curator = 'WBPerson1983';
    } elsif ($USER eq 'skd') {
        $curator = 'WBPerson51134';
    } elsif ($USER eq 'mz3') {
        $curator = 'WBPerson21950';
    }
}
open (IN, "<$filein") or die("Failed to open input file\n");
open (OUT, ">$fileout") or die("Failed to open output file\n");


while (<IN>) {
  chomp;
  my @f=split"\t";

  my $obj = $ace->fetch(Gene=>$f[);
  if (defined $obj) {
      my $wbvar_id = $obj->Public_name_for->name; 
      print ACE "Variation : \"$wbvar_id\"\n";
  }
  else {
      print ACE "Variation : \"$public_name\"\n";
  }


  print OUT "Variation \: \"$f[1]\"\nPublic_name $f[1]\nStrain $f[0]\nFlanking_sequences $f[7] $f[8]\nMapping_target CHROMOSOME_$f[3]\nDeletion\nSequenced\nEngineered_allele\nSpecies \"Caenorhabditis elegans\"\nCRISPR_Cas9\nLive\nGene $f[2]\nRemark \"Whole gene deletions made by the Rougvie, Moerman, and soon Hutter labs that replace the coding sequence with a [LoxP + myo-2p::GFP::unc-54 3Prime UTR + rps-27p::neoR::unc-54 3Prime UTR + LoxP] cassette.\"\nMethod Engineered_allele\n\n";



#my ($,$TYPE,$From,$TO,$Flank1,$Flank6,$ID2,$CHROM,$Start,$Stop);

#    if (/(pas\d+)\s+(Substitution)\s+(\S)\s+(\S)\s+(\S+)\s+(\S+)/) {
#      $ID = $1;
#      $TYPE = $2;
#      $From = $3;
#      $TO = $4;
#      $Flank1 = $5;
#      $Flank6 = $6;
#    }
}
print "Diaskeda same Poli\n"; #we had alot of fun#
close (IN);
close (OUT);
exit(0);
__END__

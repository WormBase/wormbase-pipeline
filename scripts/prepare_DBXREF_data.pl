#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Storable;

#################################################
# script usage
#
# cd ~/DATABASES
#
# prepare_DBXREF_data.pl ~wormpub/BUILD/autoace/Elegans.store WS220 Wormbase_WS220_mapping.txt
#
#################################################

my $store = shift;#'/nfs/disk100/wormpub/BUILD/autoace/wormbase.store';
my $wormbase = retrieve( $store );
my $database = shift;
my $output = "./DB_xref_files/".shift;
my $error_file = $output."_err";

if (-e $output) {$wormbase->run_command("rm $output");
print "Removing old output files and creating new ones.\n";
}

open (ERR, ">$error_file") or die "Cannot open output file $error_file\n";
open (OUT, ">$output") or die "Cannot open output file $output\n";

$database ||= $wormbase->autoace;#$ENV{'CURRENT'}; if the database isn't specified then autoace is used

print OUT "//----------------------------------------------------------------------\n";
print OUT "// Updated mappings between the ENA entries and WormBase release $database\n";
print OUT "// Mappings for coding, non-coding and pseudogenic loci. \"na\" is used when a protein ID is not available.\n";
print OUT "// File Format: <Parent Accession> Protein_ID WBGene_ID Standard_name\n//\n";
print OUT "// Examples:\n//\n//\tCDS:\t\tAF025472	AAY43980	WBGene00022688	ZK250.5b\n//\n//\tncRNA:\t\tAL022288        na      WBGene00199582  ZK1025.16\n//\n//\tPseudogene:\tAC199241	na	WBGene00044951	2RSSE.3\n//----------------------------------------------------------------------\n\n";

my %clone2acc;
$wormbase->FetchData('clone2accession',\%clone2acc,"$database/COMMON_DATA");

my $db = Ace->connect( -path=> $database) or die Ace->error;

my @genes = $db->fetch(-query => 'Find worm_genes'); 
# testing only
#my @genes = $db->fetch(-query => 'Find worm_genes where method != curated');
my $gene_count = scalar grep { defined $_ } @genes;
my $count = "0";

foreach my $gene ( @genes ) {
  $count ++;
  print "$gene_count/$count\n";

  #Pseudogene and Transcript objects.
  if ($gene->class eq "Pseudogene" or $gene->class eq "Transcript") {
    #Parent
    my $parent = $gene->S_parent->right->name;
    my $EMBL_acc = $clone2acc{$parent};
    my $pid = "na";
    my $WBG_id = $gene->Gene->name;
    my $seconday = $gene->name;
    if ($EMBL_acc and $pid and $WBG_id and $seconday) {
      print OUT "$EMBL_acc\t$pid\t$WBG_id\t$seconday\n";
    }
  }

  #some CDSs cross boundaries so are in multiple clone  EMBL files. 
  #CDS objects.....not 100% sure why they are processed like this, but preserved for now.
  else {
    my @pid_lines = $gene->Protein_id;
    foreach my $parent (@pid_lines) {
      my $EMBL_acc = $clone2acc{$parent->name};
      my $pid = "na";
      if ($gene->Protein_id and  $gene->Protein_id->right and $gene->Protein_id->right->right ) {
	$pid =$gene->at(DB_info)->Protein_id->right->right->name
      }
      my $WBG_id = $gene->Gene->name;
      my $seconday = $gene->name;
      if ($EMBL_acc and $pid and $WBG_id and $seconday) {
	print OUT "$EMBL_acc\t$pid\t$WBG_id\t$seconday\n";
      }
      else {
	print ERR "ERROR data missing for $gene : $WBG_id\n";
      }
    }
  }
}
close OUT;
close ERR;
$db->close;
print "\nYour output file is here $output\nFinished\n";
exit(0);

__END__

=pod

=head2 EMBL XREF generation

  This script will generate the file that is needed for the EMBL release XREFS

  An email will be sent to wormbase-help requesting this data.   As of Aug 2006 it is from Alastair Baldwin <abaldwin@ebi.ac.uk>

  The script runs against autoace so make sure that the protein_ids have been loaded (make_wormpep.pl -final)

  Dont gzip the file to send it as EBI autorejects files with a .gz extension

=cut

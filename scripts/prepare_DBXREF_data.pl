#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Storable;
use Getopt::Long;

#################################################
# script usage
#
# cd ~/DATABASES
#
#
# prepare_DBXREF_data.pl -database WS227 -out Wormbase_WS220_mapping.txt
#
#################################################


my ($database,$store,$out);

GetOptions ("database:s" => \$database,
	    "store:s" => \$store,
	    "out:s" => \$out,
	    "build" => \$build
	   );

#Minimum requirements
if ((!defined $database) && (!defined $build)) {die "ERROR: You must specify either -database or -build as a mandatory option\n";}

unless (defined $store) {
  if (defined $build) {
    print "Using BUILD database storable\n";
    $store = "/nfs/users/nfs_w/wormpub/BUILD/autoace/Elegans.store";
  }
  elsif (-e "/nfs/users/nfs_w/wormpub/DATABASES/".$database."Elegans.store") {
    print "Using $database storable\n";
    $store = "~wormpub/DATABASES/".$database."Elegans.store";
  }
  else {
    print "No storable specified and cannot retrieve $database store therefore defaulting to using current_DB storable\n";
    $store = "/nfs/users/nfs_w/wormpub/DATABASES/current_DB/Elegans.store";
  }
}

my $wormbase = retrieve($store);
my $version;

unless (defined $out) {
  if (defined $build) {
    $version = $wormbase->get_wormbase_version;
    $out = "Wormbase_WS".$version."_mapping.txt";
  }
  if (defined $database) {
    $out = "Wormbase_".$database."_mapping.txt";
    $version = $database;
  }
}



my $base = $wormbase->wormpub;
my $output = $base."/DATABASES/DB_xref_files/".$out;
my $error_file = $output."_err";
my $dbpath;

$dbpath = $base."/DATABASES/$database" unless (defined $build);
$dbpath = $wormbase->autoace unless (defined $database);

if (-e $output) {$wormbase->run_command("rm $output");
print "Removing old output files and creating new ones.\n";
}

open (ERR, ">$error_file") or die "Cannot open output file $error_file\n";
open (OUT, ">$output") or die "Cannot open output file $output\n";

print OUT "//----------------------------------------------------------------------\n";
print OUT "// Updated mappings between the ENA entries and WormBase release $version\n";
print OUT "// Mappings for coding, non-coding and pseudogenic loci. \"na\" is used when a protein ID is not available.\n";
print OUT "// File Format: <Parent Accession> Protein_ID WBGene_ID Standard_name\n//\n";
print OUT "// Examples:\n//\n//\tCDS:\t\tAF025472	AAY43980	WBGene00022688	ZK250.5b\n//\n//\tncRNA:\t\tAL022288        na      WBGene00199582  ZK1025.16\n//\n//\tPseudogene:\tAC199241	na	WBGene00044951	2RSSE.3\n//----------------------------------------------------------------------\n\n";

my %clone2acc;
$wormbase->FetchData('clone2accession',\%clone2acc,"$dbpath/COMMON_DATA");

my $db = Ace->connect( -path=> $dbpath) or die Ace->error;

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

  An email will be sent to wormbase-help requesting this data.   As of Aug 2011 it is from Datasubs <datasubs@ebi.ac.uk>

  The script runs against the specified database so make sure that the protein_ids have been loaded (make_wormpep.pl -final)
  
  And there are other limitations like a storable object and COMMON_DATA

  Dont gzip the file to send it as EBI.

  

=head2 Mandatory arguments

 -database WS227 // This option allows you to specify a release under ther DATABASES dir.

 or

 -build // this grabs the xref data from the ongoing build, useful if the build is in the last phases.

=head2 Optional arguments

 -out Wormbase_WS220_mapping.txt // allows you to specify a different output file name.

 -store ~wormpub/DATABASES/current_DB/Elegans.store //allows you to give a storable object for retrieval.

=cut

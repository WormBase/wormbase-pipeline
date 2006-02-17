#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Storable;

my $store = '/nfs/disk100/wormpub/BUILD/autoace/wormbase.store';
my $wormbase = retrieve( $store );
my $database = $ENV{'CURRENT'};
my %clone2acc = $wormbase->FetchData('clone2accession');

my $db = Ace->connect( -path=> $database) or die Ace->error;

my @genes = $db->fetch(-query => 'Find worm_genes');
foreach my $gene ( @genes ) {
  my $EMBL_acc = $clone2acc{$gene->Sequence->name};
  my $pid = "na";
  if($gene->Protein_id and  $gene->Protein_id->right and $gene->Protein_id->right->right ) {
    $pid =$gene->Protein_id->right->right;
  }

  my $WBG_id = $gene->Gene->name;
  my $seconday = $gene->name;

  print "$EMBL_acc\t$pid\t$WBG_id\t$seconday\n" if ($EMBL_acc and $pid and $WBG_id and $seconday);
}

$db->close;


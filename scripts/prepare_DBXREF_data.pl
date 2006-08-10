#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Storable;

my $store = '/nfs/disk100/wormpub/BUILD/autoace/wormbase.store';
my $wormbase = retrieve( $store );
my $database = $wormbase->autoace;#$ENV{'CURRENT'};
my %clone2acc = $wormbase->FetchData('clone2accession');

my $db = Ace->connect( -path=> $database) or die Ace->error;

my @genes = $db->fetch(-query => 'Find worm_genes');
foreach my $gene ( @genes ) {

  #some CDSs cross boundaries so are in multiple clone  EMBL files.
  my @pid_lines = $gene->Protein_id;
  foreach my $parent (@pid_lines) {
    my $EMBL_acc = $clone2acc{$parent->name};
    my $pid = "na";
    if ($gene->Protein_id and  $gene->Protein_id->right and $gene->Protein_id->right->right ) {
      $pid =$gene->Protein_id->right->right;
    }

    my $WBG_id = $gene->Gene->name;
    my $seconday = $gene->name;

    print "$EMBL_acc\t$pid\t$WBG_id\t$seconday\n" if ($EMBL_acc and $pid and $WBG_id and $seconday);
  }
}

$db->close;


exit(0);

__END__

=pod

=head2 EMBL XREF generation

  This script will generate the file that is needed for the EMBL release XREFS

  An email will be sent to wormbase-help requesting this data.   As of Aug 2006 it is from Alastair Baldwin <abaldwin@ebi.ac.uk>

  The script runs against autoace so make sure that the protein_ids have been loaded (make_wormpep.pl -final)

  Dont gzip the file to send it as EBI autorejects files with a .gz extension

=cut

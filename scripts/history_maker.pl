#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Tk;
require Tk::Dialog;

my $version =shift;
die "no version\n" unless $version;
my $database = glob("~wormpub/camace_orig");
my $session_file = "/tmp/history_session.$version";

my $form_cds;

# set up Tk interface
my $main_gui = MainWindow->new();
$main_gui->configure(-title => "History maker for WS${version}",
		     -background => 'blue'
		     );
$main_gui->geometry("400x100");

my $db_lbl = $main_gui->Label( -text => "$database",
				-background => 'blue',
				-foreground => 'white'
			      )->pack( -pady => '3'
				     );

my $cds_lbl = $main_gui->Label( -text => 'CDS',
				-background => 'blue',
				-foreground => 'white'
			      )->pack(-pady => '6',
				      -padx => '6',
				      -side => 'left',
				     );
my $cds_val = $main_gui->Entry( -width => '10',
				-background => 'white',
				-textvariable=> \$form_cds,
			      )->pack(-side => 'left',
				      -pady => '5',
				      -padx => '5'
				     );

my $make = $main_gui->Button( -text => "Make History",
			      -command => [\&make_history, \$form_cds]
			      )->pack(-side => 'right',
				      -pady => '2',
				      -padx => '6',
				      -anchor => "w"
				     );

my $db = Ace->connect(-path => $database);

MainLoop(); # create GUI

$db->close;

exit(0);

sub make_history 
  {
    #print "enter CDS to make history object \n";
    my $cds_form = shift;
    $cds = $$cds_form;
    $output = $session_file.$cds;
    open (HIS,">$output") or die "cant open $output\n";
    last if( $cds eq "end" );
    my $obj = $db->fetch(CDS => "$cds");

    return &nocds_message("$$cds_form") unless $obj;
    if( $obj->Method->name ne "curated" ) {
      print STDERR "I only do curated histories!\n";
      next;
    }

    my $species = $obj->Species->name;
    my $gene = $obj->Gene->name;
    my $lab = $obj->From_laboratory->name;
    my $seq = $obj->Sequence->name;

    my $protein = $obj->Corresponding_protein;
    $protein = $protein->name if $protein;

    # parent clone coords
    my $clone = $obj->Sequence;
    my @clone_CDSs = $clone->CDS_child;
    my $start;
    my $end;
    foreach my $CDS ( @clone_CDSs ) {
      next unless ($CDS->name eq "$cds");
    $start = $CDS->right->name;
      $end = $CDS->right->right->name;
    last;
    }

    #print ace format
    print HIS "Sequence : $seq\n";
    print HIS "CDS_child \"$cds:wp$version\" $start $end\n";
    print HIS "\nCDS : $cds:wp$version\n";

    foreach ($obj->Source_exons) {
      my ($start,$end) = $_->row(0);
      print HIS "Source_exons ",$start->name," ",$end->name,"\n";
    }

    print HIS "Sequence $seq\n";
    print HIS "CDS\n";
    print HIS "From_laboratory $lab\n";
    print HIS "Gene_history $gene\n" if $gene;
    print HIS "Method history\n";

    close HIS;
    system("xremote -remote 'parse $output'");

    &confirm_message("$cds:wp$version");
  }

sub nocds_message
  {
    my $cds = shift;
    my $D = $main_gui->Dialog(
			      -title => 'Invalid CDS',
			      -text  => "$cds is not a valid CDS name",
			      -default_button => 'ok',
			      -buttons        => ['ok']
			     );
    $D->configure(-background => 'red');
    my $choice = $D->Show;

  }

sub confirm_message
  {
    my $cds = shift;
    my $D = $main_gui->Dialog(
			      -title => 'History made',
			      -text  => "$cds histroy object made",
			      -default_button => 'ok',
			      -buttons        => ['ok']
			     );
    $D->configure(-background => 'cyan',
		  -foreground => 'black');

    my $choice = $D->Show;

  }

sub make {
  print "rubbish\n";
}

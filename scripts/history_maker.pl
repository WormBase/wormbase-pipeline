#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Ace;
use Tk;
require Tk::Dialog;


# This is the database used a reference for generating histories NOT the one that will be changed
my $database = glob("~wormpub/camace_orig");

# pass path to latest version of wormbase
my $version = &get_history_version("/nfs/disk100/wormpub/DATABASES/current_DB");

# file where temp ace files written
my $session_file = "/tmp/history_session.$version";

# set up Tk interface
my $main_gui = MainWindow->new();
my $form_cds; # cds variable from form

# Main window
$main_gui->configure(-title => "History maker for WS${version}",
		     -background => 'blue'
		    );
$main_gui->geometry("400x100");

# Reference database lable
my $db_lbl = $main_gui->Label( -text => "$database",
			       -background => 'blue',
			       -foreground => 'white'
			     )->pack( -pady => '3'
				    );
# CDS entry widgets
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

#make Return and Enter submit CDS 
$cds_val->bind("<Return>",[ \&make_history]);
$cds_val->bind("<KP_Enter>",[ \&make_history]);

#Make history button
my $make = $main_gui->Button( -text => "Make History",
			      -command => [\&make_history]
			    )->pack(-side => 'right',
				    -pady => '2',
				    -padx => '6',
				    -anchor => "w"
				   );
# Clear CDS entry button
my $clear = $main_gui->Button( -text => "Clear",
			       -command => [\&clear]
			     )->pack(-side => 'left',
				     -pady => '2',
				     -padx => '6',
				     -anchor => "e"
				    );

# aceperl connection to reference database
my $db = Ace->connect(-path => $database);

# create GUI
MainLoop();

$db->close;

exit(0);

sub make_history 
  {
    #print "enter CDS to make history object \n";
    $cds = $form_cds;
    return unless $cds;
    $output = $session_file.$cds;
    open (HIS,">$output") or die "cant open $output\n";
    last if( $cds eq "end" );

    $cds = &confirm_case($cds);

    my $obj = $db->fetch(CDS => "$cds");

    return &nocds_message("$cds") unless $obj;
    if ( $obj->Method->name ne "curated" ) {
      print STDERR "I only do curated histories!\n";
      next;
    }

    if ($db->fetch(CDS => "$cds:wp$version") ) {
      &error_prexisting("$cds:wp$version");
      return;
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
    print HIS "Species \"$species\"\n" if $species;
    print HIS "Method history\n";

    close HIS;
    system("xremote -remote 'parse $output'");

    &confirm_message("$cds:wp$version");
    &clear;
  }

sub clear
  {
    $cds_val->delete(0,'end');
  }

# Dialoague boxes
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

sub error_prexisting
  {
    my $cds = shift;
    my $D = $main_gui->Dialog(
			      -title => 'History exists',
			      -text  => "$cds history object already exists in database",
			      -default_button => 'ok',
			      -buttons        => ['ok']
			     );
    $D->configure(-background => 'red',
		  -foreground => 'black');

    my $choice = $D->Show;
  }

# gets the version of the latest database - confirm this is valid for your setup
sub get_history_version
  {
    my $database = shift;
    my $WS_version = `grep "NAME WS" $database/wspec/database.wrm`;
    chomp($WS_version);
    $WS_version =~ s/.*WS//;
    die "no version\n" unless $WS_version =~ /\d+/;
    $WS_version++;
    return($WS_version); 
  }

# makes sure all CDSs are uppercase ( apart from the isoform part )
sub confirm_case
  {
    my $cds = shift;
    my ($last_char) = lc(chop $cds);
    $cds = (uc $cds) . "$last_char";
    return $cds;
  }


__END__


=pod

=head1 NAME - history_maker.pl

A perl Tk interface presents to the user a simple box with a space to enter a CDS name and a button to make a history object.  The input to the acedb database is done immediately using 'xremote' and is visible after recalculating the fmap.  The target database is determined by the operating system as the last one to have been opened using the terminal that this script is run from.

To ensure that the correct database is being used a shell script should be used to launch both simultaneously.

A reference database is hardcoded ( at the moment ) and this is used to extract the relevant info needed to make a history ie source exons, gene_id etc.

Some error checking is done so that;

=over4 

  history objects cant be created for non-existant CDSs.
  history objects with the same name as existing histories will not be made.
  input case is irrelevant - all histories will be converted to uppercase clone names.

=back

=cut

#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Ace;
use Tk;
use Getopt::Long;
require Tk::Dialog;
use Coords_converter;

my $source;
my $design;
my $chromosome;
my $blast;
my $user;

GetOptions (
	    "source:s"     => \$source,
	    "design"       => \$design,
	    "chromosome:s" => \$chromosome,
	    "blast=s"      => \$blast,
	    "user:s"       => \$user,
	    "help|h"         => sub { system("perldoc $0"); exit(0);}
	   );

die "$blast doesnt exist !\n" if ($blast and !(-e $blast));

# This is the database used a reference for generating histories NOT the one that will be changed
my $database = $source ? $source : glob("~wormpub/camace_orig");

# pass path to latest version of wormbase
my $version = &get_history_version("/nfs/disk100/wormpub/DATABASES/current_DB");

# file where temp ace files written
my $session_file = "/tmp/history_session.$version";

# set up Tk interface
my $main_gui = MainWindow->new();
my $form_cds; # cds variable from form
my $form_gene; # gene variable from form

# Main window
$main_gui->configure(-title => "Curation Tool for WS${version}",
		     -background => 'blue'
		    );

my $gui_height = 300;
$gui_height += 200 if $chromosome;
$gui_height += 200 if $blast;
$main_gui->geometry("500x$gui_height");


################################################################
#intron locator frame
if ( $chromosome ) {
  my $intron_find = $main_gui->Frame( -background => "red",
				      #				    -height     => "400",
				      -label      => "Intron locator",
				      -relief     => "raised",
				      -borderwidth => 5,
				    )->pack( -pady => "20",
					     -fill => "x"
					   );

  my $intron_list = $intron_find->Scrolled("Listbox", -scrollbars => "osoe",
					   -selectmode => "single",
					   -height     => "5",
					   -width      => "100"
					  )->pack();


  my $go_to_intron = $intron_find->Button ( -text    => "Go to intron",
					    -command => [\&goto_intron,\$intron_list]
					  )->pack ( -side => 'right',
						    -pady => '2',
						    -padx => '6',
						    -anchor => "w"
						  );


  my $coords = Coords_converter->invoke;

  my $file = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/CHROMOSOME_${chromosome}.check_intron_cam.gff";
  open (INTRONS, "<$file") or die "$file\n";
  while ( <INTRONS> ) {
    #CHROMOSOME_X    .       intron  12117028        12117072        .       -       .       Confirmed_EST OSTF209B10_1      Clone:T01H10 8265 8309  Lab:HX
    my @data = split;
    $seq = $data[0]; 
    $x = $data[3];
    $y = $data[4];
    my @clone_coords = $coords->LocateSpan("$seq","$x","$y" );
    $intron_list->insert('end',"$clone_coords[0] $clone_coords[1] $clone_coords[2]");
    #last if $count++ > 15;
  }
  close INTRONS;
}
################################################################

# history_maker

my $his_maker = $main_gui->Frame( -background => "blue",
				    -height     => "400",
				    -label      => "History Maker",
				    -relief     => "raised",
				    -borderwidth => 5,
				  )->pack( -pady => "20",
					   -fill => "x"
					 );

# Reference database lable
my $db_lbl = $his_maker->Label( -text => "$database",
			       -background => 'blue',
			       -foreground => 'white'
			     )->pack( -pady => '3'
				    );
# CDS entry widgets
my $cds_lbl = $his_maker->Label( -text => 'CDS',
				-background => 'blue',
				-foreground => 'white'
			      )->pack(-pady => '6',
				      -padx => '6',
				      -side => 'left',
				     );
my $cds_val = $his_maker->Entry( -width => '10',
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
my $make = $his_maker->Button( -text => "Make History",
			      -command => [\&make_history]
			    )->pack(-side => 'right',
				    -pady => '2',
				    -padx => '6',
				    -anchor => "w"
				   );
# Clear CDS entry button
my $clear = $his_maker->Button( -text => "Clear",
			       -command => [\&clear]
			     )->pack(-side => 'left',
				     -pady => '2',
				     -padx => '6',
				     -anchor => "e"
				    );


###########################################################
##   genefinder / twinscan blesser
my $gene_blesser = $main_gui->Frame( -background => "green",
				    -height     => "400",
				    -label      => "Genefinder / Twinscan blesser",
				    -relief     => "raised",
				    -borderwidth => 5,
				  )->pack( -pady => "20",
					   -fill => "x"
					 );
# CDS entry widgets
my $gene_lbl = $gene_blesser->Label( -text => 'prediction name',
				-background => 'green',
				-foreground => 'white'
			      )->pack(-pady => '6',
				      -padx => '6',
				      -side => 'left',
				     );
my $gene_val = $gene_blesser->Entry( -width => '10',
				-background => 'white',
				-textvariable=> \$form_gene,
			      )->pack(-side => 'left',
				      -pady => '5',
				      -padx => '5'
				     );

#make Return and Enter submit CDS 
$gene_val->bind("<Return>",[ \&bless_prediction]);
$gene_val->bind("<KP_Enter>",[ \&bless_prediction]);

#Make history button
my $bless = $gene_blesser->Button( -text => "Bless this gene",
			      -command => [\&bless_prediction]
			    )->pack(-side => 'right',
				    -pady => '2',
				    -padx => '6',
				    -anchor => "w"
				   );
# Clear CDS entry button
my $clear_gene = $gene_blesser->Button( -text => "Clear",
			       -command => [\&clear_gene]
			     )->pack(-side => 'left',
				     -pady => '2',
				     -padx => '6',
				     -anchor => "e"
				    );


###########################################################
#blast hit  locator frame
if ( $blast ) {
  my $blast_find = $main_gui->Frame( -background => "purple",
				     -label      => "Blast hit locator",
				     -relief     => "raised",
				     -borderwidth => 5,
				   )->pack( -pady => "20",
					    -fill => "x"
					  );

  my $blast_list = $blast_find->Scrolled("Listbox", -scrollbars => "osoe",
					 -selectmode => "single",
					 -height     => "5",
					 -width      => "100"
					)->pack();


  my $go_to_blast = $blast_find->Button ( -text    => "Go to blast hit",
					  -command => [\&goto_intron,\$blast_list]
					)->pack ( -side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "w"
						);



  open (BLASTS, "<$blast") or die "cant open $blast\n";
  # K09E2   1660    1860    0.577   8659300
  my @data;
  while ( <BLASTS> ) {
    my @input = split;
    next if ( $input[3] < 10);
    chomp;
    push( @data, "$input[0]  $input[1]  $input[2]   $input[3]");
  }
  close BLASTS;
  foreach (sort @data){
    $blast_list->insert('end',"$_");
  }
}

###########################################################

$user = &check_user unless $user;

# aceperl connection to reference database
my $db = Ace->connect(-path => $database) unless $design;

# create GUI
MainLoop();

$db->close unless $design;

exit(0);

##############################################################
# prediction blesser 
sub bless_prediction
  {
    my $gene = $form_gene;
    return unless $gene;
    last if( $cds eq "end" );

    $cds = &confirm_case($gene);

    my $obj = $db->fetch(CDS => "$gene");
    return &error_warning("Invalid CDS","$gene is not a valid CDS name") unless $obj;
    my $method = $obj->Method->name;
    if ( ($method ne "Genefinder") or ($method eq "Twinscan") ) {
      &error_warning("Wrong method","I only bless Genefinder or Twinscan predictions, my child");
      next;
    }

    my $new_gene = &suggest_name("$gene");
    unless ( $new_gene ){
      &error_warning("No name","Cant suggest name for gene based on $gene");
      return;
    }

    $output = $session_file.$gene;
    open (BLS,">$output") or die "cant open $output\n";

    # transfer details from prediction to new
    my $species = $obj->Species->name;
    my $parent_seq = $obj->Sequence->name;

    my $clone = $obj->Sequence;
    my $lab = $clone->From_laboratory->name;
    my @clone_CDSs = $clone->CDS_child;
    my $start;
    my $end;
    foreach my $CDS ( @clone_CDSs ) {
      next unless ($CDS->name eq "$gene");
      $start = $CDS->right->name;
      $end = $CDS->right->right->name;
      last;
    }
    
    #print ace format
    #   - parent obj
    print BLS "Sequence : $parent_seq\n";
    print BLS "CDS_child \"$new_gene\" $start $end\n";

    #   - new gene
    print BLS "\nCDS : \"$new_gene\"\n";

    foreach ($obj->Source_exons) {
      my ($start,$end) = $_->row(0);
      print BLS "Source_exons ",$start->name," ",$end->name,"\n";
    }

    print BLS "Sequence $parent_seq\n";
    print BLS "CDS\n";
    print BLS "From_laboratory $lab\n";
    print BLS "Species \"$species\"\n" if $species;
    print BLS "Method curated\n";

    #get date for remark
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    print BLS "Remark \"[$date $user] Autoconversion from $gene\"\n";

    close BLS;
    #system("xremote -remote 'parse $output'");

    &confirm_message("Made new gene","Made new CDS $new_gene from $gene");
    &mail_geneace($new_gene);
  }

sub clear_gene
  {
    $gene_val->delete(0,'end');
  }

sub suggest_name
  {
    my $name = shift;
    my ($stem) = $name =~ /(.*)\./;
    my $query = "find worm_genes $stem.*";
    my @genes = $db->fetch( -query => "$query");
    my @names = map($_->name,@genes);
    my @numbers;
    foreach (@names) {
      my ($n) = $_ =~ /\.(\d+)/;
      push ( @numbers,$n);
    }
    my $max =0;

    foreach(  @numbers ) {
      $max = $_ if ($_ > $max );
    }

    my $gene_no = $max + 1;
    &error_warning("suggested","$stem.$gene_no");
    return "$stem.$gene_no";
  }

sub mail_geneace
  {
    my $gene = shift;
    $gene = "TESTING.1";
    open (MAIL,  "|/bin/mailx -r \"$user\@sanger.ac.uk\" -s \"Gene_id required for $gene\" \"mt3\@sanger.ac.uk\"");
    print MAIL "$gene\n";
    close MAIL or warn "mail not sent for $gene\n";
    return;
  }

##############################################################
# History maker
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
    return &error_warning("Invalid CDS","$cds is not a valid CDS name") unless $obj;
    if ( $obj->Method->name ne "curated" ) {
      print STDERR "I only do curated histories!\n";
      next;
    }

    if ($db->fetch(CDS => "$cds:wp$version") ) {
      &error_warning("History exists","$cds:wp$version already exists");
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

    &confirm_message("Made history","History $cds:wp$version has been created");
    &clear;
  }

sub clear
  {
    $cds_val->delete(0,'end');
  }


#############################################################################
# intron finder ( also used by blast hit finder) 

sub goto_intron
  {
    my $intron_list = shift;
    my $selection = $$intron_list->curselection();
    my @selected = split(/\s+/, $$intron_list->get( $selection ));
    &goto_location($selected[0], $selected[1], $selected[2]);
  }

############################################################################
# generic methods 
sub goto_location
  {
    my $seqobj = shift;
    my $x      = shift;
    my $y      = shift;

    my $zoomout = 500;
    $x -= $zoomout;
    $y += $zoomout;
    
    system("xremote -remote 'gif seqget $seqobj -coords $x $y; seqdisplay'");
  }

sub confirm_message
  {
    my $title = shift;
    my $text = shift;
    my $D = $main_gui->Dialog(
			      -title => "$title",
			      -text  => "$text",
			      -default_button => 'ok',
			      -buttons        => ['ok']
			     );
    $D->configure(-background => 'cyan',
		  -foreground => 'black');

    my $choice = $D->Show;
  }

sub error_warning
  {
    my $title = shift;
    my $text = shift;
    my $D = $main_gui->Dialog(
			      -title => "$title",
			      -text  => "$text",
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

# makes sure all CDSs are uppercase ( apart from the isoform or gc/tw part )
sub confirm_case
  {
    my $cds = shift;
    my ($last_char) = lc(chop $cds);
    $cds = (uc $cds) . "$last_char";
    return $cds;
  }

sub check_user
  {
    my $user = shift;
    return if( defined $user );
    my $name = `whoami`;
    chomp $name;
    if( "$name" eq "wormpub" ){
      &error_warning("WORMPUB","Please either run this as yourself or use the -user option. How else can you be blamed for your errors!");
      exit(0);
    }
    else {
      $user = $name;
    }
  }


##################################
##  NOTES ~~~~~~~
###############################
#
#  to drive FMAP display 
#
#  xremote -remote 'gif seqget SUPERLINK_CB_II -coords 111000 112000; seqdisplay'
#
#


__END__

=pod

=head1 NAME - history_maker.pl

A perl Tk interface to aid manual gene curation.

=item options

-chromosome : The chromosome to show "confirmed introns" for

-source     : The database to use as source for history and ab initio predictions

-user       : If you are not yourself enter your user to use in autgenerated comments (when blessing genefinder etc)

-design     : Doesn't make Aceperl connection - dev tool for quicker startup

=item Intron finder

Uses the GFF check files in development_release to provide a list of introns.  Select and click to go to intron in FMAP

=item History Maker

presents to the user a simple box with a space to enter a CDS name and a button to make a history object.  

=item Genefinder / Twinscan Blesser

Enter current CDS name eg AC8.gc3 and click "Bless this gene".  This will create a new CDS with the correct name based on the "worm_genes" class.

=item BLAST hit finder

Using passed file, the user is presented with a list of genomic locations where there are strong blastp hits but no current gene.
Click to go to that location in FMAP.

The input to the acedb database is done immediately using 'xremote' and is visible after recalculating the fmap.  The target database is determined by the operating system as the last one to have been opened using the terminal that this script is run from.

To ensure that the correct database is being used a shell script should be used to launch both simultaneously.

A reference database is used to extract the relevant info needed to make a history ie source exons, gene_id etc.and this is set with the -source option

Some error checking is done so that;

=over4 

  history objects cant be created for non-existant CDSs.
  history objects with the same name as existing histories will not be made.
  input case is irrelevant - all histories will be converted to uppercase clone names.

=back

=cut

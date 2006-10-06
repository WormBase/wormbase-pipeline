#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib $ENV{'CVS_DIR'};
use Ace;
use Tk;
use Getopt::Long;
require Tk::Dialog;
use Coords_converter;
use DBI;
use Wormbase;

 
######################################
# variables and command-line options #
######################################
 
my ($help, $debug, $test, $verbose, $store, $wormbase);


my $source;
my $design;
my $chromosome;
my $blast;
my $user;
my $mail;
my $intron;
my $twinscan;
my $anomaly;
my $lab;

GetOptions (
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,

	    "source:s"     => \$source,
	    "design"       => \$design,
	    "chromosome:s" => \$chromosome,
	    "blast=s"      => \$blast,
	    "user:s"       => \$user,
	    "mail"         => \$mail,
	    "intron"       => \$intron,
	    "twinscan"     => \$twinscan,
	    "anomaly"      => \$anomaly,
	    "lab=s"        => \$lab, # RW or HX
	    "help|h"       => sub { system("perldoc $0"); exit(0);}
	   );

 
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}
       
# always in test mode
$test = 1;

# keep things quiet
$debug = $user if (! defined $debug);


#################################



die "$blast doesnt exist !\n" if ($blast and !(-e $blast));

# This is the database used a reference for generating histories NOT the one that will be changed
my $database = $source ? $source : glob("~wormpub/camace_orig");

# assume the lab is HX if not specified
if (! defined $lab || $lab eq "") {
  $lab = "HX";
} 

# pass path to latest version of wormbase
my $version = &get_history_version("/nfs/disk100/wormpub/DATABASES/current_DB");

# file where temp ace files written
my $session_file = "/tmp/history_session.$version";

my $mysql;			# mysql database handle

# set up Tk interface
my $main_gui = MainWindow->new();

my $form_cds;			# cds variable from form
my $form_gene;			# gene variable from form
my $anomaly_clone;		# clone variable from anomaly form

# Main window
$main_gui->configure(-title => "Curation Tool for WS${version}",
		     -background => 'beige' # was blue
		    );

my $gui_width = 500;
$gui_width += 300 if $anomaly;
my $gui_height = 200;
$gui_height += 200 if $chromosome;
$gui_height += 200 if $blast;
$gui_height += 200 if $twinscan;
$gui_height += 300 if $anomaly;
$main_gui->geometry("${gui_width}x$gui_height");

# this is the value used to construct the window IDs in the script 'find_anomalies.pl'
my $WINDOW_SIZE = 10000;

if (defined $chromosome) {
  $chromosome =~ s/CHROMOSOME_//;
}

my $results; # the globally-accesible anomaly details that are summarised in the anomalies detail window
my $check;			# state of the check button for the weights

################################################################
#intron locator frame
if ( $intron && $chromosome ) {
  my $intron_find = $main_gui->Frame( -background => "wheat", # was red
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


  my $coords = Coords_converter->invoke(undef, undef, $wormbase);

  my $file = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/CHROMOSOME_${chromosome}.check_intron_cam.gff";
  open (INTRONS, "<$file") or die "$file\n";
  while ( <INTRONS> ) {
    #CHROMOSOME_X    .       intron  12117028        12117072        .       -       .       Confirmed_EST OSTF209B10_1      Clone:T01H10 8265 8309  Lab:HX
    my @data = split;
    my $seq = $data[0]; 
    my $x = $data[3];
    my $y = $data[4];
    my @clone_coords = $coords->LocateSpan("$seq","$x","$y" );
    $intron_list->insert('end',"$clone_coords[0] $clone_coords[1] $clone_coords[2]");
    #last if $count++ > 15;
  }
  close INTRONS;
}

################################################################

# history_maker

my $his_maker = $main_gui->Frame( -background => "lightcyan", # was blue
				  -height     => "400",
				  -label      => "History Maker",
				  -relief     => "raised",
				  -borderwidth => 5,
				)->pack( -pady => "20",
					 -fill => "x"
				       );

# Reference database lable
my $db_lbl = $his_maker->Label( -text => "$database",
				-background => 'lightcyan', # was blue
				-foreground => 'black' # was white
			      )->pack( -pady => '3'
				     );
# CDS entry widgets
my $cds_lbl = $his_maker->Label( -text => 'CDS',
				 -background => 'lightcyan',	# was blue
				 -foreground => 'black'	# was white
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
if ($twinscan) {
  my $gene_blesser = $main_gui->Frame( -background => "PaleTurquoise", # was green
				       -height     => "400",
				       -label      => "Genefinder / Twinscan blesser",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "20",
						-fill => "x"
						);
  # CDS entry widgets
  my $gene_lbl = $gene_blesser->Label( -text => 'prediction name',
				       -background => 'PaleTurquoise', # was green
				       -foreground => 'black'
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

  # make Return and Enter submit CDS 
  $gene_val->bind("<Return>",[ \&bless_prediction]);
  $gene_val->bind("<KP_Enter>",[ \&bless_prediction]);
  
  # Make history button
  my $bless = $gene_blesser->Button( -text => "Bless this gene",
				     -command => [\&bless_prediction]
				     )->pack(-side => 'right',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_gene = $gene_blesser->Button( -text => "Clear",
					  -command => [\&clear_gene, \$gene_val]
					  )->pack(-side => 'left',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}

###########################################################
#blast hit locator frame
if ( $blast ) {
  my $blast_find = $main_gui->Frame( -background => "Plum", # was purple
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
  foreach (sort @data) {
    $blast_list->insert('end',"$_");
  }
}

###########################################################
# anomalies database locator frame
if ( $anomaly ) {
  if (! defined $chromosome || $chromosome eq "") {die "Must specify -chromosome with -anomaly\n";}
  if (! defined $user || $user eq "") {die "Must specify -user with -anomaly\n";}

  my $anomaly_detail_list;
  my $zero_weight_label;

  my $coords = Coords_converter->invoke(undef, undef, $wormbase);


  my $anomaly_find = $main_gui->Frame( -background => "wheat", # was cyan
				       -label      => "Anomaly locator",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "20",
						-fill => "x"
						);


  my $anomaly_list = $anomaly_find->Scrolled("Listbox", -scrollbars => "osoe",
					     -selectmode => "single",
					     -height     => "5",
					     -width      => "100"
					     )->pack();


  my $go_to_window = $anomaly_find->Button ( -text    => "Go to anomaly region (or clone)",
					     -background => "aquamarine3",
					     -command => [\&goto_anomaly_window, \$anomaly_list, \$anomaly_detail_list, \$coords]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);


# clone entry widgets
  my $anomaly_clone_lbl = $anomaly_find->Label( -text => 'Clone',
					     -background => 'wheat',
					     -foreground => 'black'
					     )->pack(-pady => '6',
						     -padx => '6',
						     -side => 'left',
						     );

  my $clone_val = $anomaly_find->Entry( -width => '10',
				   -background => 'white',
				   -textvariable=> \$anomaly_clone,
				   )->pack(-side => 'left',
					   -pady => '5',
					   -padx => '5'
					   );

  # start the re-weighting menu stuff
  my $zero_weight;		# variable set by pop-up menu selection

  my $check_button = $anomaly_find->Checkbutton( -text => '',
						 -background => 'wheat',
						 -variable => \$check,
						 -command =>  [\&toggle_weighting, \$zero_weight_label],
						 )->pack(-pady => '6',
							 -padx => '6',
							 -side => 'left',
							 );

  $zero_weight_label = $anomaly_find->Label( -text => 'Not this anomaly',
					     -background => 'wheat',
					     -foreground => 'black'
					     )->pack(-pady => '6',
						     -padx => '0',
						     -side => 'left',
						     );

# test for popup menu for Anthony's anomaly weighting
# see: http://rcswww.urz.tu-dresden.de/CS/perl/modules/Tk/Tk-BrowseEntry.htm
  require Tk::BrowseEntry;
  my $zero_weight_list = $anomaly_find->BrowseEntry(-label => '',
						    -background => 'wheat',
						    -listwidth => '350',
						    -browsecmd => [\&reweight_anomalies, \$mysql, \$lab, \$chromosome, \$anomaly_list, \$anomaly_detail_list, \$zero_weight, \$check],
						    -variable => \$zero_weight
						    )->pack(-pady => '6',
							    -padx => '2',
							    -side => 'left',
							    );

  # create button to run Progress display
  my $display_graphs = 1;
  my $progress_button = $anomaly_find->Button ( -text    => "Progress",
					     -background => "bisque",
					     -command => [\&progress, \$display_graphs]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);

  


# connect to the mysql database
  $mysql = &connect_to_database($lab);

# populate the list of types that we may zero-weight
# and create the weight 'view' table
  &populate_zero_weight_list($mysql, $lab, $chromosome, \$zero_weight_list);


# populate $anomaly_list here from database
  &populate_anomaly_window_list($mysql, $lab, $chromosome, \$anomaly_list);



  my $anomaly_details = $main_gui->Frame( -background => "wheat", # was magenta
					  -label      => "Anomaly details in the selected region",
					  -relief     => "raised",
					  -borderwidth => 5,
					  )->pack( -pady => "20",
						   -fill => "x"
						   );

  $anomaly_detail_list = $anomaly_details->Scrolled("Listbox", -scrollbars => "osoe",
						    -selectmode => "extended",
						    -height     => "12",
						    -width      => "100"
						    )->pack();


  my $ignore_this_anomaly = $anomaly_details->Button ( -text    => "Ignore the selected anomalies",
						       -background => "orange",
						       -command => [\&ignore_anomaly, \$anomaly_detail_list]
					   )->pack ( -side => 'left',
						     -pady => '2',
						     -padx => '6',
						     -anchor => "w"
						     );

## I found that this button was too dangerous - it is too easy to click
## this by mistake and wipe out a whole window of anomalies before you
## have looked at them properly
#  my $ignore_all_window = $anomaly_details->Button ( -text    => "Ignore ALL these anomalies!",
#						     -background => "red",
#						     -command => [\&ignore_anomaly_window, \$anomaly_list, \$anomaly_detail_list]
#					      )->pack ( -side => 'left',
#							-pady => '2',
#							-padx => '6',
#							-anchor => "center"
#							);

  my $go_to_anomaly = $anomaly_details->Button ( -text    => "Go to this anomaly",
						 -background => "aquamarine3",
						 -command => [\&goto_anomaly, \$anomaly_detail_list]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);


  # calculate the progress so far for all chromosomes
  # don't display any graphs
  my $dont_display_graphs = 0;
  &progress(\$dont_display_graphs);

}

##############################################################
#
# end of setting up frames
#
##############################################################

$user = &check_user unless $user;

# aceperl connection to reference database
my $db = Ace->connect(-path => $database) unless $design;

# create GUI
MainLoop();

$db->close unless $design;

if ($anomaly) {
  $mysql->disconnect || die "error disconnecting from database", $DBI::errstr;
}

exit(0);

##############################################################
#
# Subroutines
#
##############################################################

# prediction blesser 
sub bless_prediction
  {
    my $gene = $form_gene;
    return unless $gene;
    #last if( $cds eq "end" );

    #$cds = &confirm_case($gene);
    
    my $obj = $db->fetch(CDS => "$gene");
    return &error_warning("Invalid CDS","$gene is not a valid CDS name") unless $obj;
    my $method = $obj->Method->name;
    my $stem = $obj->Sequence->name;
    my $exceptions = $obj->Sequence->Method->name;
    if (! ( ($method eq "Genefinder") || ($method eq "twinscan") ) ) {
      #if ($method ne "twinscan") {
      &error_warning("Wrong method","I only bless Genefinder or Twinscan predictions, my child");
      next;
    } elsif ($exceptions eq "Link") {
      &error_warning("Warning","This Prediction lies over clone boundaries, my child");
      next;
    }

    my $new_gene = &suggest_name("$stem");
    unless ( $new_gene ){
      &error_warning("No name","Can't suggest name for gene based on $gene");
      return;
    }

    my $output = $session_file.$gene;
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
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    print BLS "Remark \"[$date $user] Autoconversion from $gene\"\n";

    close BLS;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Made new gene","Made new CDS $new_gene from $gene");

      #pop-up window if a -mail option has not been set.
      &error_warning("Info","The gene $new_gene has been created but you have not spread the word, my child") if (!($mail));

      my @IDs = ('1983','1847','1846');
      my $person = "";
      if ($user eq "pad") {
	$person = ($IDs[0]);
      } elsif ($user eq "ar2") {
	$person = ($IDs[1]);
      } elsif ($user eq "dl1") {
	$person = ($IDs[2]);
      }

      #&suggest_name("$stem");

      &mail_geneace("$new_gene", "$person", $date) if ($mail);
    }
  }

##############################################################
sub clear_gene
  {
    my ($gene_val) = @_;
    $gene_val->delete(0,'end');
  }


##############################################################
sub suggest_name
  {
    my $stem = shift;
    my $query = "find worm_genes $stem.*";
    my @genes = $db->fetch( -query => "$query");
    my @names = map($_->name,@genes);
    my @numbers;
    foreach (@names) {
      my ($n) = $_ =~ /\.(\d+)/;
      push ( @numbers,$n);
    }
    my $max =0;

    foreach (  @numbers ) {
      $max = $_ if ($_ > $max );
    }

    my $gene_no = $max + 1;
    &error_warning("suggested","$stem.$gene_no");
    return "$stem.$gene_no";
  }


##############################################################
sub mail_geneace {
  my $gene = shift;
  my $person = shift;
  my $date = shift;

  open (MAIL,  "|/bin/mailx -r \"$user\@sanger.ac.uk\" -s \"Gene_id required for $gene\" \"mt3\@sanger.ac.uk\"");
  print MAIL "============================\n";
  print MAIL "history_maker.pl\n";
  print MAIL "Date:$date\n";
  print MAIL "============================\n\n";
  print MAIL "This is an automated request from $user.\n\n";
  print MAIL "Please create a new gene ID for $gene.\n\n";
  print MAIL "newgene.pl -seq $gene -who $person -email\n\n";
  print MAIL "Thank you for your cooperation.\n";
  close MAIL or warn "mail not sent for $gene\n";
  return;
}

##############################################################
# History maker
sub make_history 
  {
    #print "enter CDS to make history object \n";
    my $cds = $form_cds;
    return unless $cds;
    my $output = $session_file.$cds;
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
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Made history","History $cds:wp$version has been created");
      &clear;
    }
  }

sub clear
  {
    $cds_val->delete(0,'end');
  }


#############################################################################
# connect to the mysql database and return handle
# Args: lab = "HX" or "RW"
sub connect_to_database {
  my ($lab) = @_;

  my $dbsn;
  my $dbuser;
  my $dbpass;

  if ($lab eq "HX") {

    $dbsn = "DBI:mysql:database=worm_anomaly;host=ia64b";
    $dbuser = "wormadmin";
    $dbpass = "worms";

  } else {			# for St. Louis

    $dbsn = "DBI:mysql:database=worm_anomaly;host=XXXXX";
    $dbuser = "XXXXXX";
    $dbpass = "XXXXXX";

  }

  my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {AutoCommit => 1, RaiseError => 1})
      || die "cannot connect to database, $DBI::errstr";

  return $mysql;

}

#############################################################################
#  &populate_zero_weight_list($mysql, $lab, $chromosome, \$zero_weight_list);
# get the list of anomaly types that we can later use to reweight the types
# sets up the VIEW table with initial weights

sub populate_zero_weight_list {

  my ($mysql, $lab, $chromosome, $zero_weight_list_ref) = @_;

  # get the available anomaly types
  my $query = qq{ SELECT type, COUNT(*) FROM anomaly WHERE chromosome = "$chromosome" AND centre = "$lab" AND active = 1 GROUP BY type; };
  my $db_query = $mysql->prepare ( $query );
  $db_query->execute();
  $results = $db_query->fetchall_arrayref;

  # put the types in the zero-weight list with their count
  foreach my $result_row (@$results) {
    # type count
    $$zero_weight_list_ref->insert('end', $result_row->[0] . "          " . $result_row->[1]);
  }

  # set up the weight view table in the mysql database
  # we have MYSQL version 4 which doesn't support VIEWs
  # so have horrible cludge of tables with names formed from
  # user-names which we have to explicitly drop before creating anew
  my $view = "weight_$user";
  # first test to see if this view table exists already and so should be deleted
  $query = qq{ SHOW TABLES; };
  $db_query = $mysql->prepare ( $query );
  $db_query->execute();
  $results = $db_query->fetchall_arrayref;
  foreach my $result_row (@$results) {
    if ($result_row->[0] eq $view) {
      $mysql->do("DROP TABLE $view");
      #print "dropped old view\n";
    }
  }
  # now set up a fresh weight table with initial weights of '1' for each anomaly type
  $mysql->do("CREATE TABLE $view (type varchar(32), weight int(1)) AS SELECT type FROM anomaly GROUP BY type;");
  $mysql->do("UPDATE $view SET weight = 1;");



}
#############################################################################
# get the information on the available (still active) windows sorted
# by the sum of scores of anomalies in them
# populate_anomaly_window_list($mysql, $lab, $chromosome, \$anomaly_list);

sub populate_anomaly_window_list {
  my ($mysql, $lab, $chromosome, $anomaly_list_ref) = @_;

  # name of the temporary 'view' weighting table
  my $view = "weight_$user";

  # get the highest scoring 10 kb windows of one sense or the other sorted by score descending
  # group by window, sense order by 2 desc - means make the
  # SUM(thing_score) sum up rows within distinct window and distinct
  # sense and sort the output by descending SUM(thing_score)

# there is no difference in speed between these two variants of this command:
#
#  my $query = qq{ SELECT a.window, SUM(a.thing_score), a.sense, a.clone FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type    WHERE a.chromosome = "$chromosome" AND a.centre = "$lab" AND a.active = 1 AND w.weight = 1 GROUP BY window, sense ORDER BY 2 DESC; };
#
  my $query = qq{ SELECT a.window, SUM(a.thing_score * w.weight), a.sense, a.clone FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type     WHERE a.chromosome = "$chromosome" AND a.centre = "$lab" AND a.active = 1 GROUP BY window, sense ORDER BY 2 DESC; };

  print "\n";
  print "chromosome=$chromosome\n";
  print "lab=$lab\n";
  #print "query=$query\n";
  #print "\n";

  my $db_query = $mysql->prepare ( $query );

  $db_query->execute();

  $results = $db_query->fetchall_arrayref;

  my $over_10 = 0;
  my $over_5 = 0;
  my $over_2 = 0;
  my $over_1 = 0;
  my $over_half = 0;
  my $over_quarter = 0;
  my $under_quarter = 0;
  foreach my $result_row (@$results) {
    # clone, sense, score, window
    $$anomaly_list_ref->insert('end', $chromosome . " " . $result_row->[3] . " Sense: " . $result_row->[2] . " Score: " . $result_row->[1] . " ID: " . $result_row->[0] );

    # count the numbers in various score bins
    if ($result_row->[1] > 10) {
      $over_10++;
    } elsif ($result_row->[1] > 5) {
      $over_5++;
    } elsif ($result_row->[1] > 2.5) {
      $over_2++;
    } elsif ($result_row->[1] > 1) {
      $over_1++;
    } elsif ($result_row->[1] > 0.5) {
      $over_half++;
    } elsif ($result_row->[1] > 0.25) {
      $over_quarter++;
    } elsif ($result_row->[1] <= 0.25) {
      $under_quarter++;
    }
  }

  print "Number of 10Kb windows: ". scalar @$results . "\n";
  print "\twindows with score over 10:        $over_10\n";
  print "\twindows with score 5 to 10:        $over_5\n";
  print "\twindows with score 2.5 to 5:       $over_2\n";
  print "\twindows with score 1 to 2.5:       $over_1\n";
  print "\twindows with score 0.5 to 1:       $over_half\n";
  print "\twindows with score 0.25 to 0.5:    $over_quarter\n";
  print "\twindows with score less than 0.25: $under_quarter\n";

}

#############################################################################
# this is called when the check button beside the re-weighting popup menu is pressed
# it changes the label on the popup menu

sub toggle_weighting {

  my ($zero_weight_list) = @_;

  #print "In toggle_weighting\n";

  if ($check) {
    $$zero_weight_list->configure(-text => 'Only this anomaly');
    #print "check on\n";
  } else {
    $$zero_weight_list->configure(-text => 'Not this anomaly');
    #print "check off\n";

  }

}

#############################################################################
# this is called by the popup menu of anomalies to reweight to zero
# it resets the list of anomalies in anomaly_list
# it resets the list of anomalies in anomaly_detail_list
# it sets the weight of the selected anomaly to be zero
# it re-populates the anomaly_window_list


sub reweight_anomalies {
  my ($mysql_ref, $lab_ref, $chromosome_ref, $anomaly_list_ref, $anomaly_detail_list_ref, $zero_weight_ref, $check_ref) = @_;

  #print "In reweight_anomalies, zero_weight = $$zero_weight_ref check = $$check_ref\n";

  my ($type, $count) = split(/\s+/, $$zero_weight_ref);

  # name of the temporary 'view' weighting table
  my $view = "weight_$user";

  # delete whatever was in the anomaly list before
  $$anomaly_list_ref->delete( 0, "end" );

  # delete whatever was in the anomaly details list and results list before
  $$anomaly_detail_list_ref->delete( 0, "end" );
  @$results = ();

  if ($$check_ref) {		# only this anomaly
    # set the weight of all anomalies to be zero
    $$mysql_ref->do( qq{ UPDATE $view SET weight = 0; } );
    # reset the weight of the selected anomaly
    $$mysql_ref->do( qq{ UPDATE $view SET weight = 1  WHERE type = "$type"; } );
  } else {			# not this anomaly
    # set the weight of the selected anomaly to be zero
    $$mysql_ref->do( qq{ UPDATE $view SET weight = 0  WHERE type = "$type"; } );
  }

  # re-populate the anomaly_window_list
  &populate_anomaly_window_list($$mysql_ref, $$lab_ref, $$chromosome_ref, $anomaly_list_ref);



}

#############################################################################
# mark all anomalies in this 10 kb window as being inactive - they will not be displayed again

sub ignore_anomaly_window {

    my $anomaly_list = shift;
    my $anomaly_detail_list = shift;

    # if there is nothing in the list of displayed details, then return
    if ($$anomaly_detail_list->size == 0) {
      print "Nothing is in the details display, so no details will be ignored\n";
      return;
    }

    # set the status active=0 for all anomalies in this window
    my $command = "UPDATE anomaly SET active = 0 WHERE anomaly_id = ?;";
    my $sth = $mysql->prepare($command)
                or die "Couldn't prepare statement: " . $mysql->errstr;

    # go through the results list of anomaly details
    foreach my $anomaly (@{$results}) {
      my ($type, $clone, $clone_start, $clone_end, $chromosome_start, $chromosome_end, $sense, $thing_id, $thing_score, $explanation, $anomaly_id) = (@{$anomaly});

      # mark them as inactive 
      $sth->execute($anomaly_id);
      #print "Ignoring $type $clone:$clone_start..$clone_end anomaly_id $anomaly_id\n";
    }
    $sth->finish();

    # remove this from the displayed windows list
    my $selection = $$anomaly_list->curselection();
    if (defined $selection ) {
      my $details = $$anomaly_list->get( $selection );
      $$anomaly_list->delete($selection);
    } else {
      print "No current anomaly region is selected - but that's OK\n";
    }
    
    # remove all from the displayed details list
    $$anomaly_detail_list->delete(0, 'end');
    @$results = ();

}


#############################################################################
# anomaly finder

sub goto_anomaly_window {
  my $anomaly_window_list = shift;
  my $anomaly_detail_list = shift;
  my $coords_ref = shift;

  my $query;

  # delete whatever was in the anomaly details list before
  #print "delete what was in the anomaly details list before\n";
  $$anomaly_detail_list->delete( 0, "end" );
  @$results = ();

  my $selection = $$anomaly_window_list->curselection();
  my $clone = $anomaly_clone;

  # name of the temporary 'view' weighting table
  my $view = "weight_$user";

  # see if there is a clone specified, or whether we are picking a selection from the listbox
  if (defined $clone && $clone ne "" && defined $selection) {

    &confirm_message("AMBIGUOUS SELECTION", "You have selected both the clone $clone and a line from the list");
    return;

  } elsif (defined $clone && $clone ne "") {

    # display the whole clone
    &goto_location($clone, 1, 200000, '+', 0);

    # get and display the individual anomalies found in this clone
    # pull out all anomalies in this clone except those marked as active = 0 and those with zero-weighted anomaly types
    $query = qq{ SELECT a.type, a.clone, a.clone_start, a.clone_end, a.chromosome_start, a.chromosome_end, a.sense, a.thing_id, a.thing_score, a.explanation, a.anomaly_id FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type   WHERE a.clone = "$clone" AND a.active = 1 AND w.weight = 1 ORDER BY chromosome_start; };


  } elsif (defined $selection) {

    #print "get current selection\n";
    my $selected_value = $$anomaly_window_list->get( $selection );
    my @selected = split(/\s+/, $selected_value);

    # use the window * 10 K to get the chromosomal location to go to
    my $chromosome = $selected[0];
    my $window = $selected[7];
    my $window_start = $selected[7] * $WINDOW_SIZE;
    my $window_end = ($selected[7] * $WINDOW_SIZE) + $WINDOW_SIZE - 1;
    my $sense = $selected[3];

    # camace (and stlace) doesn't have the CHROMOSOME_* Sequence objects in them
    # so we now have to convert to clone coords
    my @clone_coords = ${$coords_ref}->LocateSpan("CHROMOSOME_$chromosome", $window_start, $window_end );
    &goto_location($clone_coords[0], $clone_coords[1], $clone_coords[2], $sense, 0);

    # insert the text '[Seen in this session]' at the end of the current select line (if it is not already there)
    if ($selected_value !~ /\[Seen/) {
      $$anomaly_window_list->delete($selection);
      $$anomaly_window_list->insert($selection, "$selected_value [Seen in this session]");
    }

    # get and display the individual anomalies found in this anomalies window
    # extract the new details from the database
    # pull out all anomalies in this window except those marked as active = 0 and those with zero-weighted anomaly types
    $query = qq{ SELECT a.type, a.clone, a.clone_start, a.clone_end, a.chromosome_start, a.chromosome_end, a.sense, a.thing_id, a.thing_score, a.explanation, a.anomaly_id FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type   WHERE a.chromosome = "$chromosome" AND a.window = $window AND a.sense = "$sense" and a.active = 1 AND w.weight = 1 ORDER BY chromosome_start };

  } else {

    &confirm_message("NO SELECTION", "No current selection\n"); 
    return;

  }				# end of test to see if looking at clone or list selection


  # get the query and present the details on the lower listbox
  my $db_query = $mysql->prepare ( $query );

  $db_query->execute();
  $results = $db_query->fetchall_arrayref;

#  my $len = $#{$results}+1;
#  print "Have $len rows for this window\n";
#  foreach my $result_row (@$results) {
#    print $result_row->[0] . "\t" .$result_row->[1] . "\t" .$result_row->[2] . "\t" .$result_row->[3] . "\t" .$result_row->[4] . "\t" .$result_row->[5] . "\t" .$result_row->[6] . "\t" .$result_row->[7] . "\t" .$result_row->[8]  . "\t" .$result_row->[9] . "\t" .$result_row->[10] . "\n";
#  }

  # merge similar types at the same overlapping region into a single line for better presentation
  my @boxes = &box_merge($results);

  foreach my $box (@boxes) {
    if (exists $box->{'deleted'}) {next;} # check if this box has not been deleted
    # insert the new details into the details list
    $$anomaly_detail_list->insert('end', 
				  $box->{'type'} . " score: " . 
				  $box->{'total_score'} . " " .
				  $box->{'clone'} .  ":" .
				  $box->{'clone_start'} . ".." .
				  $box->{'clone_end'} . " " .
				  $box->{'sense'} . " " .
				  $box->{'explanation'} .  " evidence: " .
				  $box->{'ID'});
    # print
    #  $box->{'type'} . " score: " . 
    #  $box->{'total_score'} . " " .
    #  $box->{'clone'} .  ":" .
    #  $box->{'clone_start'} . ".." .
    #  $box->{'clone_end'} . " " .
    #  $box->{'sense'} . " " .
    #  $box->{'explanation'} .  " evidence: " .
    #  $box->{'ID'} . "\n";
  }
}

#############################################################################

# merge the resulting regions of anomaly together to form boxed regions with scores
# derived from the sums of the alignment in the merged boxes
# returns a list of hash of details of the merged anomalies
# this stops the user from having to go through huge lists of nearly
# identical anomalies at about the same place

sub box_merge {
  my ($ref_anomalies) = @_;
 
  my @boxes = ();

  my $last_sense = '';

  foreach my $anomaly (@{$ref_anomalies}) {
    my ($type, $clone, $clone_start, $clone_end, $chromosome_start, $chromosome_end, $sense, $thing_id, $thing_score, $explanation) = (@{$anomaly});
    $last_sense = $sense;

    my $got_a_match = 0;

    foreach my $box (@boxes) {
      # find an existing box to merge to
      # if same type and has substantial overlap
      #print "$thing_id Thing $clone_start $clone_end $clone\n";
      #print "$thing_id Box   " . $box->{'clone_start'} . " " . $box->{'clone_end'} . " " . $box->{'clone'} . "\n";
      if ($box->{'clone'} eq $clone &&
	  $box->{'clone_start'} <= $clone_end && $box->{'clone_end'} >= $clone_start &&
	  $box->{'type'} eq $type &&
	  $box->{'explanation'} eq $explanation &&
	  ! exists $box->{'deleted'}) {
                                                                                                                                                       
                                                                                                                                                       
	#print "*** Box $box->{'type'} $box->{'ID'} $box->{'clone_start'}..$box->{'clone_end'}, matches $thing_id, $clone_start, $clone_end, $thing_score\n";
                                                                                                                                                       
	# add to the count of matches in this box
	$box->{'count'}++;
                                                                                                                                                       
	# add to the sum of the alignment scores
	$box->{'total_score'} += $thing_score;
                                                                                                                                                       
	# add to the protein IDs
	$box->{'ID'} .= " $thing_id";
                                                                                                                                                       
	# update start/end
	if ($box->{'clone_start'} > $clone_start) {$box->{'clone_start'} = $clone_start;}
	if ($box->{'clone_end'} < $clone_end) {$box->{'clone_end'} = $clone_end;}
	
	$got_a_match = 1;

	# check to see if we now need to merge two overlapping boxes
	my $past_first_box = 0;
	foreach my $other_box (@boxes) {
	  if (! $past_first_box) {
	    # see if this is our current $box
	    if ($other_box->{'type'} eq $type &&
		$other_box->{'clone_start'} == $box->{'clone_start'} &&
		$other_box->{'clone_end'} == $box->{'clone_end'} &&
		$other_box->{'clone'} eq $clone &&
		! exists $other_box->{'deleted'}
		) {
	      $past_first_box = 1;
	    }
	    next;
	  };

	  # now start checks for overlaps of boxes
	  if ($other_box->{'clone'} eq $box->{'clone'} &&
	      $other_box->{'clone_start'} <= $box->{'clone_end'} && $other_box->{'clone_end'} >= $box->{'clone_start'} &&
	      $other_box->{'type'} eq $box->{'type'} && 
	      $other_box->{'explanation'} eq $box->{'explanation'} &&
	      ! exists $other_box->{'deleted'}
	      ) {
	    $box->{'count'} += $other_box->{'count'};
	    $box->{'total_score'} += $other_box->{'total_score'};
	    $box->{'ID'} .= $other_box->{'ID'};
	    if ($box->{'clone_start'} > $other_box->{'clone_start'}) {$box->{'clone_start'} = $other_box->{'clone_start'}};
	    if ($box->{'clone_end'} < $other_box->{'clone_end'}) {$box->{'clone_end'} = $other_box->{'clone_end'}};
	    # delete the other box
	    #print "deleted $other_box->{'type'} $other_box->{'clone_start'} $other_box->{'clone_end'} \n";
	    $other_box->{'deleted'} = 1; # mark this box as deleted
	  }
	}
      }
    }				
                                                                                                                                                       
    # no existing boxes found, so start a new one
    if (! $got_a_match) {
      #print "New box for $thing_id, $clone_start, $clone_end, $thing_score\n";
      # want to store: $type, $clone, $clone_start, $clone_end, $chromosome_start, $chromosome_end, $sense, $thing_id, $thing_score, $explanation
      my $new_box = {};             # reference to hash
      $new_box->{'type'} = $type;
      $new_box->{'chromosome_start'} = $chromosome_start;
      $new_box->{'chromosome_end'} = $chromosome_end;
      $new_box->{'sense'} = $sense;
      $new_box->{'explanation'} = $explanation;
      $new_box->{'clone_start'} = $clone_start;
      $new_box->{'clone_end'} = $clone_end;
      $new_box->{'clone'} = $clone;
      $new_box->{'total_score'} = $thing_score;
      $new_box->{'ID'} = $thing_id;
      $new_box->{'count'} = 1;
      push @boxes, $new_box;
    }
  }

#  # now do a quick pass through to remove any instances of a single
#  # TREMBL protein mismatch because these occur all too often and we
#  # really get sick of seeing them and marking them as inactive
#  foreach my $box (@boxes) {
#    if ($box->{'count'}  == 1 && 
#	$box->{'type'} eq 'UNMATCHED_PROTEIN' && 
#	$box->{'ID'} =~ /TR:/) {
#      $box->{'deleted'} = 1;
#    }
#  }


  # sort by the start position
  if ($last_sense eq '+') {
    @boxes = sort { $a->{'chromosome_start'} <=> $b->{'chromosome_start'} } @boxes;
  } else {
    @boxes = sort { $b->{'chromosome_start'} <=> $a->{'chromosome_start'} } @boxes;
  }
                                                                                                                                            
  return @boxes;
}

#############################################################################
# mark this anomaly as being inactive - it will not be displayed again

sub ignore_anomaly {

  my $anomaly_list = shift;

  my @selection = $$anomaly_list->curselection();

  if (! @selection ) {
    &confirm_message("NO SELECTION", "No current selection\n"); 
    return;
  }

  foreach my $selection (@selection) {
    #print "Next selection = $selection\n";
    my @selected = split(/\s+/, $$anomaly_list->get( $selection ));
    #print "split line: @selected\n";
    # set the status.active=0 for this anomaly
    my $command = "UPDATE anomaly SET active = 0 WHERE anomaly_id = ?;";
    my $sth = $mysql->prepare($command)
	or die "Couldn't prepare statement: " . $mysql->errstr;

    # find the result to mark inactive in the results list of anomaly details
    # get the start and end coords for this box of anomaly details
    my $box_type = $selected[0];
    my $location = $selected[3];
    my ($box_clone, $start, $end) = ($location =~ /(\S+)\:(\d+)\.\.(\d+)/);

    #print "Box: $box_type $box_clone $start, $end\n";

    # go through the results list of anomaly details looking for ones in this box
    foreach my $anomaly (@{$results}) {
      my ($type, $clone, $clone_start, $clone_end, $chromosome_start, $chromosome_end, $sense, $thing_id, $thing_score, $explanation, $anomaly_id) = (@{$anomaly});
      #print "Detail: $type $clone $clone_start, $clone_end\n";
      if ($type eq $box_type && $box_clone eq $clone && $clone_start >= $start && $clone_end <= $end) {
	# mark this one as inactive 
	$sth->execute($anomaly_id);
	#print "Ignoring $type $clone $location anomaly_id $anomaly_id\n";
      }
    }
    $sth->finish();

    # remove this from the displayed list
    $$anomaly_list->delete($selection);

    # AND NOW BECAUSE WE HAVE REMOVED A SELECTION FROM THE LISTBOX
    # WE HAVE TO ADJUST THE VALUES OF OUR LIST OF SELECTIONS TO WORK ON 
    # IF THEY HAVE BEEN SHIFTED UP BY THE DELETION
    for (my $i = 0; $i <= $#selection; $i++) {
      if ($selection[$i] > $selection) {
	$selection[$i]--;
      }
    }
  }
}

#############################################################################
# anomaly finder - go to location with 10 bases extra around it

sub goto_anomaly {

  my $anomaly_list = shift;

  my $selection = $$anomaly_list->curselection();
  if (! defined $selection ) {
    &confirm_message("NO SELECTION", "No current selection\n"); 
    return;
  }
  my @selected = split(/\s+/, $$anomaly_list->get( $selection ));
  my $sense = $selected[4];
  my ($clone, $clone_start, $clone_end) = ($selected[3] =~ /(\S+):(\d+)\.\.(\d+)/); # split location as:F25H5:20321..21345
  &goto_location($clone, $clone_start, $clone_end, $sense, 10);

}

#############################################################################
# intron finder ( also used by blast hit finder) 

sub goto_intron
  {
    my $intron_list = shift;
    my $selection = $$intron_list->curselection();
    if (! defined $selection ) {
      &confirm_message("NO SELECTION", "No current selection\n"); 
      return;
    }
    my @selected = split(/\s+/, $$intron_list->get( $selection ));
    &goto_location($selected[0], $selected[1], $selected[2], '+', 500);
  }

############################################################################
# generic methods 
sub goto_location
  {
    my $seqobj  = shift;
    my $x       = shift;
    my $y       = shift;
    my $sense   = shift;
    my $zoomout = shift;

    if (!defined $zoomout) {
      $zoomout = 500;
    }
    $x -= $zoomout;
    $y += $zoomout;

    my $reversed_command = "";
    if (defined $sense && $sense eq '-') {
      $reversed_command = "; seqactions -rev_comp";
    }
    my $command = "xremote -remote 'gif seqget $seqobj -coords $x $y; seqdisplay $reversed_command'";

    print "command=$command\n";
    my $return_status = system($command);
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    }

  }

############################################################################
# look at the progress made in working through the scores of the windows

sub progress {

  my ($graphs_ref) = @_;

  # update the 'progress' table there is guarantee that this table
  # alredy exists, so check and create it if necessary.

  my $got_progress_table = 0;

  # first test to see if this progress table exists already
  my $query = qq{ SHOW TABLES; };
  my $db_query = $mysql->prepare ( $query );
  $db_query->execute();
  my $results = $db_query->fetchall_arrayref;
  foreach my $result_row (@$results) {
    if ($result_row->[0] eq 'progress') {
      $got_progress_table = 1;
    }
  }

  # if we have no progress table then create it
  if (! $got_progress_table) {
    $mysql->do("CREATE TABLE progress (date timestamp(8), chromosome varchar(3), centre varchar(2), over_10 int, over_5 int, over_2 int, over_1 int, over_half int, over_quarter int, under_quarter int);");
  }

  # update the progress table
  my @chromosomes  = qw(I II III IV V X);
  my @labs = qw(RW HX);

  foreach my $chr (@chromosomes) {
    foreach my $lab (@labs) {

      # assume all weights are 1 and get the scores of the windows
      my $query = qq{ SELECT SUM(a.thing_score) FROM anomaly AS a WHERE a.chromosome = "$chr" AND a.centre = "$lab" AND a.active = 1 GROUP BY window, sense; };
      my $db_query = $mysql->prepare ( $query );
      $db_query->execute();
      $results = $db_query->fetchall_arrayref;

      my $over_10 = 0;
      my $over_5 = 0;
      my $over_2 = 0;
      my $over_1 = 0;
      my $over_half = 0;
      my $over_quarter = 0;
      my $under_quarter = 0;

      foreach my $result_row (@$results) {

	# count the numbers in various score bins
	if ($result_row->[0] > 10) {
	  $over_10++;
	} elsif ($result_row->[0] > 5) {
	  $over_5++;
	} elsif ($result_row->[0] > 2.5) {
	  $over_2++;
	} elsif ($result_row->[0] > 1) {
	  $over_1++;
	} elsif ($result_row->[0] > 0.5) {
	  $over_half++;
	} elsif ($result_row->[0] > 0.25) {
	  $over_quarter++;
	} elsif ($result_row->[0] <= 0.25) {
	  $under_quarter++;
	}
      }

      #  print "\tinserting for $chr, $lab  over 10:        $over_10\n";
      #  print "\tinserting for $chr, $lab  5 to 10:        $over_5\n";
      #  print "\tinserting for $chr, $lab  2.5 to 5:       $over_2\n";
      #  print "\tinserting for $chr, $lab  1 to 2.5:       $over_1\n";
      #  print "\tinserting for $chr, $lab  0.5 to 1:       $over_half\n";
      #  print "\tinserting for $chr, $lab  0.25 to 0.5:    $over_quarter\n";
      #  print "\tinserting for $chr, $lab  less than 0.25: $under_quarter\n";

      # insert the values for this part of the chromosome into the 'progress' table
      $mysql->do(qq{ INSERT INTO progress (chromosome, centre, over_10, over_5, over_2, over_1, over_half, over_quarter, under_quarter) VALUES("$chr", "$lab", $over_10, $over_5, $over_2, $over_1, $over_half, $over_quarter, $under_quarter); });

    }
  }

  # draw graphs of progress
  if ($$graphs_ref == 1) {
    use Tk::Graph;

    # set up the window to hold the graphs
    my $D = MainWindow->new;
    $D->optionAdd('*BorderWidth' => 1); # make the style better with a thinner border

    my $frame1 = $D->Frame(-background =>'white')->pack();
    my $close_button = $frame1->Button(-text => "Close progress window",
				       -background => "bisque",
				       -command => [$D => 'destroy'],
				       )->pack ( -side => 'left',
						 -pady => '2',
						 -padx => '6',
						 -anchor => "w"
						 );


    # this frame wraps up all of the chromosomes and is scrollable
    use Tk::Pane;
    my $top = $D->Scrolled("Frame", 
			   -scrollbars => "osoe",
			   -height => 800,
			   -width  => 800,
			   )->pack(-side => 'top',
				   -fill => 'x',
				   );


    # get the data for each of the graphs
    foreach my $chr (@chromosomes) {

      my $chrom_canvas = $top->Canvas(-background => 'gray',
				  )->pack(-side => 'top',
							-fill => 'x'); 

      foreach my $lab (@labs) {

	# get the data with the count of days since 2 Oct 2006
	$query = qq{ SELECT DATEDIFF(p.date,'2006-10-02'), over_10, over_5, over_2, over_1, over_half, over_quarter, under_quarter FROM progress AS p WHERE p.chromosome = "$chr" AND p.centre = "$lab" ORDER BY 1 };
	my $db_query = $mysql->prepare ( $query );
	$db_query->execute();
	$results = $db_query->fetchall_arrayref;

	my @over_10=();
	my @over_5=();
	my @over_2=();
	my @over_1=();
	my @over_half=();
	my @over_quarter=();
	my @under_quarter=();
	my @dates=();

	foreach my $result_row (@$results) {
	  push @dates, $result_row->[0]; # get the date (days since 2 Oct 2006)
	  push @over_10, $result_row->[1];
	  push @over_5, $result_row->[2];
	  push @over_2, $result_row->[3];
	  push @over_1, $result_row->[4];
	  push @over_half, $result_row->[5];
	  push @over_quarter, $result_row->[6];
	  push @under_quarter, $result_row->[7];
	}
	
	# convert count value arrays to have one value per day
	@over_10 = &date_convert(\@over_10, \@dates);
	@over_5	= &date_convert(\@over_5, \@dates);
	@over_2 = &date_convert(\@over_2, \@dates);
	@over_1 = &date_convert(\@over_1, \@dates);
	@over_half = &date_convert(\@over_half, \@dates);
	@over_quarter = &date_convert(\@over_quarter, \@dates);
	@under_quarter = &date_convert(\@under_quarter, \@dates);
	
	# get data to draw in graphs	
	my $to_register = {
	  '> 10'  => [@over_10],
	  ' > 5'  => [@over_5],	# these spaces are to force the sort order for the key to be in this order
	  '  > 2.5'  => [@over_2],
	  '   > 1'  => [@over_1],
	  '    > 0.5'  => [@over_half],
	  '     > 0.25'  => [@over_quarter],
	  '      < 0.25'  => [@under_quarter],
	};

        # this seems to just set the last y data value to go to - weird!
	my $data = {
	  '> 10'  => $over_10[-1], # get the last value
	  ' > 5'  => $over_5[-1],  # these spaces are to force the sort order for the key to be in this order
	  '  > 2'  => $over_2[-1],
	  '   > 1'  => $over_1[-1],
	  '    > 0.5'  => $over_half[-1],
	  '     > 0.25'  => $over_quarter[-1],
	  '      < 0.25'  => $under_quarter[-1],
	};

        # place RW graphs to the left and HX to the right
	my $side = "left";
	if ($lab eq "HX") {$side = "right";}

	my $graph = $chrom_canvas->Graph(
					 -type	=> 'Line',
					 -max 	=> 700,
					 -wire     => 0, # no wire grid in background
					 -title    => "chromosome $chr $lab",
					 -linewidth => 2,
					 )->pack(
						 -side => $side,
						 );

	$graph->register($to_register);
	$graph->variable($data);

      }
    }
    MainLoop;
    


  } # end of graph stuff


}
############################################################################
# The progress values are stored in the database at irregular dates
# the Tk::Graph stuff doesn't seem to be able to deal with this sensibly
# therefore we have to create a repeated value for each day in between 
# the noted dates so that we get a graph that sensibly displays the
# data changing smoothly over time.
#
# @under_quarter = &date_convert(\@under_quarter, \@dates);


sub date_convert {
  my ($values_aref, $dates_aref) = @_;

  my @values = @{$values_aref};
  my @dates = @{$dates_aref};
  my @new_values;

  my $prev_day = -1;
  my $prev_value;
  my $value;
  foreach my $day (@dates) {
    $value = shift @values;
    #print "in day $day value $value\n";

    if ($prev_day == $day) {next;}
    if ($prev_day != -1) {
      for (my $i = $prev_day+1; $i < $day; $i++) {
	push @new_values, $prev_value;
	#print "out day $i value $prev_value\n";
      }
    }

    push @new_values, $value;
    #print "out day $day value $value\n";

    $prev_day = $day;
    $prev_value = $value;
  }

  return @new_values;

}
############################################################################

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
    $D->configure(-background => 'PaleGreen', # was cyan
		  -foreground => 'black');

    my $choice = $D->Show;
  }

############################################################################

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

############################################################################
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

############################################################################
# makes sure all CDSs are uppercase ( apart from the isoform or gc/tw part )
sub confirm_case
  {
    my $cds = shift;
    my ($last_char) = lc(chop $cds);
    $cds = (uc $cds) . "$last_char";
    return $cds;
  }

############################################################################
sub check_user
  {
    my $user = shift;
    return if( defined $user );
    my $name = `whoami`;
    chomp $name;
    if ( "$name" eq "wormpub" ) {
      &error_warning("WORMPUB","Please either run this as yourself or use the -user option. How else can you be blamed for your errors!");
      exit(0);
    } else {
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

-design     : Does not make Aceperl connection - dev tool for quicker startup

-mail       : This option emails blessed predictions to mt3 requesting a new gene ID

=item Intron finder

Uses the GFF check files in development_release to provide a list of introns.  Select and click to go to intron in FMAP

=item History Maker

presents to the user a simple box with a space to enter a CDS name and a button to make a history object.  

=item Genefinder / Twinscan Blesser

Enter current CDS name eg AC8.gc3 and click "Bless this gene".  This will create a new CDS with the correct name based on the "worm_genes" class. If you have specified the -mail option a new gene ID will automatically be requested.

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

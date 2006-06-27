#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib $ENV{'CVS_DIR'};
use Ace;
use Tk;
use Getopt::Long;
require Tk::Dialog;
use Coords_converter;
use DBI;


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


my $results;		# the globally-accesible anomaly details that are summarised in the anomalies detail window

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


  my $coords = Coords_converter->invoke;

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
if ( $anomaly && $chromosome ) {

  my $anomaly_detail_list;


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
					     -command => [\&goto_anomaly_window, \$anomaly_list, \$anomaly_detail_list]
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

# connect to the mysql database
  $mysql = &connect_to_database($lab);

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

  my $ignore_all_window = $anomaly_details->Button ( -text    => "Ignore ALL these anomalies!",
						     -background => "red",
						     -command => [\&ignore_anomaly_window, \$anomaly_list, \$anomaly_detail_list]
					      )->pack ( -side => 'left',
							-pady => '2',
							-padx => '6',
							-anchor => "center"
							);

  my $go_to_anomaly = $anomaly_details->Button ( -text    => "Go to this anomaly",
						 -background => "aquamarine3",
						 -command => [\&goto_anomaly, \$anomaly_detail_list]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);

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

    $dbsn = "DBI:mysql:database=worm_anomaly;host=ecs1f";
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
# get the information on the available (still active) windows sorted
# by the sum of scores of anomalies in them
# populate_anomaly_window_list($mysql, \$anomaly_list);

sub populate_anomaly_window_list {
  my ($mysql, $lab, $chromosome, $anomaly_list_ref) = @_;
  
  # get the highest scoring 10 kb windows of one sense or the other sorted by score descending
  # group by window, sense order by 2 desc - means make the
  # SUM(thing_score) sum up rows within distinct window and distinct
  # sense and sort the output by descending SUM(thing_score)
  my $query = qq{ SELECT a.window, SUM(a.thing_score), a.sense, a.clone FROM anomaly AS a WHERE a.chromosome = "$chromosome" AND a.centre = "$lab" AND a.active = 1 GROUP BY window, sense ORDER BY 2 DESC };

  #print "chromosome=$chromosome\n";
  #print "lab=$lab\n";
  #print "query=$query\n";
  #print "\n";

  my $db_query = $mysql->prepare ( $query );

  $db_query->execute();

  $results = $db_query->fetchall_arrayref;

  foreach my $result_row (@$results) {
    # clone, sense, score, window
    $$anomaly_list_ref->insert('end', $chromosome . " " . $result_row->[3] . " Sense: " . $result_row->[2] . " Score: " . $result_row->[1] . " ID: " . $result_row->[0] );
  }

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

  my $query;

  # delete whatever was in the anomaly details list before
  #print "delete what was in the anomaly details list before\n";
  $$anomaly_detail_list->delete( 0, "end" );
  @$results = ();

  my $selection = $$anomaly_window_list->curselection();
  my $clone = $anomaly_clone;

  # see if there is a clone specified, or whether we are picking a selection from the listbox
  if (defined $clone && $clone ne "" && defined $selection) {

    &confirm_message("AMBIGUOUS SELECTION", "You have selected both the clone $clone and a line from the list");
    return;

  } elsif (defined $clone && $clone ne "") {

    # display the whole clone
    &goto_location($clone, 1, 200000, '+', 0);

    # get and display the individual anomalies found in this clone
    # pull out all anomalies in this clone except those marked as active = 0
    $query = qq{ SELECT type, clone, clone_start, clone_end, chromosome_start, chromosome_end, sense, thing_id, thing_score, explanation, anomaly_id FROM anomaly WHERE clone = "$clone" AND active = 1 ORDER BY chromosome_start };

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
    &goto_location("CHROMOSOME_$chromosome", $window_start, $window_end, $sense, 0);

    # insert the text '[Seen in this session]' at the end of the current select line (if it is not already there)
    if ($selected_value !~ /\[Seen/) {
      $$anomaly_window_list->delete($selection);
      $$anomaly_window_list->insert($selection, "$selected_value [Seen in this session]");
    }

    # get and display the individual anomalies found in this anomalies window
    # extract the new details from the database
    # pull out all anomalies in this window except those marked as active = 0
    $query = qq{ SELECT type, clone, clone_start, clone_end, chromosome_start, chromosome_end, sense, thing_id, thing_score, explanation, anomaly_id FROM anomaly WHERE chromosome = "$chromosome" AND window = $window AND sense = "$sense" and active = 1 ORDER BY chromosome_start };

  } else {

    &confirm_message("NO SELECTION", "No current selection\n"); 
    return;

  }				# end of test to see if looking at clone or listbox selection


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

  foreach my $anomaly (@{$ref_anomalies}) {
    my ($type, $clone, $clone_start, $clone_end, $chromosome_start, $chromosome_end, $sense, $thing_id, $thing_score, $explanation) = (@{$anomaly});

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

  # sort by the reverse total_score
  @boxes = sort { $a->{'chromosome_start'} <=> $b->{'chromosome_start'} } @boxes;
                                                                                                                                                       
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

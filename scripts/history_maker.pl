#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib $ENV{'CVS_DIR'};
use Ace;
use Tk;
use Tk::SplitFrame;
use Getopt::Long;
require Tk::Dialog;
use Coords_converter;
use DBI;
use Wormbase;

 
######################################
# variables and command-line options #
######################################
 
my ($help, $debug, $test, $verbose, $store, $wormbase);
my $reference;
my $cdatabase;
my $design;
my $chromosome;
my $blast;
my $user;
my $mail;
my $intron;
my $blesser;
my $anomaly;
my $lab;
my $species;
my $display_by_clones;
my $addevidence;
my $addlastreviewed;
my $addrem;
my $ncRNA;
my $nonRNA;
my $pseudo;
my $cds_standard;
my $brugia;
my $cleangene;
my $clone;
my $v;

GetOptions (
	    "addevidence"       => \$addevidence,
	    "addlastreviewed"   => \$addlastreviewed,
	    "remark"            => \$addrem,
	    "ncrna"             => \$ncRNA,
	    "non"               => \$nonRNA,
	    "pseudo"            => \$pseudo,
	    "anomaly"           => \$anomaly,
	    "blast=s"           => \$blast,
	    "blesser"           => \$blesser,
            "brugia:s"          => \$brugia, # stores version number of brugia beta database - hard defines useful options.
	    "version:s"         => \$v, 
	    "chromosome:s"      => \$chromosome,
            "cleangene"         => \$cleangene,
	    "clone"             => \$clone,
            "debug=s"           => \$debug,
	    "design"            => \$design,
	    "display_by_clones" => \$display_by_clones,
	    "intron"            => \$intron,
	    "lab=s"             => \$lab, # RW or HX
	    "mail"              => \$mail,
	    "source:s"          => \$reference,
	    "curationdb:s"      => \$cdatabase,
	    "species:s"         => \$species,
            "store:s"           => \$store,
            "test"              => \$test,
	    "user:s"            => \$user,
            "verbose"           => \$verbose,
	    "help|h"            => sub { system("perldoc $0"); exit(0);}
	   );

 

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			   );
}
       
# keep things quiet
$debug = $user if (! defined $debug);

# Species is important now, please specify a species with -species

if (!defined $species) {
  print "Warning: Defaulting to species elegans do you want to continue? (y or n)\n";
  my $answer=<STDIN>;
  if ($answer eq "y\n") {
  }
  if ($answer eq "n\n") {
    die "Failed to launch as you need a species, please use -species <species> next time!\n";
  }
}


#################################


# you need to specify a curationdb
die "You need to specify a curation database with the -curationdb option\n" if (!defined $cdatabase);

# If using blast option you need the blast file to exist
die "$blast doesnt exist !\n" if ($blast and !(-e $blast));

# This is the database used a reference for generating histories NOT the one that will be changed
my $rdatabase = $reference ? $reference : glob("/.automount/evs-users2/root/wormpub/camace_orig");

# assume the lab is HX if not specified
if (! defined $lab || $lab eq "") {
  $lab = "HX";
} 
my $version;

# pass path to latest version of wormbase and set up brugia specifics if selected
if (defined $brugia) {
  $version = $brugia;
  $blesser = "1";
  $addevidence = "1";
  $addlastreviewed = "1";
  $clone = "1";
  $cleangene = "1";
  $addrem = "1";
}
elsif (defined $v) {
  $version = $v;
  my $cdbver = &get_history_version($wormbase->database('current'));
  $cdbver++;
  if ($v ne $cdbver) {print "\n\nWarning -version specified is not equal to $cdbver as expected is this correct?\n";}
}
else {
  $version = &get_history_version($wormbase->database('current'));
}

#If user is supplied grab a WBPersonID
if (! defined $user || $user eq "") { 
  # file where temp ace files written
  $user = "wormpub";
  # give a user warning
  print "No user was supplied so added WBPerson info cannot be included in any .ace output.\n";
}

my $person;
if (defined $user){ 
  if ($user eq 'pad') {
    $person = 'WBPerson1983';
  } elsif ($user eq 'gw3') {
    $person = 'WBPerson4025';
  } elsif ($user eq 'mh6') {
    $person = 'WBPerson4055';
  } elsif ($user eq 'jl16') {
    $person = 'WBPerson28994';
  }
}


#history nomenclature temp fix
my $wormpep_prefix  = $wormbase->wormpep_prefix;
$wormpep_prefix= lc($wormpep_prefix);

# japonica uses 'jp' not 'ja' to name the history objects - blame Phil for this!
if ($wormpep_prefix eq 'ja') {$wormpep_prefix = 'jp'}
if ($wormpep_prefix eq 'cn') {$wormpep_prefix = 'np'}

my $form_cds;			# cds variable from form
my $form_gene;			# gene variable from blesser form
my $form_gene2;                 # gene variable from blesser form
my $form_gene3;                 # gene variable from evidence form
my $form_gene4;                 # gene variable from clean_gene form
my $form_gene_rem;              # gene variable from remark form
my $form_gene_last_rev;         # gene variable from Last_reviewed form
my $form_rna;                   # gene variable from ncRNAwork form
my $form_nonRNA;                # gene variable from non_coding_transcript form
my $form_pseudo;                # pseudogene variable from Pseudogene form
my $clone_form;                 # cds variable for clone form
my $clone_form_sug;             # cds variable for clone form
my $anomaly_clone;		# clone variable from anomaly form


my $tidy_dir = glob("$cdatabase/.history_maker"); # changed from storing temp files under ~ to $cdatabase as Ele was having trouble making these files.
mkdir $tidy_dir, 0777;
system("rm -f $tidy_dir/*"); # clean out old files
my $session_file = "$tidy_dir/history_session.$version";

my $mysql;			# mysql database handle

# set up Tk interface
my $main_gui = MainWindow->new();
    my $SplitFrame = $main_gui->SplitFrame ('-orientation' => 'vertical',
        '-trimcolor'  => '#c7c7c7',
        '-background'  => 'white',
        '-sliderposition' => 250,
        '-borderwidth' => 2,
        '-sliderwidth' => 7,
        '-relief'      => 'sunken',
        '-height'      => 110,
        '-width'       => 100,
        '-padbefore'   => 0,
        '-padafter'    => 0
);

    # Values shown above are defaults.

    my $LeftLabel = $SplitFrame->Label ('-text' => '');

    my $RightLabel = $SplitFrame->Label ('-text' => '');

    $SplitFrame->pack (-expand => 'true', -fill => 'both');

    $SplitFrame->configure ('-sliderposition' => 400);




# Main window
$main_gui->configure(
		     -title => "Curation Tool for ${version}",
		     -background => 'darkgrey' # was beige
		    );
my $gui_width = 500;
$gui_width += 300 if ($anomaly or $blesser or $clone);
my $gui_height = 50; #modified 200
$gui_height += 200 if ($chromosome or $anomaly);
$gui_height += 200 if $blast;
$gui_height += 110 if $blesser;
$gui_height += 110 if $clone;
$gui_height += 130 if $addevidence;
$gui_height += 130 if $addlastreviewed;
$gui_height += 50 if ($cleangene or $ncRNA);
$gui_height += 130 if ($pseudo);
$gui_height += 80 if ($nonRNA or $addrem);
$gui_height += 300 if $anomaly;
$main_gui->geometry("${gui_width}x$gui_height");

# this is the value used to construct the window IDs in the script 'find_anomalies.pl'
my $WINDOW_SIZE = 10000;

my $chromosome_list;
if ($chromosome) {
  my @chroms = split /,/, $chromosome;
  foreach my $chrom (@chroms) {
    $chrom =~ s/CHROMOSOME_//;
    $chrom = "'".$chrom."'";
  }
  $chromosome_list = join(',', @chroms);
}

my $results; # the globally-accesible anomaly details that are summarised in the anomalies detail window
my $check;			# state of the check button for the weights


################################################################

# history_maker

my $his_maker = $LeftLabel->Frame( -background => "black", # was lightcyan
				  -height     => "400",
				  -label      => "History Maker ($wormpep_prefix)",
				  -relief     => "raised",
				  -foreground => 'whitesmoke', # new
				  -borderwidth => 5,
				)->pack( -pady => "5", #modified
					 -fill => "x"
				       );

# Reference database label
my $db_lbl = $his_maker->Label( -text => "Source: $rdatabase",
				-background => 'black', # was lightcyan
				-foreground => 'lightgrey' # was black
			      )->pack( -pady => '3'
				     );
# CDS entry widgets
my $cds_lbl = $his_maker->Label( -text => ' CDS  ID',
				 -background => 'black',	# was lightcyan
				 -foreground => 'whitesmoke'	# was black
			       )->pack(-pady => '6',
				       -padx => '6',
				       -side => 'left',
				      );
my $cds_val = $his_maker->Entry( -width => '10',
				 -background => 'whitesmoke',
				 -textvariable=> \$form_cds,
			       )->pack(-side => 'left',
				       -pady => '5',
				       -padx => '5'
				      );

#make Return and Enter submit CDS 
$cds_val->bind("<Return>",[ \&make_history]);
$cds_val->bind("<KP_Enter>",[ \&make_history]);

# Clear CDS entry button
my $clear = $his_maker->Button( -text => "Clear",
				-command => [\&clear]
			      )->pack(-side => 'right',
				      -pady => '2',
				      -padx => '6',
				      -anchor => "e"
				     );

#Make history button
my $make = $his_maker->Button( -text => "Make History",
			       -command => [\&make_history]
			     )->pack(-side => 'left',
				     -pady => '2',
				     -padx => '6',
				     -anchor => "w"
				    );



###########################################################
##   Prediction blesser
## works with genefinder, twinscan and the new Aggregates.

my $gene_val;
my $proposed_name;
if ($blesser) {

  my $gene_blesser = $main_gui->Frame( -background => "LightSteelBlue2", # was PaleTurquoise green
				       -height     => "400",
				       -label      => "Prediction blesser",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);

  # Reference database label
  my $db_lbl = $gene_blesser->Label( -text => "Source: $rdatabase",
				     -background => 'LightSteelBlue2',
				     -foreground => 'black' # was black
				   )->pack( -pady => '3'
					  );

  # CDS entry widgets
  my $gene_lbl = $gene_blesser->Label( -text => 'PRED ID',
				       -background => 'LightSteelBlue2', # was PaleTurquoise green
				       -foreground => 'black'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $gene_val = $gene_blesser->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_gene,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit CDS 
  $gene_val->bind("<Return>",[ \&bless_prediction]);
  $gene_val->bind("<KP_Enter>",[ \&bless_prediction]);
  

  # Bless button
  my $bless = $gene_blesser->Button( -text => "Bless This CDS",
				     -command => [\&bless_prediction]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "e"
					     );
  # Clear CDS entry button
  my $clear_gene = $gene_blesser->Button( -text => "Clear",
					  -command => [\&clear_gene]
					  )->pack(-side => 'left',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );



  # Supplying a chosen name for blessing Isoforms.
  my $gene_pro = $gene_blesser->Label( -text => 'Proposed Name(optional)',
					-background => 'LightSteelBlue2',
					-foreground => 'black'
				      )->pack(-side => 'left',
					      -pady => '2',
					      -padx => '6',
					      -anchor => "w"
					     );

  $proposed_name = $gene_blesser->Entry( -width => '10',
					 -background => 'whitesmoke',
					 -textvariable=> \$form_gene2,
				       )->pack(-side => 'left',
					       -pady => '2',
					       -padx => '6',
					       -anchor => "w"
					      );
  
  # Clear CDS entry button
  my $clear_proposed = $gene_blesser->Button( -text => "Clear",
					      -command => [\&clear_proposed]
					    )->pack(-side => 'right',
						    -pady => '2',
						    -padx => '6',
						    -anchor => "w"
						   );
  
##############################################


}




###########################################################
#blast hit locator frame
if ( $blast ) {
  my $blast_find = $main_gui->Frame( -background => "whitesmoke", # was purple plum
				     -label      => "Blast hit locator",
				     -relief     => "raised",
				     -borderwidth => 5,
				   )->pack( -pady => "5",
					    -fill => "x"
					  );
  # Reference database label
  my $db_lbl = $blast_find->Label( -text => "Source: $blast",
				     -background => 'whitesmoke',
				     -foreground => 'black'
				   )->pack( -pady => '3'
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

# addevidence to a CDS
###########################################################

my $cdswork;
if ($addevidence) {
  my $gene_evidence = $RightLabel->Frame( -background => "IndianRed4", #was LightGreen
				       -height     => "400",
				       -width      => "600",
				       -label      => "Populate Top Evidence hash",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $gene_evidence->Label( -text => "Source: $cdatabase",
				      -background => 'IndianRed4',
				      -foreground => 'black'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $CDS_lbl = $gene_evidence->Label( -text => ' MOD  ID',
				       -background => 'IndianRed4', #was LightGreen
				       -foreground => 'black'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $cdswork = $gene_evidence->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_gene3,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit CDS 
  $cdswork->bind("<Return>",[ \&add_evidence]);
  $cdswork->bind("<KP_Enter>",[ \&bless_evidence]);
  

  # Add Evidence button
  my $evi = $gene_evidence->Button( -text => " Tag This Model",
				     -command => [\&add_evidence]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_evi = $gene_evidence->Button( -text => "Clear",
					  -command => [\&clear_evi]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}
######### end add evidence


# add Last_reviewed to a CDS
###########################################################

my $lastrev;
if ($addlastreviewed) {
  my $last_reviewed = $RightLabel->Frame( -background => "Sandy Brown",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Populate Last_reviewed tag",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $last_reviewed->Label( -text => "Source: $cdatabase",
				      -background => 'Sandy Brown',
				      -foreground => 'black'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $CDS_lbl = $last_reviewed->Label( -text => ' MOD  ID',
				       -background => 'Sandy Brown',
				       -foreground => 'black'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $lastrev = $last_reviewed->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_gene_last_rev,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit CDS 
  $lastrev->bind("<Return>",[ \&add_last_rev]);
  $lastrev->bind("<KP_Enter>",[ \&add_last_rev]);
  

  # Add Evidence button
  my $rev = $last_reviewed->Button( -text => " Reviewed This Model",
				     -command => [\&add_last_rev]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_last_rev = $last_reviewed->Button( -text => "Clear",
					  -command => [\&clear_last_rev]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}
######### end add Last_reviewed


# addremark stub to a CDS
###########################################################

my $remwork;
if ($addrem) {
  my $gene_rem = $LeftLabel->Frame( -background => "Darkgreen",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Create a Remark stub",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $gene_rem->Label( -text => "Source: $cdatabase",
				      -background => 'Darkgreen',
				      -foreground => 'whitesmoke'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $CDS_lbl = $gene_rem->Label( -text => ' MOD  ID',
				       -background => 'Darkgreen', #was LightGreen
				       -foreground => 'whitesmoke'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $remwork = $gene_rem->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_gene_rem,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit CDS 
  $remwork->bind("<Return>",[ \&add_rem]);
  $remwork->bind("<KP_Enter>",[ \&bless_rem]);
  

  # Add Remark button
  my $remark = $gene_rem->Button( -text => " Add Remark",
				     -command => [\&add_rem]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_rem = $gene_rem->Button( -text => "Clear",
					  -command => [\&clear_rem]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}
######### end add remark stub

################ start add ncRNA tags ##############

my $rnawork;
if ($ncRNA) {
  my $ncrna_data = $RightLabel->Frame( -background => "#8470ff",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Populate standard ncRNA tags",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5",
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $ncrna_data->Label( -text => "Source: $cdatabase",
				      -background => '#8470ff',
				      -foreground => 'whitesmoke'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $ncRNA_lbl = $ncrna_data->Label( -text => ' RNA ID',
				       -background => '#8470ff', #was LightGreen
				       -foreground => 'whitesmoke'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $rnawork = $ncrna_data->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_rna,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit 
  $rnawork->bind("<Return>",[ \&add_ncrna_data]);
  $rnawork->bind("<KP_Enter>",[ \&add_ncrna_data]);
  

  # Add Evidence button
  my $data = $ncrna_data->Button( -text => "Populate Tags",
				     -command => [\&add_ncrna_data]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_data = $ncrna_data->Button( -text => "Clear",
					  -command => [\&clear_data]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}

######### end add ncRNA tags ##############

################ start add Pseudogene tags ##############

my $pseudowork;
if ($pseudo) {
  my $pseudo_data = $LeftLabel->Frame( -background => "#ffcc00",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Populate standard Pseudogene tags",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5",
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $pseudo_data->Label( -text => "Source: $cdatabase",
				      -background => '#ffcc00',
				      -foreground => 'whitesmoke'
				    )->pack( -pady => '3'
					   );

  # Pseudogene entry widgets
  my $pseudo_lbl = $pseudo_data->Label( -text => ' Pseudo ID',
				       -background => '#ffcc00', #was LightGreen
				       -foreground => 'whitesmoke'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $pseudowork = $pseudo_data->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_pseudo,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit 
  $pseudowork->bind("<Return>",[ \&add_pseudo_data]);
  $pseudowork->bind("<KP_Enter>",[ \&add_pseudo_data]);
  

  # Add Evidence button
  my $data = $pseudo_data->Button( -text => "Populate Tags",
				     -command => [\&add_pseudo_data]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear Pseudogene entry button
  my $clear_data = $pseudo_data->Button( -text => "Clear",
					  -command => [\&clear_pseudo_data]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}

######### end add Pseudogene tags ##############








################ start add non_coding_transcript tags ##############

my $nonrnawork;
if ($nonRNA) {
  my $nonrna_data = $RightLabel->Frame( -background => "#DF013A",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Populate standard non_coding_transcript  tags",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5",
						-fill => "x"
						);
  # Reference database label
  my $db_lbl = $nonrna_data->Label( -text => "Source: $cdatabase",
				      -background => '#DF013A',
				      -foreground => 'whitesmoke'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $nonRNA_lbl = $nonrna_data->Label( -text => ' RNA ID',
				       -background => '#DF013A', #was LightGreen
				       -foreground => 'whitesmoke'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $nonrnawork = $nonrna_data->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_nonRNA,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit 
  $nonrnawork->bind("<Return>",[ \&add_nonrna_data]);
  $nonrnawork->bind("<KP_Enter>",[ \&add_nonrna_data]);
  

  # Add Evidence button
  my $data = $nonrna_data->Button( -text => "Populate Tags",
				     -command => [\&add_nonrna_data]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_nonRNAdata = $nonrna_data->Button( -text => "Clear",
					  -command => [\&clear_nonRNAdata]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}

######### end add non_coding_transcript tags ##############



###### clean gene
###########################################################

my $genework;
if ($cleangene) {
  my $gene_clean = $LeftLabel->Frame( -background => "LightGreen",
				       -height     => "400",
				       -width      => "600",
				       -label      => "Clean Unwanted CDSs from loci",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);

  # Reference database label
  my $db_lbl = $gene_clean->Label( -text => "Source: $cdatabase",
				      -background => 'LightGreen',
				      -foreground => 'black'
				    )->pack( -pady => '3'
					   );

  # CDS entry widgets
  my $TAG_lbl = $gene_clean->Label( -text => 'GENE ID',
				       -background => 'LightGreen',
				       -foreground => 'black'
				       )->pack(-pady => '6',
					       -padx => '6',
					       -side => 'left',
					       );

  $genework = $gene_clean->Entry( -width => '10',
				       -background => 'whitesmoke',
				       -textvariable=> \$form_gene4,
				       )->pack(-side => 'left',
					       -pady => '5',
					       -padx => '5'
					       );

  # make Return and Enter submit CDS 
  $genework->bind("<Return>",[ \&clean_gene]);
#  $genework->bind("<KP_Enter>",[ \&bless_gene]);
  

  # Add Gene button
  my $evi = $gene_clean->Button( -text => "Clean This Gene",
				     -command => [\&clean_gene]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					     );
  # Clear CDS entry button
  my $clear_gen = $gene_clean->Button( -text => "Clear",
					  -command => [\&clear_gen]
					  )->pack(-side => 'right',
						  -pady => '2',
						  -padx => '6',
						  -anchor => "e"
						  );
}
######### end clean gene



###### clone gene
###########################################################


my ($clonework, $cloneprop);
if ($clone) {
  my $gene_clone = $main_gui->Frame( -background => "grey89",
				     -height     => "400",
				     -width      => "600",
				     -label      => "Clone a CDS/Transcript/Pseudogene from the curation database\nSource - $rdatabase",
				     -relief     => "raised",
				     -borderwidth => 5,
				   )->pack( -pady => "5", #modified
					    -fill => "x"
					  );
  # CDS entry widgets
  my $TAG_lbl = $gene_clone->Label( -text => '     ID     ',
				    -background => 'grey89',
				    -foreground => 'black'
				  )->pack(-pady => '6',
					  -padx => '6',
					  -side => 'left',
					 );
  
  $clonework = $gene_clone->Entry( -width => '10',
				   -background => 'whitesmoke',
				   -textvariable=> \$clone_form,
				 )->pack(-side => 'left',
					 -pady => '5',
					 -padx => '5'
					);




  # make Return and Enter submit CDS 
  $clonework->bind("<Return>",[ \&clone_gene]);
  # Add Gene button
  my $evi = $gene_clone->Button( -text => "Clone This Model",
				 -command => [\&clone_gene]
			       )->pack(-side => 'left',
				       -pady => '2',
				       -padx => '6',
				       -anchor => "w"
				      );
  # Clear CDS entry button
  my $clone_gen = $gene_clone->Button( -text => "Clear",
				       -command => [\&clear_cloned]
				     )->pack(-side => 'left',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "e"
					    );

############## proposed Name optional

  my $new_iso_lbl = $gene_clone->Label( -text => 'Proposed Name(optional)',
					-background => 'grey89',
					-foreground => 'black'
				      )->pack(-side => 'left',
					      -pady => '2',
					      -padx => '6',
					      -anchor => "w"
					     );

  $cloneprop = $gene_clone->Entry( -width => '10',
				   -background => 'whitesmoke',
				   -textvariable=> \$clone_form_sug,
				 )->pack(-side => 'left',
					 -pady => '2',
					 -padx => '6'
					);

  # Clear proposed name button
  my $clone_prop = $gene_clone->Button( -text => "Clear",
				       -command => [\&clear_cloned_prop]
				     )->pack(-side => 'right',
					     -pady => '2',
					     -padx => '6',
					     -anchor => "w"
					    );
}
######### end clone gene




###########################################################
# anomalies database locator frame
if ( $anomaly ) {

#  if (! defined $chromosome || $chromosome eq "") {die "Must specify -chromosome with -anomaly\n";}
  if (! defined $user || $user eq "") {die "Must specify -user with -anomaly\n";}

  my $anomaly_detail_list;
  my $zero_weight_label;

  my $coords = Coords_converter->invoke(undef, undef, $wormbase);


  my $anomaly_find = $main_gui->Frame( -background => "whitesmoke", # was cyan wheat

				       -label      => "Anomaly Region",
				       -relief     => "raised",
				       -borderwidth => 5,
				       )->pack( -pady => "5", #modified
						-fill => "x"
						);

  my $anomaly_list = $anomaly_find->Scrolled("Listbox", -scrollbars => "osoe",
					     -selectmode => "single",
					     -height     => "5",
					     -width      => "100"
					     )->pack();


  my $go_to_window = $anomaly_find->Button ( -text    => "Go to anomaly region (or clone)",
					     -background => "SlateGray2", #was aquamarine3
					     -command => [\&goto_anomaly_window, \$anomaly_list, \$anomaly_detail_list, \$coords]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);


# clone entry widgets
  my $anomaly_clone_lbl = $anomaly_find->Label( -text => 'Clone',
					     -background => 'whitesmoke', #was wheat
					     -foreground => 'black'
					     )->pack(-pady => '6',
						     -padx => '6',
						     -side => 'left',
						     );

  my $clone_val = $anomaly_find->Entry( -width => '10',
				   -background => 'whitesmoke',
				   -textvariable=> \$anomaly_clone,
				   )->pack(-side => 'left',
					   -pady => '5',
					   -padx => '5'
					   );

  # start the re-weighting menu stuff
  my $zero_weight;		# variable set by pop-up menu selection

  my $check_button = $anomaly_find->Checkbutton( -text => '',
						 -background => 'whitesmoke', #was wheat
						 -variable => \$check,
						 -command =>  [\&toggle_weighting, \$zero_weight_label],
						 )->pack(-pady => '6',
							 -padx => '6',
							 -side => 'left',
							 );

  $zero_weight_label = $anomaly_find->Label( -text => 'Not this anomaly',
					     -background => 'whitesmoke', #was wheat
					     -foreground => 'black'
					     )->pack(-pady => '6',
						     -padx => '0',
						     -side => 'left',
						     );

# test for popup menu for Anthony's anomaly weighting
# see: http://rcswww.urz.tu-dresden.de/CS/perl/modules/Tk/Tk-BrowseEntry.htm
  require Tk::BrowseEntry;
  my $zero_weight_list = $anomaly_find->BrowseEntry(-label => '',
						    -background => 'whitesmoke', #was wheat
						    -listwidth => '350',
						    -browsecmd => [\&reweight_anomalies, \$mysql, \$lab, \$chromosome_list, \$anomaly_list, \$anomaly_detail_list, \$zero_weight, \$check],
						    -variable => \$zero_weight
						    )->pack(-pady => '6',
							    -padx => '2',
							    -side => 'left',
							    );

  # create button to run Progress display
  my $display_graphs = 1;
  my $progress_button = $anomaly_find->Button ( -text    => "Progress",
					     -background => "IndianRed4", #was bisque
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
  &populate_zero_weight_list($mysql, $lab, $chromosome_list, \$zero_weight_list);


# populate $anomaly_list here from database
  &populate_anomaly_window_list($mysql, $lab, $chromosome_list, \$anomaly_list);



  my $anomaly_details = $main_gui->Frame( -background => "whitesmoke", # was magenta wheat

					  -label      => "Anomaly details in the selected region",
					  -relief     => "raised",
					  -borderwidth => 5,
					  )->pack( -pady => "5", #modified
						   -fill => "x"
						   );

  $anomaly_detail_list = $anomaly_details->Scrolled("Listbox", -scrollbars => "osoe",
						    -selectmode => "extended",
						    -height     => "12",
						    -width      => "100"
						    )->pack();


  my $ignore_this_anomaly = $anomaly_details->Button ( -text    => "Ignore the selected anomalies",
						       -background => "IndianRed4", #was orange
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
#						     -background => "IndianRed4", #was red
#						     -command => [\&ignore_anomaly_window, \$anomaly_list, \$anomaly_detail_list]
#					      )->pack ( -side => 'left',
#							-pady => '2',
#							-padx => '6',
#							-anchor => "center"
#							);

  my $go_to_anomaly = $anomaly_details->Button ( -text    => "Go to this anomaly",
						 -background => "SlateGray2", #was aquamarine3
						 -command => [\&goto_anomaly, \$anomaly_detail_list]
					      )->pack ( -side => 'right',
							-pady => '2',
							-padx => '6',
							-anchor => "e"
							);


  # calculate the progress so far for all chromosomes
  # don't display any graph
  my $dont_display_graphs = 0;
  &progress(\$dont_display_graphs);

}

##############################################################
#
# end of setting up frames
#
##############################################################

$user = &check_user unless $user;

# initial aceperl connection to databases
my $db = Ace->connect(-path => $rdatabase) unless $design;
my $cdb;
if ($cleangene || $ncRNA || $nonRNA){$cdb = Ace->connect(-path => $cdatabase) unless $design;}

# create GUI
MainLoop();

$db->close unless $design;
$cdb->close unless $design;

if ($anomaly) {
  $mysql->disconnect || die "error disconnecting from database", $DBI::errstr;
}

exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub add_last_rev
  {
    my $refgene = $form_gene_last_rev;
    return unless $refgene;

    if (!defined $person) {return &error_warning("You have not specified -user on the command line.")}

    my $class;

    &generate_message("CHECK", "Please save your database before clicking OK");
    $cdb->close();
    $cdb = Ace->connect(-path => $cdatabase);

    
    # if $refgene ends with a letter, then it is a specific isoform to
    # be tagged, else it is a sequence name and we want to find all
    # CDS/Pseudogene/Transcript objects hanging off this Gene and tag
    # them.

    my $output = $session_file.$refgene;
      
    if ($refgene =~ /[a-z]$/) {

      my $obj = $cdb->fetch(CDS       => "$refgene");
      $obj ||= $cdb->fetch(Transcript => "$refgene");
      $obj ||= $cdb->fetch(Pseudogene => "$refgene");
      return &error_warning("Invalid ID","$refgene is not a valid model name") unless $obj;
      $class = $obj->class;
      my $method = $obj->Method->name;
      
      if ($class eq 'CDS' && $method ne "curated") {return &error_warning("Only curated CDS models can have the Last_reviewed added in this way.")}
      
      # tag the objects
      open (LR,">$output") or die "cant open $output\n";
      print LR "\n$class : \"$refgene\"\n";
      print LR "Last_reviewed now $person\n\n";
      close LR;

    } else {

      open (LR,">$output") or die "cant open $output\n";

      my @cdses = $cdb->fetch(-query => "find CDS ${refgene}* WHERE Method = \"curated\"");
      foreach my $cds (@cdses){
	if ($cds eq $refgene || $cds =~ /${refgene}[a-z]$/) {
	  print LR "\nCDS : \"$cds\"\n";
	  print LR "Last_reviewed now $person\n\n";
	}
      }
    
      my @pseuds = $cdb->fetch(-query=>"find Pseudogene ${refgene}*");
      foreach my $pseud (@pseuds){
	if ($pseud eq $refgene || $pseud =~ /${refgene}[a-z]$/) {
	  print LR "\nPseudogene : \"$pseud\"\n";
	  print LR "Last_reviewed now $person\n\n";	
	}
      }
      
      my @trans = $cdb->fetch(-query=>"find Transcript ${refgene}*");
      foreach my $trans (@trans){
	if ($trans eq $refgene || $trans =~ /${refgene}[a-z]$/) {
	  print LR "\nTranscript : \"$trans\"\n";
	  print LR "Last_reviewed now $person\n\n";
	}
      }
      close LR;
    }

    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Success", "Added Last_reviewed tag to $form_gene_last_rev");
      &clear_last_rev;
    }
  }

################

sub add_evidence 
  {
    my $CLASS;
    my $refgene = $form_gene3;
    return unless $refgene;
    my $obj = $cdb->fetch(CDS       => "$refgene");
    $obj ||= $cdb->fetch(Transcript => "$refgene");
    $obj ||= $cdb->fetch(Pseudogene => "$refgene");
    return &error_warning("Invalid ID","$refgene is not a valid model name") unless $obj;
    if (defined $obj) {$CLASS = "${\$obj->class}"}
    my $method = $obj->Method->name;

    #if (!($method eq "curated")) {
    #  &error_warning("Only curated gene models can have top level evidence added in this way.");
    #  next;
    #}
    
    my $output = $session_file.$refgene;
    
    open (EVI,">$output") or die "cant open $output\n";
    # add the Evidence to the CDS


    print EVI "\n$CLASS : \"$refgene\"\n";
    if (defined $person){ print EVI "Evidence Curator_confirmed $person\n";}
    else {print EVI "Evidence\nRemark \"[No Curator specified]";}
    
    close EVI;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } 
    else {
      &confirm_message("Success", "Added evidence to $form_gene3");
      &clear_evi;
    }
  }

################

sub add_rem 
  {
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    my $CLASS;
    my $refgene = $form_gene_rem;
    return unless $refgene;
    my $obj = $cdb->fetch(CDS       => "$refgene");
    $obj ||= $cdb->fetch(Transcript => "$refgene");
    $obj ||= $cdb->fetch(Pseudogene => "$refgene");
    return &error_warning("Invalid ID","$refgene is not a valid model name") unless $obj;
    if (defined $obj) {$CLASS = "${\$obj->class}"}
    my $method = $obj->Method->name;

    my $output = $session_file.$refgene;
    
    open (REM,">$output") or die "cant open $output\n";
    # add the Evidence to the CDS


    print REM "\n$CLASS : \"$refgene\"\n";
    if (defined $person){ print REM "Remark \"[$date $user]\" Curator_confirmed $person\n";}
    
    close REM;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } 
    else {
      &confirm_message("Success", "Added Remark stub to $form_gene_rem");
      &clear_rem;
    }
  }

###################
sub add_ncrna_data
  {
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    my $refgene = $form_rna;
    return unless $refgene;
    &generate_message("CHECK", "Please save your database before clicking OK\n");
    $cdb->close();
    $cdb = Ace->connect(-path => $cdatabase);
    my $obj = $cdb->fetch(Transcript => "$refgene");
    return &error_warning("Invalid Transcript","$refgene is not a valid Transcript name") unless ($obj);
    my $method = $obj->Method->name;
    my $output = $session_file."ncRNA".$refgene;

    open (RNA,">$output") or die "cant open $output\n";
    # add the Evidence to the CDS
    print RNA "\nTranscript : \"$refgene\"\nTranscript ncRNA\nBrief_identification \"Possible non-coding RNA gene\"\nDB_remark \"C. elegans probable non-coding RNA\"\nMethod ncRNA\nRemark \"[$date $user] Created as this locus is potentially a ncRNA (with little evidence of translation) as it displays evidence of transcription and/or splicing\" Curator_confirmed $person\n\n";
    close RNA;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } 
    else {
      &confirm_message("Success", "Added Tags to $form_rna");
      &clear_data;
    }
  }

###################

###################
sub add_pseudo_data
  {
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    my $refgene = $form_pseudo;
    return unless $refgene;
    &generate_message("CHECK", "Please save your database before clicking OK\n");
    $cdb->close();
    $cdb = Ace->connect(-path => $cdatabase);
    my $obj = $cdb->fetch(Pseudogene => "$refgene");
    return &error_warning("Invalid Pseudogene","$refgene is not a valid Pseudogene name") unless ($obj);
    my $method = $obj->Method->name;
    my $output = $session_file."Pseudogene".$refgene;

    open (PSG,">$output") or die "cant open $output\n";
    # add the Evidence to the CDS
    print PSG "\nPseudogene : \"$refgene\"\nCoding_pseudogene\nDB_remark \"C. elegans predicted pseudogene\"\nMethod Pseudogene\nRemark \"[$date $user] This locus is defined as a Pseudogene as there are multiple gene modelling anomalies preventing a well structured coding sequence being annotated.\" Curator_confirmed $person\n\n";
    close PSG;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } 
    else {
      &confirm_message("Success", "Added Tags to $form_pseudo");
      &clear_pseudo_data;
    }
  }

###################



sub add_nonrna_data
  {
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
    my $refgene = $form_nonRNA;
    return unless $refgene;
    &generate_message("CHECK", "Please save your database before clicking OK\n");
    $cdb->close();
    $cdb = Ace->connect(-path => $cdatabase);
    my $obj = $cdb->fetch(Transcript => "$refgene");
    return &error_warning("Invalid Transcript","$refgene is not a valid Transcript name") unless ($obj);
    my $method = $obj->Method->name;
    my $output = $session_file."non_coding_RNA".$refgene;

    open (RNA,">$output") or die "cant open $output\n";
    # add the Evidence to the CDS
    print RNA "\nTranscript : \"$refgene\"\nTranscript ncRNA\nBrief_identification \"non-coding Transcript Isoform\"\nDB_remark \"C. elegans probable non-coding RNA\"\nMethod non_coding_transcript\nIsoform Curator_confirmed $person\nRemark \"[$date $user] Created this non_coding_transcript Isoform based on data being available that cannot be incorporated into a coding Isoform.\" Curator_confirmed $person\n\n";
    close RNA;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } 
    else {
      &confirm_message("Success", "Added Tags to $form_nonRNA");
      &clear_nonRNAdata;
    }
  }



################
  sub clean_gene
    {

      #get date for remark
      my ($day, $mon, $yr)  = (localtime)[3,4,5];
      my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
      my $refgene=$form_gene4;
      return unless $refgene;
      &generate_message("CHECK", "Have you saved your database??\n !!Failure will result in all CDS loss!!");
      $cdb->close();
      $cdb = Ace->connect(-path => $cdatabase);

      my @cdses = $cdb->fetch(-query=>"find gene $refgene ; follow Corresponding_CDS ; ! Evidence");
      return &error_warning("Invalid Gene","$refgene is not a valid Gene name") unless @cdses;
      
      # ace output for loading
      my $output = $session_file.$refgene;
      open (CLN,">$output") or die "cant open $output\n";
      
      foreach my $cds(@cdses){

	if (!defined $person){ 
	  print CLN "CDS : $cds\nEvidence\nMethod history\nGene_history $refgene\nRemark \"[$date $user] Removed automatically due to lack of evidence - based on lack of curation tags\" Curator_confirmed $person\n\nCDS : $cds\n-D Gene\n\n-R CDS $cds $cds:${wormpep_prefix}$version\n\n";
	}
	else {
	  print CLN "CDS : $cds\nEvidence Curator_confirmed $person\nMethod history\nGene_history $refgene\nRemark \"[$date $user] Removed automatically due to lack of evidence - based on lack of curation tags\" Curator_confirmed $person\n\nCDS : $cds\n-D Gene\n\n-R CDS $cds $cds:${wormpep_prefix}$version\n\n"
	}
      }
      
      close CLN;
      #load the data back to the database
      my $return_status = system("xremote -remote 'parse $output'");
      if ( ( $return_status >> 8 ) != 0 ) {
        &error_warning("WARNING", "X11 connection appears to be lost");
      } 
      else {
        &confirm_message("Success", "Cleaned CDSs from $form_gene4");
        &clear_gen;
      }
    }

##### end clear_gen






###############################################################

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
    if (! ( ($method eq "Genefinder") || ($method eq "twinscan") || ($method eq "RNASEQ.Hillier") || ($method eq "RNASEQ.Hillier.Aggregate") || ($method eq "jigsaw") || ($method =~ /cufflinks/) || ($method =~ "curated") || ($method eq "mGene")  || ($method eq "history") ) ) {

      &error_warning("Wrong method","Only Selected Predictions can be blessed");
      next;
    } elsif ($exceptions eq "Link") {
      &error_warning("Warning","This Prediction lies over clone boundaries.");
      next;
    }

    my $new_gene;
    if ($form_gene2) {
      $new_gene = $form_gene2;
    }
    else {
      $new_gene = &suggest_name("$stem");
      unless ( $new_gene ){
	&error_warning("No name","Can't suggest name for gene based on $gene");
	return;
      }
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

    #remove the exons of the target gene if it already exists.
    if (defined$form_gene2) {
      print BLS "\nCDS : \"$new_gene\"\n";
      print BLS "-D Source_exons\n\n";
    }

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
    if ($species =~ "elegans") {
      print BLS "Remark \"[$date $user] Autoconversion from $gene\" From_analysis RNASeq_Hillier_elegans\n";
    }
    else {
      print BLS "Remark \"[$date $user] Autoconversion from $gene\"";
    }
    if (defined $person){ print BLS " Remark \"[$date $user] Autoconversion from $gene\" Curator_confirmed $person\nEvidence Curator_confirmed $person\n";}
    else {print BLS "\nEvidence\n";}
    
    close BLS;
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Created", "Converted $gene into $new_gene");

      #pop-up window if a -mail option has not been set.
      &error_warning("Info","Remember to include WBGene ID info for $new_gene") if (!($mail));

      #clear the text boxes in the gui.
      &clear_gene;
      &clear_proposed;

      #&suggest_name("$stem");

      &mail_geneace("$new_gene", "$person", $date) if ($mail);
    }
  }

##########################################
# Clone a gene from the curation database - majority of code borrowed from history sub.
##########################################


sub clone_gene {
    my $cds_identifier;
    my $cds = $clone_form;
    my $suggested = $clone_form_sug;
    return unless $cds;
    my $output = $session_file."cloned".$cds;
    open (CLO,">$output") or die "cant open $output\n";
    last if( $cds eq "end" );
    my $TYPE;

    my $cloneobj = $db->fetch(CDS      => $cds);
    $cloneobj ||= $db->fetch(Transcript=> $cds);
    $cloneobj ||= $db->fetch(Pseudogene=> $cds);

    if (defined $cloneobj) {$TYPE = "${\$cloneobj->class}"}
    return &error_warning("Invalid Identifier","$cds is not a valid CDS/Transcript/Pseudogene name or does not exist in the linked database") unless $cloneobj;
    my $species = $cloneobj->Species->name;
    my $gene = $cloneobj->Gene->name;
    my $lab = $cloneobj->From_laboratory->name;
    my $seq = $cloneobj->Sequence->name;

    # capture parent clone info
    my $clone = $cloneobj->Sequence;
    my @clone_CDSs;
    if ($TYPE eq "CDS") {
      @clone_CDSs = $clone->CDS_child;
    }
    else {
      @clone_CDSs = $clone->$TYPE;
    }
    my $start;
    my $end;
    foreach my $CDS ( @clone_CDSs ) {
      next unless ($CDS->name eq "$cds");
      $start = $CDS->right->name;
      $end = $CDS->right->right->name;
      last;
    }

    if (defined $suggested){
      $cds_identifier = $suggested;
    }
    else {
      $cds_identifier = $cds.":cloned";
    }
    # cloning a Transcript_object - Y73B6BL.269b
 
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
      
    #Print Output file
    print CLO "Sequence : $seq\n";
    if ($TYPE eq 'CDS') {
      print CLO "CDS_child \"$cds_identifier\" $start $end\n";
    }
    else {
      print CLO "$TYPE \"$cds_identifier\" $start $end\n";
    }
    print CLO "\n$TYPE : $cds_identifier\n";
    if ($TYPE eq "Transcript"){
      print CLO "Brief_identification \"$species ncRNA gene\"\n";
    }
    foreach ($cloneobj->Source_exons) {
      my ($start,$end) = $_->row(0);
      print CLO "Source_exons ",$start->name," ",$end->name,"\n";
    }
    print CLO "Sequence $seq\n";
    print CLO "From_laboratory $lab\n";
    print CLO "Isoform Curator_confirmed $person\n" if $person;
    print CLO "Gene $gene\n" if $gene;
    print CLO "Species \"$species\"\n" if $species;
    if (defined $person){ 
      print CLO "Evidence Curator_confirmed $person\n" if $person;
      print CLO "Remark \"[$date $user] Cloned from the template $TYPE $clone_form\"Curator_confirmed $person\n";
    }
    else {
      print CLO "Evidence";
      print CLO "Remark \"[$date] Cloned from the template $TYPE $clone_form\"\n";
    }
    if ($TYPE eq "CDS") {
      print CLO "Method curated\n";
    }
    elsif ($TYPE eq "Transcript") {
      print CLO "Method ncRNA\n";
    }
    else {
      print CLO "Method $TYPE\n";
    }
    close CLO;
    
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Cloned Gene","Cloned $cds_identifier has been created");
      &clear_cloned;
      &clear_cloned_prop;
    }
  }


########## END Clone gene#################

##############################################################
# clear various fields
##############################################################

sub clear
  {
    $cds_val->delete(0,'end');
  }

sub clear_gene
  {
    $gene_val->delete(0,'end');
  }

sub clear_proposed
  {
    $proposed_name->delete(0,'end');
  }
sub clear_evi
  {
    $cdswork->delete(0,'end');
  }
sub clear_last_rev
  {
    $lastrev->delete(0,'end');
  }
sub clear_rem
  {
    $remwork->delete(0,'end');
  }
sub clear_data
  {
    $rnawork->delete(0,'end');
  }
sub clear_pseudo_data
  {
    $pseudowork->delete(0,'end');
  }
sub clear_nonRNAdata
  {
    $nonrnawork->delete(0,'end');
  }
sub clear_gen
  {
    $genework->delete(0,'end');
  }
sub clear_cloned
  {
    $clonework->delete(0,'end');
  }

sub clear_cloned_prop
  {
    $cloneprop->delete(0,'end');
  }

##############################################################
sub suggest_name
  {
    my $stem = shift;
    my $query = "find worm_genes $stem.*";
    my @genes = $cdb->fetch( -query => "$query");
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

    my $biotype;
    my $clone_tag; # where it lives under the Sequence object tags
    my $new_method; # the new method to give it
    my $transcript_type;
    my $transcript_type_value1;
    my $transcript_type_value2;
    my $pseudogene_type;
    my ($day, $mon, $yr)  = (localtime)[3,4,5];
    my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);

    #print "enter CDS to make history object \n";
    my $cds = $form_cds;
    return unless $cds;
    my $output = $session_file.$cds;

    last if( $cds eq "end" );

    #$cds = &confirm_case($cds);

    my $obj = $db->fetch(CDS => "$cds");
    if ($obj) {
      $biotype = 'CDS';
      $clone_tag = 'CDS_child';
      $new_method = 'history';
      if ($obj->Method->name ne "curated" ) {
        &error_warning("NOT CURATED", "I only do curated CDS histories!\n");
        return;
      }
    } elsif ($obj = $db->fetch(Pseudogene => "$cds")) {
      $biotype = 'Pseudogene';
      $clone_tag = 'Pseudogene';
      $new_method = 'history_pseudogene';
      $pseudogene_type = $obj->Type->name;
    } elsif ($obj = $db->fetch(Transcript => "$cds")) {
      $biotype = 'Transcript';
      $clone_tag = 'Transcript';
      $new_method = 'history_transcript';
      if ($obj->Method->name eq "Coding_transcript" ) {
        &error_warning("WARNING", "I don't do Coding_transcript Transcript histories!\n");
        return;
      }
      ($transcript_type, $transcript_type_value1, $transcript_type_value2) = $obj->Transcript->row(0);
    }

    if (!defined $obj) {
      &error_warning("Invalid CDS", "$cds is not a valid CDS/Pseudogene/Transcript name");
      return;
    }

    if ($db->fetch($biotype => "$cds:${wormpep_prefix}$version") ) {
      &error_warning("WARNING", "This history object already exists\n");
      return;	
    }
    
    print "Found object $cds with biotype $biotype to make history.\n";

    my $species = $obj->Species->name;
    my $gene = $obj->Gene->name;
    my $lab = $obj->From_laboratory->name;
    my $seq = $obj->Sequence->name;
    my $brief = $obj->Brief_identification;
    my $brief_identification = $brief->name if ($brief);

    # parent clone coords
    my $clone = $obj->Sequence;
    my @clone_entry;
    if ($biotype eq 'CDS') {
      @clone_entry = $clone->CDS_child;
    } elsif ($biotype eq 'Pseudogene') {
      @clone_entry = $clone->Pseudogene;
    } elsif ($biotype eq 'Transcript') {
      @clone_entry = $clone->Transcript;
    }
    my $start;
    my $end;
    foreach my $CDS ( @clone_entry ) {
      next unless ($CDS->name eq "$cds");
      $start = $CDS->right->name;
      $end = $CDS->right->right->name;
      last;
    }
    print "Found $cds in $seq start $start end $end\n";

    #print ace format
    open (HIS,">$output") or die "cant open $output\n";

    print HIS "Sequence : $seq\n";
    print HIS "$clone_tag \"$cds:${wormpep_prefix}$version\" $start $end\n";

    print HIS "\n$biotype : $cds:${wormpep_prefix}$version\n";
    print "written $biotype : $cds\n";

    foreach ($obj->Source_exons) {
      my ($start,$end) = $_->row(0);
      print HIS "Source_exons ",$start->name," ",$end->name,"\n";
    }
    print "written Source_exons\n";

    foreach ($obj->Remark) {
      my ($remark, $evidence, $evidence_value1, $evidence_value2) = $_->row(0);
      print HIS "Remark \"", $remark->name ,"\"";
      print HIS " \"", $evidence->name,"\"" if ($evidence);
      print HIS " \"", $evidence_value1->name,"\"" if ($evidence_value1);
      print HIS " \"", $evidence_value2->name,"\"" if ($evidence_value2);
      print HIS "\n";
    }
    print "written Remark\n";

    print HIS "Sequence $seq\n";
    print HIS "From_laboratory $lab\n";
    print HIS "Gene_history $gene\n" if $gene;
    print HIS "Species \"$species\"\n" if $species;
    print HIS "Evidence Curator_confirmed $person\n" if $person;
    print HIS "Evidence\n" if (!defined $person);
    print HIS "Remark \"[$date $user] Automatically generated history object\"\n";
    print HIS "Method $new_method\n";
    print HIS "Brief_identification \"$brief_identification\"\n" if ($brief_identification);
    print "written Sequence, From_laboratory, Gene_history etc.\n";

    print HIS "CDS\n" if ($biotype eq 'CDS');
    print HIS "CDS_predicted_by ${\$obj->CDS_predicted_by}\n" if $obj->CDS_predicted_by;    
 
    print HIS "$pseudogene_type" if ($pseudogene_type);

    print HIS "$transcript_type" if ($transcript_type);
    print HIS " \"",$transcript_type_value1->name,"\"" if ($transcript_type_value1);
    print HIS " \"",$transcript_type_value2->name,"\"" if ($transcript_type_value2);
    print HIS "\n" if ($transcript_type);
    print "finished writing, about to close file\n";

    close HIS;
    print "file $output closed. Size: ", (-s $output), "\n";
    system("cat $output"); # debug statement to show contents of file
    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    } else {
      &confirm_message("Made history","History $cds:${wormpep_prefix}$version has been created");
      &clear;
    }
  }




#############################################################################
# connect to the mysql database and return handle
# Args: lab = "HX" or "RW"
sub connect_to_database {
  my ($lab) = @_;

  my $dbsn;
  my $dbuser;
  my $dbpass;

  my $sqldb = "worm_anomaly";

  # if the species is anything other than 'elegans' (the default) then
  # we use the database that has a name formed from the 'base' database
  # name and that species' name

  if ($species && $species ne 'elegans') {$sqldb = $sqldb . "_" . lc $species;}

  if ($lab eq "HX") {

    $dbsn = "DBI:mysql:database=$sqldb;host=farmdb1";
    $dbuser = "wormadmin";
    $dbpass = "worms";

  } else {			# for St. Louis

    $dbsn = "DBI:mysql:database=$sqldb;host=XXXXX";
    $dbuser = "XXXXXX";
    $dbpass = "XXXXXX";

  }

  my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {AutoCommit => 1, RaiseError => 1})
      || die "cannot connect to database $sqldb, $DBI::errstr";

  return $mysql;

}

#############################################################################
#  &populate_zero_weight_list($mysql, $lab, $chromosome, \$zero_weight_list);
# get the list of anomaly types that we can later use to reweight the types
# sets up the VIEW table with initial weights

sub populate_zero_weight_list {

  my ($mysql, $lab, $chromosome_list, $zero_weight_list_ref) = @_;

  # get the available anomaly types
  my $chromosome_in = '';
  if (defined $chromosome_list) {
    $chromosome_in = "chromosome IN ( $chromosome_list ) AND ";
  }
  my $query = qq{ SELECT type, COUNT(*) FROM anomaly WHERE $chromosome_in centre = "$lab" AND active = 1 GROUP BY type; };
  my $db_query = $mysql->prepare ( $query );
  $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
  $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
  my ($mysql, $lab, $chromosome_list, $anomaly_list_ref) = @_;

  my $chromosome_in = '';
  if (defined $chromosome_list) {
    $chromosome_in = "a.chromosome IN ( $chromosome_list ) AND ";
  }

  # name of the temporary 'view' weighting table
  my $view = "weight_$user";

  # get the highest scoring 10 kb windows of one sense or the other sorted by score descending
  # group by window, sense order by 2 desc - means make the
  # SUM(thing_score) sum up rows within distinct window and distinct
  # sense and sort the output by descending SUM(thing_score)

  my $query;
  if ($display_by_clones) { # St. Louis people like to have the anomalies from a clone lumped together

    $query = qq{ SELECT a.window, SUM(a.thing_score * w.weight), a.sense, a.clone, a.chromosome FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type     WHERE $chromosome_in a.centre = "$lab" AND a.active = 1 GROUP BY clone, sense ORDER BY 2 DESC; };


  } else { # normal display by 10Kb window

# there is no difference in speed between these two variants of this command:
#
#  $query = qq{ SELECT a.window, SUM(a.thing_score), a.sense, a.clone FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type    WHERE a.chromosome = "$chromosome" AND a.centre = "$lab" AND a.active = 1 AND w.weight = 1 GROUP BY window, sense ORDER BY 2 DESC; };
#
    $query = qq{ SELECT a.window, SUM(a.thing_score * w.weight), a.sense, a.clone, a.chromosome FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type     WHERE $chromosome_in a.centre = "$lab" AND a.active = 1 GROUP BY window, sense ORDER BY 2 DESC; };
  }

  #print "\n";
  #print "chromosome=$chromosome_list\n";
  #print "lab=$lab\n";
  #print "query=$query\n";
  #print "\n";

  my $db_query = $mysql->prepare ( $query );

  $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");

  $results = $db_query->fetchall_arrayref;

  my $over_10 = 0;
  my $over_5 = 0;
  my $over_2 = 0;
  my $over_1 = 0;
  my $over_half = 0;
  my $over_quarter = 0;
  my $under_quarter = 0;
  foreach my $result_row (@$results) {
    if ($display_by_clones) { # St. Louis people like to have the anomalies from a clone lumped together
      $$anomaly_list_ref->insert('end', $result_row->[4] . " " . $result_row->[3] . " Sense: " . $result_row->[2] . " Score: " . $result_row->[1] );      
    } else { # normal display by 10Kb window
      # clone, sense, score, window
      $$anomaly_list_ref->insert('end', $result_row->[4] . " " . $result_row->[3] . " Sense: " . $result_row->[2] . " Score: " . $result_row->[1] . " ID: " . $result_row->[0] );
    }

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
  my ($mysql_ref, $lab_ref, $chromosome_list_ref, $anomaly_list_ref, $anomaly_detail_list_ref, $zero_weight_ref, $check_ref) = @_;


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
  &populate_anomaly_window_list($$mysql_ref, $$lab_ref, $$chromosome_list_ref, $anomaly_list_ref);



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
      $sth->execute($anomaly_id) or &error_warning("WARNING", "MySQL server appears to be down");
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

  } elsif (defined $clone && $clone ne "" ) {

    # display the whole clone
    &goto_location($clone, 1, 200000, '+', 0);

    # get and display the individual anomalies found in this clone
    # pull out all anomalies in this clone except those marked as active = 0 and those with zero-weighted anomaly types
    $query = qq{ SELECT a.type, a.clone, a.clone_start, a.clone_end, a.chromosome_start, a.chromosome_end, a.sense, a.thing_id, a.thing_score, a.explanation, a.anomaly_id FROM anomaly AS a INNER JOIN $view AS w ON a.type = w.type   WHERE a.clone = "$clone" AND a.active = 1 AND w.weight = 1 ORDER BY chromosome_start; };


  } elsif ($display_by_clones && defined $selection) {  # St. Louis people like to have the anomalies from a clone lumped together


    #print "get current selection\n";
    my $selected_value = $$anomaly_window_list->get( $selection );
    my @selected = split(/\s+/, $selected_value);

    # use the selected clone
    my $clone = $selected[1];

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
    my $window_start = $selected[7] * $WINDOW_SIZE + 1;
    my $window_end = ($selected[7] * $WINDOW_SIZE) + $WINDOW_SIZE - 1;
    my $sense = $selected[3];

    # camace (and stlace) doesn't have the CHROMOSOME_* Sequence objects in them
    # so we now have to convert to clone coords
    if ($species eq 'elegans') {
      $chromosome = "CHROMOSOME_$chromosome";
    }

    print "chromosome= >${chromosome}< $window_start, $window_end";
    my @clone_coords = ${$coords_ref}->LocateSpan($chromosome, $window_start, $window_end );
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

  $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
	$sth->execute($anomaly_id) or &error_warning("WARNING", "MySQL server appears to be down");
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
  $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
      $db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
    #use Tk::Graph;

    # set up the window to hold the graphs
    my $D = MainWindow->new;
    $D->optionAdd('*BorderWidth' => 1); # make the style better with a thinner border

    my $frame1 = $D->Frame(-background =>'whitesmoke')->pack();
    my $close_button = $frame1->Button(-text => "Close progress window",
				       -background => "IndianRed4", #was bisque
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
	$db_query->execute() or &error_warning("WARNING", "MySQL server appears to be down");
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
	  '  > 2.5'  => $over_2[-1],
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
					 -title    => "chromosome $chr $lab",
					 -linewidth => 3,
					 -headroom => 5,
					 #-legend => 0,
					 -look => 100, # display the last 100 days (maybe try making this 0 or undef to get all of the data)
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
    $D->configure(-background => 'DarkSeaGreen', # was cyan PaleGreen
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
    $D->configure(-background => 'IndianRed4',
		  -foreground => 'black');

    my $choice = $D->Show;
  }


############################################################################

  sub generate_message
    {
      my $title = shift;
      my $text = shift;
      my $D = $main_gui->Dialog(
                                -title => "$title",
                                -text  => "$text",
                                -default_button => 'ok',
                                -buttons        => ['ok']
                               );
      $D->configure(-background => 'Red4',
                    -foreground => 'black');
      
      my $choice = $D->Show;
    }

############################################################################
# gets the version of the currentDB database
sub get_history_version
  {
    my $rdatabase = shift;
    my $WS_version = `grep "NAME WS" $rdatabase/wspec/database.wrm`;
    chomp($WS_version);
    $WS_version =~ s/.*WS//;
    die "no version\n" unless $WS_version =~ /\d+/;

    # check to see if the Build has start - if so then use the Build's version
    my $build_file = $wormbase->autoace . "/wspec/database.wrm";
    if (-e $build_file) {
      my $Build_WS_version = `grep "NAME WS" $build_file`;
      chomp($Build_WS_version);
      $Build_WS_version =~ s/.*WS//;
      if (defined $Build_WS_version ) {$WS_version = $Build_WS_version }
    }
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
#  xremote -remote 'gif seqget CHROMOSOME_II -coords 111000 112000; seqdisplay'
#
#


__END__

=pod

=head1 NAME - history_maker.pl

A perl Tk interface to aid manual gene curation.

=over 1

=back

=head2 Options

-chromosome  : The chromosome that the curator wants to work on in the anomalies database.

-source      : The static database connection to use as source for history and ab initio predictions, clone widget etc.

-curationdb  : The active curation database that the curator is working on (cleangene option only at present).

-user        : If you are not yourself enter your user to use in autgenerated comments (when blessing genefinder etc)

-design      : Does not make Aceperl connection - dev tool for quicker startup

-mail        : This option emails blessed predictions to mt3 requesting a new gene ID

-addevidence : This option allows the curator to automatically add the top level evidence to a CDS object

-addlastreviewed : This option allows the curator to automatically add the "Last_reviewed" tag to a CDS object. It takes the sequence name to add the tag to all isoforms, or a specific isoform name to add the tag only to that isoform.

-cleangene   : This option removes all unwanted Isoforms from a given loci. Isoforms are preserved by the presence of the top level Evidence tag.

-brugia       : This automatically selects options that are useful for brugia curation and stores a version number. Version number hard sets the CDS history prefix to be :bm<var>.

-clone       : This option loads a new widget that allows a CDS model to be cloned back into the current FMAP from the reference database.

-anomaly     : This option loads a widget that pulls anomalies from the mysql anomaly database for a given -chromosome(s)

-blast       : This option takes a file of blast alignment data and allows the curator to navigate through the list easily.

-blesser     : This option loads the blesser widget that allows a curator to convert a gene model into a curated gene structure with an optional name.

-lab         : This option is used for correctly assigning the lab in C. elegans curation.

=over 1


=back

=head2 Usage Examples

  C. elegans curation

      history_maker.pl -user pad -source /nfs/wormpub/camace_orig -curationdb /nfs/wormpub/camace_pad -blesser -clone -addevidence -addlastreviewed -anomaly -clean -chromosome I,II,III,X


  Brugia curation

      history_maker.pl -user pad -lab HX -curationdb brugia_curation_v9 -source brugia_refv9 -brugia 9 -cleangene -clone

=over 1

=back

=head2 Widget Descriptions

=over 1

=item Intron finder

Uses the GFF check files in development_release to provide a list of introns.  Select and click to go to intron in FMAP

=item History Maker

presents to the user a simple box with a space to enter a CDS name and a button to make a history object.  

=item Add Evidence

This option will populate the top evidence hash with the Curator_confirmed Evidence. This is perticularly imoprtant and useful for Tier II aut removal of Isoforms via the Clean Gene option.

=item Clean Gene

This option requires that the -database is the curation database that the user is working on. It takes a Gene ID and finds any CDS_child that does not have a top level Ecvidence hash populated. These are then auto converted to history objects. A popup window exists that reminds you to save your database so that not all Isoforms are lost from a loci.

=item Clone CDS

This option/widget takes a CDS name as well as an optional target name. The function is to clone the reference state of the CDS and embed it alongside the current version. Useful if creating a new Isoform at a locus or attempting to roll back an annotation.

=item Prediction Blesser

Enter current CDS name eg AC8.gc3 and click "Bless this gene".  This will create a new CDS with the correct name based on the "worm_genes" class. If you have specified the -mail option a new gene ID will automatically be requested.

=item BLAST hit finder

Using passed file, the user is presented with a list of genomic locations where there are strong blastp hits but no current gene.
Click to go to that location in FMAP.

The input to the acedb database is done immediately using 'xremote' and is visible after recalculating the fmap.  The target database is determined by the operating system as the last one to have been opened using the terminal that this script is run from.

To ensure that the correct database is being used a shell script should be used to launch both simultaneously.

A reference database is used to extract the relevant info needed to make a history ie source exons, gene_id etc.and this is set with the -source option

Some error checking is done so that;

=over 4
  history objects cant be created for non-existant CDSs.
  history objects with the same name as existing histories will not be made.
  input case is irrelevant - all histories will be converted to uppercase clone names.

=back

=cut

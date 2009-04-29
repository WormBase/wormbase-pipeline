#!/software/bin/perl -w
#
# merge_split_camaces.pl
# 
# A script to make multiple copies of camace for curation, and merge them back again
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2009-04-29 10:58:04 $
#
# Persisting errors.
#running csh -c "reformat_acediff file 1 file2"
#Use of uninitialized value in concatenation (.) or string at /reformat_acediff line 121.
#====================

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

my $all;                   # All
my $pad;                   # Use Paul's split
my $gw3;                   # Use Gary's split
my $ar2;                   # Use Anthony's split
my $merge;                 # Merging databases
my $split;                 # Splitting databases
my $update;                # Update current database
my $debug;                 # Debug option
my $help;                  # Help menu
my $WS_version;            # Removes final wormsrv2 dependancy.
my $store;                 # Storable not needed as this is not a build script!
my $test;
my $wormbase;
my $extra;                 # remove the GeneIDupdater call as this is run outside this script.
my $email;                 # Option for child scripts that can tace a user email option.
my $nodump;                # don't dump from split camaces.
my $nosplit;               # use this option with the split if you only remove the mass_spec and tiling array data.
my $remove_only;           # just removes curation data from $database specified.
my $sdatabase;
my $load_only;

  GetOptions (
	      "all"        => \$all,
	      "pad"        => \$pad,
	      "gw3"        => \$gw3,
	      "ar2"        => \$ar2,
	      "merge"      => \$merge,
	      "split"      => \$split,
	      "update"     => \$update,
	      "help"       => \$help,
	      "debug:s"    => \$debug,
	      "version:s"  => \$WS_version,
	      "store"      => \$store,
	      "test"       => \$test,
	      "extra"      => \$extra,
	      "email:s"    => \$email,
	      "nodump"     => \$nodump,
	      "nosplit"    => \$nosplit,
	      "database:s" => \$sdatabase,
	      "remove_only" => \$remove_only,
	      "load_only" => \$load_only,
	     );



# Help pod if needed
&usage("Help") if ($help);

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug,
			     -test => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace;

if (! defined$WS_version) {
  $WS_version = $wormbase->get_wormbase_version;
}

my $next_build = $WS_version + 1;

$log->write_to ("WS_version : $WS_version\tWS_next : $next_build\n") if ($debug);

# load @databases array with user database names.
my @databases; #array to store what splits are to be merged.
push(@databases,"orig");
push(@databases,"pad") if ($pad || $all);
push(@databases,"gw3") if ($gw3 || $all);

# directory paths
my $wormpub = $wormbase->wormpub;
our $canonical = $wormbase->database('camace');
our $orig = $wormpub."/camace_orig";
our $directory   = $orig."/WS${WS_version}-WS${next_build}";
$log->write_to ("OUTPUT_DIR: ".$orig."/WS${WS_version}-WS${next_build}\n\n");

# what classes of data do we want to dump?
my @classes = ('Transposon', 'Transcript', 'CDS', 'Sequence', 'Feature', 'Feature_data', 'Pseudogene', 'dna');


#################################################
## (1) Merge split databases #1 - do the diffs ##
#################################################

if ($merge) {
  if (!-e $directory) {
    $log->write_to ("New directory : '$directory'\n\n") if ($debug);
    mkdir ($directory) or die "Failed to create ${directory}\n";
  }
  else {
    $log->write_to ("Using existing directory : '$directory'\n\n") if ($debug);
  }

  $log->write_to ("You are merging Data from " . (join '-',@databases) ."\n\n");

  # dumps the Pseudogene, Transcript, Feature and Sequence classes from the database
  &dump_camace unless ($nodump);

  shift (@databases); #remove 1st element of array as there aren't going to be any camace_orig updates
  # run acediff on the files tidy up and reformat the diff files ready to be loaded
  my $script_path = glob("~pad/wormbase/scripts");
  foreach my $database (@databases) {
    foreach my $class (@classes) {
      my $path_new = $directory . "/${class}_${database}.ace";
      my $path_ref = $directory . "/${class}_orig.ace";
      $wormbase->run_command("csh -c \"/software/worm/bin/acedb/acediff $path_ref $path_new >! $directory/${class}_diff_${database}.ace\"", $log) && die "acediff Failed for ${path_new}\n";
      $wormbase->run_script("reformat_acediff -file $directory/${class}_diff_${database}.ace -fileout $directory/update_${class}_${database}.ace", $log) && die "reformat Failed for ${class}_diff_${database}.ace\n";
    }
  }
  $log->write_to ("Phase 1 finished and all files can be found in $directory\n");
}

#################################################################
## (2) synchronises Canonical Database with the split versions ##
#################################################################

if ($update) {
  shift (@databases); 
  &update_camace;
  $log->write_to ("Phase 2 finished $canonical is now updated\n");
}

############################################################################
## (3) TransferDB calls to move Canonical Database to the split databases ##
############################################################################

if ($split) {
  $log->write_to ("Removing old split databases and Copying $canonical database to the split camaces\n");
  &split_databases unless ($nosplit);
  
  # Remove data we don't want in the build that may have got into the canonical erroneously.
  $log->write_to ("Removing any curation data that may have erroneously been loaded into the canonical db.\n");
  print "Removing any curation data that may have erroneously been loaded into the canonical db.\n" if ($debug);
  &remove_data("camace",);
  
  # Load DB_remarks into camace as these are required for EMBL dumping
  $log->write_to ("loading DB_remarks into camace\n");
  print "loading DB_remarks into camace\n" if ($debug);
  &load_curation_data("camace"); #this loads just the DB_remarks into camace for EMBL dumping.
  

  ## Check Canonical Database for errors one last time. ##
  $log->write_to ("\nChecking camace for the final time.\n-----------------------------------\n");
  print "\nChecking camace for the final time.\n" if ($debug);
  $wormbase->run_script("NAMEDB/camace_nameDB_comm.pl", $log) && die "Failed to run camace_nameDB_comm.pl\n";
  $log->write_to ("camace_nameDB_comm.pl Finished, check the build log email for errors.\n\n");
  print "camace_nameDB_comm.pl Finished, check the build log email for errors.\n\n" if ($debug);
  $wormbase->run_script("camcheck.pl -database ". $wormbase->database('camace'), $log) && die "Failed to run camcheck.pl\n";
  $log->write_to ("camcheck.pl Finished, check the build log email for errors.\n");
  print "camcheck.pl Finished, check the build log email for errors.\n" if ($debug);

  $log->write_to ("\nPhase 3 finished. All ~wormpub split camaces can now be used\n\nCheck all TransferDB log files for \"ended SUCCESSFULLY\"\n");
  print "Phase 3 finished. All ~wormpub split camaces can now be used\n\nCheck all TransferDB log files for \"ended SUCCESSFULLY\"\n" if ($debug);
}


####################################################
# Additional functionality to remove curation data #
# or load curation data into a nominated database. #
####################################################

if ($remove_only) {
  if (!defined $sdatabase) {
    die "ERROR:You have not specified a database for data to be removed!!\n";
  }
  else {
    $log->write_to ("-remove_only specified!\nRemoving curation data from $sdatabase.\n");
    print "-remove_only specified!\nRemoving curation data from $sdatabase.\n";
    &remove_data($sdatabase);
  }
}
if ($load_only) {
  if (!defined $sdatabase) {
    die "ERROR:You have not specified a database for data to be loaded!!\n";
  }
  else {
    $log->write_to ("-load_only specified!\nLoading curation data from $sdatabase.\n");
    print "-load_only specified!\nLoading curation data from $sdatabase.\n";
    &load_curation_data($sdatabase);
  }
}

# end and tidy up.
print "Diaskeda same Poli\n"; #we had alot of fun#
$log->mail();
exit(0);



##################################################################################################################################################
                                                                #################
                                                                #  Subroutines  #
                                                                #################
##################################################################################################################################################

#####################################
#(1a) Dump files from camace splits #
#####################################
sub dump_camace {
  #dumps out subset of classes from camace splits and processes the files to be loaded back to Canonical Database
  #array of classes to be dumped
  my $camace_path;
  my $path;

  foreach my $database (@databases) {

    $camace_path = $wormpub."/camace_${database}";
    $ENV{'ACEDB'} = $camace_path;
    
    foreach my $class (@classes) {
      $log->write_to ("dumping $class class from camace_${database}\n");
      $path = "$directory/" . "${class}_${database}.ace";
      &dumpace("$class",$path);
      $log->write_to ("dumped $class class from camace_${database}\n\n");
    }
  }
}

############################
#(1b) data dumping routine #
############################
sub dumpace {
  my $class    = shift;
  my $filepath = shift;
  
  my $command = "nosave\nquery find $class\nshow -a -f $filepath\nquit\n";
  
  # dump out from ACEDB
  $log->write_to ("\nFilename: $filepath\n");
  open (TACE,"| $tace") or die "Failed to open database connection\n";
  print TACE $command;
  close TACE;
}

##########################
#(2a) Loading subrouting #
##########################
sub loadace {
  my $filepath = shift;
  my $tsuser   = shift;


#####################could do with user input prompt not die!!!!!



  if ($tsuser =~ /\S+\s/) {
    &testuser($tsuser);
  }
  my $command = "pparse $filepath\nsave\nquit\n";

  # dump out from ACEDB
  $log->write_to ("\nFilename: $filepath\n");
  open (TACE,"| $tace -tsuser $tsuser") or die "Failed to open database connection\n";
  print TACE $command;
  close TACE;
}

sub testuser {
  my $subtsuser = shift;
  print "ERROR:  $subtsuser contains white space!!\nDo you want to assign a new tsuser? (y or n)\n";
  my $answer=<STDIN>;
  if ($answer eq "y\n") {
    print "Please input a new tsuser containing alphanumeric characters only:";
    my $answer2=<STDIN>;
    &testuser if ($answer2 =~ /\S+\s/);
  }
  if ($answer eq "n\n") {
    die "Failed to open database connection because of tsuser :$subtsuser\n"
  }
}

########################################
#(2b) upload data to Canonical Database#
########################################
sub update_camace {
  # upload processed diff files into Canonical Database.
  $log->write_to ("Upload diff files to $canonical");
  $ENV{'ACEDB'} = $canonical;
  
  ##  Upload the split camace update files you have just created  ##
  foreach my $database (@databases) {
    foreach my $class (@classes) {
      &loadace("$directory/update_${class}_${database}.ace", '${database}_${class}');
    }
  }

  ## Update Protein ID's and check embl sequence_versions                         ##
  $log->write_to ("\n\nRunning Gene_ID_update.pl to refresh Protein_ID\'s and clone sequence version\'s\n\n");
  if ($debug) {$wormbase->run_script("GeneID_updater.pl -proteinID -sv -version $WS_version -debug $debug -update", $log) && die "Failed to run Gene_ID_updater.pl\n";}
  else {$wormbase->run_script("GeneID_updater.pl -proteinID -sv -version $WS_version -update", $log) && die "Failed to run Gene_ID_updater.pl\n";}
  $log->write_to ("Updated Protein_ID\'s and clone Sequence_versions\n");
  
  ## Check if camace is in sync with the name server. ##
  $log->write_to ("\n\nRunning camace_nameDB_comm.pl.\n\n");
  $wormbase->run_script("NAMEDB/camace_nameDB_comm.pl", $log) && die "Failed to run camace_nameDB_comm.pl\n";
  $log->write_to ("camace_nameDB_comm.pl Finished, check the build log email for errors.\n");


  ## Check Canonical Database for errors. ##
  if ($email) {
    $wormbase->run_script("camcheck.pl -database ". $wormbase->database('camace'). "-e $email", $log) && die "Failed to run camcheck.pl\n";
  }
  else {  
    $wormbase->run_script("camcheck.pl -database ". $wormbase->database('camace'), $log) && die "Failed to run camcheck.pl\n";
  }
}


####################
#(3a) Data dispersion#
####################

sub split_databases {

#reinitialise the camace_splits
  foreach my $database (@databases) {
    $log->write_to ("Destroying $database\n");
    print "Destroying $database\n" if ($debug);
    my $split_db = $wormpub."/camace_${database}";
    $wormbase->run_command("rm -rf $split_db/database/ACEDB.wrm", $log) && die "Failed to remove $split_db/database/ACEDB.wrm\n";
    
    my $command = "y\nquit\n";
    open (TACE,"| $tace -tsuser merge_split $split_db") or die "Failed to open database connection\n";
    print TACE $command;
    close TACE;
    $log->write_to ("Destroyed $database\n");
    print "Destroyed $database\n" if ($debug);
  }

  # Transfer camace -> camace_orig
  $log->write_to ("Transfering $canonical to $orig\n");
  print "Transfering $canonical to $orig\n" if ($debug);
  $wormbase->run_script("TransferDB.pl -start $canonical -end $orig -split -database -wspec", $log) && die "Failed to run TransferDB.pl for camace_orig\n";

  # work on camace_orig to get it populated with curation data
  $log->write_to ("Refreshing curation data in $orig\n-------------------------------------\nRemoving old curation data from camace_orig\n");
  print "Refreshing curation data in $orig\n-------------------------------------\nRemoving old curation data from camace_orig\n" if ($debug);
  &remove_data($orig); #remove any stale data. 
  $log->write_to ("loading curation data into camace_orig\n");
  print "loading curation data into camace_orig\n" if ($debug);
  &load_curation_data($orig); #this loads all curation data into camace_orig (also removes bogus gene predictions)
  $log->write_to ("loading DB_remarks into camace\n");
  # Create the rest of the curation databases.
  shift @databases; #remove orig from database array
  $log->write_to ("@databases\n") if ($debug);
  foreach my $database (@databases) {
    $log->write_to ("Transfering $orig to $wormpub/camace_${database}\n");
    print "Transfering $orig to $wormpub/camace_${database}\n" if ($debug);
    $wormbase->run_script("TransferDB.pl -start  $orig -end $wormpub/camace_${database} -split -database -wspec", $log) && $log->write_to ("Failed to run TransferDB.pl for camace_${database}\n");
  }
  $log->write_to ("CAMACE SPLITS UPDATED\n");
  print "CAMACE SPLITS UPDATED\n" if ($debug);
}

###################################################################
#(3b) Load data into split database that is required to aid curation. #
###################################################################
sub load_curation_data {
  my $sub_database = shift;
  my $database_path;
  if ($sub_database eq "camace") {
    $database_path = $wormbase->database('camace');
  }
  else {
    $database_path = $sub_database;
  }
  
  $log->write_to ("Loading curation data from $database_path\n");
  $ENV{'ACEDB'} = $database_path;
  my $file;
  my @files;
  my $acefiles;
  if (-e $wormbase->acefiles."/sorted_exons.ace") {
    $acefiles = $wormbase->acefiles;
    $log->write_to ("The BUIlD is still in place, using autoace/acefiles\n");
  }
  else {
    $acefiles = $wormbase->database('current')."/acefiles";
    $log->write_to ("The BUIlD has finished...using currentDB/acefiles\n");
  }
  
  push (@files,"$wormpub/CURATION_DATA/anomalies_elegans.ace.WS${WS_version}",
	"$wormpub/CURATION_DATA/Tiling_array_data/tiling_array.ace",
	"$wormpub/CURATION_DATA/assign_orientation.WS${WS_version}.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_TEC_RED_homol.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_21urna_homol.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_twinscan.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_jigsaw.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_genefinder.ace",
	"$acefiles/misc_DB_remark.ace",
	"$acefiles/elegans_blastx.ace",
	"$acefiles/elegans_blastp.ace",
	"$wormpub/BUILD_DATA/MISC_DYNAMIC/waba.ace",
	"$acefiles/inverted_repeats.ace",
	"$acefiles/repeat_homologies.ace",
	"$acefiles/operon_coords.ace",
	"$acefiles/feature_SL1.ace",
	"$acefiles/feature_SL2.ace",
	"$acefiles/feature_polyA_signal_sequence.ace",
	"$acefiles/feature_polyA_site.ace",
	"$acefiles/feature_binding_site.ace",
	"$acefiles/feature_binding_site_region.ace",
	"$wormpub/CURATION_DATA/PAD_DATA/elegans.public_names_ws196.ace", #needs scripting to update
	"$wormpub/CURATION_DATA/PAD_DATA/genomic_signals-no-splicing.ace", #needs remapping
  	"$wormpub/BUILD_DATA/MISC_DYNAMIC/misc_mass_spec_GenniferMerrihew.ace",
	"$acefiles/mass-spec-data.ace",
       ) unless ($sub_database eq "camace");
  push (@files,"$acefiles/misc_DB_remark.ace",
        "$wormpub/CURATION_DATA/assign_orientation.WS${WS_version}.ace",
       ) if ($sub_database eq "camace");
  
  foreach $file (@files) {
    $log->write_to ("Looking for $file......................\n");
    if (-e $file) {
      &loadace("$file", '${version}_curation_data_update') or die "Failed to load $file\n";
      $log->write_to ("SUCCESS! Loaded $file into $database_path\n\n");
    }
    else {
      $log->write_to ("!!WARNING!! File: $file does not exist\n\n");
    }
  }

  # REMOVE BOGUS GENE PREDICTIONS that don't have a method from any database called.
  $log->write_to ("Removing remove bogus gene models from $database_path\n");
  my $command;
  #remove bogus gene models
  $command  = "query find ALL_genes where !method\n";
  $command  .= "kill\n";
  $command  .= "y\n";
  $command  .= "save";
  $command  .= "quit";
  open (TACE,"| $tace") or die "Failed to open database connection\n";
  print TACE $command;
  close TACE;
  $log->write_to ("Removed bogus gene models\n");
  
  # upload BLAT results to database #
  unless (($sub_database eq "camace") or ($load_only)) {
    $log->write_to ("\n\nUpdate BLAT results in $database_path\n");
    $wormbase->run_script("load_blat2db.pl -all -dbdir $database_path", $log) && die "Failed to run load_blat2db.pl\n";
  }
}

#############################
#(3c) remove curation data. #
#############################
sub remove_data {
  my $sub_database = shift;
  my $database_path;
  if ($sub_database eq "camace") {
    $database_path = $wormbase->database('camace');
  }
  else {
    $database_path = $sub_database;
  }
  $log->write_to ("Removing curation data from $database_path\n");
  $ENV{'ACEDB'} = $database_path;
  my $command;
  #remove DB_remarks
  $command  = "query find CDS\n";
  $command  .= "Edit -D DB_remark\n";
  #remove jigsaw/genefinder/twinscan predictions.....basically all CDSs that aren't curated, history or Transposon_CDSs
  $log->write_to ("Warning: If CDSs are expected to have additional methods to curated, history and Transposon_CDS, this script needs to be modified OK!\n\n");
  $command .= "query method != curated AND method != Transposon_CDS AND method != history\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  #remove most of the homol data.
  $command .= "query find Homol_data where *BLAT_*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *signal_peptide\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *inverted\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *TEC_RED\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *waba\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *curation_anomaly\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *wublastx*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *tiling*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *RepeatMasker*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *Expr*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Homol_data where *expr_anomaly*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  #remove blastp data
  $command .= "query find Protein *WP*\n";
  $command .= "Edit -D Pep_homol\n";
  $command .= "clear\n";
  #remove sequence tiling data.
  $command .= "query find Sequence *tiling*\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  #remove anomoly motif data.
  $command .= "query find Motif *curation_anomaly\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  #remove Feature mappings
  $command .= "query find Genome_sequence\n";
  $command .= "Edit -D Feature_object\n";
  $command .= "clear\n";
  #remove operon mappings
  $command .= "query find Sequence CHROM*\n";
  $command .= "Edit -D Operon\n";
  $command .= "clear\n";
  #remove public names
  $command .= "query find Gene\n";
  $command .= "Edit -D Public_name\n";
  $command .= "clear\n";
  #remove repeat data.
  $command .= "query find Feature_data *TRF\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  $command .= "query find Feature_data *inverted\n";
  $command .= "kill\n";
  $command .= "y\n";
  $command .= "clear\n";
  #remove mass_spec peptide data if canonical or 
  if (($sub_database eq "camace") or ($remove_only)) {
    $command .= "query find Homol_data where *Mass-spec\n";
    $command .= "kill\n";
    $command .= "y\n";
    $command .= "clear\n";
    $command .= "query find Protein \"MSP:*\"\n";
    $command .= "kill\n";
    $command .= "y\n";
    $command .= "clear\n";
    $command .= "query find Protein\n";
    $command .= "Edit -D Pep_homol \"MSP*\"\n";
    $command .= "clear\n";
    $command .= "query find Mass_spec_peptide MSP*\n";
    $command .= "kill\n";
    $command .= "y\n";
    $command .= "clear\n";
  }
  #Save the database.
  $command .= "save\n";
  $command .= "quit\n";
  $log->write_to ("Removing curation data data.\n");
  open (TACE,"| $tace") or die "Failed to open database connection\n";
  print TACE $command;
  close TACE;
}

####################
# Usage subroutine #
####################
sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}
__END__

=pod

=head2 NAME - merge_split_camaces

=head1 USAGE:

=over 4

=item merge_split_camaces [-options]

=back

merge_split_camaces is a wrapper with options to automate the merging of the
working split copies of camace, load the updates into the current version
of camace (with some additional housekeeping), and finally running the
TransferDB.pl jobs to (re)generate the working split copies of camace.

merge_split_camaces mandatory arguments:

=over 4

=item none

=back

merge_split_camaces optional arguments:

=over 4

=item -merge, Generate diff files from split camace databases
 
=item -update, Upload diff files to ~wormpub/DATABASES/camace and add BLAT, Locus data

=item -split, Transfer ~wormpub/DATABASES/camace into split camace databases

=item -help, Help page

=item -debug, Verbose/Debug mode

=back

=head1 RUN REQUIREMENTS:

=back

merge_split_camaces has been divorced from the /wormsrv2 disk

=head1 RUN OUTPUT:

=back

=head1 EXAMPLES:

=over 4

=item merge_split_camaces -merge -all 

=back

Dumps the relevant classes from camace_orig and the split databases, runs an
acediff to calculate the changes. This acediff file is then reformated to take
account of the acediff bug.

=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut

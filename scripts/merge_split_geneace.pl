#!/software/bin/perl -w
#
# merge_split_geneaces.pl
# 
# A script to make multiple copies of core species database for curation, and merge them back again
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2015-07-02 10:23:29 $
#====================
#perl ~/wormbase/scripts/merge_split_geneaces.pl -update -mz3 -pad -species elegans -test -version $MVERSION > /nfs/wormpub/geneace_orig/WSXXX -WSXXY/load_data.txt

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
my $mz3;                   # Use MZ3's split
my $skd;                   # Use SKD's split
my $curation;              # create a single curation database.
my $merge;                 # Merging databases
my $split;                 # Splitting databases
my $update;                # Update current database
my $debug;                 # Debug option
my $help;                  # Help menu
my $WS_version;            # Removes final wormsrv2 dependancy.
my $store;                 # Storable not needed as this is not a build script!
my $species;
my $test;
my $wormbase;
my $nodump;                # don't dump from split geneaces.
my $nosplit;               # use this option with the split if you only remove the mass_spec and tiling array dataata.
my $nochecks;              # dont run camcheck as this can take ages.
my $verbose;               # Additional printout
my $logfile;
  GetOptions (
	      "all"           => \$all,
	      "pad"           => \$pad,
	      "mz3"           => \$mz3,
	      "skd"           => \$skd,
              "curation"      => \$curation,
	      "merge"         => \$merge,
	      "split"         => \$split,
	      "update"        => \$update,
	      "help"          => \$help,
	      "debug:s"       => \$debug,
              "logfile:s"     => \$logfile,
	      "version:s"     => \$WS_version,
	      "store"         => \$store,
	      "species:s"     => \$species,
	      "test"          => \$test,
	      "nodump"        => \$nodump,
	      "nosplit"       => \$nosplit,
	      "nochecks"      => \$nochecks,
	      "verbose"       => \$verbose,
	     );



# Help pod if needed
&usage("Help") if ($help);

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug,
			     -test => $test,
			     -organism => $species
			   );
}

my $log = $logfile ? Log_files->make_log($logfile, $debug) : Log_files->make_build_associated_log($wormbase);
my $tace = $wormbase->tace;

if (!defined($WS_version)) {
  $WS_version = $wormbase->get_wormbase_version;
}

my $next_build = $WS_version + 1;

$log->write_to ("WS_version : $WS_version\tWS_next : $next_build\n") if ($debug);

# load @databases array with user database names.
my @databases; #array to store what splits are to be merged.


#Are you trying to do some work?
unless ($split || $merge || $update ) {
    $log->log_and_die("ERROR you havent specified a main function split/merge/update!!!\n");
  }

# directory paths
my $wormpub = $wormbase->wormpub;
our $canonical;
our $orig;
our $root;


#Impose a test db 
if ($test) {
    push(@databases,"testorig");
    push(@databases,"testcuration") if ($curation || $all);
}
else {
    push(@databases,"orig");
    push(@databases,"pad") if ($pad || $all);
    push(@databases,"mz3") if ($mz3 || $all);
    push(@databases,"skd") if ($skd || $all);
    push(@databases,"curation") if ($curation);
}



if ($test){
    $orig = $wormpub."/geneace_testorig";
    $canonical = $wormpub."/DATABASES/TEST_DATABASES/geneace";
}
else  {
    $orig = $wormpub."/geneace_orig";
    $canonical = $wormbase->database('geneace');

}
$root = 'geneace';

our $directory   = $orig."/WS${WS_version}-WS${next_build}";
$log->write_to ("OUTPUT_DIR: ".$orig."/WS${WS_version}-WS${next_build}\n\n");

# what classes of data do we want to dump?
my @classes = ('Gene', 'Feature', 'Variation', 'Strain');


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

  # dumps the Gene Feature Variation classes from the database
  &dump_sdb unless ($nodump);
  #remove the 1st element off the @databases array as this will be the orig reference database itself.
  shift(@databases);
  foreach my $database (@databases) {
    foreach my $class (@classes) {
      my $path_new = $directory . "/${class}_${database}.ace";
      my $path_ref;
      if ($test) {
	  $path_ref = $directory . "/${class}_testorig.ace";
      }
      else {
	  $path_ref = $directory . "/${class}_orig.ace";
      }
      if ($debug) {
	  $wormbase->run_script("acediff.pl -debug $debug -test -reference $path_ref -new $path_new -output $directory/update_${class}_${database}.ace -debug pad", $log) && die "acediff.pl Failed for ${path_new}\n";
      }
      else {
	  $wormbase->run_script("acediff.pl -reference $path_ref -new $path_new -output $directory/update_${class}_${database}.ace -debug pad", $log) && die "acediff.pl Failed for ${path_new}\n";
      }
    }
  }
  $log->write_to ("Phase 1 finished and all files can be found in $directory\n");
}

#################################################################
## (2) synchronises Canonical Database with the split versions ##
#################################################################

if ($update) {
  shift (@databases); 
  &update_canonical;
  $log->write_to ("Phase 2 finished $canonical is now updated\n");
}

#######################################################################################################################
## (3) Dump data from the canonical database and load into a stripped down curation database amd reference databases ##
#######################################################################################################################

if ($split) {
  $log->write_to ("Removing old split databases and Copying $canonical database to the split database locations\n");
  &create_dump_files;
  &split_databases unless ($nosplit);
  $log->write_to ("\nPhase 3 finished. The geneace curation database has been updated\n");
  print "Phase 3 finished. The geneace curation database has been updated\n" if ($debug);
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

###################################################
#(1a) Dump files from curation and orig databases #
###################################################
sub dump_sdb {
  #dumps out subset of classes from curation splits and processes the files to be loaded back to Canonical Database
  #array of classes to be dumped
  my $database_path;
  my $path;

  foreach my $database (@databases) {

    $database_path = $wormpub."/${root}_${database}";
    #$ENV{'ACEDB'} = $database_path;
    
    foreach my $class (@classes) {
      $log->write_to ("dumping $class class from ${root}_${database}\n");
      $path = "$directory/" . "${class}_${database}.ace";
      system("touch $path");#needed as tace doesn't produce empty files.
      &dumpace("$class",$path,$database_path);
      print "dumped $class class from ${root}_${database}\n\n" if $debug;
      $log->write_to ("dumped $class class from ${root}_${database}\n\n");
    }
  }
}

#####################################################
#(1b) data dumping routine used by dump_sdb routing #
#####################################################
sub dumpace {
  my $class    = shift;
  my $filepath = shift;
  my $ddb = shift;
  my $command = "nosave\nquery find $class\nshow -a -f $filepath\nquit\n";
  $log->write_to ("\nFilename: $filepath -> ${ddb}\n");
  print "Opening $ddb for edits\n" if ($verbose);
  open (TACE, "echo '$command' | $tace $ddb |" ) || die "Failed to open database connection to $ddb\n" ;
  while (<TACE>) {
    #keeps pipe open until command quit executes.
  }
  close TACE;
}

##########################
#(2a) Loading subrouting #
##########################
sub loadace {
  my $filepath = shift;
  my $tsuser   = shift;
  my $ldb = shift;

#####################could do with user input prompt not die!!!!!

  if ($tsuser =~ /\S+\s/) {
    &testuser($tsuser);
  }
  my $command = "pparse $filepath\nsave\nquit\n";
  
  # dump out from ACEDB
  $log->write_to ("\nFilename: $filepath -> ${ldb}\n");
print "Opening $ldb for edits\n" if ($verbose);
  open (TACE, "| $tace $ldb -tsuser $tsuser") || die "Couldn't open $ldb\n";
  print TACE $command;
  close TACE;
  $log->write_to ("SUCCESS! Loaded $filepath into ${ldb}\n\n");
  return 1;
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
sub update_canonical {
  # upload processed diff files into Canonical Database.
  $log->write_to ("Upload diff files to $canonical");
#  $ENV{'ACEDB'} = $canonical;
  
  ##  Upload the split update files you have just created  ##
  foreach my $database (@databases) {
    foreach my $class (@classes) {
      &loadace("$directory/update_${class}_${database}.ace", "${database}_${class}", "$canonical");
    }
  }


  unless ($nochecks) {
    # Gene structures
    $log->write_to ("\nRunning geneace checks\n");
    print "\nRunning geneace checks\n" if ($debug);
	$wormbase->run_script("GENEACE/geneace_check.pl -class strain -database $canonical -debug pad", $log) && die "Failed to run Strain checks\n";
	$wormbase->run_script("GENEACE/geneace_check.pl -class Allele -skipmethod Million_mutation -skipmethod SNP -skipmethod NemaGENETAG_consortium_allele -skipmethod KO_consortium_allele -skipmethod NBP_knockout_allele -skipmethod WGS_Hobert -skipmethod WGS_Rose -skipmethod WGS_Jarriault -skipmethod Mos_insertion -skipmethod Transposon_insertion -skipmethod WGS_McGrath -skipmethod CGH_allele -skipmethod WGS_Flibotte  -debug pad", $log) && die "Failed to run Allele checks\n";
    $wormbase->run_script("GENEACE/geneace_check.pl -class gene -debug pad", $log) && die "Failed to run Gene checks\n";
    $log->write_to ("Checks Finished, check the build log email for errors.\n");
  }
}


#######################
#(3a) Data dispersion needs a rewrite to generate the geneace_orig from partial dumps # 
#######################

sub split_databases {
#reinitialise the $species_splits
  foreach my $database (@databases) {
    my $split_db = $wormpub."/${root}_${database}";
    my $tsuser = "Initialise";
    if (-e $split_db."/database/ACEDB.wrm") {
      $log->write_to ("Destroying $split_db\n");
      print "Destroying $split_db\n" if ($debug);
      $wormbase->run_command("rm -rf $split_db/database/ACEDB.wrm", $log) && die "Failed to remove $split_db/database/ACEDB.wrm\n";
    }
    else {
      $log->write_to ("Databases doesn't exist so creating $split_db\n");
      print "Databases doesn't exist so creating $split_db\n";
      $wormbase->run_command("mkdir $split_db", $log) && die "Failed to create $split_db dir\n" unless (-e $split_db);
      $wormbase->run_command("mkdir $split_db/database", $log) && die "Failed to create a database dir\n" unless (-e $split_db."/database");
      $wormbase->run_command("cp -r ${wormpub}/wormbase-pipeline/wspec ${split_db}/", $log) && die "Failed to copy the wspec dir\n" unless (-e $split_db."/wspec");
      $wormbase->run_command("cp -rf ${wormpub}/wormbase-pipeline/wspec/models.wrm ${split_db}/wspec/", $log) && die "Failed to force copy the models file\n" unless (-e $split_db."/wspec/models.wrm");;
      $wormbase->run_command("chmod g+w ${split_db}/wspec/*", $log);
      $wormbase->run_command("cp -r $wormpub/wormbase-pipeline/wgf $split_db/", $log) && die "Failed to copy the wgf dir\n" unless (-e $split_db."/wgf");
    }
    my $command = "y\nsave\nquit\n";
    print "Opening $split_db for edits\n" if ($verbose);
    open (TACE, "| $tace $split_db -tsuser $tsuser") || die "Couldn't open $split_db\n";
    print TACE $command;
    close TACE;
    $log->write_to ("(Re)generated $split_db\n");
    print "(Re)generated $database\n" if ($debug);
 
    &loadace("$orig/${WS_version}_gene.ace", "Gene", "$split_db");
    &loadace("$orig/${WS_version}_var.ace", "Variation", "$split_db");
    &loadace("$orig/${WS_version}_strain.ace", "Strain", "$split_db");
    &loadace("$orig/${WS_version}_feature.ace", "Strain", "$split_db");
    if ($split_db =~ /orig/) {
        #protect the reference copy of the database now it has been updated and copied over to the curation splits.                                                                                             
	$wormbase->run_command ("touch ${split_db}/database/lock.wrm", $log)
    }
  }
  $log->write_to ("@databases SPLIT(S) UPDATED\n");
  print "@databases SPLIT(S) UPDATED\n" if ($debug);
}

##################################
# create_dump_files from geneace #
##################################
sub create_dump_files {
  my $Gene_data = "$orig/${WS_version}_gene.ace";
  my $Var_data = "$orig/${WS_version}_var.ace";
  my $Strain_data = "$orig/${WS_version}_strain.ace";
  my $Feature_data = "$orig/${WS_version}_feature.ace";

  my $refdb = "${\$wormbase->database('geneace')}";
  
  if (-e($Gene_data)){
    $log->write_to ("\nNo need to recreate Gene data for ${WS_version}\n");
    print "\nUsing $Gene_data as source of gene data\n";
  }
  else {
    my $gcommand = "\nnosave\nquery find Gene\nshow -a -f $Gene_data\nquit\n";
    # dump out from ACEDB
    $log->write_to ("\nDumped gene data\n");
    print "Opening $refdb for edits\n" if ($verbose);
    open (TACE, "| $tace $refdb -tsuser test") || die "Couldn't open $refdb\n";
    print TACE $gcommand;
    close TACE;

  }
  if (-e($Var_data)){
    $log->write_to ("\nNo need to recreate Variation data for ${WS_version}\n");
    print "\nUsing $Var_data\n";
  }
  else {
    my $vcommand = "\nnosave\nquery Find Variation WHERE method != \"WGS*\" AND method != \"SNP*\" AND method != \"Million_mutation\" AND method != \"NBP_knockout_allele\" AND method != \"KO_consortium_allele\" AND method != \"NemaGENETAG_consortium_allele\"\nshow -a -f $Var_data\nquit\n";
    # dump out from ACEDB                                                                                                                                                          
    $log->write_to ("\nDumped Var data\n");
    print "Opening $refdb for edits\n" if ($verbose);
    open (TACE, "| $tace $refdb -tsuser test") || die "Couldn't open $refdb\n";
    print TACE $vcommand;
    close TACE;
  }
    if (-e($Strain_data)){
	$log->write_to ("\nNo need to recreate Strain data for ${WS_version}\n");
	print "\nUsing $Strain_data\n";
    }
    else {
	my $scommand = "\nnosave\nquery find Strain\nshow -a -f $Strain_data\nquit\n";
    # dump out from ACEDB                                                                                                                                                          
	$log->write_to ("\nDumped strain data\n");
	print "Opening $refdb for edits\n" if ($verbose);
	open (TACE, "| $tace $refdb -tsuser test") || die "Couldn't open $refdb\n";
	print TACE $scommand;
	close TACE;
    }

    if (-e($Feature_data)){
        $log->write_to ("\nNo need to recreate Feature data for ${WS_version}\n");
        print "\nUsing $Feature_data\n";
    }
    else {
        my $fcommand = "\nnosave\nquery find Feature\nshow -a -f $Feature_data\nquit\n";
    # dump out from ACEDB                                                                                                                                                         \
                                                                                                                                                                                   
        $log->write_to ("\nDumped strain data\n");
        print "Opening $refdb for edits\n" if ($verbose);
        open (TACE, "| $tace $refdb -tsuser test") || die "Couldn't open $refdb\n";
        print TACE $fcommand;
        close TACE;
    }
  
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

=head2 NAME - merge_split_geneaces.pl

=head1 USAGE:

=over 4

=item merge_split_geneaces.pl [-options]

=back

merge_split_geneaces.pl is a wrapper with options to automate the merging of the
working split/curation copie(s), load the updates into the canonical version
of geneace (with some additional housekeeping), and finally dump mininal data to create the curation database.

merge_split_geneaces.pl mandatory arguments:

=over 4

=item none

=back

merge_split_geneaces.pl optional arguments:

=over 4

=item -merge, Generate diff files from split databases
 
=item -update, Upload diff files to ~wormpub/DATABASES/geneace

=item -split, Dumps data and regenerated the reference and curation databases ~wormpub/geneace_*

=item -help, Help page

=item -debug, Verbose/Debug mode

=back

=head1 RUN REQUIREMENTS:

=back

merge_split_geneaces.pl has been divorced from the /wormsrv2 disk

=head1 RUN OUTPUT:

=back

=head1 EXAMPLES:

=over 4

=item merge_split_geneaces.pl -merge -all 

=back

Dumps the relevant classes from orig and the split databases, runs an
acediff to calculate the changes. This acediff file is then reformated to take
account of the acediff bug.

=head1 AUTHOR - Paul Davis (Reusing subroutined written by Daniel Lawson).

Email pad@ebi.ac.uk

=cut

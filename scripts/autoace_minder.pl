#!/usr/local/bin/perl5.8.0 -w
#
# autoace_minder
# 
# Originally by Dan Lawson with many, many modifications by the Wormbase crew.
#
# Usage : autoace_minder.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2004-09-27 11:11:04 $



#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use vars;
use File::Copy;
use Coords_converter;

# is this script being run as user wormpub???
&test_user_wormpub;

##############################
# command-line options       #
##############################

my $initial;		# Start the build process 
my $unpack;		# unpack primaries	
my $gffdump;		# dump gff files
my $gffsplit;           # split gff files
my $buildpep;		# Build wormpep
my $buildrna;		# Build wormrna
my $prepare_blat;	# Prepare for blat, copy autoace before blatting, run blat_them_all.pl -dump
my $blat_est;           # run blat for ests
my $blat_mrna;          # run blat for mrnas
my $blat_ncrna;         # run blat for ncrnas
my $blat_ost;           # run blat for osts
my $blat_embl;          # run blat for non-Wormbase CDS genes in EMBL
my $blat_tc1;           # run blat for TC1 insertion sequences
my $blat_nematode;      # run blat for non-C. elegans nematode ESTs
my $blat_all;           # run all five types of blat jobs
my $addblat;		# parse (all) blat files
my $addbriggsae;        # add briggsae ace files
my $addhomol;           # parse similarity data from /ensembl_dump
my $utrs;               # generate UTR dataset
my $map;		# map PCR and RNAi
my $acefile;		# Write .acefiles
my $build;		# Build autoace 
my $builddb;		# Build autoace : DB only
my $buildchrom;		# Build autoace : CHROMOSOMES directory	
my $buildrelease;	# Build autoace : Release directory
my $buildtest;          # Test autoace build
my $test;               # In test build mode, will write to mirrored copy of autoace in ~wormpub/AUTOACE_TEST
my $quicktest;          # Same as -test but will only work with just chromosome III (where appropriate)
my $agp;		# make/check agp files
my $confirm;            # Confirm gene models (EST|mRNA)
my $ftp;		# Public release to FTP site
my $debug;		# debug mode
my $verbose;            # verbose mode - more output to screen
my $help;		# Help/Usage page
my $dbcomp;		# runs dbcomp script
my $operon;             # generate operon data
my $load;               # generic file loading routine
my $tsuser;             # tsuser setting to go with -load
my $am_option;          # track which option has been used (for logging purposes)
my $errors = 0;         # keep track of errors in each step (from bad system calls), use in subject line of email 
my $remarks;

GetOptions (
	    "agp"	     => \$agp,
	    "initial"        => \$initial,
	    "unpack"         => \$unpack,
	    "gffdump"        => \$gffdump,
	    "gffsplit"       => \$gffsplit,
	    "buildpep"       => \$buildpep,
	    "buildtest"      => \$buildtest,
	    "buildrna"       => \$buildrna,
	    "prepare_blat"   => \$prepare_blat,
	    "blat_est"       => \$blat_est,
	    "blat_ost"       => \$blat_ost,
	    "blat_mrna"      => \$blat_mrna,
	    "blat_ncrna"     => \$blat_ncrna,
	    "blat_embl"      => \$blat_embl,
	    "blat_tc1"       => \$blat_tc1,
	    "blat_nematode"  => \$blat_nematode,
	    "blat_all"       => \$blat_all,
	    "addblat"        => \$addblat,
	    "addbriggsae"    => \$addbriggsae,
	    "addhomol"       => \$addhomol,
	    "utrs"           => \$utrs,
	    "map"            => \$map,
	    "acefile"        => \$acefile,
	    "build"          => \$build,
	    "builddb"        => \$builddb,
	    "buildchrom"     => \$buildchrom,
	    "buildrelease"   => \$buildrelease,
	    "confirm"        => \$confirm,
	    "dbcomp"	     => \$dbcomp,
	    "debug:s"        => \$debug,
	    "verbose"        => \$verbose,
	    "operon"         => \$operon,
	    "load:s"         => \$load,
	    "tsuser:s"       => \$tsuser,
	    "help"           => \$help,
	    "test"           => \$test,
	    "quicktest"      => \$quicktest,
	    "remarks"        => \$remarks,
);

# Help pod if needed
&usage(0) if ($help);


# check that -test and -quicktest haven't both been set.  Also...
# if -quicktest is specified, still need to make -test true, so that test mode runs 
# for those steps where -quicktest is meaningless (can't run on only one chromosome)

if($test && $quicktest){
  &usage(21);
}
($test = 1) if ($quicktest);


##############################
# Script variables (run)     #
##############################

# double check if -tsuser not set
$tsuser = "assumed_wormpub" if (!$tsuser);

# who will receive log file?
my $maintainers = "All";

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 
my $db_path   = "$basedir/autoace";
my $scriptdir = "$basedir/scripts";               


# build flag path
our %flag = (
	     'A1'          => 'A1:Build_in_progress',
	     'A2'          => 'A2:Updated_WS_version',
	     'A3'          => 'A3:Unpack_FTP_databases',
	     'A4'          => 'A4:Primary_databases_on_wormsrv2',
	     'A5'          => 'A5:Wrote_acefiles_to_wormbase', 
	     'B1'          => 'B1:Made_autoace_database',
	     'B1:ERROR'    => 'B1:ERROR_in_loaded_acefiles',
	     'B2'          => 'B2:Dump_DNA_files',
	     'B2:ERROR'    => 'B2:ERROR_in_dumping_DNA',
 	     'B3'          => 'B3:Made_agp_files',
	     'B3:ERROR'    => 'B3:ERROR_in_agp_files',
	     'B6'          => 'B6:BLAT_analysis',
	     'B6_est'      => 'B6a:BLAT_analysis_EST',
	     'B6_mrna'     => 'B6b:BLAT_analysis_mRNA',
	     'B6_ncrna'    => 'B6f:BLAT_analysis_ncRNA',
	     'B6_embl'     => 'B6c:BLAT_analysis_EMBL',
	     'B6_ost'      => 'B6d:BLAT_analysis_OST',
	     'B6_nematode' => 'B6e:BLAT_analysis_NEMATODE',
	     'B7'          => 'B7:Upload_BLAT_data',
	     'B8'          => 'B8:Upload_wublast_data',
	     'B9'          => 'B9:Upload_briggsae_data',
	     'B10'         => 'B10:Generate_UTR_data',
	     'B11'         => 'B11:Generate_operon_data',
	     'C1'          => 'C1:Dumped_GFF_files',
	     'C1:ERROR'    => 'C1:ERROR_in_dumping_GFF_files',
	     'C2'          => 'C2:Split_GFF_files',
	     'C3'          => 'C3:Map_PCR_products',
	     'C4'          => 'C4:Confirm_gene_models',
	     'D1A'         => 'D1:Build_wormpep_initial',
	     'D1'          => 'D1:Build_wormpep_final',
	     'D2'          => 'D2:Update_pepace',
	     'D3'          => 'D3:Add_protein_data_to_autoace',
	     'D4'          => 'D4:Build_wormrna',
	     'Z1'          => 'Z1:Make_autoace_complete'
	     );

our @chrom   = ('I','II','III','IV','V','X');
    @chrom   = ('III') if ($quicktest);
# logdir for lock files
our $logdir  = "$db_path/logs";

our $WSver_file   = "$db_path/wspec/database.wrm";
our $tace         =  &tace;
our ($WS_version,$WS_previous) = &get_WS_version;

# Set up logfile
my $rundate    = &rundate;
our $log = "$basedir/logs/autoace_minder.WS${WS_version}.${rundate}.$$";
open (LOG,">$log");
LOG->autoflush();
# print logfile header
&logfile_details;


#__ PREPARE SECTION __#

# A1:Build_in_progress & A2:Updated_WS_version
# Requires: 
&initiate_build    if ($initial);


# A3:Unpack_FTP_databases & A4:Primary_databases_on_wormsrv2
# Requires: A1   
&prepare_primaries if ($unpack);


# A5:Wrote_acefiles_to_wormbase 
# Requires: A1,A4
&make_acefiles     if ($acefile);

#__ BUILD SECTION __#

# B1:Make_autoace_database & B2:Dump_DNA_files
# Requires: A1,A4,A5
&make_autoace           if ($build || $builddb || $buildchrom || $buildrelease);
&check_make_autoace     if ($buildtest);


# B3:Make_agp_files
# Requires: A1,A4,A5,B1
&make_agp		if ($agp);


# Run blat_them_all.pl -dump
# Requires: A1,A4,A5,B1
&prepare_for_blat       if ($prepare_blat);


# B6:Blat_analysis 
# Requires: A1,A4,A5,B1
&blat_jobs              if ($blat_est || $blat_ost || $blat_mrna || $blat_ncrna || $blat_embl || $blat_tc1 || $blat_nematode || $blat_all);


# B7:Upload_BLAT_data
# requires B6
if ($addblat){
  # Check that blats were actually run
  &usage(19) unless (-e "$logdir/$flag{'B6'}");  
  &load_blat_results("all");
}


# B8:Upload_wublast_data & B9:Upload_briggsae_data
# Requires: A1,A4,A5,B1
&parse_homol_data       if ($addhomol);
&parse_briggsae_data    if ($addbriggsae);


# B10:Generate_UTR_data
# Requires: A1,A4,A5,B1,B6,B7
&generate_utrs          if ($utrs);


# B11:Generate operon data
# Requires: B10
&make_operons           if ($operon);




#__ PROCESS SECTION __#

# C1:Dumped_GFF_files  
# Requires: A1,A4,A5,B1
&dump_GFFs         if ($gffdump);



# C2:Split_GFF_files   
# Requires: A1,A4,A5,B1
&split_GFFs        if ($gffsplit);


# C3:Map_PCR_products  
# Requires: A1,A4,A5,B1
&map_features      if ($map);


# C4:Confirm_gene_models
# Requires: A1,A4,A5,B1
&confirm_gene_models   if ($confirm);


#__ ANCILLIARY DATA SECTION __#

# D1:Build_wormpep                    
&make_wormpep      if ($buildpep);

# add DB_remarks
&run_command( "$scriptdir/get_pfam.pl --database /wormsrv2/autoace")  if ( $remarks );

# D4:Build_wormrna                    
&make_wormrna      if ($buildrna);


#__ CHECK SECTION __#
&dbcomp		  if ($dbcomp);


#__ RELEASE SECTION __#


#__ LOAD DATA FILE __#
&load($load,$tsuser) if ($load);


##############################################
# close log file and mail $maintainer report #
##############################################

$rundate = &rundate;
print LOG "\n# autoace_minder finished at: $rundate ",&runtime,"\n";
close LOG;

my $subject_line;
$subject_line = "BUILD REPORT: $am_option";
$subject_line = "TEST BUILD REPORT: $am_option" if ($test);

# warn about errors in subject line if there were any
if($errors == 1){
  $subject_line .= " : $errors ERROR!";
}
elsif($errors > 1){
  $subject_line .= " : $errors ERRORS!!!";
}
&mail_maintainer("$subject_line",$maintainers,$log);

##############################
# hasta luego                #
#############################

exit(0);



#################################################################################
#                                                                               #
#                                                                               # 
#                T  H  E     S  U  B  R  O  U  T  I  N  E  S                    #
#                                                                               #
#                                                                               #
#################################################################################




#################################################################################
# initiate autoace build                                                        #
#################################################################################

sub initiate_build {
  $am_option = "-initial";

  local (*FLAG);
  my $cvs_file = "$db_path/wspec/database.wrm";

  # exit if build_in_progress flag is present
  &usage(1) if (-e "$logdir/$flag{'A1'}");
    
  # get old build version number, exit if no WS version is returned
  # add 1 for new build number, but just use '666' if in test mode

  my $WS_new_name;
  my $WS_version;

  if ($test){
    $WS_version = "666"; 
    $WS_new_name = "666";
  }
  else{
    $WS_version = &get_wormbase_version;
    $WS_new_name = $WS_version +1;
  }
  &usage("No_WormBase_release_number") if (!defined($WS_version));


  # make new build_in_process flag
  system("touch $logdir/$flag{'A1'}");
  open (FLAG, ">>$logdir/$flag{'A1'}"); 
  print FLAG "WS$WS_new_name\n";
  close FLAG;
    
  # make sure that the database.wrm file is younger
  sleep 10;

  # update database.wrm using cvs
  &run_command("sed 's/WS${WS_version}/WS${WS_new_name}/' < $cvs_file > ${cvs_file}.new");
  my $status = move("$db_path/wspec/database.wrm.new", "$cvs_file") or print LOG "ERROR: renaming file: $!";
  print LOG "ERROR: Couldn't move file: $!\n" if ($status == 0);

  # make a log file in $basedir/autoace/logs
  system("touch $logdir/$flag{'A2'}");

  # check that top doesn't reveal strange processes running
  my $top = `top`;

  # add lines to the logfile
  print LOG "Updated WormBase version number to WS$WS_new_name\n";
  print LOG "You are ready to build another WormBase release\n";
  print LOG "Please tell camace and geneace curators to update their database to use the new models!!!\n\n";
  print LOG "Please also check following 'top' output to see if there are stray processes that should\n";
  print LOG "be removed:\n$top\n\n";

  # this will force a refresh of the coordinate files.
  my $coords = Coords_converter->invoke($db_path,1);
    
}
#__ end initiate_build __#

#################################################################################
# get_WS_version                                                                #
#################################################################################
# Requirements : [01] - Presence of the database.wrm file in wspec
#              : [02] - WormBase.pm modules
#
# Checks       : [01] - Fail if no WS_version is returned.

# Does         : [01] - checks the build_in_progess flag is older than the version
#              : [02] - assigns last WS version to WS_version - 1

sub get_WS_version {

  # force WS release to be 666 if in test mode
  if($test){
    $WS_version = "666";
  }
  else{
    $WS_version = &get_wormbase_version;
  }

  # exit if no WS version is returned
  &usage("No_WormBase_release_number") if (!defined($WS_version));
  
  print "$logdir/$flag{'A1'} : " . (-M "$logdir/$flag{'A1'}") . "\n" if ($verbose);
  print "$WSver_file : "         . (-M "$WSver_file") . "\n" if ($verbose);
  
  
  
  # exit if WS version (database.wrm) is older than build_in_process flag
  # i.e. the WS version has been changed since the start of the build
  # ignore if the A1 flag doesn't exist
  if (-e "$logdir/$flag{'A1'}") {
    &usage(11) if (-M "$logdir/$flag{'A1'}" < -M "$WSver_file");
  }
  
  # manipulate to assign last WS release version
  $WS_previous = $WS_version -1;
  return($WS_version,$WS_previous);
}
#__ end get_WS_version __#

#################################################################################
# prepare primary databases                                                     #
#################################################################################
#
# Requirements : [01] - Presence of the Primary_databases_used_in_build file
#
# Checks       : [01] - Fail if the build_in_progess flag is absent.
#              : [02] - Fail if the Primary_databases_used_in_build file is absent
#
# Does         : [01] - checks the Primary_database_used_in_build data
#              : [03] - writes log file A2:Update_WS_version to $basedir/autoace/logs

sub prepare_primaries {
  $am_option = "-unpack";
  # exit unless build_in_progress flag is present
  &usage(12) unless (-e "$logdir/$flag{'A1'}");

  # exit if the Primary_databases_used_in_build is absent
  &usage(13) unless (-e "$logdir/Primary_databases_used_in_build");
 
  local (*LAST_VER);
  my ($stlace_date,$brigdb_date,$citace_date,$cshace_date) = &FTP_versions;
  my ($stlace_last,$brigdb_last,$citace_last,$cshace_last) = &last_versions;
  my $options = "";

  # use test mode if autoace_minder -test was specified
  $options .= " -test" if ($test);

  # stlace
  print "\nstlace : $stlace_date last_build $stlace_last";
  unless ($stlace_last eq $stlace_date) {
    $options .= " -stlace $stlace_date";
    print "  => Update stlace";
  }
  
  # brigdb
  print "\nbrigdb : $brigdb_date last_build $brigdb_last";
  unless ($brigdb_last eq $brigdb_date) {
    $options .= " -brigace $brigdb_date";
    print "  => Update brigdb";
  }
  
  # citace
  print "\ncitace : $citace_date last_build $citace_last";
  unless ($citace_last eq $citace_date) {
    $options .= " -citace $citace_date";
    print "  => Update citace";
  }
  
  # cshace
  print "\ncshace : $cshace_date last_build $cshace_last";
  unless ($cshace_last eq $cshace_date) {
    $options .= " -cshace $cshace_date";
    print "  => Update cshace";
  }
  
  print "\n\nrunning unpack_db.pl $options\n";
  
  # confirm unpack_db details and execute
  unless ($options eq "") {
    print "Do you want to unpack these databases ?\n";
    my $answer=<STDIN>;
    &usage(2) if ($answer ne "y\n");
    
    &run_command("$scriptdir/unpack_db.pl $options");
 }
  
  # make a unpack_db.pl log file in /logs
  system("touch $logdir/$flag{'A3'}");

  if($test){
    print LOG "WARNING: Can't transfer geneace and camace from /wormsrv1.  You will have to do that by hand!\n";  
  }
  else{
    # transfer /wormsrv1/camace to $basedir/camace 
    &run_command("$scriptdir/TransferDB.pl -start /wormsrv1/camace -end $basedir/camace -database");
    # transfer /wormsrv1/geneace to $basedir/geneace 
    &run_command("$scriptdir/TransferDB.pl -start /wormsrv1/geneace -end $basedir/geneace -database");
  }

  
  #################################################
  # Check that the database have unpack correctly #
  #################################################
    
  # rewrite Primary_databases_used_in_build
  open (LAST_VER, ">$logdir/Primary_databases_used_in_build");
  print LAST_VER "stlace : $stlace_date\n"; 
  print LAST_VER "brigdb : $brigdb_date\n"; 
  print LAST_VER "citace : $citace_date\n"; 
  print LAST_VER "cshace : $cshace_date\n"; 
  close LAST_VER;
  
  # make a unpack_db.pl log file in /logs
  system("touch $logdir/$flag{'A4'}");
  
}
#__ end of prepare_primaries __#

##################
# FTP_versions   #
##################

sub FTP_versions {

  local (*STLACE_FTP,*BRIGDB_FTP,*CITACE_FTP,*CSHACE_FTP);
  
  my $stlace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl/stlace_*";
  my $brigdb_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl/brigdb_*";
  my $citace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/caltech/citace_*";
  my $cshace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/csh/cshl_*";
  my ($stlace_date,$brigdb_date,$citace_date,$cshace_date);
  
  # stlace
  open (STLACE_FTP, "/bin/ls -t $stlace_FTP |")  || die "cannot open $stlace_FTP\n";
  while (<STLACE_FTP>) { chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $stlace_date = substr($1,-2).$2.$3; last; }
  close STLACE_FTP; 
  
  # brigdb
  open (BRIGDB_FTP, "/bin/ls -t $brigdb_FTP |") || die "cannot open $brigdb_FTP\n";
  while (<BRIGDB_FTP>) { chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $brigdb_date = substr($1,-2).$2.$3; last; }
  close BRIGDB_FTP; 
  
  # citace
  open (CITACE_FTP, "/bin/ls -t $citace_FTP |") || die "cannot open $citace_FTP\n";
  while (<CITACE_FTP>) { chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $citace_date = substr($1,-2).$2.$3; last; }
  close CITACE_FTP; 
  
  # cshace
  open (CSHACE_FTP, "/bin/ls -t $cshace_FTP |") || die "cannot open $cshace_FTP\n";
  while (<CSHACE_FTP>) { chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $cshace_date = substr($1,-2).$2.$3; last; }
  close CSHACE_FTP; 
  
  # return current dates as 6-figure string
  return($stlace_date,$brigdb_date,$citace_date,$cshace_date);
}

#####################################################################################################################

sub last_versions {
    
    local (*LAST_VER);
    my ($stlace_last,$brigdb_last,$citace_last,$cshace_last);

    open (LAST_VER, "<$logdir/Primary_databases_used_in_build") || usage("Primary_databases_file_error");
    while (<LAST_VER>) {
	$stlace_last = $1 if /^stlace \: (\d+)$/;
	$brigdb_last = $1 if /^brigdb \: (\d+)$/;
	$citace_last = $1 if /^citace \: (\d+)$/;
	$cshace_last = $1 if /^cshace \: (\d+)$/;
    }
    close LAST_VER;

    usage("Absent_stlace_database") if (!defined $stlace_last);
    usage("Absent_brigdb_database") if (!defined $brigdb_last);
    usage("Absent_citace_database") if (!defined $citace_last);
    usage("Absent_cshace_database") if (!defined $cshace_last);

    # return last version dates as 6-figure string
    return($stlace_last,$brigdb_last,$citace_last,$cshace_last);

}
#__ end prepare_primaries __#


#################################################################################
# make_acefiles                                                               
#################################################################################

sub make_acefiles {
  $am_option = "-acefile";

  # exit unless build_in_progress flag is present
  &usage("Build_in_progress_absent") unless (-e "$logdir/$flag{'A1'}");
  
  # exit unless A4:Primary_databases_on_wormsrv2
  &usage("Build_in_progress_absent") unless (-e "$logdir/$flag{'A4'}");

  my $command = "$scriptdir/make_acefiles.pl";
  $command .= " -test" if ($test);
  &run_command($command);

  # make a make_acefiles log file in /logs
  system("touch $logdir/$flag{'A5'}");
  
}
#__ end make_acefiles __#


#################################################################################
# make_autoace                                                                  #
#################################################################################
# Requirements : [01] acefiles in $basedir/wormbase directories
#                [02] config file to drive loading of data (autoace.config)

# Checks       : [01] - Fail if the build_in_progess flag is absent.
#              : [02] - Fail if the Primary_databases_used_in_build file is absent

# Does         : [01] - checks the Primary_database_used_in_build data
#              : [03] - writes log file A1:Update_WS_version to $basedir/autoace/logs

sub make_autoace {
  $am_option = "-build";
  # quit if make_acefiles has not been run
  &usage(8) unless (-e "$logdir/$flag{'A5'}");
  
  if ($build || $builddb) { 

    open (EMAIL,  "|/bin/mailx -s \"WormBase build reminder\" \"wormbase\@sanger.ac.uk\" ");
    print EMAIL "Dear builder,\n\n";
    print EMAIL "You have just run autoace_minder.pl -build.  This will probably take 5-6 hours\n";
    print EMAIL "to run.  You should therefore start work on the blast pipeline. So put down that\n";
    print EMAIL "coffee and do some work.\n\n";
    print EMAIL "Yours sincerely,\nOtto\n";
    close (EMAIL);

    my $command = "$scriptdir/make_autoace.pl --database $db_path --buildautoace";
    $command .= " --test" if ($test);
    &run_command($command);
    
    # test the build for loading errors
    my $builderrors = &check_make_autoace;
    
    # errors in the make_autoace log file
    if ($builderrors > 1) {
      system("touch $logdir/$flag{'B1:ERROR'}");
      &usage("Errors_in_loaded_acefiles");
    }
    
    # make a make_autoace log file in /logs
    system("touch $logdir/$flag{'B1'}");

    # Update Common_data clone2accession info, genes2lab, and worm_genes2cgc (uses geneace)
    $command = "$scriptdir/update_Common_data.pl --build --clone2acc --genes2lab --worm_gene2cgc --worm_gene2class";
    $command .= " --test" if ($test);
    &run_command($command);
  }
  
  if ($build || $buildchrom) {
    
    # quit if you haven't built a database
    &usage(14) unless (-e "$logdir/$flag{'B1'}");
    
    # quit if you have errors in the build
    &usage("Errors_in_loaded_acefiles") if (-e "$logdir/$flag{'B1:ERROR'}");

    my $command = "$scriptdir/chromosome_dump.pl --dna --composition";

    # run in test mode?
    if($quicktest){
      $command .= " --quicktest";
    }
    elsif($test){
      $command .= " --test";
    }

    my $status = &run_command($command);
    if ($status != 0){
      system("touch $logdir/$flag{'B2:ERROR'}");
    }    
    
    # make a make_autoace log file in /logs
    system("touch $logdir/$flag{'B2'}");
    
    # make the allcmid file needed for the farm
    
    $command = "query find genome_sequence\nDNA -f /wormsrv2/autoace/allcmid\nquit\n";
    open (WRITEDB, "| $tace $basedir/autoace ") || die "Couldn't open pipe to autoace\n";
    print WRITEDB $command;
    close (WRITEDB);
    
  }
  
  if ($buildrelease) {
    
    # quit if the build is not complete
    &usage(13) unless (-e "$logdir/$flag{'B1'}");
    
    # quit if you have errors in the build
    &usage("Errors_in_loaded_acefiles") if (-e "$logdir/$flag{'B1:ERROR'}");
    
    local (*MD5SUM_IN,*MD5SUM_OUT);
    
    &run_command("$scriptdir/make_autoace.pl -database $basedir/autoace --buildrelease"); 

    
    # make a make_autoace log file in /logs
    system("touch $logdir/$flag{'D1'}");
    
    # modify the md5sum output file to remove the Sanger specific path
    open (MD5SUM_OUT, ">$basedir/autoace/release/md5sum.temp")            || die "Couldn't open md5sum file out\n";
    open (MD5SUM_IN, "<$basedir/autoace/release/md5sum.WS${WS_version}")  || die "Couldn't open md5sum file in\n";
    while (<MD5SUM_IN>) {
      s/\$basedir\/autoace\/release\///g;
      print MD5SUM_OUT $_;
    }
    close MD5SUM_IN;
    close MD5SUM_OUT;
    my $status = move("$basedir/autoace/release/md5sum.temp", "$basedir/autoace/release/md5sum.WS${WS_version}");
    print LOG "ERROR: Couldn't move file: $!\n" if ($status == 0);
  }
}
#__ end make_autoace __#

#########################################################################################################################

sub check_make_autoace {
  local (*BUILDLOOK,*BUILDLOG);
  my $log;

  # set am_option if this is not being run as part of -build
  $am_option = "-buildtest" if ($buildtest);

  print LOG &runtime, ": Entering check_make_autoace subroutine\n";

  print "Looking at log file: $basedir/logs/make_autoace.WS${WS_version}*\n" if ($verbose);
  
  open (BUILDLOOK, "ls $basedir/logs/make_autoace.WS${WS_version}* |") || die "Couldn't list logfile out\n";
  while (<BUILDLOOK>) {
    chomp;
    $log = $_;
  }
  close BUILDLOOK;
  
  print "Open log file $log\n" if ($verbose);
  
  my ($parsefile,$parsefilename);
  my $builderrors = 0;

  open (BUILDLOG, "<$log") || die "Couldn't open logfile out\n";
  while (<BUILDLOG>) {
    if (/^\* Reinitdb: started parsing (\S+)/) {
      $parsefile = $1;
    }
    if ((/^\/\/ objects processed: (\d+) found, (\d+) parsed ok, (\d+) parse failed/) && ($parsefile ne "")) {
      my $object_count = $1;
      my $error_count  = $3;
      (printf "%6s parse failures of %6s objects from file: $parsefile\n", $error_count,$object_count) if $verbose;
      if ($error_count > 0) {
	$parsefilename = $parsefile;
	$parsefilename =~ s/$basedir//;
	printf LOG "%6s parse failures of %6s objects from file: $parsefilename\n", $error_count,$object_count;
	$builderrors++;
      }
      ($parsefile,$parsefilename) = "";
    }
  }
  close BUILDLOG;


  # look for objects with no Gene tag
  print LOG "\n", &runtime, ": Looking for CDSs, Transcripts, Pseudogenes with no Gene tag\n";
  my $db = Ace->connect(-path=>$db_path, -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

  my @genes= $db->fetch(-query=>'find worm_genes NOT Gene');
  if(@genes){
    foreach (@genes){
      print LOG "ERROR: $_ has no Gene tag, please add valid Gene ID from geneace\n";
      $builderrors++;
    }
  }
  $db->close;

  print LOG &runtime, ": Finished subroutine\n\n";

  return ($builderrors);
}



#################################################################################
# make_agp                                                                      #
#################################################################################

sub make_agp {
  $am_option = "-agp";

  # This is about checking the DNA and then making agp files.  This step
  # now gets run twice during build, the first time might have errors due to 
  # EMBL synchronisation problems but the second time should be run at the end
  # of the build and errors should have gone away

  # Were there errors with generating DNA or composition, i.e. B2:ERROR file?
  &usage(17) if(-e "$logdir/$flag{'B2:ERROR'}");
  
  # have you run GFFsplitter?
  &dump_GFFs unless ((-e "$logdir/$flag{'C2'}") && (-e "$logdir/$flag{'B3'}"));


  # run three agp related scripts
  my @scripts_to_run = ("check_DNA.pl", "make_agp_file.pl", "agp2dna.pl");

  foreach my $script (@scripts_to_run){
    my $command = "$scriptdir/$script";
    # run in test mode?
    if($quicktest){
      $command .= " --quicktest";
    }
    elsif($test){
      $command .= " --test";
    }
    &run_command("$command");
  }


  # make a B3 log file if this is first run of -agp (i.e. no B3 log file there)
  system("touch $logdir/$flag{'B3'}") unless ((-e "$logdir/$flag{'B3'}"));

  # check for errors in the agp file 
  local (*AGP);
  my $agp_errors = 0;
  
  foreach my $chrom (@chrom) {
    open (AGP, "<$basedir/autoace/yellow_brick_road/CHROMOSOME_${chrom}.agp_seq.log") or die "Couldn't open agp file : $!";
    while (<AGP>) {
      $agp_errors++ if (/ERROR/);
    }
    close(AGP);
  }
    
  # Errors?
  # make 'agp_files_errors' log file in /logs this will halt the process if you are running blat
  system("touch $logdir/$flag{'B3:ERROR'}") if ($agp_errors > 1);
        
}

#__ end make_agp __#


#################################################################################
# run dbcomp.pl                                                                 #
#################################################################################

sub dbcomp{
  $am_option = "-dbcomp";	

  # need to perform class by class comparison against previous release
  &run_command("$scriptdir/dbcomp.pl");	      
}
#__ end dbcomp __#


#################################################################################
# prepare for blat jobs                                                         #
#################################################################################

sub prepare_for_blat{
  $am_option = "-prepare_blat";
  &usage(15) if (-e "$logdir/$flag{'B3:ERROR'}");
  
  # transcriptmasker run to mask ?Feature_data from raw sequences
  # note to krb. This needs bradnamisation to allow a -all flag.

  &run_command("$scriptdir/transcriptmasker.pl -all");

  # Now make blat target database using autoace (will be needed for all possible blat jobs)
  # This also makes a backup copy of the old psl files (in case you need them to refer to)

  &run_command("$scriptdir/blat_them_all.pl -dump");


} 

#################################################################################

sub blat_jobs{

  $am_option = "-blat";
  # Should only be here if there are new sequences to blat with, or genome sequence has changed.

  my $blat_dir = "$basedir/autoace/BLAT";

  # what blat jobs should I run? Do everything if blat_all selected
  # nematode should always be last job to tackle

  my @blat_jobs;
  push(@blat_jobs,"est")      if ( ($blat_est)      || ($blat_all) );
  push(@blat_jobs,"ost")      if ( ($blat_ost)      || ($blat_all) );
  push(@blat_jobs,"mrna")     if ( ($blat_mrna)     || ($blat_all) );
  push(@blat_jobs,"ncrna")    if ( ($blat_ncrna)    || ($blat_all) );
  push(@blat_jobs,"embl")     if ( ($blat_embl)     || ($blat_all) );
  push(@blat_jobs,"tc1")      if ( ($blat_tc1)      || ($blat_all) );      
  push(@blat_jobs,"nematode") if ( ($blat_nematode) || ($blat_all) );

  my $status;
  my $nematode_flag = 0; # have nematode blats been run?

  # run each blat job in turn 
  foreach my $job(@blat_jobs){

    # If all other blat jobs are finished can load all non-nematode blat data into autoace this allows 
    # gff dump and gff split to be runin later steps, this is because nematode ESTs take a long time to blat/load
    if ($job eq "nematode"){
      &load_blat_results("est","mrna","ncrna","embl","ost","tc1");
      $nematode_flag = 1;
    }
    
    # run the main blat job
    &run_command("$scriptdir/blat_them_all.pl -blat -process -$job");
    
    # also create virtual objects
    &run_command("$scriptdir/blat_them_all.pl -virtual -$job");
    
    # Run aceprocess to make cleaner files
    # slight difference for nematode files as there is no best/other distinction or intron files
    if($job eq "nematode"){
      &run_command("$scriptdir/acecompress.pl -homol ${blat_dir}/autoace.$job.ace > ${blat_dir}/autoace.blat.${job}lite.ace");
      my $status = move("${blat_dir}/autoace.blat.${job}lite.ace", "${blat_dir}/autoace.blat.$job.ace");
      print "ERROR: Couldn't move file: $!\n" if ($status == 0);      
    }
    else{
      &run_command("$scriptdir/acecompress.pl -homol ${blat_dir}/autoace.blat.$job.ace > ${blat_dir}/autoace.blat.${job}lite.ace");
      my $status = move("${blat_dir}/autoace.blat.${job}lite.ace", "${blat_dir}/autoace.blat.$job.ace");
      print "ERROR: Couldn't move file: $!\n" if ($status == 0);
      &run_command("$scriptdir/acecompress.pl -feature ${blat_dir}/autoace.good_introns.$job.ace > ${blat_dir}/autoace.good_introns.${job}lite.ace");
      $status = move("${blat_dir}/autoace.good_introns.${job}lite.ace", "${blat_dir}/autoace.good_introns.$job.ace");
      print "ERROR: Couldn't move file: $!\n" if ($status == 0);

    }

    print LOG "Finishing acecompress.pl at ",&runtime,"\n\n";
      
    # make blat job specific lock file
    system("touch $logdir/$flag{'B6_${job}'}");
    
  }
  # generic lock file
  system("touch $logdir/$flag{'B6'}");  

  # now load blat results into autoace
  # if blat_nematode was selected then only need to load just those results as other results 
  # would have been loaded above. otherwise load everything

  if ($nematode_flag == 1) {
      &load_blat_results("nematode");
  }
  else {
      &load_blat_results("all");
  } 
} 

################################################################################

# load_blat_results
# generic subroutine for loading blat data into autoace
# will load all types of blat result if 'all' is passed to the subroutine

sub load_blat_results {

    $am_option .= "-addblat";
    
    my $first_blat_type = $_[0];
    my @blat_types = @_;
    @blat_types    = ("est","mrna","ncrna","ost","embl","nematode","tc1") if ($first_blat_type eq "all");
  
    foreach my $type (@blat_types){    
	print LOG "Adding BLAT $type data to autoace at ",&runtime,"\n";
	my $file =  "$basedir/autoace/BLAT/virtual_objects.autoace.blat.$type.ace";
	&load($file,"virtual_objects_$type");

	# Don't need to add confirmed introns from nematode data (because there are none!)
	unless ( ($type eq "nematode") || ($type eq "tc1") || ($type eq "embl")|| ($type eq "ncrna") ) {
	    $file = "$basedir/autoace/BLAT/virtual_objects.autoace.ci.$type.ace"; 
	    &load($file,"blat_confirmed_introns_$type");
	    
	    $file = "$basedir/autoace/BLAT/autoace.good_introns.$type.ace";
	    &load($file,"blat_good_introns_$type");
	} 
	
	$file = "$basedir/autoace/BLAT/autoace.blat.$type.ace";           
	&load($file,"blat_${type}_data");
	
    }
    system("touch $logdir/$flag{'B7'}"); 
}

#__ end load_blat_results __#

####################################################################################


##################################
# load BLASTX data
#################################
sub parse_homol_data {
  $am_option = "-addhomol";
  my $command;

  my @files2Load = (
		    #BLAST data
		    "worm_pep_blastp.ace",
		    "worm_brigpep_blastp.ace",
		    "worm_dna_blastx.ace",
		    #motif info
		    "worm_pep_motif_info.ace",
		    "worm_brigpep_motif_info.ace",
		    #protein info
		    "ipi_hits.ace", 
		    "flybase.ace",
		    "yeast.ace",
		    "swissproteins.ace",
		    "tremblproteins.ace",
		    "brigpep.ace",
		    #other data
		    "repeat_homologies.ace",
		    "waba.ace",
		    "TRF.ace",
		   );

  foreach my $file ( @files2Load ) {

    my $newfile = "$basedir/wormbase/ensembl_dumps/$file";
    &load($newfile,$tsuser);
  }
 
  # upload_homol data log file in /logs
  system("touch $logdir/$flag{'B8'}");

}
#__ end parse_homol_data __#




##################################
# load briggsae data
#################################

sub parse_briggsae_data {
  $am_option = "-addbriggsae";

  my $file;

  # load raw briggsae data files
  my @Brigdata = ("DNA", "fosmid", "agplink", "sequence");
  foreach my $Brigdata (@Brigdata) {
    $file = "$basedir/wormbase/briggsae/briggsae_cb25.agp8_${Brigdata}.ace";
    &load($file,"briggsae_${Brigdata}");
  }
		   
  # load gene prediction sets
  my @genepredicters =("genefinder", "fgenesh", "twinscan", "ensembl", "hybrid", "rna");
  foreach my $genepredicters (@genepredicters) {
    $file = "$basedir/wormbase/briggsae/briggsae_cb25.agp8_${genepredicters}.ace"; 
    &load($file,"briggsae_${genepredicters}_genes");
  }

  # briggsae BAC end data
  $file = "/wormsrv1/briggsae/BAC_ENDS/briggsae_BAC_ends.fasta";
  &load($file,"briggsae_BAC");
                           
  $file = "/wormsrv1/briggsae/BAC_ENDS/briggsae_homol_data.ace";
  &load($file,"briggsae_BAC");
                           
  $file = "/wormsrv1/briggsae/BAC_ENDS/briggsae_BAC_ends_data.ace";
  &load($file,"briggsae_BAC");
                           
  $file = "/wormsrv1/briggsae/BAC_ENDS/bac_ends_unique.ace";
  &load($file,"briggsae_BAC");

  $file = "/wormsrv1/briggsae/BAC_ENDS/briggsae_bac_clone_ends.ace";
  &load($file,"briggsae_BAC");

  $file = "/wormsrv2/wormbase/misc_static/ortholog_WS131.ace";
  &load($file,"briggsae_orthologues");

  # upload_homol data log file in /logs
  system("touch $logdir/$flag{'B9'}");
  
}
#__ end parse_briggsae_data __#

#################################################################################
# generate_utrs                                                                 #
#################################################################################

sub generate_utrs {
  $am_option = "-utrs";
  
  # run find_utrs.pl to generate data
  my $command = "$scriptdir/find_utrs.pl -database autoace -output_dir $basedir/autoace/UTR";
  ($command .= " -test") if ($test);
  &run_command($command);

  # load data
  my $file = "$basedir/autoace/UTR/UTRs.ace"; 
  &load($file,"utrs");

  # make a log file in /autoace/logs
  system("touch $logdir/B10:Generate_UTR_data");
}
#__ end generate_utrs __#



#################################################################################################################


#################################################################################
# make_operons                                                                  #
#################################################################################

sub make_operons {
  $am_option = "-operon";
  
  # needs split gff files for this step, so should only run if UTRs have been
  # generated (which will make split GFF files
  # &usage(20) unless (-e "$logdir/$flag{'B10'}");

  # run find_utrs.pl to generate data
  &run_command("$scriptdir/map_operons.pl");

  # make a log file in /autoace/logs
  system("touch $logdir/B11:Generate_operon_data");
}
#__ end generate_utrs __#

##########################################################################################################

sub make_wormpep {
  $am_option = "-buildpep";
  # make wormpep database but also perform all the other related protein steps in the build
  unless( -e "$logdir/D1A:Build_wormpep_initial" ) {
    # make wormpep -initial
    my $command = "$scriptdir/make_wormpep.pl -initial";
    $command .= " -test" if ($test);
    &run_command("$command");
   
    #generate file to ad new peptides to mySQL database.
    $command = "$scriptdir/new_wormpep_entries.pl";
    $command .= " -test" if ($test);
    &run_command("$command");

    # update common data
    $command = "$scriptdir/update_Common_data.pl --build --cds2wormpep";
    $command .= " -test" if ($test);
    &run_command("$command");

    
    system("touch $logdir/D1A:Build_wormpep_initial");
  }
  else {
    # get Pfam domains (this step loads resulting ace file to autoace)
    &run_command("$scriptdir/GetPFAM_motifs.pl -load");

    # get interpro domains (this step loads resulting ace file to autoace)
    &run_command("$scriptdir/GetInterPro_motifs.pl -load");

    # make interpro2go connections (to be used by getProteinID)
    &run_command("$scriptdir/make_Interpro2GO_mapping.pl");

    # Get protein IDs (this step writes to ~wormpub/analysis/SWALL and loads wormpep info to autoace) 
    &run_command("$scriptdir/getProteinID.pl -load");
    
    # make wormpep -final
    my $command = "$scriptdir/make_wormpep.pl -final";
    $command .= " -test" if ($test);
    &run_command("$command");
   
    #
    # make acefile of peptides etc to add to autoace (replacement for pepace)
    &run_command("$scriptdir/build_pepace.pl");

    &run_command("$scriptdir/update_Common_data.pl -build -cds2pid");

    
    # make a make_autoace log file in /logs
    system("touch $logdir/D1:Build_wormpep_final");
  }
}
#__ end make_wormpep  __#

#################################################################################
# make_wormrna                                                                  #
#################################################################################

sub make_wormrna {
  $am_option = "-buidrna";
  &run_command("$scriptdir/make_wormrna.pl -release $WS_version");
}
#__ end make_wormrna __#


#################################################################################
# dump GFF files                                                                #
#################################################################################

sub dump_GFFs {
  $am_option .= " -gffdump";

  my $command = "$scriptdir/chromosome_dump.pl --gff";
  
  # run in test mode?
  if($quicktest){
    $command .= " --quicktest";
  }
  elsif($test){
    $command .= " --test";
  }

  &run_command($command);

  &usage(18) if(-e "$logdir/$flag{'C1:ERROR'}");

  # make dumped_GFF_file in /logs
  system("touch $logdir/C1:Dumped_GFF_files");
}
#__ end dump_GFFs __#



####################################################################################

sub split_GFFs {
  $am_option .= " -gffsplit";
  &run_command("$scriptdir/GFFsplitter.pl");

  # make GFF_splitter file in /logs
  system("touch $logdir/C2:Split_GFF_files");
}
#__ end split_GFFs __#


#################################################################################
# map things to the genome                                                      #
#################################################################################

# set of calls to map a variety of things

sub map_features {
    
  $am_option = "-map";

  my $file;

  # features
  &run_command("$scriptdir/map_features.pl -all");

  # these should be folded into a loop or into the script itself [dl]
  #Loads feature data into autoace with tsuser set to feature type.
  my @features = ("SL1", "SL2", "polyA_site", "polyA_signal",);
  foreach my $feature (@features) {
    $file = "$basedir/autoace/FEATURES/WS${WS_version}_feature_${feature}.ace";
    &load($file,"feature_${feature}");
  }
  
  # PCR products
  &run_command("$scriptdir/map_PCR_products.pl");
  
  #Oligo_sets
  &run_command ("$scriptdir/map_Oligo_set.pl");
  
  # RNAi experiments
  &run_command("$scriptdir/map_RNAi.pl");

  # alleles
  &run_command("$scriptdir/map_Alleles.pl");

  # Y2H objects
  &run_command("$scriptdir/map_Y2H.pl");
  $file = "$basedir/wormbase/misc_dynamic/misc_Y2H_connections.ace";
  &load($file,"Y2H_connections");

  # microarray connections
  &run_command("$scriptdir/map_microarray.pl");

  $file = "$basedir/wormbase/misc_dynamic/misc_microarrays.ace";
  &load($file,"microarray_connections");

}
#__ end map_features __#


#################################################################################
# confirm_gene models                                                           #
#################################################################################

sub confirm_gene_models {
  $am_option = "-confirm";

  # GFF files need to be dumped again (will now have blat data)
  # also needs GFF files to be split again
  #&dump_GFFs;
  #&split_GFFs;

  # confirm_genes from EST&OST (-est) and mRNA (-mrna) data sets
  &run_command("$scriptdir/confirm_genes.pl --est --mrna");

  # load files
  my $file = "$basedir/wormbase/misc_dynamic/misc_confirmed_by_EST.ace";
  &load($file,"genes_confirmed_by_EST_and_OST");

  $file = "$basedir/wormbase/misc_dynamic/misc_confirmed_by_mRNA.ace";
  &load($file,"genes_confirmed_by_mrna");
  
  # make dumped_GFF_file in /logs
  system("touch $logdir/C4:Confirm_gene_models");

  &run_command("$scriptdir/update_Common_data.pl -build -CDS_list");

}
#__ end confirm_gene_models __#


#################################################################################
# load  - generic subroutine for loading data into autoace                      #
#################################################################################

sub load {

  my $file = shift;

  if ($load){
    my $temp = $file;
    # remove trailing path of filename
    $temp =~ s/.*\///;
    $am_option = "-load $temp" ;
  }
  
  # tsuser is optional but if set, should replace any dots with underscores just in case 
  my $tsuser = shift; 
  $tsuser =~ s/\./_/g;

  my $flag = 0;
  unless (-e "$file"){
    $errors++;
    print LOG "ERROR: Couldn't find file named: $file\n";    
    # set flag to skip loading step
    $flag = 1;
  }

  unless($flag == 1){
    my $command = "pparse $file\nsave\nquit\n";
    
    print LOG &runtime,": adding $file info to autoace\n\n";  
    open (WRITEDB, "| $tace -tsuser $tsuser $basedir/autoace ") || die "Couldn't open pipe to autoace\n";
    print WRITEDB $command;
    close (WRITEDB);
  }

}

#################################################################################
# Open logfile                                                                  #
#################################################################################

sub logfile_details {
  $rundate = &rundate;
  print LOG "# autoace_minder.pl started at: $rundate ",&runtime,"\n";
  print LOG "# TEST MODE!!! Using $basedir/autoace\n" if ($test);
  print LOG "# WormBase/Wormpep version: WS${WS_version}\n\n";  
  print LOG "#  -initial      : Prepare for a new build, update WSnn version number\n"                 if ($initial);
  print LOG "#  -unpack       : Unpack databases from FTP site and copy Sanger dbs\n"                  if ($unpack);
  print LOG "#  -acefile      : Write .acefiles from WormBase copies of the databases\n"               if ($acefile);
  print LOG "#  -build        : Build autoace\n"                                                       if ($build);
  print LOG "#  -builddb      : Build autoace : DB only\n"                                             if ($builddb);
  print LOG "#  -buildchrom   : Build autoace : DNA data\n"                                            if ($buildchrom);

  print LOG "#  -buildrelease : Build autoace : Release directory\n"                                   if ($buildrelease);
  print LOG "#  -agp          : Make and check agp files\n"		                               if ($agp);
  print LOG "#  -dbcomp       : Check DB consistency and diffs from previous version\n"                if ($dbcomp);
  print LOG "#  -buildpep     : Build wormpep database\n"                                              if ($buildpep);
  print LOG "#  -buildrna     : Build wormrna database\n"                                              if ($buildrna);
  print LOG "#  -blat_est     : perform blat analysis on ESTs\n"                                       if ($blat_est);
  print LOG "#  -blat_ost     : perform blat analysis on OSTs\n"                                       if ($blat_ost);
  print LOG "#  -blat_mrna    : perform blat analysis on mRNAs\n"                                      if ($blat_mrna);
  print LOG "#  -blat_ncrna   : perform blat analysis on ncRNAs\n"                                     if ($blat_ncrna);
  print LOG "#  -blat_embl    : perform blat analysis on non-WormBase CDSs from EMBL\n"                if ($blat_embl);
  print LOG "#  -blat_tc1     : perform blat analysis on TC1 insertions\n"                             if ($blat_tc1);
  print LOG "#  -blat_nematode: perform blat analysis on non-C. elegans ESTs\n"                        if ($blat_nematode);
  print LOG "#  -blat_all     : perform blat analysis on everything\n"                                 if ($blat_all);
  print LOG "#  -addblat      : Load blat data into autoace\n"                                         if ($addblat);
  print LOG "#  -addhomol     : Load blast data into autoace\n"                                        if ($addhomol);
  print LOG "#  -addbriggsae  : Load briggsae data into autoace\n"                                     if ($addbriggsae);
  print LOG "#  -utrs         : Generates and load UTR data into autoace\n"                            if ($utrs);
  print LOG "#  -operon       : Generates operon dataset\n"                                            if ($operon);
  print LOG "#  -debug        : Debug mode\n"                                                          if ($debug);
  print LOG "#  -verbose      : Verbose mode\n"                                                        if ($verbose);
  print LOG "#  -gffdump      : Dump GFF files\n"                                                      if ($gffdump);
  print LOG "#  -gffsplit     : Split GFF files\n"                                                     if ($gffsplit);
  print LOG "#  -map          : map PCR and RNAi\n"                                                    if ($map);
  print LOG "#  -test         : running in test mode\n"                                                if ($test);
  print LOG "======================================================================\n\n";

}

##################################################################################
#
# Simple routine which will run commands via system calls but also check the 
# return status of a system call and complain if non-zero, increments error check 
# count, and prints a log file error
#
##################################################################################

sub run_command{
  my $command = shift;
  print LOG &runtime, ": Running $command\n";
  my $status = system($command);
  if(($status >>8) != 0){
    $errors++;
    print LOG "ERROR: $command failed. \$\? = $status\n";
  }
  print LOG &runtime, ": Finished.\n\n";

  # for optional further testing by calling subroutine
  return($status);
}


#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub usage {
    my $error = shift;

    if ($error eq "No_WormBase_release_number") {
	# No WormBase release number file
	print "The WormBase release number cannot be parsed\n";
	print "Check File: $basedir/autoace/wspec/database.wrm\n\n";
	exit(0);
    }
    elsif ($error == 1) {
	# Build_in_progress flag already exists when build initiation attempted
	print "\nautoace build aborted:\n";
	print "The 'build_in_progress' flag already exists\n";
	print "You cannot start a new build until this flag is removed\n\n";
	exit(0);
    }
    elsif ($error eq "Failed_to_update_cvs_repository") {
	# Failed to update cvs repository with new WormBase release number
	print "The 'database.wrm' file in the cvs repository was not updated\n";
	print "Try to do this by hand.\n\n";
	exit(0);
    }
    elsif ($error == 2) {
	# Abort unpack_db.pl script
	print "Abort unpack_db run.\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error eq "Primary_databases_file_error") {
	# No Primary_databases_used_in_build file
	print "The 'Primary_databases_used_in_build' is absent or unreadable.\nAbort build.\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error eq "Absent_stlace_database") {
	# No last_version date for stlace
	print "Abort unpack_db run. stlace\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error eq "Absent_brigdb_database") {
	# No last_version date for brigdb
	print "Abort unpack_db run. brigdb\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error eq "Absent_citace_database") {
	# No last_version date for citace
	print "Abort unpack_db run. citace\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error eq "Absent_cshace_database") {
	# No last_version date for cshace
	print "Abort unpack_db run. cshace\n";
	print "\n\n";
	exit(0);
    }
    elsif ($error == 7) {
	# Check build prior to building
	print "You haven't built the database yet!\n";
	exit(0);
    }
    elsif ($error == 8) {
	# Check acefile dump prior to building
	print "You haven't written up-to-date .acefiles yet!\n";
	exit(0);
    }
    elsif ($error eq "Errors_in_loaded_acefiles") {
	# Errors in loaded acefiles, check log file
	print "There were errors in the loaded acefiles.\n";
	print "Check the logfile and correct before continuing.\n\n";
	exit(0);
    }
    elsif ($error == 9) {
	# Check blat run
	print "You haven't run blat_them_all.\n";
	exit(0);
    }
    elsif ($error == 11) {
	# Build_in_progress flag is newer than WormBase version number
	print "\nautoace build aborted:\n";
	print "The WS version predates the 'build_in_progress' flag \n";
	print "Check that the WS version number is correct\n\n";
	exit(0);
    }
    elsif ($error eq "Build_in_progress_absent") {
	# Build_in_progress flag is absent when preparing primaries
	print "\nautoace build aborted:\n";
	print "The 'build_in_progress' flag is absent. \n";
	print "You can't overwrite the primary databases on wormsrv2 prior to a new build.\n\n";
	exit(0);
    }
    elsif ($error == 13) {
	# B2:Build_autoace_database absent when making release .tar.gz files
	print "\nautoace build aborted:\n";
	print "The 'Build_autoace_database' flag is absent. \n";
	print "You can't write the release .tar.gz files prior to a completing the  build.\n\n";
	exit(0);
    }
    elsif ($error == 14) {
	# Dump_dna before the database is finished
	print "\nautoace build aborted:\n";
	print "The '$flag{'B1'} flag is absent. \n";
	print "You can't write the chromosome DNA files prior to a completing the build.\n\n";
	exit(0);
    }
    elsif ($error == 15) {
	# attempted BLAT analysis with agp errors 
	print "\nautoace build aborted:\n";
	print "The '$flag{'B3:ERROR'} flag is set indicating errors in the agp file. \n";
	print "You must remove this file before you can run the BLAT analysis.\n\n";
	exit(0);
    }
    elsif ($error == 16) {
      # DEPRECATED - dont copt to midway anymore

    }
    elsif ($error == 17) {
	# atempted agp creation without clearing DNA/composition error file
	print "\nautoace build aborted:\n";
	print "The '$flag{'B2:ERROR'} flag is present indicating that the chromosome DNA and/or composition\n";
	print "files were not generated properly.  You cannot continue until this lock file is removed.\n";
	exit(0);
    }
    elsif ($error == 18) {
	# gff dumping failed
	print "\nautoace build aborted:\n";
	print "The '$flag{'C2:ERROR'} flag is present indicating that the chromosome GFF dumps did not work.\n";
	print "You cannot continue until this lock file is removed.\n";
	exit(0);
    }
    elsif ($error == 19) {
      # tried loading blat results without first runnning blat
      print "\nautoace build aborted:\n";
      print "The '$flag{'B6'} flag is not present indicating that BLAT has not been run. \n";
      print "If you didn't need to run BLAT then you must create this file before you can upload BLAT results.\n\n";
      exit(0);
    }
    elsif ($error == 20) {
      # tried loading blat results without first runnning blat
      print "\nautoace build aborted:\n";
      print "The '$flag{'B10'} flag is not present indicating that UTRs have not been made, hence GFF files. \n";
      print "have not been split.  You need split GFF files to make operons.\n\n";
      exit(0);
    }
    elsif ($error == 21) {
      # -test and -quicktest specified, only need one
      print "\nautoace build aborted:\n";
      print "You have specified -test (full-test mode) AND -quicktest (only test one chromosome where\n";
      print "appropriate).  You only need to specify one of these options.\n\n";
      exit(0);
    }

    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }


}


__END__

=pod

=head2   NAME - autoace_minder.pl

=head1 USAGE

=over 4

=item autoace_minder.pl [-options]

=back

autoace_minder.pl is a wrapper to drive the various scripts utilised in the
build of a C.elegans WS database release.

autoace_minder.pl mandatory arguments:

=over 4

=item none, (but it won\'t do anything)

=back

autoace_minder.pl OPTIONAL arguments:

=over 4

=item -acefile, Write .acefiles from WormBase copies of the databases

=item -build, Build autoace : (performs options 1 & 2 below)

=item -builddb, Build autoace : Database only

=item -buildchrom, Build autoace : Dump DNA

=item -buildrelease, Build autoace : release directory

=item -dbcomp, Check DB consistency and diffs from previous version

=item -buildpep, Build wormpep database

=item -ftp, Move WS release to the external FTP site (Full public release)

=item -debug <user>, Debug mode - log file only goes to named user

=item -verbose, Verbose mode toggle on extra command line output

=item -agp, creates and checks agp files

=item -gffdump, dump GFF files

=item -gffsplit, split GFF files

=item -map, map PCR and RNAi

=item -blat_est, map all blat EST similarities, load into autoace

=item -blat_ost, map all blat OST similarities, load into autoace

=item -blat_mrna, map all blat similarities, load into autoace

=item -blat_embl, map all blat EMBL gene similarities, load into autoace

=item -blat_nematode, map all blat other nematode ESTs similarities, load into autoace

=item -blat_all, map all blat similarities (BLAT jobs)

=item -addblat, parse all BLAT files (ESTs, OSTs, mRNAs, EMBL genes, nematode ESTs) into autoace

=item -addhomol, parse BLASTX and BLASTP data from pre-build

=item -addbriggsae, parse briggsae assembly data and briggsae gene predictions

=item -utrs, generate UTR datatset and add to autoace

=item -operon, make operons
 
=item -load <file>, load filename into autoace...filename should be valid acefile

=item -tsuser <name>, to be used in conjunction with -load (optional)

=item -test, set build to be in test mode, uses ~wormpub/AUTOACE_TEST and ~wormpub/scripts

=item -quicktest, same as -test but where appropriate will only run tests/checks on one chromosome
(chromosome III), this is meaningless for steps of the build which don't loop through chromosomes,
and for these steps -quicktest is the same as -test


=back

=cut

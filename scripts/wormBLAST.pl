#!/usr/local/bin/perl5.6.1 -w
#
# wormBLAST.pl
# 
# written by Anthony Rogers
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2004-02-19 10:01:57 $

use DBI;
use strict;
my $wormpipe_dir = glob("~wormpipe");
use lib "/nfs/acari/wormpipe/scripts";
use lib "/nfs/acari/wormpipe/scripts/BLAST_scripts";
use Wormbase;
use Getopt::Long;

#######################################
# command-line options                #
#######################################

my $chromosomes;
my $wormpep;
my $update_databases;
my $update_mySQL;
my $setup_mySQL;
my $run_pipeline;
my $run_brig;
my $dont_SQL;
my $dump_data;
my $mail;
my $test_pipeline;
my $WPver;   #  Wormpep version is passed as command line option
my $blastx;
my $blastp;
my $prep_dump;
my $cleanup;
my $errors = 0; # for tracking global error - needs to be initialised to 0
my $debug;
my $test;
my $log;
my $WS_version;
my $maintainers = "All"; # email recipients

GetOptions("chromosomes" => \$chromosomes,
	   "wormpep"     => \$wormpep,
	   "databases"   => \$update_databases,
	   "updatemysql" => \$update_mySQL,
	   "setup"       => \$setup_mySQL,
	   "run"         => \$run_pipeline,
	   "run_brig"    => \$run_brig,
	   "nosql"       => \$dont_SQL,
	   "prep_dump"   => \$prep_dump,
	   "dump"        => \$dump_data,
	   "mail"        => \$mail,
	   "testpipe"    => \$test_pipeline,
	   "version=s"   => \$WS_version,
	   "blastp"      => \$blastp,
	   "blastx"      => \$blastx,
	   "cleanup"     => \$cleanup,
	   "debug=s"     => \$debug,
	   "test"        => \$test,
	  );


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


&create_log_files;


# you can do either or both blast anaylses
if( $run_pipeline ) {
  $blastx = 1; $blastp =1;
}

&check_wormsrv2_conflicts;

die "please give a build version number ie  wormBLAST -version 114\n" unless $WS_version;
my $WS_old = $WS_version - 1;
my $scripts_dir = "$wormpipe_dir/scripts/BLAST_scripts";
#process Ids

#|         18 | gadfly3.pep         |
#|         19 | ensembl7.29a.2.pep  |
#|         20 | yeast2.pep          |
#|         23 | wormpep87.pep       |
#|         24 | slimswissprot40.pep |
#|         25 | slimtrembl21.pep    |

my %worm_dna_processIDs = ( 
			wormpep       => 2,
			brigpep       => 3, 
			ipi_human     => 4,
			yeast         => 5,
			gadfly        => 6,
			slimswissprot => 7,
			slimtrembl_1  => 8,
			slimtrembl_2  => 9,
		       );


my %wormprotprocessIds = ( 
			  wormpep       => 2,
			  brigpep       => 3, 
			  ipi_human     => 4,
			  yeast         => 5,
			  gadfly        => 6,
			  slimswissprot => 7,
			  slimtrembl_1  => 8,
			  slimtrembl_2  => 9,
			 );

#get new chromosomes
&run_command("perl5.6.1 $scripts_dir/copy_files_to_acari.pl -c") if ($chromosomes);


#get new wormpep
&run_command("perl5.6.1 $scripts_dir/copy_files_to_acari.pl -w") if ($wormpep);


my %currentDBs;   #ALSO used in setup_mySQL 
my @updated_DBs;  #used by get_updated_database_list sub - when run this array is filled with databases that have been updated since the prev build
#load in databases used in previous build
my $last_build_DBs = "$wormpipe_dir/BlastDB/databases_used_WS$WS_old";
my $database_to_use = "$wormpipe_dir/BlastDB/databases_used_WS$WS_version";


open (OLD_DB,"<$last_build_DBs") or die "cant find $last_build_DBs";
while (<OLD_DB>) {
  chomp;
  if( /(gadfly|yeast|slimswissprot|slimtrembl_1|slimtrembl_2|wormpep|ipi_human|brigpep)/ ) {
    $currentDBs{$1} = $_;
  }
}
close OLD_DB;

#check for updated Databases
if ( $update_databases ){
  print "Updating databases \n";
  open (DIR,"ls -l $wormpipe_dir/BlastDB/*.pep |") or die "readir\n";
  while (<DIR>) { 
    #  print;
    chomp;
    if( /\/(gadfly|yeast|slimswissprot|slimtrembl_1|slimtrembl_2|wormpep|ipi_human|brigpep)/ ) {
      my $whole_file = "$1"."$'";  #match + stuff after match.
      if( $1 eq "wormpep" ) {
	print "updating wormpep to version $WS_version anyway - make sure the data is there !\nCopying over will take care of setdb 'ing it\n";
	$whole_file = "wormpep".$WS_version.".pep";
	$currentDBs{$1} = "$whole_file";
	next;
      }
      if( "$whole_file" ne "$currentDBs{$1}" ) {
	#make blastable database
	print "\tmaking blastable database for $1\n";
	&run_command("setdb $wormpipe_dir/BlastDB/$whole_file");
	push( @updated_DBs,$1 );
	#change hash entry ready to rewrite external_dbs
	$currentDBs{$1} = "$whole_file";
      }
      else {
	print "\t$1 database unchanged $whole_file\n";
      }
    }
  }
  close DIR;
  
  open (NEW_DB,">$database_to_use") or die "cant write updated $database_to_use";
  foreach (keys %currentDBs){
    print NEW_DB "$currentDBs{$_}\n";
  }
  close NEW_DB;
}



#@updated_DBs = qw(ensembl gadfly yeast slimswissprot slimtrembl);
if( $mail )
  {
    &get_updated_database_list;
    
    #generate distribution request based on updated databases
    my $letter = "$wormpipe_dir/distribute_on_farm_mail";
    open (LETTER,">$letter") or die "cant open $letter to write new distribution request\n";
    print LETTER "This is a script generated email from the wormpipe Blast analysis pipeline.\nAny problems should be addessed to worm\@sanger.ac.uk.\n
=====================================================
\n";
    print LETTER "The following can be removed from /data/blastdb/Worms.\n\n";
    foreach (@updated_DBs){
      print LETTER "$_*\n";
    }
    #print LETTER "wormpep*\n";
    print LETTER "CHROMOSOME_*.dna\n";
    
    print LETTER "-------------------------------------------------------\n\nand replaced with the following files from ~wormpipe/BlastDB/\n\n";
    foreach (@updated_DBs){
      print LETTER "$currentDBs{$_}\n$currentDBs{$_}.ahd\n$currentDBs{$_}.atb\n$currentDBs{$_}.bsq\n\n";
    }
    
    #print LETTER "wormpep$WS_version.pep\nwormpep$WS_version.pep.ahd\nwormpep$WS_version.pep.atb\nwormpep$WS_version.pep.bsq\n\n";
    print LETTER "CHROMOSOME_I.dna\nCHROMOSOME_II.dna\nCHROMOSOME_III.dna\nCHROMOSOME_IV.dna\nCHROMOSOME_V.dna\nCHROMOSOME_X.dna\n\n";
    
    print LETTER "-------------------------------------------------------\n\n";
    print LETTER "Thanks\n\n ________ END ________";
    close LETTER;
  
    my $name = "Wormpipe database distribution request";
    my $maintainer = "ssg-isg\@sanger.ac.uk, sanger\@wormbase.org";     # mail to ssg and all Sanger WormBase
    print "mailing distibution request to $maintainer\n";
    &mail_maintainer($name,$maintainer,$letter);
  }



# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "worm_dna";
my $dbpass = "worms";

my @results;
my $query = "";
my $worm_dna;     #worm_dna Db handle
my $worm_pep;   #wormprot Db handle

# update mySQL database
if( $update_mySQL )
  {
    print "Updating mysql databases\n";
    #make worm_dna connection
    $worm_dna = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #internal_id number of the last clone in the worm_dna
    #in case this routine is being run twice the clone id is stored and checked
    my $last_clone;
    my $update;
    open (NEW_DB,"<$database_to_use") or die "cant read updated $database_to_use during update_mySQL";
    
    #this bits logic seems wrong - but it works so I'll fix it later - ar2
    while(<NEW_DB>){
      if( /last_clone\s(\w+)/ ){
	$last_clone = $1;
	last;
      }
      else {
	$query = "select * from clone order by clone_id desc limit 1";
	@results = &single_line_query( $query, $worm_dna );
	$last_clone = $results[0];
	$update = 1;
      }
    }
    close NEW_DB;
    
    if( $update ){   # first time writing for this build
      open (NEW_DB,">>$database_to_use") or die "cant read updated $database_to_use during update_mySQL - last_clone";
      print NEW_DB "last_clone $last_clone\n";
    }
    print "last_clone = $last_clone\n";
    
    #Make a concatenation of all six agp files from the last release to ~/Elegans  e.g.
    print "\tconcatenating agp files\n";
    &run_command("cat /wormsrv2/autoace/CHROMOSOMES/*.agp > $wormpipe_dir/Elegans/WS$WS_version.agp");
    
    #load information about any new clones
    print "\tloading information about any new clones in to $dbname\n";
    &run_command("perl5.6.1 $scripts_dir/agp2ensembl.pl -dbname worm_dna -dbhost ecs1f -dbuser wormadmin -dbpass worms -agp $wormpipe_dir/Elegans/WS$WS_version.agp -write -v -strict");
    
    #check that the number of clones in the clone table equals the number of contigs and dna objects
    my ($clone_count, $contig_count, $dna_count);
    $query = "select count(*) from clone";
    @results = &single_line_query( $query, $worm_dna );
    $clone_count = $results[0];
    
    $query = "select count(*) from contig";
    @results = &single_line_query( $query, $worm_dna );
    $contig_count = $results[0];
    
    $query = "select count(*) from dna";
    @results = &single_line_query( $query, $worm_dna );
    $dna_count = $results[0];
    
    print "checking clone contig and dna counts . . .";
    if( ($clone_count != $contig_count) or ($contig_count != $dna_count ) ){
      print "\nthe number of clones, contigs and DNAs is inconsistant\n
clones = $clone_count\ncontigs = $contig_count\ndna = $dna_count\n";
      exit(0);
    }
    else {
      print "OK\n";
    }
    
    $query = "select * from clone order by clone_id desc limit 1";
    @results = &single_line_query( $query, $worm_dna );
    my $new_last_clone = $results[0];
    print "\tnew_last_clone = $new_last_clone\n";
    
    if( $last_clone != $new_last_clone )
      {
	$query = "select contig_id from contig where contig_id > $last_clone into outfile '$wormpipe_dir/Elegans/ids.txt'";
	print &update_database( $query, $worm_dna );
	
	&run_command("perl5.6.1 $scripts_dir/InputIdManager.pl -dbname worm_dna -dbhost ecs1f -dbuser wormadmin -dbpass worms -insert -analysis SubmitContig -class contig -file $wormpipe_dir/Elegans/ids.txt");
	
      }
    $worm_dna->disconnect;

    print "\tchecking for duplicate clones\n";
    &run_command("perl5.6.1 $scripts_dir/find_duplicate_clones.pl");


    #add new peptides to MySQL database
    print "\n\nAdding new peptides to $dbname\n";
    #make protein database connection
    $dbname = "worm_pep";
    $worm_pep = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

    $query = "select * from protein order by proteinId desc limit 1";
    @results = &single_line_query( $query, $worm_pep );
    my $old_topCE = $results[0];
    if (-e "/wormsrv2/WORMPEP/wormpep$WS_version/new_entries.WS$WS_version"){
      &run_command("perl5.6.1 $scripts_dir/worm_pipeline.pl -fasta /wormsrv2/WORMPEP/wormpep$WS_version/new_entries.WS$WS_version");
    }
    else {
      die "new_entries.WS$WS_version does not exist! \nThis should have been made in autoace_minder -buildpep\n";
    }

    #check for updated ids
    @results = &single_line_query( $query, $worm_pep );
    my $new_topCE = $results[0];
    if( "$old_topCE" eq "$new_topCE" ) {
      print "\tNO new peptides were added to the $dbname mysql database\n";
    }
    else {
      print "\tnew highest proteinId is $new_topCE (old was $old_topCE )\n";
    }

    $worm_pep->disconnect;
  }



$dbuser = "wormadmin";
$dbpass = "worms";
if( $setup_mySQL )
  {  
    print "Setting up mysql ready for Blast run\n";
    #make wormprot connection
    $dbname = "worm_pep";
    my $worm_pep =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

    #make worm_brigprot connection
    $dbname = "worm_brigpep";
    my $worm_brigpep =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

    #make worm_dna connection
    $dbname = "worm_dna";
    $worm_dna = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

    &get_updated_database_list;

    # if the user passes WPversion greater than that in the current file update it anyway
    # (this means you can update the database before wormpep is made - ie during autoace_minder -build
    $currentDBs{$1} =~ /wormpep(\d+)/;
    if( $1 and ( $1<$WS_version ) ) {
      push (@updated_DBs,"wormpep$WS_version\.pep");
    }

    #update mysql with which databases need to be run against
    foreach my $database (@updated_DBs)
      {
	my $analysis = $worm_dna_processIDs{$database};
	my $db_file = $currentDBs{$database};
	print "________________________________________________________________________________\n";
	print "doing worm_dna updates . . . \n";
	$query = "update analysis set db = \"$db_file\" where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_dna );
	
	$query = "update analysis set db_file = \"/data/blastdb/Worms/$db_file\" where analysis_id = $analysis";
	print $query,"\n\n";
	&update_database( $query, $worm_dna );
	
	#delete entries so they get rerun
	$query = "delete from input_id_analysis where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_dna );
	
	$query = "delete from dna_align_feature where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_dna );
	
	
	
	print "doing worm_pep updates . . . \n";
	$analysis = $wormprotprocessIds{$database};
	$query = "update analysis set db = \"$db_file\" where analysis_id = $analysis";
	print $query,"\n";	
	&update_database( $query, $worm_pep );
	&update_database( $query, $worm_brigpep );

	$query = "update analysis set db_file = \"/data/blastdb/Worms/$db_file\" where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_pep );
	&update_database( $query, $worm_brigpep );
	
	#delete entries so they get rerun
	$query = "delete from input_id_analysis where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_pep );
	&update_database( $query, $worm_brigpep );
	
	$query = "delete from protein_feature where analysis_id = $analysis";
	print $query,"\n";
	&update_database( $query, $worm_pep );
	&update_database( $query, $worm_brigpep );
      }

    $worm_dna->disconnect;
    $worm_pep->disconnect;
    $worm_brigpep->disconnect;

  }
my $bdir = "/nfs/farm/Worms/Ensembl/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline";

if( $blastx ) {
  die "can't run pipeline whilst wormsrv2 is mounted - please exit and try again\n" if (-e "/wormsrv2");

  &run_RuleManager('worm_dna','dna');
}

if( $blastp ){

  die "can't run pipeline whilst wormsrv2 is mounted - please exit and try again\n" if (-e "/wormsrv2");

  &run_RuleManager('worm_pep','pep');
  &run_RuleManager('worm_brigpep','pep');

}



if( $prep_dump ) {
  # prepare helper files
  
  if( -e "/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_X.gff") {
    &run_command("cat /wormsrv2/autoace/CHROMOSOMES/*.gff | $scripts_dir/gff2cds.pl > /nfs/acari/wormpipe/Elegans/cds$WS_version.gff");
    &run_command("cat /wormsrv2/autoace/CHROMOSOMES/*.gff | $scripts_dir/gff2cos.pl > /nfs/acari/wormpipe/Elegans/cos$WS_version.gff");
    &run_command("$scripts_dir/prepare_dump_blastx.pl > $wormpipe_dir/dumps/accession2clone.list");
    &run_command("cp /wormsrv2/WORMPEP/wormpep$WS_version/wormpep.diff$WS_version $wormpipe_dir/dumps/");
    &run_command("cp /wormsrv2/WORMPEP/wormpep$WS_version/new_entries.WS$WS_version $wormpipe_dir/dumps/");
    &run_command("cp /wormsrv2/autoace/COMMON_DATA/CE2gene.dat $wormpipe_dir/dumps/");
    &run_command("cp /wormsrv2/autoace/COMMON_DATA/gene2CE.dat $wormpipe_dir/dumps/");
    
    system("touch $wormpipe_dir/DUMP_PREP_RUN");
  }
    else {
      print " cant find GFF files at /wormsrv2/autoace/CHROMOSOMES/ \n ";
      exit(1);
    }
  }

if( $dump_data )
  {
    unless ( -e "$wormpipe_dir/DUMP_PREP_RUN" ) {
      print "Please run wormBLAST.pl -prep_dump version $WS_version    before dumping\n\nTo dump you CAN NOT have wormsrv2 mounted\n\n";
      exit(0);
    }

    # need this to dump new databases in full
    &get_updated_database_list;
    my $anal_list = join(',',@updated_DBs);

    # Dump data for new peptides - into separate file to append
    print "Dumping new peptides for worm_pep\n";
    &run_command("perl5.6.1 $wormpipe_dir/scripts/Dump_blastp.pl -version $WS_version -matches -database worm_pep -new_peps $wormpipe_dir/dumps/new_entries.WS$WS_version");

    # updated databases (eg gadfly) need to be redumped
    print "Dumping all peptides for updated databases for worm_pep\n";
    &run_command("perl5.6.1 $wormpipe_dir/scripts/Dump_blastp.pl -version $WS_version -matches -database worm_pep -all -analysis $anal_list");

    # . . and brigpep
    print "Dumping all peptides for updated databases for worm_brig\n";
    &run_command("perl5.6.1 $wormpipe_dir/scripts/Dump_blastp.pl -version $WS_version -matches -database worm_brigpep -all -analysis $anal_list");


    # dump blastx
    print "Dumping blastx for analysis $anal_list\n";
    &run_command("perl5.6.1 $scripts_dir/dump_blastx_new.pl -version $WS_version -analysis $anal_list");

    # dump motifs for elegans and brig
    print "Dumping motifs\n";
    &run_command("perl5.6.1 $scripts_dir/dump_motif.pl");
    &run_command("perl5.6.1 $scripts_dir/dump_motif.pl -database worm_brigpep");

    # Dump extra info for SWALL proteins that have matches. Info retrieved from the dbm databases on /acari/work2a/wormpipe/
    print "Creating acefile of SWALL proteins with homologies\n";
    &run_command("perl5.6.1 $scripts_dir/write.swiss_trembl.pl -swiss -trembl");

    print "Creating acefile of matched IPI proteins\n";
    &run_command("perl5.6.1 $scripts_dir/write_ipi_info.pl");
  }


if( $cleanup ) {
  print "clearing up files generated in this build\n";
# files to move to ~wormpub/last-build/
#   /acari/work2a/wormpipe/dumps/
#    ipi_hits_list
#    trembllist.txt
#    swisslist.txt
#    best_blastp_hits
#    best_blastp_hits_brigprot

#  ~wormpipe/Elegans
#    WS99.agp
#    cds99.gff
#    cos99.gff
#    ids.txt

# to delete
#   /acari/work2a/wormpipe/dumps/
#      *.ace
#      *.log
  my $clear_dump = "/acari/work2a/wormpipe/dumps";
  print "Removing . . . \n";
  print "\t$clear_dump/*.ace\n";  system("rm -f $clear_dump/*.ace") && warn "cant remove ace files from $clear_dump";
  print "\t$clear_dump/*.log\n";  system("rm -f $clear_dump/*.log") && warn "cant remove log files from $clear_dump";

  print "Removing files currently in $wormpipe_dir/last_build/n";
  system(" rm -f $wormpipe_dir/last_build/*.gff");
  system(" rm -f $wormpipe_dir/last_build/*.agp");

  print "\nmoving the following to ~wormpipe/last_build . . \n";
  print "\t$clear_dump/*.txt\n"; system("mv -f $clear_dump/*.txt $wormpipe_dir/last_build/") && warn "cant move $clear_dump/*.txt\n";
  print "\t$clear_dump/ipi*\n"; system("mv -f $clear_dump/ipi* $wormpipe_dir/last_build/") && warn "cant move $clear_dump/ipi*\n";
  print "\t$clear_dump/best_blastp\n"; system("mv -f $clear_dump/best_blastp* $wormpipe_dir/last_build/") && warn "cant move $clear_dump/best_blast*\n";
  print "\t$wormpipe_dir/Elegans/*\n"; system("mv -f $wormpipe_dir/Elegans/* $wormpipe_dir/last_build/") && warn "cant move $wormpipe_dir/Elegans/*\n";

  print "\nRemoving the $wormpipe_dir/DUMP_PREP_RUN lock file\n"; system("rm -f $wormpipe_dir/DUMP_PREP_RUN") && warn "cant remove $wormpipe_dir/DUMP_PREP_RUN\n";

  print "\nRemoving farm output and error files from /acari/scatch5/ensembl/Worms/*\n"; system("rm -rf /acari/scratch5/ensembl/Worms/*") && warn "cant clear up /acari/scatch5/ensembl/Worms/*\n";

  print "\n\nCLEAN UP COMPLETED\n\n";
}


&wait_for_pipeline_to_finish if $test_pipeline; # debug stuff


##############################################
# close log file and mail $maintainer report #
##############################################

my $rundate    = `date +%y%m%d`; chomp $rundate;
print LOG "\n# wormBLAST.pl finished at: $rundate ",&runtime,"\n";
close LOG;

# warn about errors in subject line if there were any
if($errors == 0){
  &mail_maintainer("BUILD BLAST REPORT: ",$maintainers,$log);
}
elsif ($errors ==1){
  &mail_maintainer("BUILD BLAST REPORT: $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer("BUILD BLAST REPORT: $errors ERRORS!!!",$maintainers,$log);
}


exit(0);



###############################################################################################
#
#
#                          T  H  E     S  U  B  R  O  U  T  I  N  E  S
#
#
#
################################################################################################


sub wait_for_pipeline_to_finish
  {
    my $finished = 0;
    while( $finished == 0 ) {
      my $jobsleft = `cjobs`;
      chomp $jobsleft;
      if( $jobsleft == 0 ){
	$finished = 1;
	print "pipeline finished\n" ;
      }
      else {
	print "$jobsleft jobsleft (Im going to sleep for that long! )\n";
	
	sleep $jobsleft;
      }
    }
    print "Pipeline finished - waiting 60 secs to make sure everything is through\n";
    sleep 60;
    print "DONE\n\n";
    return;
  }

sub check_wormsrv2_conflicts
  {
    if( ( -e "/wormsrv2" ) && ($blastx || $blastp || $dump_data) )
      {
	print "no can do - to run the blast pipeline wormsrv2 can NOT be mounted.  Your options conflict with this.\nthe following options REQUIRE wormsrv2 to not be mounted \n\t-blastx\t-blastp\t-dump\t-run\t-testpipe\t-run_brig\t-cleanup\n\n";
	exit (1);
      }
    elsif( !(-e "/wormsrv2") && ($chromosomes || $wormpep  || $update_mySQL || $prep_dump) ) {
      print "The following option need access to wormsrv2 and it aint there (rsh wormsrv2 )\n";
      print "-chromosomes\t-wormpep\t-databases\t-updatemysql\t-prep_dump\n";
      exit(1);
    }
    else {
      print "wormsrv2 requirements and command line options checked and seem OK \n";
    }
  }

sub update_database
  {
    if( $dont_SQL ){
      return;
    }
    else{
      my $query = shift;
      my $db = shift;
      my $sth = $db->prepare( "$query" );
      $sth->execute();
      return;
    }
  }

sub single_line_query
  {
    if( $dont_SQL ){
      my @bogus = qw(3 3 3 3 3 3);
      return @bogus;
    }
    else{
      my $query = shift;
      my $db = shift;
      my $sth = $db->prepare( "$query" );
      $sth->execute();
      my @results = $sth->fetchrow_array();
      $sth->finish();
      return @results;
    }
  }


sub get_updated_database_list
  {
    @updated_DBs = ();
    #process new databases
    open (OLD_DB,"<$last_build_DBs") or die "cant find $last_build_DBs";
    my %prevDBs;
   # my %currentDBs;
    
    #get database file info from databases_used_WS(xx-1) (should have been updated by script if databases changed
    while (<OLD_DB>) {
      chomp;
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl_1|slimtrembl_2|wormpep|ipi_human|brigpep)/ ) {
	$prevDBs{$1} = $_;
      }
    }  
    open (CURR_DB,"<$database_to_use") or die "cant find $database_to_use";
    while (<CURR_DB>) {
      chomp;
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl_1|slimtrembl_2|wormpep|ipi_human|brigpep)/ ) {
	$currentDBs{$1} = $_;
      }
    }
    close CURR_DB;
    
    #compare old and new database list
    foreach (keys %currentDBs){
      if( "$currentDBs{$_}" ne "$prevDBs{$_}" ){
	push( @updated_DBs, "$_");
      }
    }
  }

############################################################################################
sub create_log_files{

  # create main log ile using script name for
  $log        = "/nfs/acari/wormpipe/logs/wormBLAST.${WS_version}.${rundate}.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$0\n";
  print LOG "started at ",&runtime,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}


sub run_RuleManager
  {
    my ($dbname, $moltype ) = @_;
    my $script;
    $script = "$bdir/RuleManager3.pl" if $moltype eq "dna";
    $script = "$bdir/RuleManager3Prot.pl" if $moltype eq "pep";

    die "invalid or no moltype passed to run_RuleManager : $moltype\n" unless $script;
    &run_command("perl $script -dbhost $dbhost -dbname $dbname -dbpass $dbpass -dbuser $dbuser");

  }

##########################################################
sub run_command{
  my $command = shift;
  print LOG &runtime, ": started running $command\n";
  my $status = system($command);
  if($status != 0){
    $errors++;
    print LOG "ERROR: $command failed\n";
  }
  print LOG &runtime, ": finished running\n\n";

  # for optional further testing by calling subroutine
  return($status);
}
############################################################


__END__
  
=pod

  
=head2 NAME - wormBLAST.pl

=head1 USAGE
  
=over 4
 
=item wormBLAST.pl [-chromosomes -wormpep -databases -datemysql -setup -run -nosql -dump -mail]
  
=back
  
This script:
  
I<wormBLAST.pl MANDATORY arguments:>

B<NONE>

I<wormBLAST.pl  Overview:>

This script is a collection of subroutines that automate the BLAST part of the Wormbase build process.  It keeps track of which databases are being used in the files ~wormpipe/BlastDB/databases_used_WSXX and databases_used_WSXX-1.  There are certain problems with this at the moment mounting acari and wormsrv2.  Most of this can be run from wormsrv2, but when actually running and dumping do this on acari.  The Wormbase.pm module is used by this script, so until all the Pipeline scripts are under CVS we'll have to live with it.


I<wormBLAST.pl  OPTIONAL arguments:>

B<-version XX>    Certain parts of this script dont work if wormsrv2 is mounted which makes getting the current build version via Wormbase.pm difficult!  So, now Wormbase.pm is not used at all and you have to enter the version of the database being built on the command line.  If you dont the script will stop and tell you so.

B<-chromosomes> Runs ~wormpipe/Pipeline/copy_files_to_acari.pl -c to copy the newly formed chromosomes from /wormsrv2 and cats them in to one

B<-wormpep>      Runs ~wormpipe/Pipeline/copy_files_to_acari.pl -c to take the new wormpepXX file from /wormsrv2/WORMPEP and creates a BLASTable database (setdb)

B<-databases>    Checks existing database files (slimtremblXX.pep etc) against what was used last time and updates them.  It takes the new .pep file and runs setdb.  Modifies the databases_used_WSXX file to reflect changes.  This will ensure that the update propegates to the remainder of the process.

B<-mail>           Creates and sends a mail to systems to distribute new databases over the farm.  Compares the database files used in this and the previous build to tell system which files to remove and what to replace them with.

B<-updatemysql>    Updates any new databases to be used in MySQL.  

B<-setup>          Prepares MySQL for the pipeline run.  Performs the 'delete' commands so that new data is included in the run

B<you must NOT have wormsrv2 mounted when running blast jobs>

B<-run>            Actually starts the BLAST analyses.  Does a single analysis at a time based on what new databases are being used, plus a couple of "do everything runs" to finish it all off.
The BLAST pipeline is limited so that we can only have one RuleManager running at a time.  These means that even if all the balstx jobs have been submitted we cant start the blastp run until they have all gone through.  Therefore the script allows the running of these separately by using the -blastx and -blastp options.  Both can be entered together (same as -run).  If this is done the script will monitor the progress of the blastx jobs (done 1st) and only start the blastp run when this has finished.

B<-blastx>      Submits blastx jobs based on updated databases

B<-blastp>      Submits blastp jobs based on updated databases

B<-dump>        Dumps data from MySQL after anaylsis is complete.

B<-nosql>       Debug option where SQL calls to databases are not performed. 

=back

=over 4

=head1 REQUIREMENTS

=over 4

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

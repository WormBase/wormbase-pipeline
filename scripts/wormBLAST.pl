#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
my $wormpipe_dir = glob("~wormpipe");
# no longer use wormbase.pm routines are in here. Script doesn't like it if wormsrv2 mounted !
#use libs "$wormpipe_dir/wormbase/scripts/"
#use lib "/wormsrv2/scripts/";
#use Wormbase;
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
my $dont_SQL;
my $dump_data;
my $mail;
my $test_pipeline;
my $WPver;   #  Wormpep version is passed as command line option
my $blastx;
my $blastp;
my $prep_dump;

GetOptions("chromosomes" => \$chromosomes,
	   "wormpep"     => \$wormpep,
	   "databases"   => \$update_databases,
	   "updatemysql" => \$update_mySQL,
	   "setup"       => \$setup_mySQL,
	   "run"         => \$run_pipeline,
	   "nosql"       => \$dont_SQL,
	   "prep_dump"   => \$prep_dump,
	   "dump"        => \$dump_data,
	   "mail"        => \$mail,
	   "testpipe"    => \$test_pipeline,
	   "version=s"   => \$WPver,
	   "blastp"      => \$blastp,
	   "blastx"      => \$blastx
	  );

# you can do either or both blast anaylses
if( $run_pipeline ) {
  $blastx = 1; $blastp =1;
}

&check_wormsrv2_conflicts;

#$WPver = &get_wormbase_version unless $WPver;
die "please give a build version number ie  wormBLAST -version 00\n" unless $WPver;
my $WP_old = $WPver - 1;
my $scripts_dir = "$wormpipe_dir/scripts/BLAST_scripts";
#process Ids

#|         18 | gadfly3.pep         |
#|         19 | ensembl7.29a.2.pep  |
#|         20 | yeast2.pep          |
#|         23 | wormpep87.pep       |
#|         24 | slimswissprot40.pep |
#|         25 | slimtrembl21.pep    |

my %worm01processIDs = ( wormpep => 23, 
			 yeast => 20,
			 ensembl => 19,
			 gadfly  => 18,
			 slimswissprot => 24,
			 slimtrembl =>25
		       );

#|          7 | yeast2.pep          | 
#|          8 | gadfly3.pep         |
#|          9 | ensembl7.29a.2.pep  | 
#|         11 | wormpep87.pep       | 
#|         13 | slimswissprot40.pep | 
#|         14 | slimtrembl21.pep    |

my %wormprotprocessIds = ( wormpep => 11, 
			   ensembl => 9,
			   gadfly  => 8,
			   yeast => 7,
			   slimswissprot => 13,
			   slimtrembl =>14
			 );

if( $chromosomes ) {
  #get new chromosomes
  print "Updating chromosomes\n";
  `$scripts_dir/copy_files_to_acari.pl -c`;
}

if( $wormpep ) {
  #get new wormpep
  print "Updating wormpep . . \n";
  `$scripts_dir/copy_files_to_acari.pl -w`;
}

my %currentDBs;   #ALSO used in setup_mySQL 
my @updated_DBs;  #used by get_updated_database_list sub - when run this array is filled with databases that have been updated since the prev build
#load in databases used in previous build
my $last_build_DBs = "$wormpipe_dir/BlastDB/databases_used_WS$WP_old";
my $database_to_use = "$wormpipe_dir/BlastDB/databases_used_WS$WPver";


open (OLD_DB,"<$last_build_DBs") or die "cant find $last_build_DBs";
while (<OLD_DB>) {
  chomp;
  if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl|wormpep)/ ) {
    $currentDBs{$1} = $_;
  }
}
close OLD_DB;

#check for updated Databases
if ( $update_databases )
  {
    print "Updating databases \n";
    open (DIR,"ls -l $wormpipe_dir/BlastDB/*.pep |") or die "readir\n";
    while (<DIR>) { 
      #  print;
      chomp;
      if( /\/(ensembl|gadfly|yeast|slimswissprot|slimtrembl|wormpep)/ ) {
	my $whole_file = "$1"."$'";  #match + stuff after match.
	if( "$whole_file" ne "$currentDBs{$1}" ) {
	  #make blastable database
	  print "\tmaking blastable database for $1\n";
	   `/usr/local/pubseq/bin/setdb $wormpipe_dir/BlastDB/$whole_file`;
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
    open (LETTER,">$letter");
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
    
    #print LETTER "wormpep$WPver.pep\nwormpep$WPver.pep.ahd\nwormpep$WPver.pep.atb\nwormpep$WPver.pep.bsq\n\n";
    print LETTER "CHROMOSOME_I.dna\nCHROMOSOME_II.dna\nCHROMOSOME_III.dna\nCHROMOSOME_IV.dna\nCHROMOSOME_V.dna\nCHROMOSOME_X.dna\n\n";
    
    print LETTER "-------------------------------------------------------\n\n";
    print LETTER "Thanks\n\n ________ END ________";
    close LETTER;
  
    my $name = "Wormpipe database distribution request";
    my $maintainer = "ar2\@sanger.ac.uk";
    print "mailing distibution request to $maintainer\n";
    &mail_maintainer($name,$maintainer,$letter);
  }



# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
#my $dbuser = "wormro";
my $dbname = "worm01";
my $dbpass = "worms";

my @results;
my $query = "";
my $worm01;  #worm01 Db handle
my $wormprot;#wormprot Db handle
#update mySQL database
if( $update_mySQL )
  {
    print "Updating mysql databases\n";
    #make worm01 connection
    $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #internal_id number of the last clone in the worm01
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
	$query = "select * from clone order by internal_id desc limit 1";
	@results = &single_line_query( $query, $worm01 );
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
    `cat /wormsrv2/autoace/CHROMOSOMES/*.agp > $wormpipe_dir/Elegans/WS$WPver.agp`;
    
    #load information about any new clones
    print "\tloading information about any new clones\n";
    `$scripts_dir/agp2ensembl.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -agp $wormpipe_dir/Elegans/WS$WPver.agp -write -v -strict`;
    
    #check that the number of clones in the clone table equals the number of contigs and dna objects
    my ($clone_count, $contig_count, $dna_count);
    $query = "select count(*) from clone";
    @results = &single_line_query( $query, $worm01 );
    $clone_count = $results[0];
    
    $query = "select count(*) from contig";
    @results = &single_line_query( $query, $worm01 );
    $contig_count = $results[0];
    
    $query = "select count(*) from dna";
    @results = &single_line_query( $query, $worm01);
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
    
    $query = "select * from clone order by internal_id desc limit 1";
    @results = &single_line_query( $query, $worm01 );
    my $new_last_clone = $results[0];
    print "\tnew_last_clone = $new_last_clone\n";
    
    if( $last_clone != $new_last_clone )
      {
	$query = "select id from contig where internal_id > $last_clone into outfile '$wormpipe_dir/Elegans/ids.txt'";
	print &update_database( $query, $worm01 );
	
	`$scripts_dir/InputIdManager.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -insert -analysis SubmitContig -class contig -file $wormpipe_dir/Elegans/ids.txt`;
	
      }
    $worm01->disconnect;
    
    print "\tchecking for duplicate clones\n";
    `$scripts_dir/find_duplicate_clones.pl`;
    
    
    #add new peptides to MySQL database
    print "\n\nAdding new peptides to wormprot\n";
    #make wormprot connection
    $dbname = "wormprot";
    $wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    $query = "select * from protein order by proteinId desc limit 1";
    @results = &single_line_query( $query, $wormprot );
    my $old_topCE = $results[0];
    
    if (-e "/wormsrv2/WORMPEP/wormpep$WPver/new_entries.WS$WPver"){
      `$scripts_dir/worm_pipeline.pl -f /wormsrv2/WORMPEP/wormpep$WPver/new_entries.WS$WPver`;
    }
    else {
      die "new_entries.WS$WPver does not exist! \nThis should have been made in autoace_minder -buildpep\n";
    }
    
    
    #check for updated ids
    @results = &single_line_query( $query, $wormprot );
    my $new_topCE = $results[0];
    if( "$old_topCE" eq "$new_topCE" ) {
      print "\tNO new peptides were added to the wormprot mysql database\n";
    }
    else {
      print "\tnew highest proteinId is $new_topCE (old was $old_topCE )\n";
    }
    
    $wormprot->disconnect;
  }



$dbuser = "wormadmin";
$dbpass = "worms";
if( $setup_mySQL )
  {  
    print "Setting up mysql ready for Blast run\n";
    #make wormprot connection
    $dbname = "wormprot";
    my $wormprot =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    $dbname = "worm01";
    #make worm01 connection
    $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    &get_updated_database_list;
    
    # if the user passes WPversion greater than that in the current file update it anyway
    # (this means you can update the database before wormpep is made - ie during autoace_minder -build
    $currentDBs{$1} =~ /wormpep(\d+)/;
    if( $1 and ( $1<$WPver ) ) {
      push (@updated_DBs,"wormpep$WPver\.pep");
    }
    #   mysql -h ecs1f -u wormadmin -p worm01
    #   delete from InputIdAnalysis where analysisId = 23; (DNA clones BLASTX'd against wormpep)
    #   delete from feature where analysis = 23;
    #   use wormprot;
    #   delete from InputIdAnalysis where analysisId = 11 (wormpep proteins BLASTP'd against other proteins)
    #   delete from protein_feature where analysis = 11
    
    #update mysql with which databases need to be run against
    foreach my $database (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$database};
	my $db_file = $currentDBs{$database};
	print "________________________________________________________________________________\n";
	print "doing worm01 updates . . . \n";
	$query = "update analysisprocess set db = \"$db_file\" where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $worm01 );
	
	$query = "update analysisprocess set db_file = \"/data/blastdb/Worms/$db_file\" where analysisId = $analysis";
	print $query,"\n\n";
	&update_database( $query, $worm01 );
	
	#delete entries so they get rerun
	$query = "delete from InputIdAnalysis where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $worm01 );
	
	$query = "delete from feature where analysis = $analysis";
	print $query,"\n";
	&update_database( $query, $worm01 );
	
	
	
	print "doing wormprot updates . . . \n";
	$analysis = $wormprotprocessIds{$database};
	$query = "update analysisprocess set db = \"$db_file\" where analysisId = $analysis";
	print $query,"\n";
	
	&update_database( $query, $wormprot );
	$query = "update analysisprocess set db_file = \"/data/blastdb/Worms/$db_file\" where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $wormprot );
	
	#delete entries so they get rerun
	$query = "delete from InputIdAnalysis where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $wormprot );
	
	$query = "delete from protein_feature where analysis = $analysis";
	print $query,"\n";
	&update_database( $query, $wormprot );
      }
    
    $worm01->disconnect;
    $wormprot->disconnect;
    
  }
my $bdir = "/nfs/farm/Worms/EnsEMBL/branch-ensembl-121/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline";
if( $blastx ) 
  {   
    die "can't run pipeline whilst wormsrv2 is mounted - please exit and try again\n" if (-e "/wormsrv2");
    #make sure we have the databases to work on.
    &get_updated_database_list;
    #run worm01 stuff
    
    `cp -f $bdir/pipeConf.pl.worm01 $bdir/pipeConf.pl` and die "cant copy pipeConf worm01 file\n";
    
    #any updated databases
    # worm01 stuff
    foreach (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$_};
	`perl $bdir/RuleManager3.pl -once -flushsize 5 -analysis $analysis`;
      }
    `perl $bdir/RuleManager3.pl -once -flushsize 5`;#finish off anything that didn't work
    
    &wait_for_pipeline_to_finish if $blastp;
  }

if( $blastp )
  {
    &wait_for_pipeline_to_finish;
    die "can't run pipeline whilst wormsrv2 is mounted - please exit and try again\n" if (-e "/wormsrv2");
    #make sure we have the databases to work on.
    &get_updated_database_list;
    
    #run wormpep stuff
    `cp -f $bdir/pipeConf.pl.wormprot $bdir/pipeConf.pl` and die "cant copy pipeConf wormprot file\n";   
    
    #for anything updated
    foreach (@updated_DBs)
      {
	my $analysis = $wormprotprocessIds{$_};
	`perl $bdir/RuleManager3Prot.pl -once -flushsize 5 -analysis $analysis`;
      }
    `perl $bdir/RuleManager3Prot.pl -once -flushsize 5`; # finish off anything that didn't work + PFams and low complexity, signalp, ncoils,transmembrane
  }

if( $prep_dump ) 
  {
       # prepare helper files

    if( -e "/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_X.gff") {
      print 
      `cat /wormsrv2/autoace/CHROMOSOMES/*.gff | $scripts_dir/gff2cds.pl > /nfs/acari/wormpipe/Elegans/cds$WPver.gff`;
      `cat /wormsrv2/autoace/CHROMOSOMES/*.gff | $scripts_dir/gff2cos.pl > /nfs/acari/wormpipe/Elegans/cos$WPver.gff`;
      `$scripts_dir/prepare_dump_blastx.pl > $wormpipe_dir/dumps/accession2clone.list`;
      `cp /wormsrv2/WORMPEP/wormpep$WPver/wormpep.diff$WPver $wormpipe_dir/dumps/`;
      `cp /wormsrv2/autoace/COMMOM_DATA/CE2gene.dat $wormpipe_dir/dumps/`;
      `cp /wormsrv2/autoace/COMMOM_DATA/gene2CE.dat $wormpipe_dir/dumps/`;

      `touch $wormpipe_dir/DUMP_PREP_RUN`;
    }     
    else {
      print " cant find GFF files at /wormsrv2/autoace/CHROMOSOMES/ \n ";
      exit(1);
    }
  }

if( $dump_data )
  {
    unless ( -e "$wormpipe_dir/DUMP_PREP_RUN" ) {
      print "Please run wormBLAST.pl -prep_dump version $WPver    before dumping\n\nTo dump you CAN NOT have wormsrv2 mounted\n\n";
      exit(0);
    }
    # Dump
    print "Dumping blastp\n";
    #`$wormpipe_dir/scripts/Dump_new_prot_only.pl -all -version $WPver`;
    print "Dumping blastx\n";
    `$scripts_dir/dump_blastx_new.pl -w $wormpipe_dir/BlastDB/wormpep$WPver.pep -a ~/Elegans/WS$WPver.agp -g ~/Elegans/cds$WPver.gff -c ~/Elegans/cos$WPver.gff -m`;
    print "Dumping motifs\n";
      `$scripts_dir/dump_motif.pl`;
    
    # Dump extra info for SWALL proteins that have matches. Info retrieved from the dbm databases on /acari/work2a/wormpipe/
    print "Creating acefile of SWALL proteins with homologies\n";
    `$scripts_dir/write.swiss_trembl.pl -swiss -trembl`;
  }


&wait_for_pipeline_to_finish if $test_pipeline; # debug stuff

exit(0);

sub wait_for_pipeline_to_finish
  {
    my $finished = 0;
    while( $finished == 0 ) {
      my $jobsleft = `bjobs | wc -l`;
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
    if( ($chromosomes || $wormpep  || $update_databases || $update_mySQL || $prep_dump) && ($blastx || $blastp || $dump_data) ) 
      {
	print "no can do - to run the blast pipeline wormsrv2 can NOT be mounted.  Your options conflict with this.\nthe following options REQUIRE wormsrv2 - \n\t-chromosomes\t-wormpep\t-updatemysql\n\n";
	exit (1);
      }
    else {
      print "command line options checked and seem OK (in terms of wormsrv2 requirements )\n\n";
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
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl|wormpep)/ ) {
	$prevDBs{$1} = $_;
      }
    }  
    open (CURR_DB,"<$database_to_use") or die "cant find $database_to_use";
    while (<CURR_DB>) {
      chomp;
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl|wormpep)/ ) {
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


sub get_wormbase_version {

    my $WS_version = `grep "NAME WS" /wormsrv2/autoace/wspec/database.wrm`;
    chomp($WS_version);
    $WS_version =~ s/.*WS//;    
    return($WS_version);
}

sub mail_maintainer {
    my ($name,$maintainer,$logfile) = @_;
    $maintainer = "dl1\@sanger.ac.uk, ar2\@sanger.ac.uk, ck1\@sanger.ac.uk, krb\@sanger.ac.uk" if ($maintainer eq "All");
    open (OUTLOG,  "|/bin/mailx -s \"$name\" $maintainer ");
    if ( $logfile )
      {
	open (READLOG, "<$logfile");
	while (<READLOG>) 
	  { print OUTLOG "$_";
	  }
	close READLOG;
      }
    else {
      print OUTLOG "$name";
    }
    close OUTLOG;
  }

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

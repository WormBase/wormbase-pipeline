#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
my $antdir = glob("~ar2");
#use lib "/wormsrv2/scripts/";
use lib "/nfs/team71/worm/ar2/wormbase_cvs/scripts/";
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
my $dont_SQL;
my $dump_data;

GetOptions("chromosomes" => \$chromosomes,
	   "wormpep"     => \$wormpep,
	   "databases"   => \$update_databases,
	   "updatemysql" => \$update_mySQL,
	   "setup"       => \$setup_mySQL,
	   "run"         => \$run_pipeline,
	   "nosql"       => \$dont_SQL,
	   "dump"        => \$dump_data
	  );

my $wormpipe_dir = glob("~wormpipe");
my $WPver = &get_wormbase_version;
my $WP_old = $WPver - 1;

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
  `$wormpipe_dir/Pipeline/copy_files_to_acari.pl -c`;
}

if( $wormpep ) {
  #get new wormpep
  `$wormpipe_dir/Pipeline/copy_files_to_acari.pl -w`;
}

my %currentDBs;   #ALSO used in setup_mySQL 

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

my @updated_DBs;
#check for updated Databases
if ( $update_databases )
  {
    open (DIR,"ls -l $wormpipe_dir/BlastDB/*.pep |") or die "readir\n";
    while (<DIR>) { 
      #  print;
      chomp;
      if( /\/(ensembl|gadfly|yeast|slimswissprot|slimtrembl|wormpep)/ ) {
	my $whole_file = "$1"."$'";  #match + stuff after match.
	if( "$whole_file" ne "$currentDBs{$1}" ) {
	  #make blastable database
	  `setdb $whole_file`;
	  push( @updated_DBs,$1 );
	  #change hash entry ready to rewrite external_dbs
	  $currentDBs{$1} = "$whole_file";
	 	}
	else {
	  print "$1 database unchanged $whole_file\n";
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
&mail_maintainer($name,$maintainer,$letter);




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
    #make worm01 connection
     $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #internal_id number of the last clone in the worm01
    my $last_clone;
    
    $query = "select * from clone order by internal_id desc limit 1";
    @results = &single_line_query( $query, $worm01 );
    $last_clone = $results[0];
    print "last_clone = $last_clone\n";
    
    #Make a concatenation of all six agp files from the last release to ~/Elegans  e.g.
    `cat /wormsrv2/current_DB/CHROMOSOMES/*.agp > $wormpipe_dir/Elegans/WS$WPver.agp`;
    
    #load information about any new clones
    `$wormpipe_dir/Pipeline/agp2ensembl.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -agp $wormpipe_dir/Elegans/WS$WPver.agp -write -v -strict`;
    
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
    
    if( ($clone_count != $contig_count) or ($contig_count != $dna_count ) ){
      print "the number of clones, contigs and DNAs is inconsistant\n
clones = $clone_count\ncontigs = $contig_count\ndna = $dna_count\n";
      exit(0);
    }
    
    $query = "select * from clone order by internal_id desc limit 1";
    @results = &single_line_query( $query, $worm01 );
    my $new_last_clone = $results[0];
    print "new_last_clone = $new_last_clone\n";
    $last_clone = 3487;

    if( $last_clone != $new_last_clone )
      {
	$query = "select id from contig where internal_id > $last_clone into outfile '$wormpipe_dir/Elegans/ids.txt'";
	print &update_database( $query, $worm01 );
	
	`$wormpipe_dir/Pipeline/InputIdManager.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -insert -analysis SubmitContig -class contig -file $wormpipe_dir/Elegans/ids.txt`;
	
      }
     $worm01->disconnect;
    `$wormpipe_dir/Pipeline/find_duplicate_clones.pl`;


     #add new peptides to MySQL database
     #make wormprot connection
     $dbname = "wormprot";
     $wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
       || die "cannot connect to db, $DBI::errstr";

     $query = "select * from protein order by proteinId desc limit 1";
     @results = &single_line_query( $query, $wormprot );
     my $old_topCE = $results[0];
    # `$wormpipe_dir/Pipeline/worm_pipeline.pl -f /wormsrv2/WORMPEP/wormpep$WPver/new_entries.WS$WPver`;
     
     #check for updated ids
     @results = &single_line_query( $query, $wormprot );
     my $new_topCE = $results[0];
     if( "$old_topCE" eq "$new_topCE" ) {
       print "NO new peptides were added to the wormprot mysql database\n";
     }
     else {
       print "new highest proteinId is $new_topCE (old was $old_topCE )\n";
     }

     $wormprot->disconnect;
   }



$dbuser = "wormadmin";
$dbpass = "worms";
if( $setup_mySQL )
  {    
    #make wormprot connection
    $dbname = "wormprot";
    my $wormprot =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    $dbname = "worm01";
    #make worm01 connection
    $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #process new databases
    open (OLD_DB,"<$last_build_DBs") or die "cant find $last_build_DBs";
    my %prevDBs;
    my %currentDBs;
    my @updated_DBs;
    
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
    
    #update mysql with which databases need to be run against
    foreach my $database (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$database};
	my $db_file = $currentDBs{$database};
	print "________________________________________________________________________________\n";
	$query = "update analysisprocess set db = \"$db_file\" where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $worm01 );
	$query = "update analysisprocess set db_file = \"/data/blastdb/Worms/$db_file\" where analysisId = $analysis";
	print $query,"\n\n";
	&update_database( $query, $worm01 );

	$analysis = $wormprotprocessIds{$database};
	$query = "update analysisprocess set db = \"$db_file\" where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $wormprot );
	$query = "update analysisprocess set db_file = \"/data/blastdb/Worms/$db_file\" where analysisId = $analysis";
	print $query,"\n";
	&update_database( $query, $wormprot );
      }
    
    $worm01->disconnect;
    $wormprot->disconnect;

  }
my $bdir = "/nfs/farm/Worms/EnsEMBL/branch-ensembl-121/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline";
if( $run_pipeline )
  {   
    #run worm01 stuff
   
    `cp -f $bdir/pipeConf.pl.worm01 $bdir/pipeConf.pl` and die "cant copy pipeConf worm01 file\n";

    `perl RuleManager3.pl -once -flushsize 3 -analysis 23`;# - wormpep
    #anything updated
    foreach (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$_};
	`perl RuleManager3.pl -once -flushsize 5 -analysis $analysis`;
      }
    
    #run wormpep stuff
    `cp -f $bdir/pipeConf.pl.wormprot $bdir/pipeConf.pl` and die "cant copy pipeConf wormprot file\n";
    `perl RuleManager3.pl -once -flushsize 5 -analysis $wormprotprocessIds{wormpep}`;# - wormpep
    
    `perl RuleManager3.pl -once -flushsize 5`;  #finish off anything that didn't work
    #anything updated
    foreach (@updated_DBs)
      {
	my $analysis = $wormprotprocessIds{$_};
	`perl RuleManager3Prot.pl -once -flushsize 5 -analysis $analysis`;
      }
    `perl RuleManager3prot.pl -once -flushsize 5`; # finish off anything that didn't work + PFams and low complexity, signalp, ncoils, transmembrane
  }

if( $dump_data )
  {
    #dont forget the other stuff before this
    `$wormpipe_dir/Pipeline/dump_blastp.pl -w ~/BlastDB/wormpep$WPver.pep -s`;
    `$wormpipe_dir/Pipeline/dump_blastx_new.pl -w ~/BlastDB/wormpep$WPver.pep -a ~/Elegans/WS$WPver.agp -g ~/Elegans/cds87.gff -c ~/Elegans/cos87.gff -m`;
    `$wormpipe_dir/Pipeline/dump_motif.pl`;
  }

exit(0);
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

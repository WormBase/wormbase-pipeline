#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use lib "/wormsrv2/scripts/";
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

GetOptions("chromosomes" => \$chromosomes,
	   "wormpep"     => \$wormpep,
	   "databases"   => \$update_databases,
	   "updatemysql" => \$update_mySQL,
	   "setup"       => \$setup_mySQL,
	   "run"         => \$run_pipeline,
	   "nosql"       => \$dont_SQL
	  );

my $wormpipe_dir = glob("~wormpipe");
my $WPver = &get_wormbase_version;

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
  `$wormpipe_dir/Pipeline/copy_files_to_acari.pl`;
}

if( $wormpep ) {
  #get new wormpep
  `$wormpipe_dir/Pipeline/copy_files_to_acari.pl`;
}

my %currentDBs;   #ALSO used in setup_mySQL
my @updated_DBs;
#check for updated Databases
if ( $update_databases )
  {
    #read in existing 
    my $extDB_file = "$wormpipe_dir/BlastDB/external_dbs";
    open (OLD_DB,"<$extDB_file") or die "cant find $extDB_file";
    
    while (<OLD_DB>) {
      chomp;
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl)/ ) {
	$currentDBs{$1} = $_;
      }
    }
    
    open (DIR,"ls -l $wormpipe_dir/BlastDB/*.pep |") or die "readir\n";
    while (<DIR>) { 
      #  print;
      chomp;
      if( /\/(ensembl|gadfly|yeast|slimswissprot|slimtrembl)/ ) {
	my $whole_file = "$1"."$'";  #match + suff after match.
	if( "$whole_file" ne "$currentDBs{$1}" ) {
	  #make blastable database
	  push( @updated_DBs,$1 );
	  #change hash entry ready to rewrite external_dbs
	  $currentDBs{$1} = "$&";
	 	}
	else {
	  print "$1 database unchanged $&\n";
	}
      }
    }
  }



#@updated_DBs = qw(ensembl gadfly yeast slimswissprot slimtrembl);

#generate distribution request based on updated databases
my $letter = "/wormsrv2/logs/distribute_on_farm_mail";
open (LETTER,">$letter");
print LETTER "This is a script generated email from the wormpipe Blast analysis pipeline.\nAny problems should be addessed to worm\@sanger.ac.uk.\n
=====================================================
\n";
print LETTER "The following can be removed from /data/blastdb/Worms.\n\n";
foreach (@updated_DBs){
  print LETTER "$_*\n";
}
print LETTER "wormpep*\n";
print LETTER "CHROMOSOME_*.dna\n";

print LETTER "-------------------------------------------------------\n\nand replaced with the following files from ~wormpipe/BlastDB/\n\n";
foreach (@updated_DBs){
  print LETTER "$currentDBs{$_}\n$currentDBs{$_}.ahd\n$currentDBs{$_}.atb\n$currentDBs{$_}.bsq\n\n";
}

print LETTER "wormpep$WPver.pep\nwormpep$WPver.pep.ahd\nwormpep$WPver.pep.atb\nwormpep$WPver.pep.bsq\n\n";
print LETTER "CHROMOSOME_I.dna\nCHROMOSOME_II.dna\nCHROMOSOME_III.dna\nCHROMOSOME_IV.dna\nCHROMOSOME_V.dna\nCHROMOSOME_X.dna\n\n";

print LETTER "-------------------------------------------------------\n\n";
print LETTER "Thanks\n\n ________ END ________";
close LETTER;

my $name = "Wormpipe database distribution request";
my $maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer($name,$maintainer,$letter);




# mysql database parameters
my $dbhost = "ecs1f";
#my $dbuser = "wormadmin";
my $dbuser = "wormro";
my $dbname = "worm01";
my $dbpass = "";

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
    `cat /wormsrv2/current_DB/CHROMOSOMES/*.agp > ~/Elegans/WS$WPver.agp`;
    
    #load information about any new clones
    `~/Pipeline/agp2ensembl.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -agp ~/Elegans/WS$WPver.agp -write -v -strict`;
    
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
    
    if( $last_clone != $new_last_clone )
      {
	$query = "select id from contig where internal_id > $last_clone into outfile '$wormpipe_dir/Elegans/ids.txt'";
	print &single_line_query( $query, $worm01 );
	
	`$wormpipe_dir/Pipeline/InputIdManager.pl -dbname worm01 -dbhost ecs1f -dbuser wormadmin -dbpass worms -insert -analysis SubmitContig -class contig -file ~/Elegans/ids.txt`;
	
      }
     $worm01->disconnect;
    `$wormpipe_dir/Pipeline/find_duplicate_clones.pl`;
  }


if( $setup_mySQL )
  {    
    #make wormprot connection
    $dbname = "wormprot";
    my $wormprot =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #make worm01 connection
    $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
    
    #update mySQL database with new wormpep info
    $query = "update analysisprocess set db = 'wormpep$WPver.pep' where analysisId = 11";
    &single_line_query( $query, $wormprot );
    $query = "update analysisprocess set db_file = '/data/blastdb/Worms/wormpep$WPver.pep' where analysisId = 11";
    &single_line_query( $query, $wormprot );
    $query = "update analysisprocess set db = 'wormpep$WPver.pep' where analysisId = 23";
    &single_line_query( $query, $worm01 );
    $query = "update analysisprocess set db_file = '/data/blastdb/Worms/wormpep$WPver.pep' where analysisId = 23";
    &single_line_query( $query, $worm01 );
    
    #processes to be re-run
    $query = "delete from InputIdAnalysis where analysisId = 23";
    &single_line_query( $query, $worm01 );
    $query = "delete from feature where analysis = 23";
    &single_line_query( $query, $worm01 );
    $query = "delete from InputIdAnalysis where analysisId = 11";
    &single_line_query( $query, $wormprot );
    $query = "delete from protein_feature where analysis = 11";
    &single_line_query( $query, $wormprot );
    #process new databases
    my $extDB_file = "$wormpipe_dir/BlastDB/external_dbs";
    open (OLD_DB,"<$extDB_file") or die "cant find $extDB_file";
    my %currentDBs;
    my @updated_DBs;
    
    #get database file info from external_dbs (should have been updated by script if databases changed
    while (<OLD_DB>) {
      chomp;
      if( /(ensembl|gadfly|yeast|slimswissprot|slimtrembl)/ ) {
	$currentDBs{$1} = $_;
      }
    }
    #@updated_DBs = qw(ensembl gadfly yeast slimswissprot slimtrembl);
    foreach my $database (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$database};
	my $db_file = $currentDBs{$database};
	$query = "update analysisprocess set db = '$db_file' where analysisId = $analysis";
	print $query,"\n";
	&single_line_query( $query, $worm01 );
	$query = "update analysisprocess set db_file = '/data/blastdb/Worms/$db_file' where analysisId = $analysis";
	print $query,"\n\n\n";
	&single_line_query( $query, $worm01 );

	$analysis = $wormprotprocessIds{$database};
	$query = "update analysisprocess set db = '$db_file' where analysisId = $analysis";
	print $query,"\n";
	&single_line_query( $query, $wormprot );
	$query = "update analysisprocess set db_file = '/data/blastdb/Worms/$db_file' where analysisId = $analysis";
	print $query,"\n";
	&single_line_query( $query, $wormprot );
      }
    
    $worm01->disconnect;
    $wormprot->disconnect;

  }

if( $run_pipeline )
  {   
#    #make wormprot connection
#    $dbname = "wormprot";
#    my $wormprot =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
#      || die "cannot connect to db, $DBI::errstr";
    
#    #make worm01 connection
#    $worm01 = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
#      || die "cannot connect to db, $DBI::errstr";



    #run worm01 stuff
    my $bdir = "/nfs/farm/Worms/EnsEMBL/branch-ensembl-121/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline";
    `cp -f $bdir/pipeConf.pl.worm01 $bdir/pipeConf.pl` and die "cant copy pipeConf worm01 file\n";

    `perl RuleManager3.pl -once -flushsize 3 -analysis 23`;# - wormpep
    #anything updated
    foreach (@updated_DBs)
      {
	my $analysis = $worm01processIDs{$database};
	`perl RuleManager3.pl -once -flushsize 5 -analysis $analysis`;
      }
    
    #run wormpep stuff
    `cp -f $bdir/pipeConf.pl.wormprot $bdir/pipeConf.pl` and die "cant copy pipeCon fwormprot file\n";
    `perl RuleManager3.pl -once -flushsize 5 -analysis $wormprotprocessIDs{wormpep}`;# - wormpep
    
    `perl RuleManager3.pl -once -flushsize 5`;  #finish off anything that didn't work
    #anything updated
    foreach (@updated_DBs)
      {
	my $analysis = $wormprotprocessIDs{$database};
	`perl RuleManager3Prot.pl -once -flushsize 5 -analysis $analysis`;
      }
    `perl RuleManager3prot.pl -once -flushsize 5`; # finish off anything that didn't work + PFams and low complexity, signalp, ncoils, transmembrane
    
#    $worm01->disconnect;
#    $wormprot->disconnect;
  }

exit(0);


sub single_line_query
  {
    if( $dont_SQL ){
      my @bogus = qw(3 3 3 3 3 3);
      return @bogus;
    }
    else{
      my $query = shift;
      my $db = shift;
      my $sth = $worm01->prepare( "$query" );
      $sth->execute();
      my @results = $sth->fetchrow_array();
      $sth->finish();
      return @results;
    }
  }

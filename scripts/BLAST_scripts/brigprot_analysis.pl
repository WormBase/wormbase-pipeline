#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
my $wormpipe_dir = glob("~wormpipe");
use Getopt::Long;

#######################################
# command-line options                #
#######################################


my ($debug, $WPver, $mysql, $run, $dump);
GetOptions("debug" => \$debug,
	   "version" => \$WPver,
	   "mysql"   => \$mysql,
	   "dump"    => \$dump
	  );

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "worm_brigprot";
my $dbpass = "worms";

my @results;
my $query = "";
my $worm_brigprot;


$dbuser = "wormadmin";
$dbpass = "worms";
if( $mysql )  
  {
    print "Setting up mysql ready for Blast run\n";
    #make wormprot connection
    my $worm_brigprot =  DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";
  
    #update mysql with new wormpep version
    
    my $analysis = "11";
    my $db_file = "wormpep".$WPver.".pep";
    print "________________________________________________________________________________\n";
    print "doing worm01 updates . . . \n";
    $query = "update analysisprocess set db = \"$db_file\" where analysisId = $analysis";
    print $query,"\n";
    &update_database( $query, $worm_brigprot );
    
    $query = "update analysisprocess set db_file = \"/data/blastdb/Worms/$db_file\" where analysisId = $analysis";
    print $query,"\n\n";
    &update_database( $query, $worm_brigprot );
    
    #delete entries so they get rerun
    $query = "delete from InputIdAnalysis where analysisId = $analysis";
    print $query,"\n";
    &update_database( $query, $worm_brigprot );
    
    $query = "delete from feature where analysis = $analysis";
    print $query,"\n";
    &update_database( $query, $worm_brigprot );
    
    $worm_brigprot->disconnect;
    
  }

if( $run )
  {
    die "can't run pipeline whilst wormsrv2 is mounted - please exit and try again\n" if (-e "/wormsrv2");

    my $bdir = "/nfs/farm/Worms/EnsEMBL/branch-ensembl-121/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline";

    `cp -f $bdir/pipeConf.pl.worm_brigprot $bdir/pipeConf.pl` and die "cant copy pipeConf worm_brigprot file\n";   

    `perl $bdir/RuleManager3Prot.pl -once -flushsize 5`; # just do everything
  }

if( $dump )
  {
    unless ( -e "$wormpipe_dir/DUMP_PREP_RUN" ) {
      print "Please run wormBLAST.pl -prep_dump version $WPver    before dumping\n\nTo dump you CAN NOT have wormsrv2 mounted\n\n";
      exit(0);
    }
    # Dump
    print "Dumping blastp\n";
    `perl5.6.1 $wormpipe_dir/scripts/Dump_new_prot_only.pl -all -version $WPver -brigprot -analysis 11 -matches`;
    print "Dumping motifs\n";
      `perl5.6.1 $wormpipe_dir/scripts/BLAST_scripts/dump_motif.pl -database worm_brigprot`;
  }



exit(0);


sub update_database
  {
    my $query = shift;
    my $db = shift;
    my $sth = $db->prepare( "$query" );
    $sth->execute();
    return;
  }

sub single_line_query
  {
    my $query = shift;
    my $db = shift;
    my $sth = $db->prepare( "$query" );
    $sth->execute();
    my @results = $sth->fetchrow_array();
    $sth->finish();
    return @results;
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

  
=head2 NAME - brigprot_analysis.pl

=head1 USAGE
  
=over 4
 
=item 
  
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

#!/usr/local/ensembl/bin/perl -w

use DBI;
use strict;
my $wormpipe_dir = glob("~wormpipe");
use Getopt::Long;

#######################################
# command-line options                #
#######################################


my ($debug, $WPver, $mysql, $run, $dump, $dump_all);
GetOptions("debug" => \$debug,
	   "version:s" => \$WPver,
	   "mysql"   => \$mysql,
	   "run"     => \$run,
	   "dump"    => \$dump, 
	   "dump_all"=> \$dump_all
	  );

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "worm_brigprot";
my $dbpass = "worms";

my @results;
my $query = "";
my $worm_brigprot;

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
    print "doing $dbname updates . . . \n";
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

    $query = "delete from protein_feature where analysis = $analysis";
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

if( $dump or $dump_all )
  {
    unless ( -e "$wormpipe_dir/DUMP_PREP_RUN" ) {
      print "Please run wormBLAST.pl -prep_dump version $WPver    before dumping\n\nTo dump you CAN NOT have wormsrv2 mounted\n\n";
      exit(0);
    }
    # Dump
    print "Dumping blastp\n";
    `perl $wormpipe_dir/scripts/Dump_new_prot_only.pl -all -version $WPver -brigprot -analysis 11 -matches` if $dump;
    `perl $wormpipe_dir/scripts/Dump_new_prot_only.pl -all -version $WPver -brigprot -matches` if $dump_all;
    print "Dumping motifs\n";
      `perl $wormpipe_dir/scripts/BLAST_scripts/dump_motif.pl -database worm_brigprot`;
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
  
I<brigprot_analysis.pl MANDATORY arguments:>

B<NONE>

I<wormBLAST.pl  Overview:>

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

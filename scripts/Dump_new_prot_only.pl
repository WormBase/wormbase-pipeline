#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use lib "/wormsrv2/scripts/";
use lib "/nfs/team71/worm/ar2/wormbase_cvs/scripts/";
use Wormbase;
use Getopt::Long;
use Wormpep;

#######################################
# command-line options                #
#######################################
my $test;

my $chromosomes;
my $wormpep;
my $update_databases;
my $update_mySQL;
my $setup_mySQL;
my $run_pipeline;
my $dont_SQL;
my $dump_data;


GetOptions("test"        => \$test,
	   "chromosomes" => \$chromosomes,
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
$WPver-- if $test;

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

my %processIds2prot_analysis = ( 11 => "wublastp_worm",
				 9  => "wublastp_ensembl",
				 8  => "wublastp_fly",
				 7  => "wublastp_yeast",
				 13 => "wublastp_slimswissprot",
				 14 => "wublastp_slimtrembl"
			       );




#get list of wormpeps to dump from wormpep.diffXX
my @peps2dump;
open( DIFF,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.diff$WPver") or die "cant opne diff file\n";
while (<DIFF>)
  {
    my $pep;
    if( /new/ ){
       /CE\d+/;
       $pep = $&;
    }
    elsif ( /changed/ ){
      /-->\s+(CE\d+)/;
      $pep = $1;
    }
    if( $pep ){
      push( @peps2dump, $pep );
    }
  }
close DIFF;


# mysql database parameters
my $dbhost = "ecs1f";
#my $dbuser = "wormadmin";
my $dbuser = "wormro";
my $dbname = "wormprot";
my $dbpass = "";

my @results;
my $query = "";
my $wormprot;#wormprot Db handle

# set a E-value threshold (used in the sql query)
my $e_threshold = 1.e-6;
$wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";


my $sth_p = $wormprot->prepare ( q{ SELECT proteinId, length
                                 FROM protein
                             ORDER BY proteinId
                             } );

# protein_feature table
my $sth_f = $wormprot->prepare ( q{ SELECT proteinId,analysis,
                                      start, end,
                                      hId, hstart, hend,  
                                      evalue, cigar
                                 FROM protein_feature
                                WHERE proteinId = ? and evalue < ?
                             ORDER BY hId
	  	  	     } );

open (OUT,">wublastp.ace") or die "cant open out file\n";
my $count;

my %CE2gene = &CE2gene;
my %gene2CE = &gene2CE;

my %matched_genes;

foreach my $pep (@peps2dump)
  {
#    print OUT "//\nProtein : \"WP:$pep\"\n";

    #get gene names and create new wormpep obj
    my $gene = &justGeneName($CE2gene{$pep});
    my $pepObj = Wormpep->new($pep, $gene);

    #retreive data from mysql
    $sth_f->execute($pep, $e_threshold);
    my $ref_results = $sth_f->fetchall_arrayref;
    my ($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
    foreach my $result_row (@$ref_results)
      {
	($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar) = @$result_row;
	
	
	# take -log(10) of evalues
	if ($e != 0) {
	  $e = sprintf ("%.2f", -(log($e))/log(10));
	}
	else {
	  $e = -1;
	}
	if( $e > 40 ) {        # arbitrarilly set to 40 - not sure what should be ar2.
	  my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
	  if( $analysis == 11 )   #wormpep
	    {
	      #send the CE id rather than gene name
	      my $worm_protein = $gene2CE{$homolID};
	      next if ("$worm_protein" eq "$pep" );

	      my $pepgene = &justGeneName($homolID);
	      $data[4] = $worm_protein;
	      $pepObj->add_worm_data($pepgene,@data);
	    }
	  else {
	    $pepObj->add_data(@data);
	  }
	}
	else {
	  print "irrelevant $e \n";
	}
      }
    
    #write ace info
    my @cigar = split(/:/,"$cigar");
    
    foreach (@cigar){
      print OUT "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
      my @align = split(/\,/,$_);
      print "@align\n";
      print OUT "$align[0] $align[1]\n";
      print OUT "\n";
    }
    
    $count++;    last if( $count > 4);
  }

close OUT;

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
      my $sth = $db->prepare( "$query" );
      $sth->execute();
      my @results = $sth->fetchrow_array();
      $sth->finish();
      return @results;
    }
  }
       
sub gene2CE
  {
    my @array;
    open (FH,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.table$WPver") or die "cant open wormpep.table$WPver\n";
    while(<FH>){
      my @data = split(/\s+/,$_);
      push( @array, $data[0], $data[1]);
    }
    return @array;
  }




sub CE2gene
  {
    my @array;
    open (FH,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.accession$WPver") or die "cant open wormpep.accession$WPver\n";
    while(<FH>)
      {
	chomp;
	my @data = split;
	if( $data[1] ){
	  push(@array, $data[0],$data[1]);
	}
      }
    return @array;
  }

sub justGeneName
  {
    my $test = shift;
    if( $test =~ m/(\w+\.\d+)/ ){
      return $1;
    }
    else {
      return "nonworm";
    }
  }

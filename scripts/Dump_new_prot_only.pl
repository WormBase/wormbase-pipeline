#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;

#######################################
# command-line options                #
#######################################
my $test;



GetOptions("test"        => \$test,
	  );

my $wormpipe_dir = glob("~wormpipe");
my $WPver = &get_wormbase_version;
$WPver-- if $test;


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
my $pep;
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
my $e_threshold = 40;
$wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

#to get all peptides
#my $sth_p = $wormprot->prepare ( q{ SELECT proteinId, length
#                                 FROM protein
#                             ORDER BY proteinId
#                             } );

# protein_feature table
my $sth_f = $wormprot->prepare ( q{ SELECT proteinId,analysis,
                                      start, end,
                                      hId, hstart, hend,  
                                      -log10(evalue), cigar
                                 FROM protein_feature
                                WHERE proteinId = ? and -log10(evalue) > 40
                             ORDER BY hId
	  	  	     } );

open (OUT,">wublastp.ace") or die "cant open out file\n";
my $count;

my %CE2gene = &CE2gene;
my %gene2CE = &gene2CE;

my %matched_genes;



foreach $pep (@peps2dump)
  {
    my %worm_matches;
    my %fly_matches;
    my %human_matches;
    my %yeast_matches;
    my %swiss_matches;
    my %trembl_matches;

    my $gene = &justGeneName($CE2gene{$pep});

    #retreive data from mysql
    $sth_f->execute($pep);
    my $ref_results = $sth_f->fetchall_arrayref;
    my ($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
    foreach my $result_row (@$ref_results)
      {
	($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar) = @$result_row;
	print "$analysis\n";
	my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
	if( $analysis == 11 )   #wormpep
	  {
	    #send the CE id rather than gene name
	    my $worm_protein = $gene2CE{$homolID};
	    next if ("$worm_protein" eq "$pep" );
	    
	    my $pepgene = &justGeneName($homolID);
	    $data[4] = $worm_protein;
	  }
	elsif( $analysis == $wormprotprocessIds{'gadfly'}  ) { #gadfly peptide set also has isoforms
	  &addData ( \%fly_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'human'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%human_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'yeast'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%yeast_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'slimswissprot'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%swiss_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'slimtrembl'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%trembl_matches, \@data );
	}
      }
    
    
    &dumpData ($proteinId,\%worm_matches,\%human_matches,\%fly_matches,\%yeast_matches,\%swiss_matches,\%trembl_matches);
  }

close OUT;

exit(0);

sub dumpData
  {
    my $matches;
    my $pid = shift;
    print OUT "\nProtein : $pid\n";
    while( $matches = shift) {   #pass reference to the hash to dump
      #write ace info
      foreach (keys %$matches ){
        foreach my $data (@$matches{$_}) {
	  my @cigar = split(/:/,"$$data[0]->[8]");  #split in to blocks of homology
	  
	  foreach (@cigar){
	    #print OUT "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
	    print OUT "Pep_homol ";
            print OUT "$$data[0]->[4] ";   #  homolID
            print OUT "$$data[0]->[1] ";   #  analysis
            print OUT "$$data[0]->[7] ";   #  e value
            print OUT "$$data[0]->[2] ";   #  HomolStart
	    print OUT "$$data[0]->[3] ";   #  HomolEnd
	    print OUT "$$data[0]->[5] ";   #  pepHomolStar
	    print OUT "$$data[0]->[6] ";   #  pepHomolEnd
	    print OUT "Align ";

	    my @align = split(/\,/,$_);
	    print "@align\n";
	    print OUT "$align[0] $align[1]\n";
	  }
	}
      }
    }
  }

sub addData 
  {
    my $match = shift;
    my $data = shift;
    my $homol = $$data[4];
    if ( $$match{$homol} ){
      my $existing_e = $$match{$homol}[0]->[7];  #1st array - 7th index
      my $this_e = $$data[7];
      if( $this_e > $existing_e ) { #replace what is currently stored for this homology
	$$match{$homol} = ();
	$$match{$homol}[0] = [ @$data ];
      }
      elsif( $this_e == $existing_e )  {
	push( @{$$match{$homol} },[@$data ]);
      }
    }
    else {
      $$match{$homol}[0] = [ @$data ];
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

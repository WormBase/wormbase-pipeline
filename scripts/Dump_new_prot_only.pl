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
my $all;

GetOptions("test"     => \$test,
	   "all"      => \$all
	  );
my @sample_peps = @_;

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


my $runtime    = `date +%H:%M:%S`; chomp $runtime;

#get list of wormpeps to dump from wormpep.diffXX or wormpep.tableXX depending on wether u want to dump all or just new
my @peps2dump;
my $pep;
my $pep_input = shift;
while ( $pep_input ) {
  push(@peps2dump, $pep_input);
  $pep_input = shift;
}
unless (@peps2dump)  {
  if( $all ) {
    open( DIFF,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.table$WPver") or die "cant opne diff file\n";
    while(<DIFF>) {
      chomp;
      my @tabledata = split;
      push( @peps2dump, $tabledata[1]);
    }
  }
  else {
    open( DIFF,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.diff$WPver") or die "cant opne diff file\n";
    while (<DIFF>)
      {
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
  }
  close DIFF;
}

# mysql database parameters
my $dbhost = "ecs1f";
#my $dbuser = "wormadmin";
my $dbuser = "wormro";
my $dbname = "wormprot";
my $dbpass = "";

my @results;
my $query = "";
my $wormprot;#wormprot Db handle

# set a E-value threshold ( (not actually ) used in the sql query)
my $e_threshold = 40;
$wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

# get results from mysql for each specific peptide
my $sth_f = $wormprot->prepare ( q{ SELECT proteinId,analysis,
                                      start, end,
                                      hId, hstart, hend,  
                                      -log10(evalue), cigar
                                 FROM protein_feature
                                WHERE proteinId = ? and -log10(evalue) > 40
                             ORDER BY hId
	  	  	     } );

open (OUT,">wublastp.ace") or die "cant open out file\n";
my $recip_file = "wublastp_recip.ace";
open (RECIP,">$recip_file") or die "cant open recip file\n";
print OUT "$runtime\n";
my $count;

our %CE2gene = &CE2gene;
my %gene2CE = &gene2CE;

foreach $pep (@peps2dump)
  {
    my %worm_matches;
    my %fly_matches;
    my %human_matches;
    my %yeast_matches;
    my %swiss_matches;
    my %trembl_matches;


    #retreive data from mysql
    $sth_f->execute($pep);
    my $ref_results = $sth_f->fetchall_arrayref;
    my ($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
    foreach my $result_row (@$ref_results)
      {
	($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar) = @$result_row;
	my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
	if( $analysis == 11 )   #wormpep
	  {
	   # my $gene = &justGeneName($CE2gene{$pep});
	    &addWormData ( \%worm_matches, \@data );
	  }
	elsif( $analysis == $wormprotprocessIds{'gadfly'}  ) { #gadfly peptide set also has isoforms
	  &addFlyData ( \%fly_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'ensembl'}  ) { # others dont have isoforms so let adding routine deal with them
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
    
    
    &dumpData ($pep,\%worm_matches,\%human_matches,\%fly_matches,\%yeast_matches,\%swiss_matches,\%trembl_matches) if (%worm_matches or %human_matches or %fly_matches or %yeast_matches or %swiss_matches or %trembl_matches);
  }

$runtime    = `date +%H:%M:%S`; chomp $runtime;
print OUT "$runtime\n";
close OUT;
close RECIP;

#process the recip file so that proteins are grouped

#open(TEST,"ls -l $outdir | grep $clone.dna |");

open (SORT,"sort -t \" \" -k2 $recip_file | ");
open (RS,">recip_sorted.ace") or die "rs";


#Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 694 822
#Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 773 903
#Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 870 1002

my $current_pep;
while (<SORT>) {
  chomp;
  my @info = split(/line/,$_);
  if( $current_pep ) {
    if( "$current_pep" eq "$info[0]" ){
      print RS "Pep_homol $info[1]\n";
    }
    else {
      $current_pep = $info[0];
      print RS "\n$current_pep\n";
      print RS "Pep_homol $info[1]\n";
    }
  }
  else {
    $current_pep = $info[0];
    print RS "\n$current_pep\n";
  }
}
  

exit(0);


sub dumpData
  {
    my $matches;
    my $pid = shift; 
    print OUT "\nProtein : WP:$pid\n";
    while( $matches = shift) {   #pass reference to the hash to dump
      #write ace info
      my $output_count = 0;
    HOMOLOGUES:
      foreach (sort {$$matches{$b}[0]->[7] <=> $$matches{$a}[0]->[7]} keys %$matches ){
        foreach my $data_array_ref (@$matches{$_}) {
	  last HOMOLOGUES if $output_count++ ==  10; # only output the top 10

	  foreach my $data (@$data_array_ref) {
	    my @cigar = split(/:/,"$$data[8]");  #split in to blocks of homology

	    #need to convert form gene name to CE id for worms (they are stored as genes to compare isoforms)
	    if( "$$data[1]" eq "wublastp_worm" ) {
	      my $gene = $$data[4]; 
	      $$data[4] = $gene2CE{"$gene"};
	    }
	    
	    foreach (@cigar){
	      #print OUT "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
	      print OUT "Pep_homol ";
	      print OUT "$$data[4] ";   #  homolID
	      print OUT "$$data[1] ";   #  analysis
	      print OUT "$$data[7] ";   #  e value
	      print OUT "$$data[2] ";   #  HomolStart
	      print OUT "$$data[3] ";   #  HomolEnd
	      print OUT "$$data[5] ";   #  pepHomolStar
	      print OUT "$$data[6] ";   #  pepHomolEnd
	      print OUT "Align ";
	      
	      my @align = split(/\,/,$_);
	      print OUT "$align[0] $align[1]\n";
	      
	      
	      #and print out the reciprocal homology to different file
	      print RECIP "Protein : WP:$$data[4] line "; #  matching peptide
	      print RECIP "$pid ";              #worm protein
	      print RECIP "$$data[1] ";   #  analysis
	      print RECIP "$$data[7] ";   #  e value
	      print RECIP "$$data[2] ";   #  HomolStart
	      print RECIP "$$data[3] ";   #  HomolEnd
	      print RECIP "$$data[5] ";   #  pepHomolStar
	      print RECIP "$$data[6] ";   #  pepHomolEnd
	      print RECIP "Align ";
	      print RECIP "$align[0] $align[1]\n";
	    }
	  }
	}
      }
    }
  }


sub addFlyData 
  {
    my $match = shift;   #hash to add data to 
    my $data = shift;    #array data to analyse
    my $homol = $$data[4];

    my $i = 0;
    foreach (keys %$match ) { # check against all previously matched if there is a matching protein
      next unless $$match{$_}[0]->[0]; # previously removed isoforms still have valid keys
      my $existing_gene = &justGeneName( "$_" );
      my $homol_gene = &justGeneName( "$homol" );

      if( "$homol_gene" eq "$existing_gene" ) { 
	# the result being add is an isoform of a protein already matched - check and replace if nec.
	my $existing_e = $$match{$_}[0]->[7];  #1st array - 7th index   evalue will be same in all arrays
	my $this_e = $$data[7];

	if( $this_e > $existing_e ) { #replace what is currently stored for this homology
	  undef $$match{$_}[0];
	  $$match{$homol}[0] = [ @$data ];
	}
	elsif( $this_e == $existing_e )  {
	  push( @{$$match{'$_'}},[ @$data ]) if ("$homol" eq "$_"); #if protein matches 2 isoforms with same evalue just keep 1st one.
	}
	return; # shouldn't need to check against rest
      } 
    }
    # if we get to here it must be a match to something that is not an isoform of things already matched or a self-match
    $$match{$homol}[0] = [ @$data ];
  }


sub addWormData 
  {
    my $match = shift;   #hash to add data to 
    my $data = shift;    #array data to analyse
    my $my_gene = &justGeneName( $CE2gene{ $$data[0] } ) ;
    my $homol = $$data[4];
    my $homol_gene = &justGeneName( $homol );
    return if ("$homol_gene" eq "$my_gene");    # self match or isoform of same gene

    #have we already matched an isoform of this protein
    # CE26000 | wublast_worm | start | end | Y73B6BL.34    fields or result array
    my $i = 0;
    foreach (keys %$match ) { # check against all previously matched if there is a matching protein
      next unless $$match{$_}[0]->[0]; # previously removed isoforms still have valid keys
      my $existing_gene = &justGeneName( "$_" );

      if( "$homol_gene" eq "$existing_gene" ) { 
	# the result being add is an isoform of a protein already matched - check and replace if nec.
	my $existing_e = $$match{$_}[0]->[7];  #1st array - 7th index   evalue will be same in all arrays
	my $this_e = $$data[7];

	if( $this_e > $existing_e ) { #replace what is currently stored for this homology
	  undef $$match{$_}[0];
	  $$match{$homol}[0] = [ @$data ];
	}
	elsif( $this_e == $existing_e )  {
	  push( @{$$match{"$_"}},[ @$data ]) ;#if ("$homol" eq "$_"); #if protein matches 2 isoforms with same evalue just keep 1st one.
	}
	return; # shouldn't need to check against rest
      } 
    }
    # if we get to here it must be a match to something that is not an isoform of things already matched or a self-match
    $$match{$homol}[0] = [ @$data ];
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
    
    if( $test =~ m/(\w+\.\d+)/ ){   # worm gene
      return $1;
    }
    elsif( $test =~ m/(\w+)-\p{IsUpper}{2}/ ) {
      return $1;
    }
    else {
      return "$test";
    }
  }

#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Common_data;
use Getopt::Long;

my $maintainers = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;
my $log = "/wormsrv2/logs/Dump_new_prot_only.pl.$rundate";
open( LOG, ">$log") || die "cant open $log";
print LOG "Dump_new_prot_only.pl log file $rundate ",&runtime,"\n";
print LOG "-----------------------------------------------------\n\n";

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
my $acedb_database;
if( $test ) {
  $WPver-- ;
  $acedb_database = "/wormsrv1/antace";
  $maintainers = "ar2\@sanger.ac.uk";
}
else {
  $acedb_database = "/wormsrv2/autoace";
}
  


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

our %org_prefix = ( 'wublastp_worm' => 'WP',
		    'wublastp_ensembl' => 'ENSEMBL',
		    'wublastp_fly' => 'GADFLY',
		    'wublastp_yeast' => 'SGD',
		    'wublastp_slimswissprot' => 'SW',
		    'wublastp_slimtrembl' => 'TR'
		  );

#get list of wormpeps to dump from wormpep.diffXX or wormpep.tableXX depending on wether u want to dump all or just new
my @peps2dump;
my $pep;
my $pep_input = shift;
while ( $pep_input ) {
  print LOG "Dumping blastp results for $pep_input\n";
  push(@peps2dump, $pep_input);
  $pep_input = shift;
}
unless (@peps2dump)  {
  if( $all ) {
    print LOG &runtime," : Dumping all current wormpep proteins ( Wormpep$WPver )\n";
    my %tmp_peps;
    &CE2gene(\%tmp_peps);
    foreach (keys %tmp_peps) {
      push( @peps2dump, $_ );
    }
  }
  #  open( DIFF,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep$WPver") or die "cant open Wormpep$WPver file\n";
#    while(<DIFF>) {
#      if( /^>/ ) {
#	chomp;
#	my @tabledata = split;
#	push( @peps2dump, $tabledata[1]);
#      }
  
  
  else {
    open( DIFF,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep.diff$WPver") or die "cant opne diff file\n";
    print LOG &runtime," : Dumping updated proteins ( wormpep.diff$WPver )\n";
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
    close DIFF;
  }
 # close DIFF;
}

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "wormprot";
my $dbpass = "";
$runtime    = `date +%H:%M:%S`; chomp $runtime;
print LOG "\n",&runtime," : Connecting to database : $dbname on $dbhost as $dbuser\n";

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
my $output = "/wormsrv2/wormbase/ensembl_dumps/blastp_ensembl.ace";
open (OUT,">$output") or die "cant open $output\n";
my $recip_file = "/wormsrv2/tmp/wublastp_recip.ace";
open (RECIP,">$recip_file") or die "cant open recip file\n";
my $count;

our %CE2gene;
&CE2gene(\%CE2gene);

my %gene2CE;
&gene2CE(\%gene2CE);;

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

close OUT;
close RECIP;

print LOG &runtime," : Data extraction complete\n\n";


#process the recip file so that proteins are grouped
print LOG &runtime," : processing the reciprocal data file for efficient loading.\n";
open (SORT,"sort -t \" \" -k2 $recip_file | ");
open (RS,">>$output") or die "rs"; #append this sorted data to the main homols file and load together.

print RS "\n\n"; # just to make sure we're not adding to last object.

#sample from reciprocal match file
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
  
my $wormpub = glob("~wormpub");
my $tace =  &tace;
my $command;

print LOG &runtime," : finished\n\n______END_____";

close LOG;
`rm -f $recip_file`;

&mail_maintainer("Dump_new_proteins_only.pl",$maintainers,$log);

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
	      print OUT $org_prefix{"$$data[1]"}":$$data[4] ";   #  homolID
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
	      print RECIP "Protein : ",$$data[4] line "; #  matching peptide
	      print RECIP "WP:$pid ";              #worm protein
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


__END__

=pod

=head2 NAME - Dump_new_prot_only.pl

=head1 USAGE 

=over 4

=item Dump_new_prot_only.pl

=back

This script is to extract blastp results from the mysql database.
You can choose to get all of current Wormpep proteins, just those updated in the last build or any specified on the command line.

=over 4

=item none

=back

OPTIONAL arguments:

default is just updated proteins

-all get all of current Wormpep proteins

Dump_new_prot_only.pl CE00095 - get result for only this peptide

=over 4

=item none

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use Getopt::Long;
use DB_File;

#######################################
# command-line options                #
#######################################
my ($test, $debug, $help, $all, $WPver, $analysisTOdump);
GetOptions ("debug=s"   => \$debug,
	    "test"    => \$test,
	    "help"    => \$help,
	    "all"     => \$all,
	    "analysis=s" => \$analysisTOdump,
	    "version=s" =>\$WPver
           );

my $maintainers = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $wormpipe_dir = "/acari/work2a/wormpipe";
my $log = "$wormpipe_dir/Dump_new_prot_only.pl.$rundate";
my $best_hits = "$wormpipe_dir/dumps/best_blastp_hits";
open (BEST, ">$best_hits") or die "cant open $best_hits for writing\n";

open ( LOG, ">$log") || die "cant open $log";
print LOG "Dump_new_prot_only.pl log file $rundate\n";
print LOG "-----------------------------------------------------\n\n";

# to be able to include only those proteins that have homologies we need to record those that do
# this file is then used by write_ipi_info.pl
my $ipi_file = "$wormpipe_dir/dumps/ipi_hits_list";
open (IPI_HITS,">$ipi_file") or die "cant open $ipi_file\n";


# help page
&usage("Help") if ($help);

# no debug name
print "DEBUG = \"$debug\"\n\n" if $debug;
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);


my @sample_peps = @_;
my $wormpipe = glob("~wormpipe");
my $output = "$wormpipe_dir/dumps/blastp_ensembl.ace";
my $recip_file = "$wormpipe_dir/dumps/wublastp_recip.ace";

$WPver-- if( $test );
  

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
			   slimtrembl_1 =>14,
			   slimtrembl_2 =>16,
			   ipi_human => 15,
			   brigpep   => 17
			 );

my %processIds2prot_analysis = ( 11 => "wublastp_worm",
				 9  => "wublastp_ensembl",
				 8  => "wublastp_fly",
				 7  => "wublastp_yeast",
				 13 => "wublastp_slimswissprot",
				 14 => "wublastp_slimtrembl",
				 16 => "wublastp_slimtrembl",# slimtrembl is too large so is split
				 15 => "wublastp_human",
				 17 => "wublastp_briggsae"
			       );

our %org_prefix = ( 'wublastp_worm' => 'WP',
		    'wublastp_ensembl' => 'ENSEMBL',
		    'wublastp_fly' => 'GADFLY',
		    'wublastp_yeast' => 'SGD',
		    'wublastp_slimswissprot' => 'SW',
		    'wublastp_slimtrembl' => 'TR',
		    'wublastp_briggsae'  => 'BP',
		    'wublastp_ipi_human' => 'IP'      # should never actually get included
		  );
# gene CE info from COMMON_DATA files copied to ~wormpipe/dumps in prep_dump
undef $/;
our %CE2gene;
open (C2G ,"<$wormpipe/dumps/CE2gene.dat" );
my $in_data = <C2G>;
my $VAR1;
eval $in_data;
die if $@;
close C2G;
%CE2gene = (%$VAR1);

undef $VAR1;

my %gene2CE;
open (G2C ,"<$wormpipe/dumps/gene2CE.dat" );
$in_data = <G2C>;
eval $in_data;
die if $@;
close C2G;
%gene2CE = (%$VAR1);
$/ = "\n";

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
    print LOG " : Dumping all current wormpep proteins ( Wormpep$WPver )\n";
    foreach (keys %CE2gene) {
      push( @peps2dump, $_ );
    }
  }
  
  else {
    open( DIFF,"<$wormpipe/Elegans/wormpep.diff$WPver") or die "cant opne diff file $wormpipe/Elegans/wormpep.diff$WPver\n";
    print LOG " : Dumping updated proteins ( wormpep.diff$WPver )\n";
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
$dbname .= "_dev" if $test;
my $dbpass = "";
my $runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "\n : Connecting to database : $dbname on $dbhost as $dbuser\n";

my @results;
my $query = "";
my $wormprot;#wormprot Db handle

# E-value threshold  used in the sql query is set so that -log10(evalue) >  1. This is equivalent to an
# evalue of 1e-3 (0.001)
$wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

# get results from mysql for each specific peptide
my $sth_f;
if ( $analysisTOdump ) {

  $sth_f = $wormprot->prepare ( q{ SELECT proteinId,analysis,
				     start, end,
				     hId, hstart, hend,  
				     -log10(evalue), cigar
				       FROM protein_feature
					 WHERE proteinId = ? and -log10(evalue) > 1
					   AND analysis = ?
					   ORDER BY hId
				   } );  
}
else {
  $sth_f = $wormprot->prepare ( q{ SELECT proteinId,analysis,
				     start, end,
				     hId, hstart, hend,  
				     -log10(evalue), cigar
				       FROM protein_feature
					 WHERE proteinId = ? and -log10(evalue) > 1
					   ORDER BY hId
					 } );
}
open (OUT,">$output") or die "cant open $output\n";

# reciprocals of matches ie if CE00000 matches XXXX_CAEEL the homology details need to be written for XXXX_CAEEL 
# as well.  These are put in a separate file and post processed so that all matches for XXXX_CAEEL are loaded 
# in one go for efficient loading ( cf acecompress.pl )

open (RECIP,">$recip_file") or die "cant open recip file\n";

dbmopen our %ACC2DB, "$wormpipe_dir/dumps/acc2db.dbm", 0666 or die "cannot open acc2db \n";

my $count;


foreach $pep (@peps2dump)
  {
    my %worm_matches;
    my %fly_matches;
    my %human_matches;
    my %yeast_matches;
    my %swiss_matches;
    my %trembl_matches;
    my %brig_matches;


    #retreive data from mysql
    if( $analysisTOdump ) {
      $sth_f->execute($pep, $analysisTOdump);
    }
    else {
      $sth_f->execute($pep);
    }

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
	  #&addData ( \%human_matches, \@data ); superceded by ipi_human
	}
	elsif( $analysis == $wormprotprocessIds{'yeast'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%yeast_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'slimswissprot'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%swiss_matches, \@data );
	}
	elsif( ( $analysis == $wormprotprocessIds{'slimtrembl_1'} ) || ( $analysis == $wormprotprocessIds{'slimtrembl_2'} ) ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%trembl_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'ipi_human'}  ) { # others dont have isoforms so let adding routine deal with them
	  &addData ( \%human_matches, \@data );
	}
	elsif( $analysis == $wormprotprocessIds{'brigpep'}  ) {
	  &addData ( \%brig_matches, \@data );
	}
      }
    
    
    &dumpData ($pep,\%worm_matches,\%human_matches,\%fly_matches,\%yeast_matches,\%swiss_matches,\%trembl_matches,\%brig_matches) if (%worm_matches or %human_matches or %fly_matches or %yeast_matches or %swiss_matches or %trembl_matches or %brig_matches);
  }

close OUT;
close RECIP;
close BEST;
close IPI_HITS;

print "sorting ipi_hits file . . ";
`mv $ipi_file $ipi_file._tmp`;
`sort -u $ipi_file._tmp > $ipi_file`;
`rm -f $ipi_file._tmp`;
print "DONE\n";


print LOG " : Data extraction complete\n\n";


#process the recip file so that proteins are grouped
print LOG " : processing the reciprocal data file for efficient loading.\n";
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

print LOG " : finished\n\n______END_____";

close LOG;
#`rm -f $recip_file`;


print "\nEnd of dump - moving $output to /wormsrv2/wormbase/ensembl_dumps/\n";

# Copy resulting file to wormsrv2 - leave in original place for subsequent script write.swiss_trembl
# `/usr/bin/rcp $output /wormsrv2/wormbase/ensembl_dumps/`;

dbmclose %ACC2DB;


exit(0);

sub dumpData
  {
    my $matches;
    my $pid = shift;
    my %BEST;
    print OUT "\nProtein : \"WP:$pid\"\n";
    while( $matches = shift) {   #pass reference to the hash to dump
      #write ace info
      my $output_count = 0;
      my $best = 0; #flag for best matches resets for each analysis
    HOMOLOGUES:
#Use of uninitialized value in numeric comparison (<=>) at /nfs/team71/worm/ar2/wormbase/scripts/Dump_new_prot_only.pl line 291.

      foreach (sort {$$matches{$b}[0]->[7] <=> $$matches{$a}[0]->[7]} keys %$matches ){
        foreach my $data_array_ref (@$matches{$_}) {
	  last HOMOLOGUES if $output_count++ ==  10; # only output the top 10

	  foreach my $data (@$data_array_ref) {
	    print "@$data\n" if ($debug);
	    my @cigar = split(/:/,"$$data[8]");  #split in to blocks of homology

	    #need to convert form gene name to CE id for worms (they are stored as genes to compare isoforms)
	    my $prefix = $org_prefix{"$$data[1]"};
	    if( "$$data[1]" eq "wublastp_worm" ) {
	      my $gene = $$data[4]; 
	      $$data[4] = $gene2CE{"$gene"};
	      next unless $$data[4];
	    }
	    
	    # sort out prefix - mainly for ipi_human where it can be ENS, SW, TR, LL etc
	    elsif ( "$$data[1]" eq "wublastp_human" ) {
	      $prefix = &getPrefix("$$data[4]");
	    }
	    if ($best == 0) {
	      $BEST{$$data[1]} = "$prefix:$$data[4]"."score$$data[7]";
	      $best = 1;
	    }
	    foreach (@cigar){
	      #print OUT "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
	      print OUT "Pep_homol ";
	      print OUT "\"$prefix:$$data[4]\" ";   #  homolID
	      print OUT "$$data[1] ";   #  analysis
	      print OUT "$$data[7] ";   #  e value
	      print OUT "$$data[2] ";   #  HomolStart
	      print OUT "$$data[3] ";   #  HomolEnd
	      print OUT "$$data[5] ";   #  pepHomolStar
	      print OUT "$$data[6] ";   #  pepHomolEnd
	      print OUT "Align ";
	
	      my @align = split(/\,/,$_);
	      print OUT "$align[0] $align[1]\n";
	
	      unless ("$$data[1]" eq "wublastp_worm")  #no need for WORMPEP
		{		
		  #and print out the reciprocal homology to different file
		  #prints out on single line. "line" is used to split after sorting

		  print RECIP "Protein : \"$prefix:$$data[4]\" line "; #  matching peptide
		  print RECIP "\"WP:$pid\" ";              #worm protein
		  print RECIP "$$data[1] ";   #  analysis
		  print RECIP "$$data[7] ";   #  e value
		  print RECIP "$$data[5] ";   #  HomolStart
		  print RECIP "$$data[6] ";   #  HomolEnd
		  print RECIP "$$data[2] ";   #  pepHomolStar
		  print RECIP "$$data[3] ";   #  pepHomolEnd
		  print RECIP "Align ";
		  print RECIP "$align[1] $align[0]\n";

		}
	    }
	  }
	}
      }
    }
    # output best matches
    my @to_output = qw(wublastp_worm wublastp_human wublastp_briggsae wublastp_fly wublastp_yeast wublastp_slimswissprot wublastp_slimtrembl);
    print BEST "$pid";
    foreach my $ana (@to_output) {
      if ($BEST{$ana}) {
	my($homol, $score) =  split(/score/,$BEST{$ana});
	print BEST ",$homol,$score"
      }
      else {
	print BEST ",,";
      }
    }
    print BEST "\n";
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

sub getPrefix 
  {
    my $name = shift;
    if( $ACC2DB{$name} ) {
      print IPI_HITS "$name\n";
      return $ACC2DB{$name} 
    }
    # NOTE this is only the prefix - not the method (it will look like wublastp_ipi_human ENSEMBL:ENS00342342 etc)
    if( $name =~ /ENS\w+/ ) {
      return $org_prefix{'wublastp_ensembl'};
    }
    else {
      if (length $name > 6 ) {
      return $org_prefix{'wublastp_slimswissprot'};
      }
      else {
      return $org_prefix{'wublastp_slimtrembl'};
      }
    }
  }

sub usage {
     my $error = shift;

     if ($error eq "Help") {
         # Normal help menu
         system ('perldoc',$0);
         exit (0);
     }
     elsif ($error eq "Debug") {
         # No debug bod named
         print "You haven't supplied your name\nI won't run in debug mode
         until i know who you are\n";
        exit (0);
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

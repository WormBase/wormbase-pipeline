#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use Getopt::Long;
use DB_File;

my $file = shift;
#######################################
# command-line options                #
#######################################
my ($test, $debug, $verbose, $help, $all, $WPver, $analysisTOdump, $just_matches, $matches, $list, $database);
GetOptions ("debug=s"      => \$debug,
	    "verbose"      => \$verbose,
	    "test"         => \$test,
	    "help"         => \$help,
	    "all"          => \$all,
	    "analysis=s"   => \$analysisTOdump,
	    "version=s"    => \$WPver,
	    "just_matches" => \$just_matches,
	    "matches"      => \$matches,
	    "dumplist=s"   => \$list,
	    "database=s"     => \$database
           );

die "please give me a mysql protein database eg -database worm_pep\n" unless $database;
my @sample_peps = @_;
my $dbname = $database;

my $maintainers = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $wormpipe = glob("~wormpipe");
my $wormpipe_dir = "/acari/work2a/wormpipe";

my $best_hits = "$wormpipe_dir/dumps/${database}_best_blastp_hits";
my $ipi_file = "$wormpipe_dir/dumps/${database}_ipi_hits_list";
my $output = "$wormpipe_dir/dumps/${database}_blastp.ace";
my $recip_file = "$wormpipe_dir/dumps/${database}_wublastp_recip.ace";

open (BEST, ">$best_hits") or die "cant open $best_hits for writing\n";

my $log = "$wormpipe_dir/Dump_new_prot_only.pl.$dbname.$rundate";
open ( LOG, ">$log") || die "cant open $log";
print LOG "Dump_new_prot_only.pl log file $rundate\n";
print LOG "-----------------------------------------------------\n\n";

# to be able to include only those proteins that have homologies we need to record those that do
# this file is then used by write_ipi_info.pl
open (IPI_HITS,">$ipi_file") or die "cant open $ipi_file\n";

# help page
&usage("Help") if ($help);

# no debug name
print "DEBUG = \"$debug\"\n\n" if $debug;
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);

$WPver-- if( $test );
  

#|          7 | yeast2.pep          |
#|          8 | gadfly3.pep         |
#|          9 | ensembl7.29a.2.pep  |
#|         11 | wormpep87.pep       |
#|         13 | slimswissprot40.pep |
#|         14 | slimtrembl21.pep    |

my %wormprotprocessIds = ( 'wormpep'       => '2',
			   'brigpep'       => '3',
			   'ipi_human'     => '4',
			   'yeast'         => '5',
			   'gadfly'        => '6',
			   'slimswissprot' => '7',
			   'slimtrembl'    => '8',
			   'remanei'       => '15',
			 );

my %processIds2prot_analysis = ( '2' => "wublastp_worm",
				 '3' => "wublastp_briggsae",
				 '4' => "wublastp_human",
				 '5' => "wublastp_yeast",
				 '6' => "wublastp_fly",
				 '7' => "wublastp_slimswissprot",
				 '8' => "wublastp_slimtrembl",
				 '15'=> "wublastp_remanei",
			       );

our %org_prefix = ( 'wublastp_worm'          => 'WP',
		    'wublastp_ensembl'       => 'ENSEMBL',
		    'wublastp_fly'           => 'FLYBASE',
		    'wublastp_yeast'         => 'SGD',
		    'wublastp_slimswissprot' => 'SW',
		    'wublastp_slimtrembl'    => 'TR',
		    'wublastp_briggsae'      => 'BP',
		    'wublastp_ipi_human'     => 'IP', # should never actually get included
		    'wublastp_remanei'       => 'RP',
		  );

my %database2species = ( 'worm_pep'     => 'Caenorhabditis elegans',
			 'worm_brigpep' => 'Caenorhabditis briggsae',
			 'worm_remapep' => 'Caenorhabditis remanei',
			 'species_test' => 'Caenorhabditis test',
		       );

my $QUERY_SPECIES = $database2species{"$database"};
 
#connect to DB_File databases for species determination and establish hashes

my $db_files = "/acari/work2a/wormpipe/swall_data";
my (%SWISSORG, %TREMBLORG);
dbmopen %SWISSORG, "$db_files/swissprot2org", 0666 or die "cannot open swissprot2org DBM file $db_files/swissprot2org";
unless (-s "$db_files/swissprot2des") {
  die "swissprot2des not found or empty";
}

dbmopen %TREMBLORG, "$db_files/trembl2org", 0666 or die "cannot open trembl2org DBM file";
unless (-s "$db_files/trembl2des") {
  die "trembl2des not found or empty";
}

# gene CE info from COMMON_DATA files copied to ~wormpipe/dumps in prep_dump
undef $/;
our %CE2gene;
open (C2G ,"<$wormpipe/dumps/wormpep2cds.dat" );
my $in_data = <C2G>;
my $VAR1;
eval $in_data;
die if $@;
close C2G;
%CE2gene = (%$VAR1);

undef $VAR1;

my %gene2CE;
open (G2C ,"<$wormpipe/dumps/cds2wormpep.dat" );
$in_data = <G2C>;
eval $in_data;
die if $@;
close C2G;
%gene2CE = (%$VAR1);
$/ = "\n";

my @results;

open (OUT,">$output") or die "cant open $output\n";

# reciprocals of matches ie if CE00000 matches XXXX_CAEEL the homology details need to be written for XXXX_CAEEL 
# as well.  These are put in a separate file and post processed so that all matches for XXXX_CAEEL are loaded 
# in one go for efficient loading ( cf acecompress.pl )
print "opening $recip_file";
open (RECIP,">$recip_file") or die "cant open recip file $recip_file: $!\n";

dbmopen our %ACC2DB, "$wormpipe_dir/dumps/acc2db.dbm", 0666 or warn "cannot open acc2db \n";

my $count;
my $count_limit = 10;
$count_limit = 1 if ($just_matches);

open (BLAST,"<$file") or die "file $file\n";

my $current_pep;  
my %worm_matches;
my %fly_matches;
my %human_matches;
my %yeast_matches;
my %swiss_matches;
my %trembl_matches;
my %brig_matches;
my %rem_matches;

my %type_count;

while (<BLAST>) {

  #49483856        CE00003 1       407     1       407     B0303.14        2       2114    219.327902      100     1,1,0


  #prepare data

  chomp;
  my @data_line = split;

  #assign to vars
  my ($featureId, $proteinId,  $myHomolStart, $myHomolEnd, $pepHomolStart, $pepHomolEnd, $homolID, $analysis, $score, $e, $percent_id, $cigar) = @data_line;
  #$e = -log10($e);

  # only want current proteins
  #next unless ($CE2gene{$proteinId});

  #check if next protein
  if ( $current_pep and $current_pep ne $proteinId ) {  
    &dumpData ($current_pep,\%worm_matches,\%human_matches,\%fly_matches,\%yeast_matches,\%swiss_matches,\%trembl_matches,\%brig_matches, \%rem_matches) if (%worm_matches or %human_matches or %fly_matches or %yeast_matches or %swiss_matches or %trembl_matches or %brig_matches or %rem_matches);

    #undef all hashes
    %worm_matches = ();
    %fly_matches = ();
    %human_matches = ();
    %yeast_matches = ();
    %swiss_matches = ();
    %trembl_matches = ();
    %brig_matches = ();
    %rem_matches = ();

    %type_count = ();
    #print "$proteinId\n";

  }
  $current_pep = $proteinId;

  next unless $processIds2prot_analysis{$analysis}; #only does those with valid process ids
  next if( ( defined $type_count{$analysis} ) and ( $type_count{$analysis} > 9 ) and ( $e ne "NULL") );

  $e = 1000 if( $e eq "NULL"); # mysql gives -log10(v small no) as NULL 

  my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
  next if( ($database eq "worm_pep") and !( $CE2gene{ $data[0] } ) );  # mysql database retains dead proteins but we dont want to dump them

  my $added = 0;
  if ( $analysis == $wormprotprocessIds{'wormpep'} ) #wormpep
    {
      # my $gene = &justGeneName($CE2gene{$pep});
      $added = &addWormData ( \%worm_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'gadfly'}  ) { #gadfly peptide set also has isoforms
      $added = &addFlyData ( \%fly_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'yeast'}  ) { # others dont have isoforms so let adding routine deal with them
      $added = &addData ( \%yeast_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'slimswissprot'}  ) { # others dont have isoforms so let adding routine deal with them
      $added = &addData ( \%swiss_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'slimtrembl'} ) { # others dont have isoforms so let adding routine deal with them
      $added = &addData ( \%trembl_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'ipi_human'}  ) { # others dont have isoforms so let adding routine deal with them
      $added = &addData ( \%human_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'brigpep'}  ) {
      $added = &addData ( \%brig_matches, \@data );
    } elsif ( $analysis == $wormprotprocessIds{'remanei'}  ) {
      $added = &addData ( \%rem_matches, \@data );
    }

  #this keeps track of how many hits are stored for each analysis.  Once we have 10 we can ignore the rest as the list is sorted.
  if ($added == 1) {
    if ( defined $type_count{$analysis} ) {
      $type_count{$analysis}++;
    } else {
      $type_count{$analysis} = 1;
    }
  }
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

  print RS "\n\n";		# just to make sure we're not adding to last object.

  #sample from reciprocal match file
  #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 694 822
  #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 773 903
  #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 870 1002

  undef $current_pep;
  while (<SORT>) {
    chomp;
    my @info = split(/line/,$_);
    if ( $current_pep ) {
      if ( "$current_pep" eq "$info[0]" ) {
	print RS "Pep_homol $info[1]\n";
      } else {
	$current_pep = $info[0];
	print RS "\n$current_pep\n";
	print RS "Pep_homol $info[1]\n";
      }
    } else {
      $current_pep = $info[0];
      print RS "\n$current_pep\n";
    }
  }

my $wormpub = glob("~wormpub");

print LOG " : finished\n\n______END_____";

close LOG;
#`rm -f $recip_file`;


print "\nEnd of dump/\n";

# Copy resulting file to wormsrv2 - leave in original place for subsequent script write.swiss_trembl
# `/usr/bin/rcp $output /wormsrv2/wormbase/ensembl_dumps/`;

dbmclose %ACC2DB;
dbmclose %SWISSORG;
dbmclose %TREMBLORG;


  exit(0);

  sub dumpData
    {
      my $matches;
      my $pid = shift;
      my %BEST;
      my $prot_pref = "WP";
      $prot_pref = "BP" if ($database eq "worm_brigpep" ) ; # not generic for when 3rd species arrives
      print OUT "\nProtein : \"$prot_pref:$pid\"\n";
      while ( $matches = shift) { #pass reference to the hash to dump
	#write ace info
	my $output_count = 0;
	my $best = 0;		#flag for best matches resets for each analysis

	###############################################
	#      debug loop

	#     my $stop;
	#      next unless (%$matches);
	#      foreach my $wormpep_d (keys %$matches ) {
	#	foreach my $match_d (@{$$matches{$wormpep_d}} ){
	#	  print "@{$match_d}\n";
	#	  print "no e @{$match_d}->[0]"."@{$match_d}->[4]\n" unless @{$match_d}->[7];
	#	}
	#      }

	#     end debug loop  
	###############################################

      HOMOLOGUES:
	#Use of uninitialized value in numeric comparison (<=>) at /nfs/team71/worm/ar2/wormbase/scripts/Dump_new_prot_only.pl line 291.

	foreach (sort {$$matches{$b}[0]->[7] <=> $$matches{$a}[0]->[7]} keys %$matches ) {
	  foreach my $data_array_ref (@$matches{$_}) {
	    last HOMOLOGUES if $output_count++ ==  $count_limit; # only output the top 10

	    foreach my $data (@$data_array_ref) {
	      print "@$data\n" if ($verbose);
	      my @cigar = split(/:/,"$$data[8]"); #split in to blocks of homology

	      #need to convert form gene name to CE id for worms (they are stored as genes to compare isoforms)
	      my $prefix = $org_prefix{"$$data[1]"};
	      if ( "$$data[1]" eq "wublastp_worm" ) {
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
		next HOMOLOGUES if $just_matches; # dont bother with all the rest
	      }	
	      print OUT "Pep_homol ";
	      print OUT "\"$prefix:$$data[4]\" "; #  homolID
	      print OUT "$$data[1] ";	#  analysis
	      print OUT "$$data[7] ";	#  e value
	      print OUT "$$data[2] ";	#  HomolStart
	      print OUT "$$data[3] ";	#  HomolEnd
	      print OUT "$$data[5] ";	#  pepHomolStar
	      print OUT "$$data[6] ";	#  pepHomolEnd
	      print OUT "Target_species \"",&species_lookup($$data[1], $$data[4]),"\"\n";


	      print RECIP "Protein : \"$prefix:$$data[4]\" line "; #  matching peptide
	      print RECIP "\"$prot_pref:$pid\" ";	#worm protein
	      print RECIP "$$data[1] "; #  analysis
	      print RECIP "$$data[7] "; #  e value
	      print RECIP "$$data[5] "; #  HomolStart
	      print RECIP "$$data[6] "; #  HomolEnd
	      print RECIP "$$data[2] "; #  pepHomolStar
	      print RECIP "$$data[3] "; #  pepHomolEnd
	      print RECIP "Target_species \"$QUERY_SPECIES\"\n";
	      
	      foreach (@cigar) {
		#print OUT "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
		print OUT "Pep_homol ";
		print OUT "\"$prefix:$$data[4]\" "; #  homolID
		print OUT "$$data[1] ";	#  analysis
		print OUT "$$data[7] ";	#  e value
		print OUT "$$data[2] ";	#  HomolStart
		print OUT "$$data[3] ";	#  HomolEnd
		print OUT "$$data[5] ";	#  pepHomolStar
		print OUT "$$data[6] ";	#  pepHomolEnd
		print OUT "Align ";
	
		my @align = split(/\,/,$_);
		print OUT "$align[0] $align[1]\n";
	
		unless ("$$data[1]" eq "wublastp_worm")	#no need for WORMPEP
		  {		
		    #and print out the reciprocal homology to different file
		    #prints out on single line. "line" is used to split after sorting

		    print RECIP "Protein : \"$prefix:$$data[4]\" line "; #  matching peptide
		    print RECIP "\"$prot_pref:$pid\" ";	#worm protein
		    print RECIP "$$data[1] "; #  analysis
		    print RECIP "$$data[7] "; #  e value
		    print RECIP "$$data[5] "; #  HomolStart
		    print RECIP "$$data[6] "; #  HomolEnd
		    print RECIP "$$data[2] "; #  pepHomolStar
		    print RECIP "$$data[3] "; #  pepHomolEnd
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
	  my $e = 10**(-$score);
	  print BEST ",$homol,"; printf BEST "%g",$e,"\n";
	} else {
	  print BEST ",,";
	}
      }
      print BEST "\n";
    }


  sub addFlyData 
    {
      my $match = shift;	#hash to add data to 
      my $data = shift;		#array data to analyse
      my $homol = $$data[4];

      my $i = 0;
      foreach (keys %$match ) { # check against all previously matched if there is a matching protein
	my $existing_gene = &justGeneName( "$_" );
	my $homol_gene = &justGeneName( "$homol" );

	if ( "$homol_gene" eq "$existing_gene" ) { 
	  # the result being add is an isoform of a protein already matched - check and replace if nec.
	  my $existing_e = $$match{$_}[0]->[7];	#1st array - 7th index   evalue will be same in all arrays
	  my $this_e = $$data[7];

	  if ( $this_e > $existing_e ) { #replace what is currently stored for this homology
	    delete $$match{$_};
	    $$match{$homol}[0] = [ @$data ];
	  } elsif ( $this_e == $existing_e ) {
	    push( @{$$match{'$_'}},[ @$data ]) if ("$homol" eq "$_"); #if protein matches 2 isoforms with same evalue just keep 1st one.
	  }
	  return 0;		# shouldn't need to check against rest
	} 
      }
      # if we get to here it must be a match to something that is not an isoform of things already matched or a self-match
      $$match{$homol}[0] = [ @$data ];
      return 1;
    }


  sub addWormData 
    {
      my $match = shift;	#hash to add data to 
      my $data = shift;		#array data to analyse
      my $homol = $$data[4];
      my $homol_gene = &justGeneName( $homol );

      if ( $database eq "worm_pep" ) {
	my $my_gene = &justGeneName( $CE2gene{ $$data[0] } ) ;
	return 0 if ("$homol_gene" eq "$my_gene"); # self match or isoform of same gene
      }

      #have we already matched an isoform of this protein
      # CE26000 | wublast_worm | start | end | Y73B6BL.34    fields or result array

      foreach (keys %$match ) { # check against all previously matched if there is a matching protein
	next unless $$match{$_}[0]->[0]; # previously removed isoforms still have valid keys
	my $existing_gene = &justGeneName( "$_" );

	if ( "$homol_gene" eq "$existing_gene" ) { 
	  # the result being add is an isoform of a protein already matched - check and replace if nec.
	  my $existing_e = $$match{$_}[0]->[7];	#1st array - 7th index   evalue will be same in all arrays
	  my $this_e = $$data[7];

	  if ( $this_e > $existing_e ) { #replace what is currently stored for this homology
	    delete $$match{$_};
	    $$match{$homol}[0] = [ @$data ];
	  } elsif ( $this_e == $existing_e ) {
	    push( @{$$match{"$_"}},[ @$data ]) ; #if ("$homol" eq "$_"); #if protein matches 2 isoforms with same evalue just keep 1st one.
	  }
	  return 0;		# shouldn't need to check against rest
	} 
      }
      # if we get to here it must be a match to something that is not an isoform of things already matched or a self-match
      $$match{$homol}[0] = [ @$data ];
      return 1;
    }

  sub addData 
    {
      my $match = shift;
      my $data = shift;
      my $homol = $$data[4];
      return 0 if( $homol eq "$$data[0]" );
      if ( $$match{$homol} ) {
	my $existing_e = $$match{$homol}[0]->[7]; #1st array - 7th index
	my $this_e = $$data[7];
	if ( $this_e > $existing_e ) { #replace what is currently stored for this homology
	  $$match{$homol} = ();
	  $$match{$homol}[0] = [ @$data ];
	} elsif ( $this_e == $existing_e ) {
	  push( @{$$match{$homol} },[@$data ]);
	}
      } else {
	$$match{$homol}[0] = [ @$data ];
	return 1;
      }
    }


  sub justGeneName
    {
      my $test = shift;
      unless( defined $test ) { print "NO TEST\n";return}
      if ( $test =~ m/(\w+\.\d+)/ ) { # worm gene
	return $1;
      } elsif ( $test =~ m/(\w+)-\p{IsUpper}{2}/ ) {
	return $1;
      } else {
	return "$test";
      }
    }

  sub getPrefix 
    {
      my $name = shift;
      if ( $ACC2DB{$name} ) {
	print IPI_HITS "$name\n";
	return $ACC2DB{$name} 
      }
      # NOTE this is only the prefix - not the method (it will look like wublastp_ipi_human ENSEMBL:ENS00342342 etc)
      if ( $name =~ /ENS\w+/ ) {
	return $org_prefix{'wublastp_ensembl'};
      } else {
	if (length $name > 6 ) {
	  return $org_prefix{'wublastp_slimswissprot'};
	} else {
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
    } elsif ($error eq "Debug") {
      # No debug bod named
      print "You haven't supplied your name\nI won't run in debug mode
         until i know who you are\n";
      exit (0);
    }
  }

  sub log10 {
    my $n = shift;
    return -999 if $n == 0;
    return log($n)/log(10);
  } 


sub species_lookup
    {
      my $analysis = shift;
      my $protein = shift;

      my %species = ( 'wublastp_worm'          => 'Caenorhabditis elegans',
		      'wublastp_briggsae'      => 'Caenorhabditis briggsae',
		      'wublastp_remanei'       => 'Caenorhabditis remanei',
		      'wublastp_human'         => 'Homo sapiens',
		      'wublastp_yeast'         => 'Saccharomyces cerevisiae',
		      'wublastp_fly'           => 'Drosophila melanogaster',
		      'wublastp_slimswissprot' => $SWISSORG{$protein},
		      'wublastp_slimtrembl'    => $TREMBLORG{$protein}
		      );
      my $species = $species{"$analysis"};

      $species = 'unknown' unless$species;
      return $species;
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

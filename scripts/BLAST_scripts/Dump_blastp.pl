#!/usr/local/ensembl/bin/perl -w

use DBI;
use strict;
use Getopt::Long;
use DB_File;
use File::Path;

#######################################
# command-line options                #
#######################################
my ($test, $debug, $verbose, $help, $all, $WPver, @analysisTOdump, $just_matches, $matches, $list, $database, $new_peps, $final, $hits, $specify);
GetOptions (# usual stuff
	    "debug=s"      => \$debug,
	    "verbose"      => \$verbose,
	    "test"         => \$test,
	    "help"         => \$help,
	    "version=s"    => \$WPver,

	    # analyses and required output data
	    "analysis=s"   => \@analysisTOdump,
	    "just_matches" => \$just_matches,
	    "matches"      => \$matches,

	    # select which proteins to dump
	    "dumplist=s"   => \$list,
	    "database=s"   => \$database,
	    "new_peps=s"   => \$new_peps,
	    "all"          => \$all,
	    "specify"      => \$specify,

	    "final"        => \$final,
	    "hits"         => \$hits
           );

@analysisTOdump = split(/,/,join(',',@analysisTOdump));

my @sample_peps = @_;

my $maintainers  = "All";
my $rundate      = `date +%y%m%d`; chomp $rundate;
my $wormpipe     = glob("~wormpipe");

my $wormpipe_dir = "/acari/work2a/wormpipe"; $wormpipe_dir .= "/test" if $test;
my $dbname       = "$database";
my $output_dir   = "$wormpipe_dir/dumps/$database/blastp/ACE";
my $best_hit_dir = "$wormpipe_dir/dumps/$database/HITS";
my $ipi_file     = "$wormpipe_dir/dumps/$database/blastp_ipi";

# make sure these dir exist;
my @dirs = ($output_dir, $best_hit_dir,$wormpipe_dir ) ;
foreach ( @dirs ){
 # &mkpath( $_) unless ( -e $_ );
}

# set up output file handles;
my %output;
my %recip_output;
my %best_output;

&generate_full_file if $final; # this option will just cat all the files together ready for transfer

my $log          = "$wormpipe/logs/dump_blastp.pl.$dbname.$rundate.$$"; $log.= ".test" if $test;
open ( LOG, ">$log") || die "cant open $log";
print LOG "$0 log file $rundate at, ", time,"\n";
print LOG "Dumping analysis @analysisTOdump\n";
print LOG "-----------------------------------------------------\n\n";

my $count;
my $count_limit = 10;
$count_limit = 1 if ($just_matches);


$WPver-- if( $test );
 
my %wormprotprocessIds = ( 'wormpep'       => '2',
			   'brigpep'       => '3',
			   'ipi_human'     => '4',
			   'yeast'         => '5',
			   'gadfly'        => '6',
			   'slimswissprot' => '7',
			   'slimtrembl_1'  => '8',
			   'slimtrembl_2'  => '9',
			 );

# analysis_ids to WormBase homol_data
my %processIds2prot_analysis = ( '2' => "wublastp_worm",
				 '3' => "wublastp_briggsae",
				 '4' => "wublastp_human",
				 '5' => "wublastp_yeast",
				 '6' => "wublastp_fly",
				 '7' => "wublastp_slimswissprot",
				 '8' => "wublastp_slimtrembl",
				 '9' => "wublastp_slimtrembl",# slimtrembl is too large so is split
			       );

#Wormbase specific prefix identifiers 
our %org_prefix = ( 'wublastp_worm'          => 'WP',
		    'wublastp_ensembl'       => 'ENSEMBL',
		    'wublastp_fly'           => 'GADFLY',
		    'wublastp_yeast'         => 'SGD',
		    'wublastp_slimswissprot' => 'SW',
		    'wublastp_slimtrembl'    => 'TR',
		    'wublastp_briggsae'      => 'BP',
		    'wublastp_ipi_human'     => 'IP'      # should never actually get included
		  );


my %database2prefix = ( 'worm_pep'      => 'WP',
			'worm_brigpep' => 'BP'
		      );

# choose which routine is used to add data to the store.  Elegans and Gadfly need different ones coz of their isoforms.
my %addDataSubs = ( '2' => \&addWormData,
		    '3' => \&addData,
		    '4' => \&addData,
		    '5' => \&addData,
		    '6' => \&addFlyData,
		    '7' => \&addData,
		    '8' => \&addData,
		    '9' => \&addData
		  );
		


#This is done here as it uses the above hashes.
&check_options;


# help page
&usage("Help") if ($help);

# no debug name
print "DEBUG = \"$debug\"\n\n" if $debug;
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);

# gene CE info from COMMON_DATA files copied to ~wormpipe/dumps in prep_dump
our %CE2gene;
my %gene2CE;

&load_data("$wormpipe/dumps/gene2CE.dat",\%gene2CE);
&load_data("$wormpipe/dumps/CE2gene.dat",\%CE2gene);
&generate_best_hits if $hits;


###########    Estatblish database connections    #############

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbpass = "";
my $runtime = `date +%H:%M:%S`; chomp $runtime;

print LOG "\n : Connecting to database : $dbname on $dbhost as $dbuser\n";

my @results;
my $query = "";
my $db;# Db handle

# E-value threshold  used in the sql query is set so that -log10(evalue) >  1. This is equivalent to an
# evalue of 1e-3 (0.001)
$db = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";


###########         Get proteins to dump         ############

#get list of proteins to dump from passed list file, diff file,database or %CE2gene depending on wether u want to dump all or just new
my @peps2dump;
my $pep;
my $pep_input = shift;

while ( $pep_input ) {
  print LOG "Dumping blastp results for $pep_input\n";
  push(@peps2dump, $pep_input);
  $pep_input = shift;
}
unless (@peps2dump)  {
  
  if( $database ne "worm_pep" ) {  # get all the "other species" proteins from respective database
    print LOG " : Dumping all proteins from $database\n";
    my $protein_query = $db->prepare ( q{ SELECT proteinId
					   FROM protein
					 } );

    $protein_query->execute;
    while( my $aref = $protein_query->fetchrow_arrayref ) {
      $pep = $aref->[0];
      push( @peps2dump, $pep);
    }
  }
  else { 
    if ( $all ) {
      print LOG " : Dumping all current wormpep proteins ( Wormpep$WPver )\n";
      foreach (keys %CE2gene) {
	push( @peps2dump, $_ );
      }
    }

    if ($list) {
      print LOG " : Dumping wormpep proteins from $list\n";
      open( LIST,"<$list") or die "cant open list file $list:\t$!\n";
      while (<LIST>) {
	chomp;
	push( @peps2dump, $_ );
      }
    }
    # new only really works for WORMPEP
    if ( $new_peps) {
      open( DIFF,"<$new_peps") or die "cant open file containing diff to previous database : $new_peps:\t$!\n";
      print LOG " : Dumping updated proteins ( using $new_peps )\n";
      while (<DIFF>) {
	if ( />(\w+)/ ) {
	  push( @peps2dump, $1 );
	}
      }
      close DIFF;
    }
  }
}

print LOG "got ",scalar( @peps2dump )," proteins to dump\n";

my $recip_file;
my $output;
my $best_hits;
###################     retreive data from mysql  ######################
 
dbmopen our %ACC2DB, "$wormpipe_dir/dumps/acc2db.dbm", 0666 or warn "cannot open acc2db: $!\n";
my  $sth_f = $db->prepare ( &build_SQL );  

foreach ( @peps2dump) { 
#  print LOG "\nQuery for $_ start at ",time,"\n";
  $sth_f->execute($_);  
#  print LOG "\nQuery return at ",time,"\n\n";

  my %matches;

  my $last_pep;
  my $last_analysis;

  my $ref_results = $sth_f->fetchall_arrayref;

  foreach my $result_row (@$ref_results) {
    my ($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
    ($proteinId, $analysis,  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar) = @$result_row;
    next unless $CE2gene{$proteinId};
    unless( defined $e) {
      $e = 1000;
    }

    #  # if we are moving on to a new protein - dump the last one.
    #  if ( defined ($last_pep ) and ( $last_pep ne $proteinId) ) {
    #    &dumpData ($last_pep,$last_analysis,\%matches) if (%matches);
    #    %matches = ();
    #    undef $last_pep;
    #    undef $last_analysis;
    #  }

      # if we are moving on to a new analysis - dump the last one
      if ( defined ($last_analysis) and ( $last_analysis != $analysis )) {
        &dumpData ($proteinId, $last_analysis,\%matches) if (%matches);
        %matches = ();
        undef $last_analysis;
      }

    #  # upate tracking vars
    #  $last_pep = $proteinId;
    #  $last_analysis = $analysis;

    # mysql gives -log10(v small no) as NULL 
    my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e, $cigar);
    my $add_data_sub = $addDataSubs{$analysis};
    $add_data_sub->( \%matches, \@data );

    $last_analysis = $analysis;
  }
  &dumpData ($_,$last_analysis,\%matches) if (%matches);
  %matches = ();
  undef $last_analysis;
}

#close the output file
foreach (keys %output) {
  close $output{$_};
  close $best_output{$_};
  close $recip_output{$_};
}


#&process_reciprocal_file("$recip_file");


close IPI_HITS;

print LOG " : Data extraction complete\n\n";
print LOG " : finished at ",time,"\n\n______END_____";

close LOG;
dbmclose %ACC2DB;

exit(0);

sub dumpData
  {
    my $matches;
    my $pid = shift;
    my %best;
    my $prot_pref = $database2prefix{$database};
    my $this_anal = shift;
    my $org =$processIds2prot_analysis{$this_anal} ;

    unless ($output{$org}) {
      my $mode =  defined ($all) ? ">" : ">>" ; # append or overwrite
      my $ace = "$output_dir/${org}_blastp.ace";
      open( $output{$org},     ,$mode, "$ace")               or die "cant open $ace: $!\n";
      open( $recip_output{$org},$mode, "${ace}_recip")       or die "cant open ${ace}_recip: $!\n";
      open( $best_output{$org} ,$mode, "${org}_best_blastp") or die "cant open ${org}_best_blastp: $!\n";
    }

    my $output_FH = $output{$org};
    my $best_FH = $best_output{$org};
    my $recip_FH = $recip_output{$org};
      
    print $output_FH "\nProtein : \"$prot_pref:$pid\"\n";
    while( $matches = shift) {   #pass reference to the hash to dump
      my $output_count = 0;

    HOMOLOGUES:
      foreach (sort {$$matches{$b}[0]->[7] <=> $$matches{$a}[0]->[7]} keys %$matches ){
        foreach my $data_array_ref (@$matches{$_}) {
	  last HOMOLOGUES if $output_count++ ==  $count_limit; # only output the top 10

	  foreach my $data (@$data_array_ref) {
	    print "@$data\n" if ($verbose);
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
	
	    # data is sorted on score so store the 1st one thru as best
	    unless ( %best ) {
	      $best{'id'}    = "$prefix:$$data[4]";
	      $best{'score'} = "$$data[7]";
	      next HOMOLOGUES if $just_matches; # dont bother with all the rest
	    }

	    foreach (@cigar){
	      #print $output_FH "Pep_homol \"$homolID\" $processIds2prot_analysis{$analysis} $e $myHomolStart $myHomolEnd $pepHomolStart $pepHomolEnd Align ";
	      print $output_FH "Pep_homol ";
	      print $output_FH "\"$prefix:$$data[4]\" ";   #  homolID
	      print $output_FH "$$data[1] ";   #  analysis
	      print $output_FH "$$data[7] ";   #  e value
	      print $output_FH "$$data[2] ";   #  HomolStart
	      print $output_FH "$$data[3] ";   #  HomolEnd
	      print $output_FH "$$data[5] ";   #  pepHomolStar
	      print $output_FH "$$data[6] ";   #  pepHomolEnd
	      print $output_FH "Align ";
	
	      my @align = split(/\,/,$_);
	      print $output_FH "$align[0] $align[1]\n";
	
	      unless ("$$data[1]" eq "wublastp_worm")  #no need for WORMPEP
		{		
		  #and print out the reciprocal homology to different file
		  #prints out on single line. "line" is used to split after sorting

		  print $recip_FH "Protein : \"$prefix:$$data[4]\" line "; #  matching peptide
		  print $recip_FH "\"$prot_pref:$pid\" ";              #worm protein
		  print $recip_FH "$$data[1] ";   #  analysis
		  print $recip_FH "$$data[7] ";   #  e value
		  print $recip_FH "$$data[5] ";   #  HomolStart
		  print $recip_FH "$$data[6] ";   #  HomolEnd
		  print $recip_FH "$$data[2] ";   #  pepHomolStar
		  print $recip_FH "$$data[3] ";   #  pepHomolEnd
		  print $recip_FH "Align ";
		  print $recip_FH "$align[1] $align[0]\n";

		}
	    }
	  }
	}
      }
    }
    # output best matches 

    print $best_FH "$pid,";
    if (%best) {
      my $e = 10**(-$best{'score'});
      print $best_FH $best{'id'},","; printf $best_FH "%g",$e,", \n";
    }
    else {
      print $best_FH ",,";
    }
    print$best_FH "\n";
  }


sub addFlyData 
  {
    my $match = shift;   #hash to add data to 
    my $data = shift;    #array data to analyse
    my $homol = $$data[4];

    my $i = 0;
    foreach (keys %$match ) { # check against all previously matched if there is a matching protein
      my $existing_gene = &justGeneName( "$_" );
      my $homol_gene = &justGeneName( "$homol" );

      if( "$homol_gene" eq "$existing_gene" ) { 
	# the result being add is an isoform of a protein already matched - check and replace if nec.
	my $existing_e = $$match{$_}[0]->[7];  #1st array - 7th index   evalue will be same in all arrays
	my $this_e = $$data[7];

	if( $this_e > $existing_e ) { #replace what is currently stored for this homology
	  delete $$match{$_};
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
    my $homol = $$data[4];
    my $homol_gene = &justGeneName( $homol );

    if( $database eq "worm_pep" ) {
      my $my_gene = &justGeneName( $CE2gene{ $$data[0] } ) ;
      return if ("$homol_gene" eq "$my_gene"); # self match or isoform of same gene
    }

    #have we already matched an isoform of this protein
    # CE26000 | wublast_worm | start | end | Y73B6BL.34    fields or result array

    foreach (keys %$match ) { # check against all previously matched if there is a matching protein
      next unless $$match{$_}[0]->[0]; # previously removed isoforms still have valid keys
      my $existing_gene = &justGeneName( "$_" );

      if( "$homol_gene" eq "$existing_gene" ) { 
	# the result being add is an isoform of a protein already matched - check and replace if nec.
	my $existing_e = $$match{$_}[0]->[7];  #1st array - 7th index   evalue will be same in all arrays
	my $this_e = $$data[7];

	if( $this_e > $existing_e ) { #replace what is currently stored for this homology
	  delete $$match{$_};
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

sub load_data
  {
    my $file = shift;
    my $hash = shift;
    undef $/;
    open (FH ,"<$file" );
    my $in_data = <FH>;
    my $VAR1;
    eval $in_data;
    die if $@;
    close FH;
    $/ = "\n";
    %{$hash} = (%$VAR1);
  }


sub check_options
  {
    if ($help) {
      # Normal help menu
      system ('perldoc',$0);
      exit (0);
    }
    unless( defined $database ) {
      die "You must provide a database eg \t -database worm_pep\n";
    }
    return if $hits;

    if( $new_peps and !(-e $new_peps ) ) {
      die "cant find the file of new proteins $new_peps\n";
    }

    # want to end up with integer analysis_id's on analysisTOdump but allow them to be passed as string eg yeast
    if ( @analysisTOdump ) {
      my @tmp_anals = @analysisTOdump;
      @analysisTOdump = ();
      foreach my $anal ( @tmp_anals ) {
	if(lc $anal =~ /[a .. z]/ ) { 
	  if ( $wormprotprocessIds{$anal} ) {
	    print LOG "Dumping analysis $anal\n";
	    push(@analysisTOdump,$wormprotprocessIds{$anal} );
	    
	  } 
	  else {
	    print LOG "$_ is not a valid analysis.  Try these . . \n";
	    foreach ( keys %wormprotprocessIds ) { 
	      print LOG "\t$_\n";
	    }
	    die "FAILED - see $log for details\n";
	  }
	}
	else {
	  push(@analysisTOdump,$anal );
	}
      }
    }
    else {
      print LOG "Dumping all analyses\n";
      foreach ( keys %wormprotprocessIds ) { 
	push(@analysisTOdump,$wormprotprocessIds{$_} );
      }
    }
    die "You need to tell me which proteins to dump\nUse the -all, -new_peps, -dumplist or non-elegans database\n" unless($specify || $new_peps || $list || $all || ($database ne "worm_pep") );
  }


sub process_reciprocal_file
  {
    # The *.line file is appended to each build so will include everything. This is processed to a single acefile
    my $file = shift;
    $file =~ /(.*)\.line$/;
    my $ace_file = "$1.ace";

    #process the recip file so that proteins are grouped
    print LOG " : processing the reciprocal data file $file for efficient loading.\n";
    open (SORT,"sort -t \" \" -k2 $file | ");
    open (RS,">$ace_file") or die "cant open $ace_file for sorting reciprocal hits\n"; #append this sorted data to the main homols file and load together.

    print RS "\n\n";		# just to make sure we're not adding to last object.

    #sample from reciprocal match file
    #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 694 822
    #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 773 903
    #Protein : AGO1_ARATH line CE32063 wublastp_slimswissprot 187.468521 633 882 747 1014 Align 870 1002

    my $current_pep;
    while (<SORT>) {
      chomp;
      next unless /\w/;
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
    close RS;
    close SORT;
  }

sub generate_full_file
  {
    die "Which database dump files ? \n\teg --database worm_pep";
    print "cat ing files into ${database}_ensembl_protein_info.ace\n";
    system("cat $output_dir/*.ace > $wormpipe_dir/dumps/${database}_ensembl_protein_info.ace");
    print "DONE\n";
    exit(0);
  }

sub generate_best_hits
  {
    my %best_hits;

    # read data in to hash
    foreach my $analysis ( keys %processIds2prot_analysis ) {
      my $anal_name = $processIds2prot_analysis{$analysis};
      next if $best_hits{$anal_name}; # slimtrembl has 2 analyses
      $best_hits{$anal_name}->{'file'} = "$best_hit_dir/".$anal_name."_best_hits";
      $best_hits{$anal_name}->{'id'} = $analysis;
      open( HIT_FILE,"<$best_hits{$anal_name}->{'file'}") or die "cant open $best_hits{$anal_name}->{'file'}\t$!\n";
      while( <HIT_FILE> ) {
	chomp;
	my @data = split(/\,/,$_);
	push  @{$best_hits{$anal_name}->{$data[0]}},$data[1];
	push  @{$best_hits{$anal_name}->{$data[0]}},$data[2];

      }
      close HIT_FILE;
    }

    # write the file
    open (BEST, ">$best_hit_dir/best_hits") or die "cant write best_hits $best_hit_dir/best_hits\t$!\n";
    foreach my $pep (sort keys %CE2gene ) {
      print BEST "\n$pep";
      foreach my $analysis (sort { $best_hits{$a}->{'id'}<=>$best_hits{$b}->{'id'} } keys %best_hits ) {
	next if $analysis eq "file";
	if ( $best_hits{$analysis}->{$pep} ) {
	  print BEST ",",join(',',@{$best_hits{$analysis}->{$pep}});
	}
	else {
	  print BEST ",,";
	}
      }
    }
    exit(0);
  }

sub build_SQL
  {
     my $sql = "SELECT protein_id,analysis_id, seq_start, seq_end, hit_id, hit_start, hit_end, -log10(evalue), cigar
		FROM protein_feature
		WHERE (-log10(evalue) > 1 or evalue = 0) " ;

     # add protein_ids unless we're doing all of 'em
     unless ( $all ) {
       $sql .= "AND ( ";
       my $or_protein_id = join(' OR protein_id = ',map('"'.$_.'"',@peps2dump));
       $sql .= "protein_id = $or_protein_id ) ";
     }

     # specify which analyses to do
     if ( @analysisTOdump ) {
       my  $x = join(" OR analysis_id = ",@analysisTOdump);
       $sql .= " AND ( analysis_id = $x )";
     }

     $sql .= " ORDER BY protein_id,analysis_id,hit_id";
#########################################################################
     #trying getting per protein
     $sql = "SELECT protein_id,analysis_id, seq_start, seq_end, hit_id, hit_start, hit_end, -log10(evalue), cigar
		FROM protein_feature
		WHERE (-log10(evalue) > 1 or evalue = 0) 
                AND  protein_id = ? ";

     # specify which analyses to do
     if ( @analysisTOdump ) {
       my  $x = join(" OR analysis_id = ",@analysisTOdump);
       $sql .= " AND ( analysis_id = $x )";
     }

     $sql .= " ORDER BY analysis_id, hit_id";


     return $sql;
  }

sub make_dir_struct
  {
    my $path= shift;
    my $built;
    my @dirs = split(/\//,$path);
    foreach ( @dirs ) {
      $built .= $_;
      mkdir $built unless ( -e $built );
    }
  }

__END__

=pod

=head2 NAME - Dump_blastp.pl

=head1 USAGE 

=over 4

=item Dump_blastp.pl

=back

This script is to extract blastp results from the mysql database.

To save time redumping the same data every build this script can now dump data more specifically.  If a blast database (eg gadfly ) hasnt been updated then only hits to new proteins will be dumped and appended to the existing file.

Therefore each analysis data is stored in a separate file.

If applied to a non WORMPEP mysql database we assume that the peptide set in that database ( eg brigpep ) is static and so dump all proteins for the specified analyses.

You can choose to get all of current Wormpep proteins, just those updated in the last build or any specified on the command line.
You can also select which analyses to dump

=head2 REQUIRED arguments:

B<-database>     you must specify which database to dump from. Assumed to be on ecs1f.  If not the script will need modifying

B<-version>      WS version youre building


=head2 OPTIONAL arguments:

=item CHOOSE WHICH PROTEINS

=over 4

If the mysql database selected isnt the WORMPEP one (currently worm_pep) all proteins will be dumped for which ever analyses are selected

For the WORMPEP database the following options are applicable . . 

B<-all>          get all of current Wormpep proteins

B<-new_peps>     this should be set as the wormpep.diffWSXX file generated by the wormpep building scripts.

B<-dumplist>     a file containing a list of protein names that you want to dump eg CE23435


default is just updated proteins


=item SELECT DATA TO DUMP

=over4

B<-analysis>     commma separated list of the analyses to be dumped as database analysis ids or database name 

eg    C<-analysis 3,4,gadfly>
will dump blastp data for wormpep, brigpep and gadfly for the select set of proteins

if no analysis is selected then all will be dumped.

B<-matches>     will output a "best_blastp_hits" file for the selected data.  These will need some post processing to create a full proteome file

B<-just_matches> just creates the "best_blastp_hits" file with NO blast data.

=back

=head1 EXAMPLES

=over4

I<The usual usage ( appending new WORMPEP entries to the existing output files ) will be >

Dump_blastp.pl --database worm_pep --version XXX --new_peps ~/dumps/wormpep.diffXXX --matches

When a target database ( eg Gadfly ) has changed  . .

Dump_blastp.pl--database worm_pep --version XXX  --matches -analysis gadfly

      NOTE The new proteins dump will need to be done in addition to this, so there will be some overlap.

If you just want the new protein hits against slimswissprot and yeast( analysis_id 5) (with best hits matches

Dump_blastp.pl--database worm_pep --version XXX --new_peps ~/dumps/wormpep.diffXXX --matches -analysis slimswissprot,5

Updated analysis dump for Brigpep ( for yeast )

Dump_blastp.pl--database worm_brigprot --version XXX -analysis yeast



=head1 REQUIREMENTS

=over 4

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

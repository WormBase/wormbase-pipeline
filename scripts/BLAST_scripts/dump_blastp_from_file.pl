#!/software/bin/perl -w

use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use DB_File;
use Wormbase;
use Storable;
use Log_files;
use File::Copy;


#######################################
# command-line options                #
#######################################
my ($test, $debug, , $store, $verbose, $help, $species);
my ($all, $analysisTOdump, $just_matches, $matches, $list, $database, $dumps_dir, $file);
GetOptions ("debug=s"      => \$debug,
	    "verbose"      => \$verbose,
	    "test"         => \$test,
	    "help"         => \$help,
	    "all"          => \$all,
	    "analysis=s"   => \$analysisTOdump,
	    "just_matches" => \$just_matches,
	    "matches"      => \$matches,
	    "dumplist=s"   => \$list,
	    "database=s"   => \$database,
	    "store:s"      => \$store,
	    "species:s"    => \$species,
            "dumpdir:s"    => \$dumps_dir,
	    "file:s"       => \$file,
           );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
			   );
}
$species = $wormbase->species; #in case defaulting to elegans when species not set.
my $log = Log_files->make_build_log($wormbase);

$log->log_and_die("please give me a mysql protein database eg -database worm_ensembl_elegans\n") unless $database;
$log->log_and_die("Please supply a valid location for the dumps") if not defined $dumps_dir or not -d $dumps_dir;
$log->log_and_die("Missing -file filename on command line\n") if not defined  $file;

my @sample_peps = @_;
my $dbname = $database;

my $maintainers = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;

my $best_hits = "$dumps_dir/${species}_best_blastp_hits";
my $ipi_file = "$dumps_dir/${species}_ipi_hits_list";
my $output = "$dumps_dir/${species}_blastp.ace";
my $recip_file = "$dumps_dir/${database}_wublastp_recip.ace";

open (BEST, ">$best_hits") or $log->log_and_die("cant open $best_hits for writing\n");

# $log->write_to("Dump_new_prot_only.pl log file $rundate\n"); ### don't see how that could have ever worked after ecs1f was axed
$log->write_to("-----------------------------------------------------\n\n");

# to be able to include only those proteins that have homologies we need to record those that do
# this file is then used by write_ipi_info.pl
open (IPI_HITS,">$ipi_file") or $log->log_and_die("cant open $ipi_file\n");

# help page
&usage("Help") if ($help);

# no debug name
#print "DEBUG = \"$debug\"\n\n" if $debug;
#&usage("Debug") if ((defined $debug) && ($debug eq ""));
#
# assign $maintainers if $debug set
#($maintainers = $debug . '\@sanger.ac.uk') if ($debug);


my %processIds2prot_analysis = ( 
				'wormpepp'        => 'wublastp_worm',
				'brigpepp'        => 'wublastp_briggsae',
				'ipi_humanp'      => 'wublastp_human',
				'yeastp'          => 'wublastp_yeast',
				'gadflyp'         => 'wublastp_fly',
				'slimswissprotp'  => 'wublastp_slimswissprot',
				'slimtremblp'     => 'wublastp_slimtrembl',
				'remapepp'        => 'wublastp_remanei',
				'ppapepp'         => 'wublastp_pristionchus',
				'jappepp'         => 'wublastp_japonica',
				'brepepp'         => 'wublastp_brenneri',
				'brugpepp'        => 'wublastp_brugia',
				'ovolpepp'        => 'wublastp_ovolvulus',
                                'srapepp'         => 'wublastp_sratti'
			       );


our %org_prefix = ( 
	            'wublastp_worm'          => 'WP',
		    'wublastp_ensembl'       => 'ENSEMBL',
		    'wublastp_fly'           => 'FLYBASE',
		    'wublastp_yeast'         => 'SGD',
		    'wublastp_slimswissprot' => 'SW',
		    'wublastp_slimtrembl'    => 'TR',
		    'wublastp_briggsae'      => 'BP',
		    'wublastp_ipi_human'     => 'IP', # should never actually get included
		    'wublastp_remanei'       => 'RP',
		    'wublastp_pristionchus'  => 'PP',
		    'wublastp_japonica'      => 'JA',
		    'wublastp_brenneri'      => 'CN',
		    'wublastp_brugia'        => 'BM',
		    'wublastp_ovolvulus'     => 'OV',
                    'wublastp_sratti'        => 'SRP',
		  );

my $QUERY_SPECIES = $wormbase->full_name;
 
#connect to DB_File databases for species determination and establish hashes

my $wormpipe_dir  = $ENV{'PIPELINE'};
my $db_files_dir  = "$wormpipe_dir/swall_data";
my $dbm_files_dir = "$wormpipe_dir/dumps";

my %file_mapping = ( 
	"$db_files_dir/swissprot2org" => "/tmp/${species}_dir/swissprot2org",
	"$db_files_dir/trembl2org"    => "/tmp/${species}_dir/trembl2org",
	"$dbm_files_dir/acc2db.dbm" => "/tmp/${species}_dir/acc2db.dbm",
);
mkdir "/tmp/${species}_dir", 0777; 
while (my($from,$to)=each %file_mapping ){
  unlink $to if -e $to;
  copy $from,$to or $log->log_and_die("Copy of $from to $to failed: $!");
}

my (%SWISSORG, %TREMBLORG);
tie (%SWISSORG, 'DB_File', "/tmp/${species}_dir/swissprot2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open swissprot2org DBM file /tmp/${species}_dir/swissprot2org\n");
unless (-s "$db_files_dir/swissprot2org") { $log->log_and_die("swissprot2org not found or empty");}
tie (%TREMBLORG, 'DB_File', "/tmp/${species}_dir/trembl2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open /tmp/${species}_dir/trembl2org DBM file\n");
unless (-s "$db_files_dir/trembl2org") { $log->log_and_die("trembl2org not found or empty");}

# gene CE info from COMMON_DATA files
# using flattened arrays to merge hashes ... don't try this at home ... or with big hashes
my %CE2gene;
my %gene2CE;
$wormbase->FetchData('wormpep2cds',\%CE2gene);
$wormbase->FetchData('cds2wormpep',\%gene2CE);

my %species_accessor=$wormbase->species_accessors;
foreach my $key(keys %species_accessor){
 my %tmp_hash;
 print STDERR "getting data for $key\n";
 my $commondata = $species_accessor{$key}->common_data;
 $species_accessor{$key}->FetchData('wormpep2cds',\%tmp_hash) if -e "$commondata/wormpep2cds.dat";
 %CE2gene=(%CE2gene,%tmp_hash);
 $species_accessor{$key}->FetchData('cds2wormpep',\%tmp_hash) if -e "$commondata/cds2wormpep.dat";
 %gene2CE=(%gene2CE,%tmp_hash);
}

my @results;

unlink $output if -e $output; # remove leftover blastp acefile from the last build
open (OUT,">$output") or $log->log_and_die("cant open $output\n");

# reciprocals of matches ie if CE00000 matches XXXX_CAEEL the homology details need to be written for XXXX_CAEEL 
# as well.  These are put in a separate file and post processed so that all matches for XXXX_CAEEL are loaded 
# in one go for efficient loading ( cf acecompress.pl )
print "opening $recip_file";
open (RECIP,">$recip_file") or $log->log_and_die("cant open recip file $recip_file: $!\n");

our %ACC2DB;
tie (%ACC2DB, 'DB_File', "/tmp/${species}_dir/acc2db.dbm", O_RDWR|O_CREAT, 0666, $DB_HASH) or warn "cannot open /tmp/${species}_dir/acc2db.dbm \n";


my $count;
my $count_limit = 10;
$count_limit = 1 if ($just_matches);

open (BLAST,"<$file") or $log->log_and_die("file $file\n");

my $current_pep;  
my %worm_matches;
my %fly_matches;
my %human_matches;
my %yeast_matches;
my %swiss_matches;
my %trembl_matches;
my %brig_matches;
my %rem_matches;
my %ppa_matches;
my %jap_matches;
my %bre_matches;
my %bru_matches;
my %ovo_matches;
my %sra_matches;

my %type_count;

while (<BLAST>) {

#  +--------------------+-----------------------+-----------+---------+-----------+---------+-----------+-------------+-------+---------+------------+
#  | protein_feature_id | translation_stable_id | seq_start | seq_end | hit_start | hit_end | hit_id    | logic_name  | score | evalue  | perc_ident |
#  +--------------------+-----------------------+-----------+---------+-----------+---------+-----------+-------------+-------+---------+------------+
#  |                  1 |         2L52.1        |        27 |     168 |         1 |     157 | PF07801.2 |    wormpepX | 229.9 | 5.7e-66 |          0 |
#  +--------------------+-----------------------+-----------+---------+-----------+---------+-----------+-------------+-------+---------+------------+


  #prepare data

  chomp;
  my @data_line = split;

  #assign to vars
  my ($featureId, $cdsid,  $myHomolStart, $myHomolEnd, $pepHomolStart, $pepHomolEnd, $hit_id, $analysis, $score, $e, $percent_id) = @data_line;
  # $e = -log10($e);
  my $proteinId=($gene2CE{$cdsid}||$cdsid);
  my $homolID=($gene2CE{$hit_id}||$hit_id);

  next if $proteinId eq $homolID; # self-hit removal

  # check if next protein
  if ( $current_pep and $current_pep ne $proteinId ) {  
    &dumpData ($current_pep,
               \%worm_matches,
               \%human_matches,
               \%fly_matches,
               \%yeast_matches,
               \%swiss_matches,
               \%trembl_matches,
               \%brig_matches, 
               \%rem_matches,
               \%jap_matches,
               \%bre_matches,
               \%ppa_matches,
               \%bru_matches,
               \%ovo_matches,
               \%sra_matches) 
        if (%worm_matches or 
            %human_matches or 
            %fly_matches or 
            %yeast_matches or 
            %swiss_matches or 
            %trembl_matches or 
            %brig_matches or 
            %rem_matches or 
            %jap_matches or 
            %bre_matches or 
            %ppa_matches or 
            %bru_matches or
            %ovo_matches or
            %sra_matches
        );

    #undef all hashes
    %worm_matches = ();
    %fly_matches = ();
    %human_matches = ();
    %yeast_matches = ();
    %swiss_matches = ();
    %trembl_matches = ();
    %brig_matches = ();
    %rem_matches = ();
    %ppa_matches = ();
    %jap_matches = ();
    %bre_matches = ();
    %bru_matches = ();
    %ovo_matches = ();
    %sra_matches = ();
    
    %type_count = ();

  }
  $current_pep = $proteinId;

  next if( ( defined $type_count{$analysis} ) and ( $type_count{$analysis} > 9 ) and ( $e ne "NULL") );

  $e = 1000 if( $e eq "NULL"); # mysql gives -log10(v small no) as NULL 

  my @data = ($proteinId, $processIds2prot_analysis{$analysis},  $myHomolStart, $myHomolEnd, $homolID, $pepHomolStart, $pepHomolEnd, $e);

  my $added = 0;

  if ( $analysis eq 'wormpepp' ){ # wormpep
    $added = &addWormData ( \%worm_matches, \@data );
  } elsif ( $analysis eq 'gadflyp'  ) { # gadfly peptide set also has isoforms
    $added = &addFlyData ( \%fly_matches, \@data );
    
    # others dont have isoforms so let adding routine deal with them
  } elsif ( $analysis eq 'yeastp'  ) {
    $added = &addData ( \%yeast_matches, \@data );
  } elsif ( $analysis eq 'slimswissprotp'  ) {
    $added = &addData ( \%swiss_matches, \@data );
  } elsif ( $analysis eq 'slimtremblp') {
    $added = &addData ( \%trembl_matches, \@data );
  } elsif ( $analysis eq 'ipi_humanp') {
    $added = &addData ( \%human_matches, \@data );
  } elsif ( $analysis eq 'brigpepp') {
    $added = &addWormData ( \%brig_matches, \@data );
  } elsif ( $analysis eq 'remapepp') {
    $added = &addWormData ( \%rem_matches, \@data );
  } elsif ( $analysis eq 'ppapepp') {
    $added = &addWormData ( \%ppa_matches, \@data );
  } elsif ( $analysis eq 'jappepp') {
    $added = &addWormData ( \%jap_matches, \@data);
  } elsif ( $analysis eq 'brepepp') {
    $added = &addWormData ( \%bre_matches,\@data);
  } elsif ( $analysis eq 'brugpepp') {
    $added = &addWormData ( \%bru_matches,\@data);
  } elsif ( $analysis eq 'ovolpepp') {
    $added = &addWormData ( \%ovo_matches,\@data);
  } elsif ( $analysis eq 'srapepp') {
    $added = &addWormData ( \%sra_matches,\@data);
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

&dumpData ($current_pep,
           \%worm_matches,
           \%human_matches,
           \%fly_matches,
           \%yeast_matches,
           \%swiss_matches,
           \%trembl_matches,
           \%brig_matches, 
           \%rem_matches,
           \%jap_matches,
           \%bre_matches,
           \%ppa_matches,
           \%bru_matches,
           \%ovo_matches,
           \%sra_matches) 
    if (%worm_matches or 
        %human_matches or 
        %fly_matches or 
        %yeast_matches or 
        %swiss_matches or 
        %trembl_matches or 
        %brig_matches or 
        %rem_matches or 
        %jap_matches or 
        %bre_matches or 
        %ppa_matches or 
        %bru_matches or 
        %ovo_matches or
        %sra_matches);


close OUT;
close RECIP;
close BEST;
close IPI_HITS;

print "\nsorting ipi_hits file . . ";
$wormbase->run_command("mv $ipi_file $ipi_file._tmp",        $log);
$wormbase->run_command("sort -u $ipi_file._tmp > $ipi_file", $log);
$wormbase->run_command("rm -f $ipi_file._tmp",         		$log);
print "DONE\n";
$log->write_to(" : Data extraction complete\n\n");


##########################
# reciprocal hits part   #
##########################

# process the recip file so that proteins are grouped
$log->write_to(": processing the reciprocal data file for efficient loading.\n");
open (SORT,"sort -t \" \" -k2 $recip_file | ");
open (RS,">>$output") or $log->log_and_die("Failed to open $output for appending");

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
$log->write_to(" : finished\n\n______END_____");

untie %ACC2DB;
untie %SWISSORG;
untie %TREMBLORG;

# cleanup
while (my($from,$to)=each %file_mapping ){
	 unlink $to or $log->log_and_die("delete of $to failed: $!");
}

print "\nEnd of dump/\n";
$log->mail;
exit(0);

sub dumpData {
      my $matches;
      my $pid = shift;
      my %BEST;
      my $prot_pref = $wormbase->wormpep_prefix;
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
	# Use of uninitialized value in numeric comparison (<=>) at /nfs/team71/worm/ar2/wormbase/scripts/Dump_new_prot_only.pl line 291.

	foreach (sort {$$matches{$b}[0]->[7] <=> $$matches{$a}[0]->[7]} keys %$matches ) {
	  foreach my $data_array_ref (@$matches{$_}) {
	    last HOMOLOGUES if $output_count++ ==  $count_limit; # only output the top 10

	    foreach my $data (@$data_array_ref) {
	      print "@$data\n" if ($verbose);

	      # need to convert form gene name to CE id for worms (they are stored as genes to compare isoforms)
	      my $prefix = $org_prefix{"$$data[1]"};
	      if ( "$$data[1]" eq "wublastp_worm" ) {
		      #my $gene = $$data[4]; 
		      #$$data[4] = $gene2CE{"$gene"};
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
	      print OUT "$$data[5] ";	#  pepHomolStart
	      print OUT "$$data[6] ";	#  pepHomolEnd
	      print OUT "Target_species \"",&species_lookup($$data[1], $$data[4]),"\"\n";

	      print RECIP "Protein : \"$prefix:$$data[4]\" line "; #  matching peptide
	      print RECIP "\"$prot_pref:$pid\" ";	# worm protein
	      print RECIP "$$data[1] "; #  analysis
	      print RECIP "$$data[7] "; #  e value
	      print RECIP "$$data[5] "; #  HomolStart
	      print RECIP "$$data[6] "; #  HomolEnd
	      print RECIP "$$data[2] "; #  pepHomolStart
	      print RECIP "$$data[3] "; #  pepHomolEnd
	      print RECIP "Target_species \"$QUERY_SPECIES\"\n";
	    }
	  }
	}
    }
      # output best matches
      print BEST "$pid";
      foreach my $ana (values %processIds2prot_analysis) {
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


sub addFlyData {
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


sub addWormData {
      my $match = shift;	# hash to add data to 
      my $data = shift;		# array data to analyse
      my $homol = $$data[4];
      my $homol_gene = &justGeneName( $homol );

      #if ( $database eq "worm_remanei" || $database eq 'worm_pristionchus' || $database eq 'worm_japonica' || $database eq 'worm_brenneri') {
	my @_genes=split(/\s/,$$data[0]);
	foreach my $_g(@_genes) {
	 	my $my_gene = &justGeneName( $CE2gene{ $$data[0] } ) ;
		return 0 if ("$homol_gene" eq "$my_gene");       # self match or isoform of same gene
		return 0 if (&justGeneName($CE2gene{$homol_gene}) eq $my_gene); # elegans case
	}
      #}

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

sub addData {
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


sub justGeneName {
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

sub getPrefix {
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


sub species_lookup {
      my $analysis = shift;
      my $protein = shift;

      my %species = ( 'wublastp_worm'          => 'Caenorhabditis elegans',
		      'wublastp_briggsae'      => 'Caenorhabditis briggsae',
		      'wublastp_brugia'        => 'Brugia malayi',
		      'wublastp_ovolvulus'     => 'Onchocerca volvulus',
		      'wublastp_sratti'        => 'Strongyloides ratti',
		      'wublastp_remanei'       => 'Caenorhabditis remanei',
		      'wublastp_japonica'      => 'Caenorhabditis japonica',
		      'wublastp_brenneri'      => 'Caenorhabditis brenneri',
		      'wublastp_human'         => 'Homo sapiens',
		      'wublastp_yeast'         => 'Saccharomyces cerevisiae',
		      'wublastp_fly'           => 'Drosophila melanogaster',
		      'wublastp_pristionchus'  => 'Pristionchus pacificus',
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

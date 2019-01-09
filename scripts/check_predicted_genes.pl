#!/software/bin/perl -w
#
# Check_gene.pl
# Replacement for check_predicted_genes.pl
# 
# by gw3
#


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $build, $species, $database);


my %databases = ( # names of common databases and their species
		 camace       => 'elegans',
		 autoace      => 'elegans',
		 current_DB   => 'elegans',
		 WS           => 'elegans',
		 elegans      => 'elegans',
		 briggsae     => 'briggsae',
		 brenneri     => 'brenneri',
		 brugia       => 'brugia',
		 japonica     => 'japonica',
		 remanei      => 'remanei',
		 pristionchus => 'pristionchus',
		 tmuris       => 'tmuris',
		 ovolvulus    => 'ovolvulus',
		 sratti       => 'sratti',
);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database, # database to check
	    "build"      => \$build,    # Checks specific to a full database containing genes and models and proteins.
	    "species:s"  => \$species,
	    );


if (!defined $database) {die "-database is not specified\n"}
if (!defined $species) {
  foreach my $db (keys %databases) {
    if ($database =~ /$db/) {
      $species = $databases{$db};
    }
  }
  if (!defined $species) {die "-species is not defined and it cannot be deduced from the name of the database $database.\n"}
}



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# Establish database connection etc.
$log->log_and_die("Please use -database <path> specify a valid database directory.\n\n") if (!defined($database));
# Specify which tace to use if you are using -program flag
my $tace = $wormbase->tace;
my $db = Ace->connect(-path=>$database) or  $log->log_and_die("Couldn't connect to $database\n". Ace->error);

# create separate arrays for different classes of errors (1 = most severe, 4 = least severe)
our (@error1, @error2, @error3, @error4, @error5,);


my $cds_regex = $wormbase->cds_regex;
my $cds_regex_noend = $wormbase->seq_name_regex;
my $speciesfn = $wormbase->full_name;

################################
#         Main Body            # 
################################

my %sequence_names;
my %sequence_classes;
my %sequence_structures;


# these CDS sequence do not start with the canonical AUG codon
# many of them are micro-ORF genes taken from the sORF database
# but about 10 of them are larger genes
  my %non_canonical_initiation = (
                                'W09G3.8b' => 1, # 'gtg'
                                'B0041.12' => 1, # 'gtg'
                                'B0361.14' => 1, # 'ttg'
                                'C07G2.5' => 1, # 'ttg'
                                'C08B6.18' => 1, # 'ttg'
                                'C09E8.6' => 1, # 'gtg'
                                'C13F10.8' => 1, # 'gtg'
                                'C14C10.10' => 1, # 'ctg'
                                'C26H9A.4' => 1, # 'gtg'
                                'C34D10.5' => 1, # 'ctg'
                                'C34E10.12b' => 1, # 'gtg'
                                'C37C3.8c' => 1, # 'att' well conserved across the Caenorhabditids
                                'C48A7.16' => 1, # 'gtg'
                                'C53H9.4' => 1, # 'gtg'
                                'C54A12.4' => 1, # 'ctg' described in WormBook
                                'C56G2.22' => 1, # 'ttg'
                                'F07F6.2' => 1, # 'att' looks reasonable
                                'F08B12.11' => 1, # 'ctg'
                                'F08C6.18' => 1, # 'ctg'
                                'F13G3.15' => 1, # 'gtg'
                                'F20B6.10a' => 1, # 'ttg'
                                'F20B6.10b' => 1, # 'gtg'
                                'F27D4.12' => 1, # 'ttg'
                                'F43C9.2c' => 1, # 'ttg'
                                'F44E5.17' => 1, # 'ctg'
                                'F45H11.8' => 1, # 'atc'
                                'F52G3.7' => 1, # 'ttg'
                                'F53F10.2c' => 1, # 'ttg'
                                'F53F10.2d' => 1, # 'ttg'
                                'F53F10.2e' => 1, # 'ctg'
                                'F53F10.2f' => 1, # 'gtg'
                                'F53G2.12b' => 1, # 'ctg'
                                'F53G2.12c' => 1, # 'gtg'
                                'K04G2.11' => 1, # scbp-2
                                'K11D12.14' => 1, # 'ctg'
                                'R10E9.4' => 1, # 'ctg'
                                'R12B2.13' => 1, # 'ttg'
                                'T03G6.7' => 1, # 'ctg'
                                'T04C10.11' => 1, # 'ctg'
                                'T05A8.11' => 1, # 'ctg'
                                'T06D10.8' => 1, # 'ttg'
                                'T06E8.8' => 1, # 'ctg'
                                'T09F3.7' => 1, # 'ctg'
                                'T17H7.4m' => 1, # 'gtg'
                                'T19D2.13a' => 1, # 'gtg'
                                'T19D2.13b' => 1, # 'ctg'
                                'T23F2.13' => 1, # 'gtg'
                                'T24D1.7' => 1, # 'ttg'
                                'T28F3.11' => 1, # 'gtg'
                                'W02D3.13b' => 1, # 'ctg'
                                'W09G10.9' => 1, # 'atc' looks reasonable
                                'Y32F6B.9' => 1, # 'ttg'
                                'Y39G8B.13' => 1, # 'ttg'
                                'Y44E3A.7' => 1, # 'gtg'
                                'Y55B1BL.2' => 1, # 'gtg'
                                'Y62E10A.2' => 1, # 'ttg'
                                'Y75B8A.63' => 1, # 'ctg'
                                'Y92H12BL.6' => 1, # 'gtg'
                                'ZC328.8a' => 1, # 'ctg'
			       );






# check CDS/Transcripts/Pseudogenes
&quick_tests() if !$build; # tests for the curation databases - we assume that all these problems have been fixed by the time we do the Build

# main Gene checks
&main_gene_checks(); # tests for both curation databases and the Build databases

# single query checks - tests on whole classes
&single_query_tests();  # tests for both curation databases and the Build databases


# checks for the full Build
&extra_build_checks() if $build;

# protein checks
&protein_checks() if $build;





&print_results();

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



# check on CDS/Transcripts/Pseudogenes
sub  quick_tests {



  my $qclass;
  my @qclasses = ("Pseudogene", "Transcript", "CDS");
  my %seen;
  
  # generic tests for all classes
  my %queries = (
		 "has no Method"               => "where !Method",
		 "has no Source_exons"         => "where !Source_exons",
		 "has no S_parent"             => "where !S_parent",
		 "has no Species"              => "where !Species",
		 "has an unusual Species name" => "where Species != \"$speciesfn\"",
		 "should be removed"           => 'where "*iso*"',
		);
  
  
  # tests for just CDS
  my %CDS_queries = (
		     "has no CDS tag"                          => "where !CDS",
		     "has a bad Method tag"                    => '; Method != "history"; Method != "curated"; Method != "Transposon_CDS"; Method != "Genefinder"; Method != "not_confirmed_isoformer"; Method != "RNASEQ.Hillier"; Method != "mGene"; Method != "RNASEQ.Hillier.Aggregate"; Method != "isoformer"; Method != "jigsaw"; Method != "twinscan" ; Method != "ensembl" ; Method != "Trinity_CDS" ; Method != "genBlastG"',
		     "has no Corresponding_transposon tag"     => "where Method = Transposon_CDS; !Corresponding_transposon",
		     "has no Gene tag"                         => "where Method = Transposon_CDS AND !Gene",
		     "has no Gene tag for this CDS"            => "where Method = curated AND !Gene",
		     "has no From_laboratory tag"              => "where Method = Transposon_CDS AND !From_laboratory",
		     "has no From_laboratory tag for this CDS" => "where Method = curated AND !From_laboratory",
		    );
  

  # tests for just Transcript
  my %Transcript_queries = (
			    "has a bad Method tag"                       => '; Method != "Coding_transcript"; Method != "history_transcript"; Method != "*RNA*"; Method != "non_coding_transcript"; Method != "non_coding_transcript_isoformer" ; Method != "cufflinks"; Method != "TIGR_BEST"',
                            "has no Corresponding_CDS tag"               => '; Method = Coding_transcript AND !Corresponding_CDS',
                            "has no Gene tag"                            => '; Method = Coding_transcript AND !Gene',
                            "has no Gene tag for this ncRNA"             => '; Method = "*RNA*" AND !Gene',
			    "has no Transcript tag"                      => '; Method = "*RNA*" AND !Transcript',
                            "has no Transcript tag value"                => '; Method = "*RNA*" AND Transcript AND NOT NEXT',
                            "has no Brief_identification tag"            => '; Method = "*RNA*" AND !Brief_identification',
                            "has no From_laboratory tag for this ncRNA"  => '; Method = "*RNA*" AND !From_laboratory',
#                            "has no Database tag for this ncRNA"         => '; Method = "*RNA*" AND !Database',
#                            "has no DB_remark tag for this ncRNA"        => '; Method = "*RNA*" AND !DB_remark',
			   );
  

  # tests for just Pseudogene
  my %Pseudogene_queries = (
			    "has a bad Method tag"                        => '; Method != "history_pseudogene"; Method != "Pseudogene"; Method != tRNA_pseudogene; Method != Transposon_Pseudogene; Method != rRNA_pseudogene',
 #                           "has no DB_remark tag"                        => '; !DB_remark',
                            "has no From_laboratory tag"                  => '; !From_laboratory',
                            "has no Type tag"                             => '; !Type AND Method != "history_pseudogene"',
                            "has no Gene tag"                             => '; Method = Pseudogene AND !Gene',
                            "has no Gene tag for this Pseudogene_history" => '; Method = Pseudogene_history AND !Gene',
                            "has no Brief_identification tag"             => '; !Brief_identification AND Method != "history_pseudogene"',
			   );
  
  
  my $s;
  foreach $qclass (@qclasses) {
    $log->write_to("\nQuick test of $qclass\n");
    %seen = ();
    my %ignore;
    
    # the CDS and Transcript ghost objects with only an Evidence tag set are not interesting and should be ignored
    my @evidence = $db->fetch (-query => "FIND $qclass where Evidence AND !Method AND !Species AND !Sequence AND !Source_exons AND !Gene AND !Remark AND !Gene_history");
    foreach my $g (@evidence) {
      my $gg=$g->name; 
      $ignore{$gg}=1;
    }
    
    # the CDS and Transcript objects in the merged autoace from Tier III species are not interesting
    my @tier3 = $db->fetch (-query => "FIND $qclass where Database = WormBase AND NEXT = TierIII");
    foreach my $g (@tier3) {
      my $gg=$g->name; 
      $ignore{$gg}=1;
    }
    


    foreach my $query_text (keys %queries) {
      my $query = $queries{$query_text};
      if ($query eq 'where "*iso*"' && !$build) {next}
      do_a_quick_test($qclass, $query, $query_text, \%seen, \%ignore);
    }
    
    ##########################################
    # CDS/Transcript/Pseudogene - only tests
    
    if ($qclass eq "CDS") {
      foreach my $query_text (keys %CDS_queries) {
	my $query = $CDS_queries{$query_text};
	do_a_quick_test($qclass, $query, $query_text, \%seen, \%ignore);
      }      
      
    } elsif ($qclass eq "Transcript") {
      foreach my $query_text (keys %Transcript_queries) {
	my $query = $Transcript_queries{$query_text};
	do_a_quick_test($qclass, $query, $query_text, \%seen, \%ignore);
      }      
      
    } elsif ($qclass eq "Pseudogene") {
      foreach my $query_text (keys %Pseudogene_queries) {
	my $query = $Pseudogene_queries{$query_text};
	do_a_quick_test($qclass, $query, $query_text, \%seen, \%ignore);
      }      
      
    }

  }


}


sub do_a_quick_test {
  my ($qclass, $query, $query_text, $seen, $ignore) = @_;

  my @genes = $db->fetch (-query => "FIND $qclass $query");
  foreach my $g (@genes) {
    my $s = '';
    my $gg=$g->name;
    if ($gg =~ /^temp_gene/) {next}
    if ($ignore->{$gg}) {next}
    if (exists $seen->{$gg}) {
      $s=' (seen already)';
    }
    $log->write_to("ERROR: $qclass $gg $query_text $s\n"); 
    $seen->{$gg}=1;
  }
}

##########################################
# main gene checks

sub main_gene_checks {

    
  my %checked_small_intron = (
			    'R74.3a' => 1, # 23bp intron which switches frame to give an alternate 3' end to the protein. This isoform encodes the active form of xbp-1.
			    'R74.3c' => 1, # 23bp intron which switches frame to give an alternate 3' end to the protein.
			   );



  my @Models = $db->fetch (-query => 'Find All_genes where (Species = "'.$speciesfn.'")'); # CDS, Transcripts, Pseudogenes, etc

  my $model_count = scalar @Models;
  $log->write_to("\nFound $model_count $speciesfn Gene models to check.\n\n");



 CHECK_GENE:
  
  while (my $gene_model = shift @Models) {
    my $gene_model_name = $gene_model->name;

    unless (defined $gene_model->Method) {
      push(@error1,"ERROR: $gene_model appears to be incomplete - no Method\n");
      next CHECK_GENE;
    }

    my $method_name = $gene_model->Method->name;

    unless ($gene_model_name =~ /$cds_regex/) {
      push(@error1,"WARNING: The name of $gene_model_name is invalid\n") if (($method_name !~ /history|tRNA|Transposon|pre_miRNA|miRNA_primary_transcript|non_coding_transcript_isoformer/) && ($gene_model_name !~ /\S+.t\d+/));
    }
    
    my @exon_coord1 = sort by_number ($gene_model->get('Source_exons',1));
    my @exon_coord2 = sort by_number ($gene_model->get('Source_exons',2));
    
    # check for duplicated sequence names
    if (exists $sequence_names{$gene_model_name}) {
      my $class = $gene_model->class;
      push(@error1, "ERROR: $class $gene_model is both a $method_name and a $sequence_classes{$gene_model_name} $sequence_names{$gene_model_name}\n");
    }
    
    $sequence_names{$gene_model_name} = $method_name; # save all the names and methods
    $sequence_classes{$gene_model_name} = $gene_model->class;
    
    # check for duplicated sequence structures
    if ($method_name eq 'curated' || $method_name eq 'non_coding_transcript' || ($method_name =~ /RNA/ && $method_name !~ /pseudogene/ && $gene_model_name =~ /$cds_regex/)) {
      my ($gene_name) = ($gene_model_name =~ /($cds_regex_noend)/); # get the sequence name without an isoform letter at the end
      # make a hash key out of the exon starts and ends
      my $hash_key = join(':', @exon_coord1) . ',' . join(':', @exon_coord2);
      my $sequence = $gene_model->Sequence;
      my $class = $gene_model->class;
      my @clone_locations;
      if ($class eq 'CDS') {
	unless (defined $sequence->CDS_child) {
	  push(@error1, "ERROR: $class $gene_model has S_parent issues\n");
	  next;
	}
	@clone_locations = $sequence->CDS_child;
      } elsif ($class eq 'Transcript') {
	unless (defined $sequence->Transcript) {
	  push(@error1, "ERROR: $class $gene_model has S_parent issues\n");
	  next;
	} 
	@clone_locations = $sequence->Transcript;	    
      }
      my ($target_start, $target_end);
      my $found = 0;
      foreach my $target_location ( @clone_locations ) {
	next unless ($target_location->name eq $gene_model_name);
	$found = 1;
	$target_start = $target_location->right->name;
	if (!defined $target_start) {
	  push(@error1, "ERROR: $class $gene_model_name has no start position in the Sequence object\n");
	}
	$target_end = $target_location->right->right->name;
	if (!defined $target_end) {
	  push(@error1, "ERROR: $class $gene_model_name has no end position in the Sequence object\n");
	}
	last;
      }
      if (!$found) {
	push(@error1, "ERROR: $class $gene_model_name was not found in the Sequence object\n");
      }
      $hash_key = $sequence->name . ':' . $target_start . ':' . $target_end . ':' . $hash_key; # NB appending the exons hash_key made above
      if (exists $sequence_structures{$hash_key}) {
	my $other_isoform = $sequence_structures{$hash_key};
	my $class = $gene_model->class;
	unless ($gene_model =~ (/F59A3.6/)) { # F59A3.6e and F59A3.6c are identical duplicated isoforms in different locations in this complex locus.
	  # see if the two similar structures are at the same position on the same clone
	  push(@error1, "ERROR: $class $gene_model has the same structure as $other_isoform\n");
	}
      }
      $sequence_structures{$hash_key} = $gene_model_name;
    }
    
    #Check that the source_exons and the the span are the same size.
    
    
    
    unless (($method_name eq 'Transposon') || ($method_name eq 'history_transposon')) {
      if (!defined($exon_coord2[0])) {
	push(@error1, "ERROR: $gene_model has a problem with its exon co-ordinates\n");
	next;
      }
      if (($exon_coord2[0] < "1") && ($method_name eq 'curated')){
	push(@error1, "ERROR: $gene_model has a problem with its exon co-ordinates\n");
	next;
      }
    }
    
    if (!exists $checked_small_intron{$gene_model_name}) {
      for (my $i=1; $i<@exon_coord2; $i++) {
	my $intron_size = ($exon_coord1[$i] - $exon_coord2[$i-1] -1); 
	if (($intron_size < 25) && ($method_name eq 'curated')) {
	  my $artificial_intron_Ribosomal_slippage = $gene_model->get('Ribosomal_slippage');
	  my $artificial_intron_Low_quality_sequence = $gene_model->get('Low_quality_sequence');
	  if ($intron_size < 5 && !$artificial_intron_Ribosomal_slippage && !$artificial_intron_Low_quality_sequence) {
	    push(@error4,"ERROR: $gene_model has a very small intron ($intron_size bp) and no Ribosomal_slippage or Low_quality_sequence tag\n");
	  } elsif (!$artificial_intron_Ribosomal_slippage && !$artificial_intron_Low_quality_sequence) {
	    push(@error4,"Warning: $gene_model has a small intron ($intron_size bp) and no Artificial_intron tag\n") if ($species eq 'elegans');
	  }
	}
      }
    }
    
    for (my $i=0; $i<@exon_coord1; $i++) {
      my $start = $exon_coord1[$i];
      my $end = $exon_coord2[$i];
      for (my $j = $i+1; $j<@exon_coord1; $j++) {
	if (($end > $exon_coord1[$j]) && ($start < $exon_coord2[$j])) {
	  push(@error1,"ERROR: $gene_model exon inconsistency, exons overlap\n") if ($method_name !~ /history/);
	}
      }
    }
    
    if ($method_name eq 'curated' && check_sequence_span(\@exon_coord1, \@exon_coord2, $gene_model)) {
      push(@error1,"ERROR: $gene_model exon inconsistency, not the same span as the Sequence span\n");      
    }
    
    # check that 'Start_not_found' and 'End_not_found' tags present? (CDS specific.....extended to all genes :) )
    my $start_tag = "";
    my $start_tag_val = "";
    my $end_tag = "";
    my $genetic_code = "";
    
    if ($gene_model->get('Start_not_found')) {
      $start_tag = "present";
      if ($gene_model->get('Start_not_found')->right) {
        $start_tag_val = $gene_model->get('Start_not_found')->right->name;
      } else {
	push(@error4,"Warning: $gene_model Start_not_found tag present but does not contain a value assuming should be 1\n");      
	$start_tag_val = 1;
      }
      # Throw an error only for elegans curated CDS with Start_not_found that are not in our list of known non-canonical initiation sites 
      if ($species eq 'elegans' && $method_name eq 'curated' && !exists $non_canonical_initiation{$gene_model_name}) {
	push(@error2,"ERROR: $gene_model Start_not_found tag present\n");	  
      }
    }
    
    if ($gene_model->get('End_not_found')) {
      $end_tag = "present";
      if ($species eq 'elegans' && $method_name eq 'curated') {
	push(@error2,"ERROR: $gene_model End_not_found tag present\n");
      }
    }
    
    if ($gene_model->get('Genetic_code')) {
      $genetic_code = $gene_model->get('Genetic_code')->right->name;
    }
    
    #All Isoforms should have the Isoform tag set. cds_regex
    
    if (($gene_model_name =~  (/$cds_regex/)) && ($method_name !~  /history/)) {
      my $Isoform = $gene_model->at('Properties.Isoform');
      
      if ($gene_model_name =~  (/\S+[a-z]$/)) {
	push(@error3, "ERROR: $gene_model [$gene_model_name] requires an Isoform\n") unless (defined $Isoform);
      }
      if ($gene_model_name =~  (/S+\d$/)) {
	push(@error3, "ERROR: $gene_model [$gene_model_name] has Isoform tag set but no isoform letter in sequence name\n") if (defined $Isoform);
      }
    }
    
	


    
    ###################################
    #All gene predictions should have #
    ###################################
    
    # check that 'Sequence' tag is present and if so then grab parent sequence details
    my $source;
    my $source_type;
    if (!defined($gene_model->Sequence)){
      push(@error1,"ERROR: $gene_model has no Sequence, cannot check DNA\n");
      next CHECK_GENE;
    }
    
    # check that history genes have a history method.
    if ($method_name !~ /history/ && $gene_model_name =~ /$cds_regex\:\w+/) {
      push(@error3, "ERROR: $gene_model history object doesn't have a history method.\n");
    }
    
    # check that history genes are renamed.
    if (($method_name =~ /history/ && !($gene_model_name =~ /\:/) && !($gene_model_name =~ /WBTransposon/))) {
      push(@error3, "ERROR: $gene_model needs to be renamed as it is part of history.\n");
    }
    
    if ($method_name eq "Transposon") {
      next CHECK_GENE;
    }
    
    #Gene ID checks.
    unless ($method_name eq "Transposon_CDS") {
      my $Gene_ID     = $gene_model->at('Visible.Gene.[1]');
      my $Genehist_ID = $gene_model->at('Visible.Gene_history.[1]');
      
      # Can't have both Gene and Gene_history.
      if ((defined $Genehist_ID) && (defined $Gene_ID)) {
	push(@error2, "ERROR: Gene Model $gene_model contains both a Gene and a Gene_history tag, Please fix.\n");
      }
      
      #History genes have to have a Gene_history ID of 8 digits.
      if ($gene_model_name =~ (/$cds_regex\:\w+/)) {
	if (defined $Genehist_ID) {
	  push(@error2, "ERROR: The Gene ID '$Genehist_ID' in $gene_model is invalid!\n") unless ($Genehist_ID =~ /WBGene[0-9]{8}/);
	}  else {
	  push(@error2, "ERROR: $gene_model does not have the Gene_history populated\n");
	}
	if (defined $Gene_ID) {
	  push(@error2, "ERROR: $gene_model should not contain a Gene ID under Gene\n");
	}
      }
      
      #non history
      elsif ($gene_model_name =~ (/$cds_regex/)) {
	if (defined $Gene_ID) {
	  push(@error2, "ERROR: The Gene ID '$Gene_ID' in $gene_model is invalid!\n") unless ($Gene_ID =~ /WBGene[0-9]{8}/);
	}  else {
	  push(@error2, "ERROR: $gene_model does not have a Gene ID!\n") unless ($method_name eq 'Transposon_Pseudogene');
	}
	if (defined $Genehist_ID) {
	  push(@error2, "ERROR: $gene_model should not contain a Gene ID under Gene_history\n");
	}
      }

      #####################################################################################################

      # then run misc. sequence integrity checks
      my $dna = $gene_model->asDNA();
      if (!$dna) {
	push(@error1,"ERROR: $gene_model can't find any DNA to analyse\n");
	next CHECK_GENE;
      }
      
      # feed DNA sequence to function for checking
      &test_gene_sequence_for_errors($gene_model,$start_tag,$start_tag_val,$end_tag,$dna,$method_name,$genetic_code);
    }
  }

  # now check that the isoform names are consistent
  foreach my $sequence_name (keys %sequence_names) {
    if (! defined $sequence_names{$sequence_name}) {
      #push(@error1, "ERROR: The $sequence_classes{$sequence_name} '$sequence_name' has no Method\n");
      next;
    }
    # don't want to look at history objects
    if ($sequence_names{$sequence_name} =~ /history/) {next}
    
    # does the isoform have multiple letters in
    if ($sequence_name =~ /\w+\.\d+([a-z]{2,})$/) {
      push(@error1, "ERROR: The sequence_name '$sequence_name' is invalid! Multiple letters in the isoform name.\n")
    }
    # if it is an isoform name, check for non-isoforms 
    elsif ($sequence_name =~ /(\w+\.\d+)[a-z]$/) {
      my $base = $1;
      if (exists $sequence_names{$base}) {
	if (($sequence_names{$base} eq 'miRNA_primary_transcript' || $sequence_names{$base} eq 'pre_miRNA') && (($sequence_names{"${base}a"} && $sequence_names{"${base}a"} eq 'miRNA') || ($sequence_names{"${base}b"} && $sequence_names{"${base}b"} eq 'miRNA'))) {
	  next;
	} 
	# ignore the primary and mature miRNA forms
	push(@error1, "ERROR: The $sequence_names{$base} sequence '$base' and the $sequence_names{$sequence_name} sequence '$sequence_name' both exist!\n")
      }
    }
  }


}

##########################################
# single query checks - tests on whole classes

sub single_query_tests {


  # Check for non-standard Methods in CDS class
  my @CDSfilter;
  
  if ($build){
    @CDSfilter = $db->fetch (-query => 'FIND CDS; Method != Transposon_CDS; Method != Transposon_Pseudogene; Method != curated; Method != history');
  } else {
    @CDSfilter = $db->fetch (-query => 'FIND CDS; Method != Transposon_CDS; Method != Transposon_Pseudogene; Method != curated; Method != history; Method != Genefinder; Method != twinscan; Method != jigsaw; Method != mGene; Method != RNASEQ.Hillier; Method != RNASEQ.Hillier.Aggregate; Method != cufflinks*; Method != genBlastG; Method != *isoformer; Method != ensembl');
  }

  
  foreach my $CDSfilter (@CDSfilter) {
    push(@error4, "ERROR! CDS:$CDSfilter contains an invalid Method please check\n");
  }


}

##########################################
# extra build_checks for when we have a full Build

sub extra_build_checks {

  my @genes = $db->fetch (-query => 'Find Gene where Live AND Species = "$speciesfn" AND Sequence_name AND NOT Corresponding_CDS AND NOT Corresponding_pseudogene AND NOT Corresponding_transcript');
  foreach my $gene (@genes) {
      push(@error1, "ERROR! Gene:$gene is Live but not attached to a current gene model\n\n");
  }

  my @cds = $db->fetch (-query => 'FIND CDS; Method = curated AND NOT Corresponding_CDS');
  foreach my $cds (@cds) {
      push(@error1, "ERROR! CDS:$cds is Method=curated but does not have a Corresponding_protein\n");
  }

}

##########################################
# Protein checks

sub protein_checks {

  my @proteins = $db->fetch (-query => 'Find Protein where Live AND NOT Peptide');
  foreach my $protein (@proteins) {
      push(@error1, "ERROR! Protein:$protein is Live but does not have a Peptide\n");
  }
  


}


##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################


#################################################################
# results
#################################################################

sub print_results {

  # print warnings to log file, log all category 1 errors, and then fill up.

  my $count_errors = 0;
  my @error_list = ( \@error1, \@error2, \@error3,  \@error4, \@error5);
  foreach my $list (@error_list) {
    foreach my $error (@{$list}) {
      $count_errors++;
      if ($error =~ /ERROR/i) {$log->error;}
      if ($count_errors > 500 && ! $verbose) {
	$log->write_to("More than 500 lines output... listing truncated ...\n");
	last;
      }
      $log->write_to("$count_errors $error");
    }
  }
}

#################################################################
# Subroutines
#################################################################

sub test_gene_sequence_for_errors {
  my $gene_model = shift;
  my $start_tag = shift;
  my $start_tag_val = shift;
  my $end_tag = shift;
  my $dna = shift;
  my $method_name =shift;
  my $genetic_code=shift;


  my %checked_small_genes = (
			     'B0035.22' => 1,
			     'B0041.12' => 1,
			     'B0250.18b' => 1,
			     'B0302.4' => 1,
			     'B0348.2b' => 1,
			     'B0513.2b' => 1,
			     'B0513.3b' => 1,
			     'B0546.4c' => 1,
			     'BE0003N10.6e' => 1,
			     'C01B4.7d' => 1,
			     'C01F6.9b' => 1,
			     'C01F6.9c' => 1,
			     'C04G2.1b' => 1,
			     'C08B6.18' => 1,
			     'C08H9.18' => 1,
			     'C09E8.6' => 1,
			     'C12D8.21b' => 1,
			     'C13G3.12' => 1,
			     'C16A11.5e' => 1,
			     'C18D4.12' => 1,
			     'C32B5.18' => 1,
			     'C33A12.19b' => 1,
			     'C33A12.3b' => 1,
			     'C33H5.13c' => 1,
			     'C34E10.12b ' => 1,
			     'C47A10.16' => 1,
			     'C48A7.16' => 1,
			     'C49F5.13' => 1,
			     'C51E3.10b' => 1,
			     'D1046.18' => 1,
			     'D2063.4c' => 1,
			     'F01E11.20' => 1,
			     'F08A7.1b' => 1,
			     'F08B12.11' => 1,
			     'F08G5.5c' => 1,
			     'F11A6.15' => 1,
			     'F13G3.10c' => 1,
			     'F13H10.8b' => 1,
			     'F19G12.9' => 1,
			     'F20C5.2f' => 1,
			     'F20D12.12' => 1,
			     'F21C10.17' => 1,
			     'F22F4.4b' => 1,
			     'F23B2.13b' => 1,
			     'F23C8.14b' => 1,
			     'F27B10.1c' => 1,
			     'F27D4.4b' => 1,
			     'F28H7.10b' => 1,
			     'F32F2.1g' => 1,
			     'F35C8.9' => 1,
			     'F36A4.5b' => 1,
			     'F36H9.4c' => 1,
			     'F37C4.8b' => 1,
			     'F37H8.6' => 1,
			     'F38C2.1d' => 1,
			     'F40F12.9c' => 1,
			     'F44D12.16' => 1,
			     'F47G6.1c' => 1,
			     'F52E10.4d' => 1,
			     'F52G3.7' => 1,
			     'F54D7.6' => 1,
			     'F54D7.7' => 1,
			     'F55A4.13' => 1,
			     'F55H12.7b' => 1,
			     'F56H9.2b' => 1,
			     'F59A6.1' => 1,
			     'H03A11.13 ' => 1,
			     'H12I19.115' => 1,
			     'H24K24.4e' => 1,
			     'H40L08.8' => 1,
			     'K02E2.11' => 1,
			     'K07C5.13a' => 1,
			     'K08D9.10' => 1,
			     'K08E5.5' => 1,
			     'K08F8.15' => 1,
			     'K10C9.8c' => 1,
			     'K11H12.8d' => 1,
			     'LLC1.2b' => 1,
			     'M02E1.3' => 1,
			     'M142.8b' => 1,
			     'M176.16' => 1,
			     'R07B7.12b' => 1,
			     'R11D1.12b' => 1,
			     'T03G6.7' => 1,
			     'T04H1.13' => 1,
			     'T06A1.7b' => 1,
			     'T09F5.20a' => 1,
			     'T11F9.21b' => 1,
			     'T12G3.11' => 1,
			     'T13A10.64' => 1,
			     'T21H3.4c' => 1,
			     'T23F6.3c' => 1,
			     'T24A6.1b' => 1,
			     'T24B8.3c' => 1,
			     'T24B8.4d' => 1,
			     'T26H5.11' => 1,
			     'T28F3.11' => 1,
			     'W02A2.6c' => 1,
			     'W02D3.13b' => 1,
			     'W03F8.6c' => 1,
			     'W03F8.6d' => 1,
			     'W08A12.2b' => 1,
			     'W09G3.8b' => 1,
			     'W10D9.5b' => 1,
			     'Y105C5A.1285' => 1,
			     'Y10G11A.90' => 1,
			     'Y113G7A.22' => 1,
			     'Y116A8C.28e' => 1,
			     'Y17D7B.10' => 1,
			     'Y19D10A.4d' => 1,
			     'Y39G8B.13' => 1,
			     'Y40B1A.1b' => 1,
			     'Y41D4B.13c' => 1,
			     'Y41E3.1d' => 1,
			     'Y46G5A.47' => 1,
			     'Y47D7A.17' => 1,
			     'Y50D7A.13b' => 1,
			     'Y51H4A.938a' => 1,
			     'Y54F10AM.15' => 1,
			     'Y54G2A.28c' => 1,
			     'Y57G11C.1142' => 1,
			     'Y62E10A.12b' => 1,
			     'Y62E10A.8b' => 1,
			     'Y65A5A.24' => 1,
			     'Y66H1A.8a' => 1,
			     'Y67H2A.10c' => 1,
			     'Y69A2AR.24a' => 1,
			     'Y69A2AR.50' => 1,
			     'Y69A2AR.51' => 1,
			     'Y71G12A.5' => 1,
			     'Y73C8B.9a' => 1,
			     'Y73C8B.9b' => 1,
			     'Y77E11A.6b' => 1,
			     'ZC328.8a' => 1,
			     'ZK180.8a' => 1,
			     'ZK377.15' => 1,
			     'ZK381.62' => 1,
			     'ZK792.12' => 1,
			     'ZK863.7b' => 1,
			     'ZK897.10a' => 1,
			     'ZK897.1f' => 1,
			    );


  my $gene_model_name = $gene_model->name;
  
  if ($start_tag eq "present") {
    unless (defined $start_tag_val) {
      $start_tag_val = '1';
    }
    if ($start_tag_val eq 0) {
      $start_tag_val = '1';
    }
  }
  
  if ($method_name ne 'history'){
    # trim DNA sequence to just A,T,C,G etc.
    $dna =~ s/\n//g;
    my $length_gene_name = length($gene_model)+1;
    $dna = substr($dna, $length_gene_name);
    if (!$dna){
      push(@error1, "ERROR: $gene_model has a problem with it's DNA connection.\n");
      next CHECK_GENE;
    }
    # calculate other necessary values
    my $gene_model_length = length($dna);
    my $remainder;
    
    if (($gene_model->Method eq 'curated') && ($gene_model->Start_not_found)) {
      my $extra = $gene_model->Start_not_found->name;
      my $length_calc = $gene_model_length + $extra;
      $remainder = $length_calc%3;
    } else {
      $remainder = $gene_model_length%3;
    }

    my $start_codon = substr($dna,0,3);
    my $stop_codon = substr($dna,-3,3);   
    
    # check for length errors(CDS specific)
    

    my $warning;
      
    if (!exists $checked_small_genes{$gene_model}) {
      if (($gene_model_length < 35) && ($method_name eq 'curated')) {
	$warning = "WARNING: $gene_model is very short ($gene_model_length bp),";
	if (defined($gene_model->at('Properties.Coding.Confirmed_by'))) {
	  $warning .= "gene is Confirmed\n";
	} elsif (defined($gene_model->at('Visible.Matching_cDNA'))) {
	  $warning .= "gene is Partially_confirmed\n";
	} else {
	  $warning .= "gene is Predicted\n";
	}
	push(@error3, $warning);
      } elsif (($gene_model_length < 55) && ($method_name eq 'curated')) {
	if (defined($gene_model->at('Properties.Coding.Confirmed_by'))) {
	  $warning = "WARNING: $gene_model is short ($gene_model_length bp) and is Confirmed\n";
	}
	elsif (defined($gene_model->at('Visible.Matching_cDNA'))) {
	  $warning .= "WARNING: $gene_model is short ($gene_model_length bp) and is Partially_confirmed\n";
	} else {
	  $warning .= "WARNING: $gene_model is short ($gene_model_length bp) and is Predicted\n";
	}
	push(@error5, $warning);
      }
    }

      
    # Is the gene prediction complete?
    if (($remainder != 0) && ($method_name eq 'curated')) {
      if (($end_tag ne "present") && ($start_tag ne "present")) {
	push(@error1,"ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n");
      } elsif (($start_tag ne "present") && (!exists $non_canonical_initiation{$gene_model_name})) {
	push(@error1,"ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found tag MISSING and not a known non-canonical start codon gene\n");
      }
    }
    if (defined $gene_model->Sequence) {
      unless (($gene_model_name =~ /MTCE/) || ($gene_model->Sequence->name =~ /MITO/) || ($gene_model->Sequence->name =~ /Cbre_Contig2951/) || ($gene_model->Sequence->name =~ /Cbre_Contig2389/) || ($gene_model->Sequence->name =~ /Cbre_Contig1932/)) {
	# look for incorrect stop codons (CDS specific)
	if (($stop_codon ne 'taa') && ($stop_codon ne 'tga') && ($stop_codon ne 'tag') && ($method_name eq 'curated')) {
	  if ($end_tag ne "present") {
	    push(@error1, "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n");
	  } else {
	    push(@error2,"ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag present\n");
	  }
	}
	# look for incorrect start codons(CDS specific)
	if (exists $non_canonical_initiation{$gene_model_name}) {
	  print "INFORMATION: $gene_model  intentionally utilises a novel '$start_codon' start codon......Ignoring\n" if $verbose;
	}
	if (($start_codon ne 'atg') && ($method_name eq 'curated') && ($start_tag ne "present") && (!exists $non_canonical_initiation{$gene_model_name})) {
	  if (($start_tag ne "present")) {
	    push(@error1,"ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n");
	  } else {
	    push(@error2, "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag present\n");
	  } 
	}
	
	# check for internal stop codons (CDS specific)
	my $i;
	my $j;
	if ($start_tag_val) {
	  print "start_tag_val final = $start_tag_val\n" if ($verbose);
	  if ($start_tag_val eq 3) {#if 3 +2
	    $i = '2';
	  } elsif ($start_tag_val eq 2) {#if 2 +1
	    $i = '1';
	  } else {
	    $i = '0';
	  }
	} else {
	  $i = '0';
	}
	
	for (;$i<$gene_model_length-3;$i+=3) {
	  # hold position of codon in $j
	  $j=$i+1;
	  my $codon =substr($dna,$i,3);
	  if (($codon eq "taa") || ($codon eq "tag") || ($codon eq "tga")) {      
	    my $previous_sequence = substr($dna, $j-11,10);
	    my $following_sequence = substr($dna, $j+2, 10);
	    my $offending_codon = substr($dna, $j-1, 3);
	    if (($method_name eq 'curated')) {
	      push(@error1, "ERROR: $gene_model internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n") unless ($genetic_code eq 'Selenocysteine');      
	    }
	  }
	}
	if ($species eq "elegans") {
	  # look for non-ACTG characters
	  if ($dna =~ /[^acgt]/i) {
	    $dna =~ s/[acgt]//g;
	    push(@error2, "ERROR: $gene_model DNA sequence contains the following non-ATCG characters: $dna\n"); 
	  }
	}
      } ##MTCE exclusion loop
    } else {  ##Sequence_parent loop
      push(@error1, "ERROR: $gene_model has no SParent Sequence connection.\n") unless (defined $gene_model->Transposon);
    }
  }
}
############################################################################################
# check that the span of the exons is the same size as the span of the position in the Sequence
sub check_sequence_span {
  my ($exon_coord1, $exon_coord2, $gene_model) = @_;
  
  my $result=0;
  
  my $exon_start = $exon_coord1->[0];
  my $exon_end = $exon_coord2->[$#{$exon_coord2}];
  my $exon_span = $exon_end - $exon_start + 1;
  
  # get the Sequence start-end span 
  my $target = $gene_model->name;
  my $sequence = $gene_model->Sequence;
  if (!defined $sequence) {
    push(@error1, "ERROR: $gene_model has no SParent Sequence connection.\n");
    return $result;
  }
  my @clone_locations = $sequence->CDS_child;
  my ($target_start, $target_end);
  foreach my $target_location ( @clone_locations ) {
    next unless ($target_location->name eq $target);
    $target_start = $target_location->right->name;
    $target_end = $target_location->right->right->name;
    last;
  }
  if (defined $target_start && defined $target_end) {
    my $sequence_span = abs($target_start - $target_end) + 1;
    if ($exon_span != $sequence_span) {
      $result = 1;
      return $result;
    }
  } else {
    push(@error1, "ERROR: $gene_model has a malformed SParent Sequence entry.\n");
    return $result;    
  }
  
}


############################################################################################
sub by_number{ $a <=> $b;}
############################################################################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut

#!/software/bin/perl -w
#
# check_predicted_genes.pl
#
# by Keith Bradnam
#
# Last updated on: $Date: 2014-11-25 15:03:53 $
# Last updated by: $Author: mh6 $
#
# see pod documentation at end of file for more information about this script

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use IO::Handle;
use Getopt::Long;
use Log_files;
use Storable;

my ($verbose, $db_path, $basic, $test1, $debug, $store, $test,$build,$species,$incomplete,$nogenome);

GetOptions ("verbose"    => \$verbose, # prints screen output and checks the CDS class instead of All_genes.
	    "database=s" => \$db_path, # Path to the database you want to check.
	    "basic"      => \$basic,   # Ignores some of the checks.
	    "debug:s"    => \$debug,   # turns on more printing and errorlogging
	    "test1:s"    => \$test1,   # only checks the CDSs from 1 clone (elegans) or 1 CDS in other TierII databases.
	    "store:s"    => \$store,   # 
	    "test"       => \$test,    # Test build
	    "build"      => \$build,   # Checks specific to a full database containing genes and models.
	    "species:s"   => \$species,  # used to hold briggsae/brenneri/remanei for some checks.
	    "incomplete" => \$incomplete, # used to avoid start/end not found warnings
	    "nogenome"  => \$nogenome, # for testing annotations when you don't have the genome loaded.
	   );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species
			   );
}

unless (defined $species){$species = 'elegans';}

#Establish Log
my $log = Log_files->make_build_log($wormbase);

# Establish database connection etc.
$log->log_and_die("Please use -database <path> specify a valid database directory.\n\n") if (!defined($db_path));
# Specify which tace to use if you are using -program flag
my $tace = $wormbase->tace;
my $db = Ace->connect(-path=>$db_path) or  $log->log_and_die("Couldn't connect to $db_path\n". Ace->error);

# create separate arrays for different classes of errors (1 = most severe, 4 = least severe)
our (@error1, @error2, @error3, @error4, @error5,);
my $errorcountCDS = "0";

#Permitted exceptions where the error has been checked.
#List of verified small genes or have confirmed small introns.
my @checkedgenes = ('F56H11.3b','F13E9.8','Y45F10D.7','T21C12.8','F58G11.2','F54E4.3','F43D9.3a','F43D9.3b','F36H1.3','F02C12.1','C40H5.1', 'H12D21.1', 'H12D21.12','H12D21.13', 'H12D21.14',  'H12D21.15','W06A7.5','Y50E8A.17', 'ZC412.6', 'ZC412.7', 'R74.3a');


my $cds_regex = $wormbase->cds_regex;
my $cds_regex_noend = $wormbase->seq_name_regex;
my $speciesfn = $wormbase->full_name;

################################
#         Main Body            # 
################################

my %sequence_names;
my %sequence_classes;
my %sequence_structures;

if ($build) {
  &extra_build_checks;
}
else {
  # Fetch Gene Models to be tested (All_genes) #

  my @Predictions;
  my @bad_genes;
  my $bad_genes;
  my $Bcount;
  if ($test1) {
    $log->write_to("Only checking genes on ${test1}.......\n");
    @Predictions = $db->fetch (-query => "FIND All_genes WHERE Sequence = ${test1} AND method");
    foreach my $Predictions(@Predictions) {
      print $Predictions->name."\n";
    }
    $log->write_to("Checking for Gene models with no method");
    @bad_genes = $db->fetch (-query => 'FIND All_genes ${test1}* AND !method');
    $Bcount = @bad_genes;
    if ($Bcount ne '0') {
      $log->write_to("\nError: $Bcount models have no method, please check\n");
      foreach $bad_genes(@bad_genes) {
	$log->write_to("Error: $bad_genes has no method\n") if ($verbose);
      }
    }
    else {
      $log->write_to("\nThat's ok, $Bcount models have no method :)\n\n");
    }
  }
  else {

    # Quick tests
######################################################################################
    my $qclass;
    my @qclasses = ("CDS","Pseudogene","Transcript");
    foreach $qclass (@qclasses) {
      my @tmpbad_genes = $db->fetch (-query => "FIND $qclass where !method");
      splice @bad_genes,0,0,@tmpbad_genes;
      my @no_se = $db->fetch (-query => "FIND $qclass where !Source_exons");
      splice @bad_genes,0,0,@no_se;
      my @no_Sparent_genes = $db->fetch (-query => "FIND $qclass where !S_parent");
      splice @bad_genes,0,0,@no_Sparent_genes;
      my @no_species = $db->fetch (-query => "FIND $qclass where !Species");
      splice @bad_genes,0,0,@no_species;
      #my @no_lab = $db->fetch (-query => "FIND $qclass where !From_laboratory");
      #splice @bad_genes,0,0,@no_lab;
    }
    my @no_CDS = $db->fetch (-query => "FIND CDS where !CDS");
    splice @bad_genes,0,0,@no_CDS;


    my %h = map {($_,1)} @bad_genes;
    my @bad_genesuniq = keys %h;

    $Bcount = @bad_genesuniq;
    if ($Bcount ne '0') {
      $log->write_to("Error: $Bcount models have a quicktest issue (No method/Source_exons/S_parent/Species/CDS), please check\n");
      print "Error: $Bcount models have a quicktest issue (No method/Source_exons/S_parent/Species/CDS), please check\n" if ($debug);
      foreach $bad_genes(@bad_genesuniq) {
	$log->write_to("Error: $bad_genes has one of the above issues\n"); #if ($verbose);
	print "Error: $bad_genes has one of the above issues\n" if ($debug);
      }
    }
    else {
      $log->write_to("\nGreat $Bcount models have basic issues :)\n\n");
      print "Great $Bcount models have basic issues :)\n" if ($debug);
    }
    print STDERR "Fetching all genes...\n" if $verbose;
    #@Predictions = $db->fetch (-query => 'FIND All_genes where method AND Species = "$speciesfn"');
    @Predictions = $db->fetch (-query => 'Find All_genes where (Species = "'.$wormbase->full_name.'")');
  }
  
  my $gene_model_count=@Predictions;
  $log->write_to("Checking $gene_model_count Predictions in $db_path\n\n");
  print "\nChecking $gene_model_count Predictions...\n\n";
  &main_gene_checks(\@Predictions);

  print "\nDoing single query tests\n";
  &single_query_tests;

  # print warnings to log file, log all category 1 errors, and then fill up.
  $log->write_to("WARNING: $errorcountCDS incomplete CDSs found, please check the database for CDS !method.\n\n") unless $errorcountCDS eq '0';
  my $count_errors = 0;
  my @error_list = ( \@error1, \@error2, \@error3,  \@error4, \@error5);
  foreach my $list (@error_list) {
    foreach my $error (@{$list}) {
      $count_errors++;
      if ($error =~ /ERROR/) {$log->error;}
      $log->write_to("$count_errors $error");
      last if $count_errors > 1900;
    }
  }
}

$db->close;
$log->mail;
exit(0);


#################################################################
# Main Subroutines
#################################################################

sub main_gene_checks {
  my ($pred_ref) = @_;
  
 CHECK_GENE:
  
  while (my $gene_model = shift @$pred_ref) {
    my $gene_model_name = $gene_model->name;
     print "Checking $gene_model_name\n" if ($verbose);

    unless (defined $gene_model->Method) {
      $errorcountCDS ++;
      print "$gene_model appears to be incomplete: it has no method.\n" if $verbose;
      push(@error1,"ERROR: $gene_model appears to be incomplete\n");
      next CHECK_GENE;
    }
    my $method_test = $gene_model->Method->name;
    unless ($method_test =~ /ransposon/) {
      unless (defined $gene_model->Source_exons) {
	print "$gene_model appears to be incomplete: it has no Source_exons.\n" if $verbose;
	push(@error1,"ERROR: $gene_model appears to be incomplete: it has no Source_exons.\n");
	next CHECK_GENE;
      }
    }
    unless ($gene_model_name =~ /$cds_regex/) {
      push(@error1,"warning $gene_model_name invalid\n") if (($method_test !~ /history|tRNA|Transposon|miRNA_primary_transcript/) && ($gene_model_name !~ /\S+.t\d+/));
    }
    my @exon_coord1 = sort by_number ($gene_model->get('Source_exons',1));
    my @exon_coord2 = sort by_number ($gene_model->get('Source_exons',2));
    my $i;
    my $j;

  
    # check for duplicated sequence names
    if (exists $sequence_names{$gene_model_name}) {
      my $class = $gene_model->class;
      push(@error1, "ERROR: $class $gene_model is both a $method_test and a $sequence_classes{$gene_model_name} $sequence_names{$gene_model->name}\n") unless ($method_test eq 'Coding_transcript');
      print "ERROR: $class $gene_model is both a $method_test and a $sequence_classes{$gene_model->name} $sequence_names{$gene_model->name}\n" unless ($method_test eq 'Coding_transcript');
    }
    $sequence_names{$gene_model->name} = $method_test; # save all the names and methods
    $sequence_classes{$gene_model->name} = $gene_model->class;

    # check for duplicated sequence structures
    if ($method_test eq 'curated' || $method_test eq 'non_coding_transcript') {
      my ($gene_name) = ($gene_model_name =~ /($cds_regex_noend)/);
      # make a hash key out of the exon starts and ends
      my $hash_key = join(':', @exon_coord1) . ',' . join(':', @exon_coord2);
      if (exists $sequence_structures{$gene_name}{$hash_key}) {
	my $other_isoform = $sequence_structures{$gene_name}{$hash_key};
	my $class = $gene_model->class;
	unless ($gene_model =~ (/F59A3.6/)) {
	  push(@error1, "ERROR: $class $gene_model has the same structure as $other_isoform\n");
	  print "ERROR: $class $gene_model has the same structure as $other_isoform\n";
	}
      }
      $sequence_structures{$gene_name}{$hash_key} = $gene_model->name;
    }

    #Check that the source_exons and the the span and the same size.

    

    unless (($method_test eq 'Transposon') || ($method_test eq 'history_transposon')) {
      if (!defined($exon_coord2[0])) {
	print "ERROR: $gene_model has a problem with it\'s exon co-ordinates\n";
	push(@error1, "ERROR: $gene_model has a problem with it\'s exon co-ordinates\n");
	next;
      }
      if (($exon_coord2[0] < "1") && ($method_test eq 'curated')){
	push(@error1, "ERROR: $gene_model has a problem with it\'s exon co-ordinates\n");
	print "ERROR: $gene_model has a problem with it\'s exon co-ordinates\n";
	next;
      }
    }

    for ($i=1; $i<@exon_coord2;$i++) {
      my $intron_size = ($exon_coord1[$i] - $exon_coord2[$i-1] -1); 
      my @ck = grep(/^$gene_model/, @checkedgenes);
      push(@ck,"0"); #this gets rid of undef issues.
      unless (defined $ck[0]) {
	if (($intron_size < 34) && ($method_test eq 'curated')) {
	  push(@error4,"ERROR: $gene_model has a very small intron ($intron_size bp)\n");
	}
	push(@error5,"WARNING: $gene_model has a small intron ($intron_size bp)\n") if (($intron_size > 33) && ($intron_size < 39) && (!$basic) && ($method_test eq 'curated') && ($ck[0] ne $gene_model_name));
      }
    }

    for ($i=0; $i<@exon_coord1; $i++) {
      my $start = $exon_coord1[$i];
      my $end = $exon_coord2[$i];
      for ($j=$i+1;$j<@exon_coord1;$j++) {
	if (($end > $exon_coord1[$j]) && ($start < $exon_coord2[$j])) {
	  print "ERROR: $gene_model exon inconsistency, exons overlap\n" if $verbose;
	  push(@error1,"ERROR: $gene_model exon inconsistency, exons overlap\n") if ($method_test !~ /history/);
	}
      }
    }

    if ($method_test eq 'curated' && check_sequence_span(\@exon_coord1, \@exon_coord2, $gene_model)) {
      print "ERROR: $gene_model exon inconsistency, not the same span as the Sequence span\n" if $verbose;
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
      }
      else {
	print "Warning: $gene_model Start_not_found tag present but does not contain a value assuming should be 1\n";
	#print "CDS : \"$gene_model\"\nStart_not_found 1\n\n" if $debug;
	$start_tag_val = 1;
      }
      unless (($incomplete) || ($gene_model->name eq "K04G2.11")) {
	push(@error2,"ERROR: $gene_model Start_not_found tag present\n"); 
	print "ERROR: $gene_model Start_not_found tag present\n" if ($verbose);
      }
    }

    if ($gene_model->get('End_not_found')) {
      $end_tag = "present";
      unless ($incomplete) {
	push(@error2,"ERROR: $gene_model End_not_found tag present\n");
	print "ERROR: $gene_model End_not_found tag present\n" if $verbose;
      }
    }

    if ($gene_model->get('Genetic_code')) {
      $genetic_code = $gene_model->get('Genetic_code')->right->name;
    }

    #All Isoforms should have the Isoform tag set. cds_regex

    if (($gene_model_name =~  (/$cds_regex/)) && ($method_test !~  /history/)) {
      my $Isoform = $gene_model->at('Properties.Isoform');
      
      if ($gene_model_name =~  (/\S+[a-z]$/)) {
	push(@error3, "ERROR: $gene_model [$gene_model_name] requires an Isoform\n") unless (defined $Isoform);
      }
      if ($gene_model_name =~  (/S+\d$/)) {
	push(@error3, "ERROR: $gene_model [$gene_model_name] has Isoform tag set but no isoform letter in sequence name\n") if (defined $Isoform);
      }
    }

	


    #############################
    # Pseudogene Specific       #
    #############################

    # check that Pseudogenes have a type tag.
    if ($method_test eq "Pseudogene" && $gene_model_name =~ (/$cds_regex/)) {
      my $prob_prediction = $gene_model->at('Type');
      push(@error3, "ERROR: The Pseudogene $gene_model does not have a Type Tag!\n") if (!defined($prob_prediction));
    }

    ############################
    # Transcript_specific      #
    ############################

    if ($method_test =~ (/transcript/) or ($method_test =~ (/RNA/)) && $gene_model_name =~ (/\w+\d+\.\d+\Z/)) {
      my $prob_prediction = $gene_model->at('Visible.Brief_identification');
      unless ($method_test =~ (/history_transcript/)) {push(@error3, "ERROR: The Transcript $gene_model does not have a Brief_identification and will throw an error in the build :(!\n") if (!defined($prob_prediction));
					  }
    }

    ###################################
    #All gene predictions should have #
    ###################################

    # check that 'Sequence' tag is present and if so then grab parent sequence details
    my $source;
    my $source_type;
    if (defined($gene_model->Sequence)){
      $source = $gene_model->Sequence->name;
      $source_type = "Sequence";
    }
    elsif (defined($gene_model->Transposon)){
      $source = $gene_model->Transposon->name;
      $source_type = "Transposon";
    }
    else {
      push(@error1,"ERROR: $gene_model has no Parent, cannot check DNA\n");
      print "ERROR: $gene_model has no Sequence tag, cannot check DNA\n" if $verbose;
      next CHECK_GENE;
    }

    # check species is correct
    my $species;
    ($species) = ($gene_model->get('Species'));
    push(@error3,"ERROR: $gene_model species is $species\n") if ($species ne "$speciesfn");
    print "ERROR: $gene_model species is $species\n" if ($species ne "$speciesfn" && $verbose);

    # check Method isn't 'hand_built'
    push(@error3,"ERROR: $gene_model method is hand_built\n") if ($method_test eq 'hand_built');
    print "ERROR: $gene_model method is hand_built\n" if ($method_test eq 'hand_built' && $verbose);

    # check From_laboratory tag is present.
    if (($method_test ne 'Genefinder') && ($method_test ne 'twinscan') && ($method_test ne 'jigsaw') && ($method_test ne 'RNASEQ.Hillier.Aggregate')) {
      my $laboratory = ($gene_model->From_laboratory);
      push(@error3, "ERROR: $gene_model does not have From_laboratory tag\n") if (!defined($laboratory));
      print "ERROR: $gene_model does not have From_laboratory tag\n" if (!defined($laboratory) && $verbose);
    }

    # check that history genes have a history method.
    if ($method_test !~ /history/ && $gene_model_name =~ /$cds_regex\:\w+/) {
      push(@error3, "ERROR: $gene_model history object doesn't have a history method.\n");
    }

    # check that history genes are renamed.
    if (($method_test =~ /history/ && !($gene_model_name =~ /\:/) && !($gene_model_name =~ /WBTransposon/))) {


      push(@error3, "ERROR: $gene_model needs to be renamed as it is part of history.\n");
    }

    if ($method_test eq "Transposon") {
      next CHECK_GENE;
    }

    #Gene ID checks.
    unless ($method_test eq "Transposon_CDS") {
      my $Gene_ID     = $gene_model->at('Visible.Gene.[1]');
      my $Genehist_ID = $gene_model->at('Visible.Gene_history.[1]');

      unless ($nogenome) {
	# Can't have both Gene and Gene_history.
	if ((defined $Genehist_ID) && (defined $Gene_ID)) {
	  push(@error2, "ERROR: Gene Model $gene_model contains both a Gene and a Gene_history tag, Please fix.\n");
	}
	
	#History genes have to have a Gene_history ID of 8 digits.
	if ($gene_model_name =~ (/$cds_regex\:\w+/)) {
	  if (defined $Genehist_ID) {
	    push(@error2, "ERROR: The Gene ID '$Genehist_ID' in $gene_model is invalid!\n") unless ($Genehist_ID =~ /WBGene[0-9]{8}/);
	  } 
	  else {
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
	  } 
	  else {
	    push(@error2, "ERROR: $gene_model does not have a Gene ID!\n") unless ($method_test eq 'Transposon_Pseudogene');
	  }
	  if (defined $Genehist_ID) {
	    push(@error2, "ERROR: $gene_model should not contain a Gene ID under Gene_history\n");
	  }
	}
      }

      #####################################################################################################

      # then run misc. sequence integrity checks
      unless ($nogenome) {
	my $dna = $gene_model->asDNA();
	if (!$dna) {
	  push(@error1,"ERROR: $gene_model can't find any DNA to analyse\n");
	  print "ERROR: $gene_model can't find any DNA to analyse\n" if $verbose;
	  next CHECK_GENE;
	}
	
	# feed DNA sequence to function for checking
	unless ((defined $gene_model) && (defined$start_tag) && (defined$end_tag) && (defined$dna) && (defined$method_test)) {print "$_\n";}
	&test_gene_sequence_for_errors($gene_model,$start_tag,$start_tag_val,$end_tag,$dna,$method_test,$genetic_code);
      }
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
	if ($sequence_names{$base} eq 'miRNA_primary_transcript' && (($sequence_names{"${base}a"} && $sequence_names{"${base}a"} eq 'miRNA') || ($sequence_names{"${base}b"} && $sequence_names{"${base}b"} eq 'miRNA'))) {
	  next;
	} 
	# ignore the primary and mature miRNA forms
	push(@error1, "ERROR: The $sequence_names{$base} sequence '$base' and the $sequence_names{$sequence_name} sequence '$sequence_name' both exist!\n")
      }
    }
  }
}


#####################################
# Additional Tests on whole classes #
#####################################
sub single_query_tests {
  
  {
    #Transcript checks from camcheck
    my @Transcripts= $db->fetch(-query=>'find elegans_RNA_genes where method != "history_transcript" AND NOT Transcript');
    if(@Transcripts){
      foreach (@Transcripts){
        $log->error;
        $log->write_to("ERROR: $_ has no Transcript tag, this will cause errors in the build\n");
      }
    }
    else {
      $log->write_to("\nTranscripts all have Transcript set\n");
    }
  }
  
  {
    # Transposon checks
    if ($species eq "elegans") {
      my @Transposons= $db->fetch(-query=>'find Transposon');
      my $Transposon_no = @Transposons;
      unless ($Transposon_no eq "758"){print "\nChange in Transposon_numbers required 758 actual $Transposon_no - has additional Transposon annotation been done?\n"}
    }
  }    

  # Check for non-standard methods in CDS class
  my @CDSfilter;
  if ($species eq 'brugia') {
    @CDSfilter = $db->fetch (-query => 'FIND CDS; method != Transposon_CDS; method != Transposon_Pseudogene; method != curated; method !=history; method !=Genefinder; method !=twinscan; method !=jigsaw; method !=mGene; method !=RNASEQ.Hillier; method !=RNASEQ.Hillier.Aggregate; method !=cufflinks*; method !=genBlastG');
  }
  else {
    @CDSfilter = $db->fetch (-query => 'FIND CDS; method != Transposon_CDS; method != Transposon_Pseudogene; method != curated; method !=history; method !=Genefinder; method !=twinscan; method !=jigsaw; method !=mGene; method !=RNASEQ.Hillier; method !=RNASEQ.Hillier.Aggregate');
  }
  foreach my $CDSfilter (@CDSfilter) {
    push(@error4, "ERROR! CDS:$CDSfilter contains an invalid method please check\n");
  }
}
	
##########################################
# Additional Tests on the build instance #
##########################################

sub extra_build_checks {
  my $Gene_model;
  my @Gene_models = $db->fetch (-query => 'Find Gene where Live AND Species = "$speciesfn" AND Sequence_name AND NOT Corresponding_CDS AND NOT Corresponding_pseudogene AND NOT Corresponding_transcript');
  my $count = scalar(@Gene_models);
  print "\nThere are $count genes that are Live but not attached to a current gene model\n\n";
  $log->write_to("\nThere are $count genes that are Live but not attached to a current gene model\n\n");
  foreach $Gene_model (@{Gene_models}){
    $log->write_to("$Gene_model \n");
  }
}

#################################################################
# Subroutines
#################################################################

sub test_gene_sequence_for_errors{
  my $gene_model = shift;
  my $start_tag = shift;
  my $start_tag_val = shift;
  my $end_tag = shift;
  my $dna = shift;
  my $method_test =shift;
  my $genetic_code=shift;

  if ($start_tag eq "present") {
    unless (defined $start_tag_val) {
      $start_tag_val = '1';
    }
    print "start tag $start_tag value = $start_tag_val\n" if ($verbose);
    if ($start_tag_val eq 0) {
      $start_tag_val = '1';
      print "start tag $start_tag value modified = $start_tag_val\n" if ($verbose);
    }
  }
  
  if ($method_test eq '*history'){
    print "Not checking the codon or dna usage of $gene_model as it is a history model\n" if ($verbose);
  }
  if ($method_test ne 'history'){
    # trim DNA sequence to just A,T,C,G etc.
    $dna =~ s/\n//g;
    my $length_gene_name = length($gene_model)+1;
    $dna = substr($dna, $length_gene_name);
    if (!$dna){
      push(@error1, "$gene_model has a problem with it's DNA connection.\n");
      next CHECK_GENE;
    }
    # calculate other necessary values
    my $gene_model_length = length($dna);
    my $remainder;

    if ($incomplete) {
      if (($gene_model->Method eq 'curated') && ($gene_model->Start_not_found)) {
	my $extra = $gene_model->Start_not_found->name;
	my $length_calc = $gene_model_length + $extra;
	$remainder = $length_calc%3;
      }
      else {
	$remainder = $gene_model_length%3;
      }
    }
    else {
      $remainder = $gene_model_length%3;
    }
    my $start_codon = substr($dna,0,3);
    my $stop_codon = substr($dna,-3,3);   
    my $Lab = ($gene_model->get('From_laboratory'));
    
    # check for length errors(CDS specific)
    my @ck;
    if (!$basic) {
      my $warning;
      @ck = grep(/^$gene_model/, @checkedgenes);
      push (@ck,"0");
      if (($gene_model_length < 75) && ($method_test eq 'curated') && ($ck[0] ne $gene_model->name)) {
	$warning = "WARNING: $gene_model is very short ($gene_model_length bp),";
	print "WARNING: $gene_model is very short ($gene_model_length bp), " if $verbose;
	if (defined($gene_model->at('Properties.Coding.Confirmed_by'))) {
	  $warning .= "gene is Confirmed\n";
	  print "gene is Confirmed\n" if $verbose;
	}
	elsif (defined($gene_model->at('Visible.Matching_cDNA'))) {
	  $warning .= "gene is Partially_confirmed\n";
	  print "gene is Partially_confirmed\n" if $verbose;
	}
	else {
	  $warning .= "gene is Predicted\n";
	  print "gene is predicted\n" if $verbose;
	}
	push(@error3, $warning) unless ($basic);
      }
      elsif (($gene_model_length < 100) && ($method_test eq 'curated')) {
	if (defined($gene_model->at('Properties.Coding.Confirmed_by'))) {
	  $warning = "WARNING: $gene_model is short ($gene_model_length bp) and is Confirmed\n";
	  print "WARNING: $gene_model is short ($gene_model_length bp) and is Confirmed\n" if $verbose;
	}
	elsif (defined($gene_model->at('Visible.Matching_cDNA'))) {
	  $warning .= "WARNING: $gene_model is short ($gene_model_length bp) and is Partially_confirmed\n";
	  print "WARNING: $gene_model is short ($gene_model_length bp) and is Partially_confirmed\n" if $verbose;
	}
      else {
	$warning .= "WARNING: $gene_model is short ($gene_model_length bp) and is Predicted\n";
	print "WARNING: $gene_model is short ($gene_model_length bp) and is Predicted\n" if $verbose;
      }
	push(@error5, $warning);
      }
    }
    # Is the gene prediction complete?
    if (($remainder != 0) && ($method_test eq 'curated')) {
      if (($end_tag ne "present") && ($start_tag ne "present")) {
	push(@error1,"ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n");
	print "ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n" if $verbose;
      }
      else {
	push(@error2,"ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n") unless ($incomplete);
	print "ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n" if $verbose;
      }
    }
    if (defined $gene_model->Sequence) {
      unless (($gene_model->name =~ /MTCE/) || ($gene_model->Sequence->name =~ /MITO/)) {
	# look for incorrect stop codons (CDS specific)
	if (($stop_codon ne 'taa') && ($stop_codon ne 'tga') && ($stop_codon ne 'tag') && ($method_test eq 'curated')) {
	  if ($end_tag ne "present") {
	    push(@error1, "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n");
	    print "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n" if $verbose;
	    #print "CDS : \"$gene_model\"\nEnd_not_found\n\n" if ($debug);
	  }
	  else {
	    push(@error2,"ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag present\n") unless ($incomplete);
	    print "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag present\n" if $verbose;
	  }
	}
	# look for incorrect start codons(CDS specific)
	if ($gene_model->name eq "K04G2.11") {
	  print "WARNING: $gene_model  intentionally utilises a novel '$start_codon' start codon......Ignoring\n" if $verbose;
	}
	if (($start_codon ne 'atg') && ($method_test eq 'curated') && ($start_tag ne "present") && ($gene_model->name ne "K04G2.11")) {
	  if (($start_tag ne "present")) {
	    push(@error1,"ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n");
	    print "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n" if $verbose;
	  }  
	  else {
	    push(@error2, "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag present\n") unless ($incomplete);
	    print "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag present\n" if $verbose;
	  } 
	}
	
	# check for internal stop codons (CDS specific)
	my $i;
	my $j;
	if ($start_tag_val) {
	  print "start_tag_val final = $start_tag_val\n" if ($verbose);
	  if ($start_tag_val eq 3) {#if 3 +2
	    $i = '2';
	  }
	  elsif ($start_tag_val eq 2) {#if 2 +1
	    $i = '1';
	  }
	  else {
	    $i = '0';
	  }
	}
	else {
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
	    if (($method_test eq 'curated')) {
	      push(@error1, "ERROR: $gene_model internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n") unless ($genetic_code eq 'Selenocysteine');      
	      print "ERROR: $gene_model internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n" if (($verbose) && ($genetic_code ne 'Selenocysteine'));
	    }
	  }
	}
	if ($species eq "elegans") {
	  # look for non-ACTG characters
	  if ($dna =~ /[^acgt]/i) {
	    $dna =~ s/[acgt]//g;
	    push(@error2, "ERROR: $gene_model DNA sequence contains the following non-ATCG characters: $dna\n"); 
	    print "ERROR: $gene_model DNA sequence contains the following non-ATCG characters: $dna\n" if $verbose;
	    
	  }
	}
      } ##MTCE exclusion loop
    } ##Sequence_parent loop
    else {
      push(@error1, "$gene_model has no SParent Sequence connection.\n") unless (defined $gene_model->Transposon);
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
    push(@error1, "$gene_model has no SParent Sequence connection.\n");
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
    push(@error1, "$gene_model has a malformed SParent Sequence entry.\n");
    return $result;    
  }

}


############################################################################################
sub by_number{ $a <=> $b;}
############################################################################################

__END__

=pod

=head1 NAME - check_predicted_genes.pl

=back


=head1 USAGE

=over 4

=item check_predicted_genes.pl path_to_acedb_database [log_file] 

=back

This script is designed to check the validity of predicted gene objects from any wormbase (or
other acedb) database.  The script will analyse objects in the 'All_genes' subclass.
The script emails the top 20 problems each day (sorted by severity).


=over 4

=item MANDATORY arguments: -database <database_path>

This argument must be a path to a valid acedb database, e.g.
check_predicted_genes.pl -database ~wormpub/DATABASES/camace

=back

=over 4

=item OPTIONAL arguments: -log <logfile>, -verbose, -basic

If the file specified by -log already exists, the script will append the output to that file.
Otherwise it will attempt to write to a new file by that name.  If no log file is specified,
the script will generate a log file in ~wormpub/BUILD/logs.

If -verbose is specified, output will be written to screen as well as to the log file

If -basic is specified, the script will skip some simple size checks (small introns and small genes)

=back



=head1 DOCUMENTATION

=over 4

=back

check_predicted_genes.pl performs a wide range of checks on ?CDS objects to ensure
that they are valid objects and will behave properly within the database.  The script was
written to be called from within the I<camcheck> script which performs a wider range of error
checking on the camace database, but the script can also be run as a standalone program to check
the status of other wormbase databases such as autoace, etc.

The script checks predicted genes for the following:

=over 4

=item 1.

Incorrect sequence length (i.e. length in bp is not a multiple of three)

=item 2.

Improper start codon with no 'Start_not_found' tag present


=item 3.

Improper end codon with no 'End_not_found' tag present

=item 4.

Internal stop codons

=item 5.

Non ATCG characters in the sequence

=item 6.

CDSs which dont have a 'Species = species_full_name' tag-value pair

=item 7.

CDSs which belong to superlink objects

=item 8.

CDSs which have 'Method = hand_built'

=item 9.

Presence of 'Sequence' tag

=item 10.

Inconsistencies in exon coordinates, i.e. where the coordinates of any two exons in a gene might overlap

=item 11.

Checks for presence of 'From_laboratory' tag

=head1 SEE ALSO

L<camcheck>

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

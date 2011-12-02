#!/software/bin/perl -w
#
# check_predicted_genes.pl
#
# by Keith Bradnam
#
# Last updated on: $Date: 2011-12-02 10:33:44 $
# Last updated by: $Author: klh $
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

my ($verbose, $db_path, $basic, $test1, $debug, $store, $test,$build);

GetOptions ("verbose"    => \$verbose, # prints screen output.
	    "database=s" => \$db_path, # Path to the database you want to check.
	    "basic"      => \$basic,   # Ignores some of the checks.
	    "debug:s"    => \$debug,   # turns on more printing and errorlogging
	    "test1:s"    => \$test1,   # only checks the CDSs from 1 clone in the database.
	    "store:s"    => \$store,   # 
	    "test"       => \$test,    # Test build
	    "build"      => \$build,   # Checks specific to a full database containing genes and models.
	   );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			   );
}

#Establish Log
my $log = Log_files->make_build_log($wormbase);

# Establish database connection etc.
$log->log_and_die("Please use -database <path> specify a valid database directory.\n\n") if (!defined($db_path));
# Specify which tace to use if you are using -program flag
my $tace = $wormbase->tace;
my $db = Ace->connect(-path=>$db_path) or  $log->log_and_die("Couldn't connect to $db_path\n". Ace->error);

# create separate arrays for different classes of errors (1 = most severe, 4 = least severe)
our (@error1, @error2, @error3, @error4, @error5);

#Permitted exceptions where the error has been checked.
#List of verified small genes or have confirmed small introns.
my @checkedgenes = ('F56H11.3b','F13E9.8','Y45F10D.7','T21C12.8','F58G11.2','F54E4.3','F43D9.3a','F43D9.3b','F36H1.3','F02C12.1','C40H5.1', 'H12D21.1', 'H12D21.12','H12D21.13', 'H12D21.14',  'H12D21.15','W06A7.5','Y50E8A.17', 'ZC412.6', 'ZC412.7', 'R74.3a');



################################
#         Main Body            # 
################################
my @Predictions;
my %sequence_names;
my %sequence_classes;

if ($build) {
  &extra_build_checks;
}
else {
  # Fetch Gene Models to be tested (All_genes) #
  if ($test1) {
    $log->write_to("Only checking genes on ${test1}.......\n");
    #  @Predictions = $db->fetch('All_genes','Y32B12C*');
    @Predictions = $db->fetch('All_genes',"${test1}*");
    foreach my $Predictions(@Predictions) {
      print $Predictions->name."\n";
    }
  }
  else {
    @Predictions = $db->fetch (-query => 'FIND All_genes');
  }
  my $gene_model_count=@Predictions;
  $log->write_to("Checking $gene_model_count Predictions in $db_path\n\n");
  print "\nChecking $gene_model_count Predictions in '$db_path'\n\n" if $verbose;

  &main_gene_checks;
  &single_query_tests;
  
  
  # print warnings to log file, log all category 1 errors, and then fill up.
  my $count_errors =0;
  my @error_list = ( \@error1, \@error2, \@error3,  \@error4, \@error5);
  foreach my $list (@error_list) {
    foreach my $error (@{$list}) {
      $count_errors++;
      $log->write_to("$count_errors $error");
      last if $count_errors > 190;
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
CHECK_GENE:
foreach my $gene_model ( @Predictions ) {
  #next unless ($gene_model eq "Y32B12C.2b"); #stop the script at a specified gene. debug line
  print STDOUT "$gene_model\n" if $verbose;
  my $method_test = ($gene_model->Method);
  my @exon_coord1 = sort by_number ($gene_model->get('Source_exons',1));
  my @exon_coord2 = sort by_number ($gene_model->get('Source_exons',2));
  my $i;
  my $j;

#  if (!defined($method_test)) {print "$gene_model\n";}
  
  # check for duplicated sequence names
  if (exists $sequence_names{$gene_model->name}) {
    my $class = $gene_model->class;
    push(@error1, "ERROR: $class $gene_model is both a $method_test and a $sequence_classes{$gene_model->name} $sequence_names{$gene_model->name}\n");
    print "ERROR: $class $gene_model is both a $method_test and a $sequence_classes{$gene_model->name} $sequence_names{$gene_model->name}\n";
  }
  $sequence_names{$gene_model->name} = $method_test; # save all the names and methods
  $sequence_classes{$gene_model->name} = $gene_model->class;

  if ($method_test ne 'Transposon') {
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
      push(@error5,"WARNING: $gene_model has a small intron ($intron_size bp)\n") if (($intron_size > 33) && ($intron_size < 39) && (!$basic) && ($method_test eq 'curated') && ($ck[0] ne $gene_model->name));
    }
  }
	
  for ($i=0; $i<@exon_coord1; $i++) {
    my $start = $exon_coord1[$i];
    my $end = $exon_coord2[$i];
    for ($j=$i+1;$j<@exon_coord1;$j++) {
      if (($end > $exon_coord1[$j]) && ($start < $exon_coord2[$j])) {
	print "ERROR: $gene_model exon inconsistency, exons overlap\n" if $verbose;
	push(@error1,"ERROR: $gene_model exon inconsistency, exons overlap\n") if ($method_test !~ /^history/);
      }
    }
  }

  # check that 'Start_not_found' and 'End_not_found' tags present? (CDS specific.....extended to all genes :) )
  my $start_tag = "";
  my $end_tag = "";

  if ($gene_model->get('Start_not_found')) {
    $start_tag = "present";
    push(@error2,"ERROR: $gene_model Start_not_found tag present\n");
    print "ERROR: $gene_model Start_not_found tag present\n" if $verbose;
  }

  if ($gene_model->get('End_not_found')) {
    $end_tag = "present";
    push(@error2,"ERROR: $gene_model End_not_found tag present\n");
    print "ERROR: $gene_model End_not_found tag present\n" if $verbose;
  }

  #All Isoforms should have the Isoform tag set.
  if (($gene_model->name =~ (/\w+\.\d+[a-z]$/)) && ($method_test !~  /^history/)) {
    my $Isoform = $gene_model->at('Properties.Isoform');
    push(@error3, "ERROR: $gene_model requires an Isoform\n") if !$Isoform;
  }

  #All with an Isoform tag should be named correctly
  if (($gene_model->name =~ (/\w+\.\d+$/)) && ($method_test !~ /^history/)) {
    my $Isoform = $gene_model->at('Properties.Isoform');
    push(@error3, "ERROR: $gene_model requires an Isoform\n") if $Isoform;
  }

  #Test for erroneous Isoform tags.??
  if (($gene_model->name =~ (/\b\w+\.[0-9]{1,2,3}\b/)) && ($method_test !~ /^history/)) {
    my $Isoform = $gene_model->at('Properties.Isoform');
    push(@error3, "ERROR: $gene_model contains an invalid Isoform tag\n") if $Isoform;
  }

  #############################
  # Pseudogene Specific       #
  #############################
	
  # check that Pseudogenes have a type tag.
  if ($method_test eq "Pseudogene" && $gene_model->name =~ (/\w+\d+\.\d+\Z/)) {
    my $prob_prediction = $gene_model->at('Type');
    push(@error3, "ERROR: The Pseudogene $gene_model does not have a Type Tag!\n") if (!defined($prob_prediction));
  }

  ############################
  # Transcript_specific      #
  ############################
 
  if ($method_test =~ (/transcript/) or ($method_test =~ (/RNA/)) && $gene_model->name =~ (/\w+\d+\.\d+\Z/)) {
    my $prob_prediction = $gene_model->at('Visible.Brief_identification');
    push(@error3, "ERROR: The Transcript $gene_model does not have a Brief_identification and will throw an error in the build :(!\n") if (!defined($prob_prediction));
  }
 
  ###################################
  #All gene predictions should have #
  ###################################

  # check that 'Sequence' tag is present and if so then grab parent sequence details
  my $source;
  if (defined($gene_model->Sequence)){
    $source = $gene_model->Sequence->name;
  }
  elsif (defined($gene_model->Transposon)){
    $source = $gene_model->Transposon->name;
  }
  else {
    push(@error1,"ERROR: $gene_model has no Parent, cannot check DNA\n");
    print "ERROR: $gene_model has no Sequence tag, cannot check DNA\n" if $verbose;
    next CHECK_GENE;
  }
	
  # check species is correct
  my $species;
  ($species) = ($gene_model->get('Species'));
  push(@error3,"ERROR: $gene_model species is $species\n") if ($species ne "Caenorhabditis elegans");
  print "ERROR: $gene_model species is $species\n" if ($species ne "Caenorhabditis elegans" && $verbose);
	
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
  if ($method_test !~ /^history/ && $gene_model->name =~ /\w+\.\w+\:\w+/) {
    push(@error3, "ERROR: $gene_model history object doesn't have a history method.\n");
  } 
	
  # check that history genes are renamed.
  if ($method_test =~ /^history/ && !($gene_model->name =~ /\:/)) {
    push(@error3, "ERROR: $gene_model needs to be renamed as it is part of history.\n");
  }

  if ($method_test eq "Transposon") {
    next CHECK_GENE;
  }

  #Gene ID checks.
  my $Gene_ID     = $gene_model->at('Visible.Gene.[1]');
  my $Genehist_ID = $gene_model->at('Visible.Gene_history.[1]');
	
  #curated Gene modles eg. C14C10.3  C14C10.33 and C14C10.3a have to have an 8 digit gene id.
  if ($gene_model->name =~ /\w+\.\d+\w?$/) {
    if (defined $Gene_ID) {
      push(@error2, "ERROR: The Gene ID '$Gene_ID' in $gene_model is invalid!\n") unless ($Gene_ID =~ /WBGene[0-9]{8}/);
    } else {
      push(@error2, "ERROR: $gene_model does not have a Gene ID!\n") unless (($method_test eq 'Transposon_Pseudogene') && (defined $Genehist_ID));
    }
  }
  #History genes have to have a Gene_history ID of 8 digits.
  elsif ($gene_model->name =~ (/\w+\.\w+\:\w+/)) {
    if (defined $Genehist_ID) {
      push(@error2, "ERROR: The Gene ID '$Genehist_ID' in $gene_model is invalid!\n") unless ($Genehist_ID =~ /WBGene[0-9]{8}/);
    } else {
      push(@error2, "ERROR: $gene_model does not have the Gene_history populated\n");
    }
  }
  # Can't have both Gene and Gene_history.
  if ((defined $Genehist_ID) && (defined $Gene_ID)) {
    push(@error2, "ERROR: Gene Model $gene_model contains both a Gene and a Gene_history tag, Please fix.\n");
  }
	    
  ###################################################################################################################

  # then run misc. sequence integrity checks
  my $dna = $gene_model->asDNA();
  if (!$dna) {
    push(@error1,"ERROR: $gene_model can't find any DNA to analyse\n");
    print "ERROR: $gene_model can't find any DNA to analyse\n" if $verbose;
    next CHECK_GENE;
  }
	    
  # feed DNA sequence to function for checking
  &test_gene_sequence_for_errors($gene_model,$start_tag,$end_tag,$dna,$method_test);
}

# now check that the isoform names are consistent
foreach my $sequence_name (keys %sequence_names) {
  # don't want to look at history objects
  if ($sequence_names{$sequence_name} =~ /history/) {next}

  # does the isoform have multiple letters in
  if ($sequence_name =~ /\w+\.\d+([a-z]{2,})$/) {
    push(@error1, "ERROR: The sequence_name '$sequence_name' is invalid! Multiple letters in the isoform name.\n")


  # if it is an isoform name, check for non-isoforms
  } elsif ($sequence_name =~ /(\w+\.\d+)[a-z]$/) {
    my $base = $1;
    if (exists $sequence_names{$base}) {
      if ($sequence_names{$base} eq 'miRNA_primary_transcript' && $sequence_names{"${base}a"} eq 'miRNA') {next} # ignore the primary and mature miRNA forms
      push(@error1, "ERROR: The $sequence_names{$base} sequence '$base' and the $sequence_names{$sequence_name} sequence '$sequence_name' both exist!\n")
    }
  }
}

}


#####################################
# Additional Tests on whole classes #
#####################################
sub single_query_tests {

#Transcript checks from camcheck
  my @Transcripts= $db->fetch(-query=>'find elegans_RNA_genes NOT Transcript');
  if(@Transcripts){
    foreach (@Transcripts){
      $log->write_to("ERROR: $_ has no Transcript tag, this will cause errors in the build\n");
    }
  }
  else {
    $log->write_to("\nTranscripts OK\n");
  }

# Transposon checks
  my @Transposons= $db->fetch(-query=>'find Transposon');
  my $Transposon_no = @Transposons;
  unless ($Transposon_no eq "234"){print "\nChange in Transposon_numbers\n"}
  

  # Check for non-standard methods in CDS class
  my @CDSfilter = $db->fetch (-query => 'FIND CDS; method != Transposon_CDS; method != Transposon_Pseudogene; method != curated; method !=history; method !=Genefinder; method !=twinscan; method !=jigsaw; method !=mGene; method !=RNASEQ.Hillier !=RNASEQ.Hillier.Aggregate');
  foreach my $CDSfilter (@CDSfilter) {
    push(@error4, "ERROR! CDS:$CDSfilter contains an invalid method please check\n");
  }
}
	
##########################################
# Additional Tests on the build instance #
##########################################

sub extra_build_checks {
  my $Gene_model;
  my @Gene_models = $db->fetch (-query => 'Find Gene where Live AND Species = "Caenorhabditis elegans" AND Sequence_name AND NOT Corresponding_CDS AND NOT Corresponding_pseudogene AND NOT Corresponding_transcript');
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
  my $end_tag = shift;
  my $dna = shift;
  my $method_test =shift;
	    
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
  my $remainder = $gene_model_length%3;
  my $start_codon = substr($dna,0,3);
  my $stop_codon = substr($dna,-3,3);   
  my $Lab = "";
  ($Lab) = ($gene_model->get('From_laboratory'));
	    
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
    } else {
      push(@error2,"ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n");
      print "ERROR: $gene_model length ($gene_model_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n" if $verbose;
    }
  }

  # look for incorrect stop codons (CDS specific)
  if (($stop_codon ne 'taa') && ($stop_codon ne 'tga') && ($stop_codon ne 'tag') && ($method_test eq 'curated') && ($Lab eq "HX")) {
    if ($end_tag ne "present") {
      push(@error1, "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n");
      print "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n" if $verbose;
    } else {
      push(@error2,"ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag present\n");
      print "ERROR: $gene_model '$stop_codon' is not a valid stop codon. End_not_found tag present\n" if $verbose;
    }
  }

  # look for incorrect start codons(CDS specific)
  if (($start_codon ne 'atg') && ($method_test eq 'curated') && ($Lab eq "HX")) {
    if (($start_tag ne "present")) {
      push(@error1,"ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n");
      print "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n" if $verbose;
    } else {
      push(@error2, "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag present\n");
      print "ERROR: $gene_model '$start_codon' is not a valid start codon. Start_not_found tag present\n" if $verbose;
    }
  }

  # check for internal stop codons (CDS specific)
  my $i;
  my $j;
  for ($i=0; $i<$gene_model_length-3;$i+=3) {
    # hold position of codon in $j
    $j=$i+1;
    my $codon =substr($dna,$i,3);
    if (($codon eq "taa") || ($codon eq "tag") || ($codon eq "tga")) {      
      my $previous_sequence = substr($dna, $j-11,10);
      my $following_sequence = substr($dna, $j+2, 10);
      my $offending_codon = substr($dna, $j-1, 3);
      if (($method_test eq 'curated') && ($Lab eq "HX")) {
	push(@error1, "ERROR: $gene_model internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n");      
	print "ERROR: $gene_model internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n" if $verbose;
      }
    }
  }
  
  # look for non-ACTG characters
  if ($dna =~ /[^acgt]/i) {
    $dna =~ s/[acgt]//g;
    push(@error2, "ERROR: $gene_model DNA sequence contains the following non-ATCG characters: $dna\n"); 
    print "ERROR: $gene_model DNA sequence contains the following non-ATCG characters: $dna\n" if $verbose;
  }
}


sub by_number{ $a <=> $b;}

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
the status of other wormbase databases such as autoace, stlace etc.

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

CDSs which dont have a 'Species = Caenorhabditis elegans' tag-value pair

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

#!/usr/local/bin/perl -w

# check_predicted_genes.pl
#
# by Keith Bradnam aged 12 and a half
# 15/06/01
#
# see pod documentation at end of file for more information about this script

use Ace;
#use IO::Handle;
use strict;
$|=1;


################################
#
# Establish database connection
#
################################

die "Please specify a path to the database as a command line parameter.\n\n" if ($ARGV[0] eq "");

my $db_path = $ARGV[0];

# Specify which tace to use if you are using -program flag



# RELEASE.SUPPORTED version is (yet again) buggy, rolling back to an older version
# for now

my $tace = glob("~acedb/RELEASE.2001_05_01.BUILD/bin.ALPHA_4/tace");
#my $tace = glob("~acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace");

#my $tace = glob("~acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace");
my $db = Ace->connect(-path=>$db_path, -program=>$tace) || die "Couldn't connect to $db_path\n", Ace->error;

#otherwise just use the default style of Ace database connection
#my $db = Ace->connect(-path=>$db_path) || die "Couldn't connect to $db_path\n", Ace->error;




################################
#
# Toggle amount of output
#
################################

# toggle reporting of genes with improper stop/start codons and/or genes which have length's 
# not divisible by three but have 'Start_not_found' and/or 'End_not_found' tags.
# 'ON' setting means that you will get full output, 'OFF' means that you will get restricted output

#my $verbose = 'ON';
my $verbose = 'OFF';

# set up log file for output, use specified command line argument where present or write to screen
# if log file specified and it exists, then append.  Else write to new file.

if(defined($ARGV[1])){
  my $log = "$ARGV[1]";

  if(-e $log){
    open(LOG,">>$log");  
  }
  else{
    open(LOG,">$log");
  }
  # make LOG the default location for 'print' commands
#  LOG->autoflush();
  select(LOG);
}


################################
#
# Fetch Predicted genes
#
################################

# load all predicted genes into an array
my @predicted_genes = $db->fetch (-query => 'FIND Predicted_gene');

# count genes in database
my $gene_count=@predicted_genes;
print "\nChecking $gene_count predicted genes in '$db_path'\n\n";



################################
#
# Run checks on genes
#
################################

CHECK_GENE:
foreach my $gene (@predicted_genes){
  # get gene
  my $gene_object = $db->fetch(Sequence=>$gene);
  
  # check that 'Source' tag is present and if so then grab parent sequence details
  my $source = $gene_object->Source;
  if (!defined ($source)){
    print "Gene error - $gene: has no Source tag, cannot check DNA\n";
    next CHECK_GENE;
  }
  my $parent = $db->fetch(Sequence=>$source);
  
  # check to see if any predicted gene belongs to superlink object...they shouldn't
  if ($parent =~ m/SUPERLINK/){
    print "Gene error - $gene: belongs to a superlink object ($parent)\n";
  }


  # check coordinate system exons in relation to each other
  my @exon_coord1 = $gene_object->get('Source_Exons',1);
  my @exon_coord2 = $gene_object->get('Source_Exons',2);

  my $i;
  my $j;
  for($i=0; $i<@exon_coord1; $i++){
    my $start = $exon_coord1[$i];
    my $end = $exon_coord2[$i];
    for ($j=$i+1;$j<@exon_coord1;$j++){
      if(($end > $exon_coord1[$j]) && ($start < $exon_coord2[$j])){
	print "Gene error - $gene: exon inconsistency, exons overlap\n";
      }
    }
  }

  # check coordinates system of exons in relation to parent
  #
  # haven't done this bit yet
  if ($parent !~ m/LINK/){    
    my ($parent_length) = $parent->DNA(2);
    $parent_length = 0 if (!defined($parent_length));
  }





  # check that 'Start_not_found' and 'End_not_found' tags present?
  my $start_tag = "";
  my $end_tag = "";
  $start_tag = "present" if ($gene_object->get('Start_not_found'));
  $end_tag = "present"  if ($gene_object->get('End_not_found'));


  # check species is correct
  my $species = "";
  ($species) = ($gene_object->get('Species'));
  print "Gene error - $gene: species is $species\n" if ($species ne "Caenorhabditis elegans");
  

  # check Method isn't 'hand_built'
  my $method = "";
  ($method) = ($gene_object->get('Method'));
  print "Gene error - $gene: method is hand_built\n" if ($method eq 'hand_built');

  # check From_laboratory tag is present
  my $laboratory = ($gene_object->From_laboratory);
  print "Gene error - $gene: does not have From_laboratory tag\n" if (!defined($laboratory));


  # then run misc. sequence integrity checks
  my $dna = $gene->asDNA();
  if(!$dna){
    print "Gene error - $gene: can't find any DNA to analyse\n";
    next CHECK_GENE;
  }

  # feed DNA sequence to function for checking
  &test_gene_sequence_for_errors($gene,$start_tag,$end_tag,$dna,$verbose);

}
close(LOG) if(defined($ARGV[1]));

#################################################################

sub test_gene_sequence_for_errors{

  my $gene = shift;
  my $start_tag = shift;
  my $end_tag = shift;
  my $dna = shift;
  my $verbose = shift;

  # trim DNA sequence to just A,T,C,G etc.
  $dna =~ s/\n//g;
  my $length_gene_name = length($gene)+1;
  $dna = substr($dna, $length_gene_name);

  # calculate other necessary values
  my $gene_length = length($dna);
  my $remainder = $gene_length%3;
  my $start_codon = substr($dna,0,3);
  my $stop_codon = substr($dna,-3,3);   

  # check for length errors
  if ($remainder != 0){
    if (($end_tag ne "present") && ($start_tag ne "present")){
      print "Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n";
    }
    elsif($verbose eq 'ON'){
      print "Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n";
    }   
  }   


  # look for incorrect stop codons
  if (($stop_codon ne 'taa') && ($stop_codon ne 'tga') && ($stop_codon ne 'tag')){
    if($end_tag ne "present"){
      print "Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n";
    }
    elsif($verbose eq 'ON'){
      print "Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag present\n";
    }
  }


  # look for incorrect start codons
  if ($start_codon ne 'atg'){
   if($start_tag ne "present"){
      print "Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n";
    }
    elsif($verbose eq 'ON'){
      print "Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag present\n";
    }
  }

  
  # check for internal stop codons
  my $i;
  my $j;

  for($i=0; $i<$gene_length-3;$i+=3){
    # hold position of codon in $j
    $j=$i+1;
    my $codon =substr($dna,$i,3);
    if(($codon eq "taa") || ($codon eq "tag") || ($codon eq "tga")){      
      my $previous_sequence = substr($dna, $j-11,10);
      my $following_sequence = substr($dna, $j+2, 10);
      my $offending_codon = substr($dna, $j-1, 3);
      print "Gene error - $gene: internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n";      
    }
  }			       

  # look for non-ACTG characters
  if ($dna =~ /[^acgt]/i){
    $dna =~ s/[acgt]//g;
    print "Gene error - $gene: DNA sequence contains the following non-ATCG characters: $dna\n" 
  }

}

__END__

=pod

=head1 NAME - check_predicted_genes.pl

=back


=head1 USAGE

=over 4

=item check_predicted_genes.pl path_to_acedb_database [log_file] 

=back

This script is designed to check the validity of predicted gene objects from any wormbase (or
other acedb) database.  The script will only analyse objects in the 'Predicted_gene' class.


=over 4

=item MANDATORY arguments:

First argument must be a path to a valid acedb database, e.g.
check_predicted_genes.pl /wormsrv2/camace

=back

=over 4

=item OPTIONAL arguments:

A log file can be named as a second optional argument.  If such a file already
exists, output will be appended to it.  If it doesn't exist, then such a file will be 
created.  If no log file is specified then all output will be sent to STDOUT as the 
default.  

=back



=head1 DOCUMENTATION

=over 4

=back

check_predicted_genes.pl performs a wide range of checks on predicted gene objects to ensure
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

Sequences which don't have a 'Species = Caenorhabditis elegans' tag-value pair

=item 7.

Subsequences which belong to superlink objects

=item 8.

Subsequences which have 'Method = hand_built'

=item 9.

Presence of 'Source' tag

=item 10.

Inconsistencies in exon coordinates, i.e. where the coordinates of any two exons in 
a gene might overlap

=item 11.

Checks for presence of 'From_laboratory' tag


=head1 SEE ALSO

L<camcheck>


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

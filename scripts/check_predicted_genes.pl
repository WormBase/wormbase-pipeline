#!/usr/local/bin/perl5.6.1 -w

# check_predicted_genes.pl
#
# by Keith Bradnam aged 12 and a half
#
# Last updated on: $Date: 2002-11-22 09:19:50 $
# Last updated by: $Author: krb $
#
# see pod documentation at end of file for more information about this script

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;
use IO::Handle;
use Getopt::Long;

$|=1;

our($verbose,$db_path,$log);

GetOptions ("verbose"    => \$verbose,
	    "database=s" => \$db_path,    
	    "log=s"        => \$log);

# verbose mode
# toggle reporting of genes with improper stop/start codons and/or genes which have length's 
# not divisible by three but have 'Start_not_found' and/or 'End_not_found' tags.
# 'ON' setting means that you will get full output, 'OFF' means that you will get restricted output,
# Also outputs gene names as you cycle through main loop


################################
# Establish database connection
################################

die "Please use -database <path> specify a valid database directory.\n\n" if (!defined($db_path));

# Specify which tace to use if you are using -program flag

my $tace = &tace;
my $db = Ace->connect(-path=>$db_path, -program=>$tace) || die "Couldn't connect to $db_path\n", Ace->error;


# set up log file for output, use specified command line argument where present or write to screen
# if log file specified and it exists, then append.  Else write to new file.

if($log){
  if(-e $log){
    open(LOG,">>$log");  
  }
  else{
    open(LOG,">$log");
  }
}
else{
  my $rundate    = `date +%y%m%d`; chomp $rundate;
  $log = "/wormsrv2/logs/check_predicted_genes.log.$rundate.$$";
  open(LOG,">$log") || die "cant open $log\n";
}
print LOG "\ncheck_predicted_genes.pl started at ",`date`,"\n";


LOG->autoflush();



################################
# Fetch Predicted genes
################################

my @predicted_genes = $db->fetch (-query => 'FIND Predicted_gene');
my $gene_count=@predicted_genes;
print LOG "Checking $gene_count predicted genes in '$db_path'\n\n";
print "\nChecking $gene_count predicted genes in '$db_path'\n\n" if $verbose;;


################################
# Run checks on genes
################################

# create separate arrays for different classes of errors (1 = most severe, 4 = least severe)
our @error1;
our @error2;
our @error3;
our @error4;


CHECK_GENE:
foreach my $gene (@predicted_genes){
  print STDOUT "$gene\n" if $verbose;
  # get gene
  my $gene_object = $db->fetch(Sequence=>$gene);
  
  # check that 'Source' tag is present and if so then grab parent sequence details
  my $source = $gene_object->Source;
  if (!defined ($source)){
    push(@error1,"Gene error - $gene: has no Source tag, cannot check DNA\n");
    print "Gene error - $gene: has no Source tag, cannot check DNA\n" if $verbose;
    next CHECK_GENE;
  }


  # check coordinate system exons in relation to each other
  my @exon_coord1 = sort by_number ($gene_object->get('Source_Exons',1));
  my @exon_coord2 = sort by_number ($gene_object->get('Source_Exons',2));

  my $i;
  my $j;

  for($i=1; $i<@exon_coord2;$i++){
    my $intron_size = ($exon_coord1[$i] - $exon_coord2[$i-1] -1);
    print "Gene warning - $gene has a small intron ($intron_size bp)\n" if (($intron_size < 31)  && $verbose);
    push(@error3,"Gene warning - $gene has a very small intron ($intron_size bp)\n") if ($intron_size > 20 && $intron_size < 26);
    push(@error4,"Gene warning - $gene has a small intron ($intron_size bp)\n") if ($intron_size > 25 && $intron_size < 31);

  }

  for($i=0; $i<@exon_coord1; $i++){
    my $start = $exon_coord1[$i];
    my $end = $exon_coord2[$i];
    for ($j=$i+1;$j<@exon_coord1;$j++){
      if(($end > $exon_coord1[$j]) && ($start < $exon_coord2[$j])){
	print "Gene error - $gene: exon inconsistency, exons overlap\n" if $verbose;
	push(@error1,"Gene error - $gene: exon inconsistency, exons overlap\n");
      }
    }
  }



  # check that 'Start_not_found' and 'End_not_found' tags present?
  my $start_tag = "";
  my $end_tag = "";

  if ($gene_object->get('Start_not_found')){
    $start_tag = "present";
    push(@error2,"Gene warning - $gene: Start_not_found tag present\n");
    print "Gene warning - $gene: Start_not_found tag present\n" if $verbose;
  }
  
  if ($gene_object->get('End_not_found')){
    $end_tag = "present";
    push(@error2,"Gene warning - $gene: End_not_found tag present\n");
    print "Gene warning - $gene: End_not_found tag present\n" if $verbose;
  }

  # check species is correct
  my $species = "";
  ($species) = ($gene_object->get('Species'));
  push(@error3,"Gene error - $gene: species is $species\n") if ($species ne "Caenorhabditis elegans");  
  print "Gene error - $gene: species is $species\n" if ($species ne "Caenorhabditis elegans" && $verbose);

  # check Method isn't 'hand_built'
  my $method = "";
  ($method) = ($gene_object->get('Method'));
  push(@error3,"Gene error - $gene: method is hand_built\n") if ($method eq 'hand_built');
  print "Gene error - $gene: method is hand_built\n" if ($method eq 'hand_built' && $verbose);

  # check From_laboratory tag is present
  my $laboratory = ($gene_object->From_laboratory);
  push(@error3, "Gene error - $gene: does not have From_laboratory tag\n") if (!defined($laboratory));
  print "Gene error - $gene: does not have From_laboratory tag\n" if (!defined($laboratory) && $verbose);


  # then run misc. sequence integrity checks
  my $dna = $gene->asDNA();
  if(!$dna){
    push(@error1,"Gene error - $gene: can't find any DNA to analyse\n");
    print "Gene error - $gene: can't find any DNA to analyse\n" if $verbose;
    next CHECK_GENE;
  }

  # feed DNA sequence to function for checking
  &test_gene_sequence_for_errors($gene,$start_tag,$end_tag,$dna,$gene_object);

}


# print warnings to log file, log all category 1 errors, and then fill up a top
# 20 with what is left

my $count_errors =0;
foreach my $error (@error1){
  $count_errors++;
  print LOG "$count_errors) $error";
}
while ($count_errors < 20){
  foreach my $error (@error2){
    $count_errors++;
    print LOG "$count_errors) $error";
    last if $count_errors > 20;
  }
  foreach my $error (@error3){
    $count_errors++;
    print LOG "$count_errors) $error";
    last if $count_errors > 20;
  }
}



print LOG "\ncheck_predicted_genes.pl ended at ",`date`,"\n";;
close(LOG);
$db->close;
exit(0);

#################################################################

sub test_gene_sequence_for_errors{

  my $gene = shift;
  my $start_tag = shift;
  my $end_tag = shift;
  my $dna = shift;
  my $gene_object = shift;

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
  my $warning;
  if ($gene_length < 100){
    $warning = "Gene warning - $gene is very short ($gene_length bp), ";
    print "Gene warning - $gene is very short ($gene_length bp), " if $verbose;
    if(defined($gene_object->at('Properties.Coding.Confirmed_by'))){
      $warning .= "gene is Confirmed\n";
      print "gene is Confirmed\n" if $verbose;
    }
    elsif(defined($gene_object->at('Visible.Matching_cDNA'))){
      $warning .= "gene is Partially_confirmed\n";
      print "gene is Partially_confirmed\n" if $verbose;
    }
    else{
      $warning .= "gene is Predicted\n";
      print "gene is predicted\n" if $verbose;
    }
    push(@error3, $warning);
  }
  elsif($gene_length < 150){
    if (defined($gene_object->at('Properties.Coding.Confirmed_by'))){
      $warning = "Gene warning - $gene is short ($gene_length bp) and is Confirmed\n";
      print "Gene warning - $gene is short ($gene_length bp) and is Confirmed\n" if $verbose;
    }    
    elsif(defined($gene_object->at('Visible.Matching_cDNA'))){
      $warning .= "Gene warning - $gene is short ($gene_length bp) and is Partially_confirmed\n";
      print "Gene warning - $gene is short ($gene_length bp) and is Partially_confirmed\n" if $verbose;
    }
    else{
      $warning .= "Gene warning - $gene is short ($gene_length bp) and is Predicted\n";
      print "Gene warning - $gene is short ($gene_length bp) and is Predicted\n" if $verbose;
    }
    push(@error4, $warning);
  }

  if ($remainder != 0){
    if (($end_tag ne "present") && ($start_tag ne "present")){
      push(@error1,"Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n");
      print "Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found & End_not_found tags MISSING\n" if $verbose;
    }
    else{
      push(@error2,"Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n");
      print "Gene error - $gene: length ($gene_length bp) not divisible by 3, Start_not_found and/or End_not_found tag present\n" if $verbose;
    }   
  }   


  # look for incorrect stop codons
  if (($stop_codon ne 'taa') && ($stop_codon ne 'tga') && ($stop_codon ne 'tag')){
    if($end_tag ne "present"){
      push(@error1, "Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n");
      print "Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag MISSING\n" if $verbose;
    }
    else{
      push(@error2,"Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag present\n");
      print "Gene error - $gene: '$stop_codon' is not a valid stop codon. End_not_found tag present\n" if $verbose;
    }
  }


  # look for incorrect start codons
  if ($start_codon ne 'atg'){
   if($start_tag ne "present"){
      push(@error1,"Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n");
      print "Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag MISSING\n" if $verbose;
    }
    else{
      push(@error2, "Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag present\n");
      print "Gene error - $gene: '$start_codon' is not a valid start codon. Start_not_found tag present\n" if $verbose;
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
      push(@error1, "Gene error - $gene: internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n");      
      print "Gene error - $gene: internal stop codon at position $j ...$previous_sequence $offending_codon $following_sequence...\n" if $verbose;   
    }
  }			       

  # look for non-ACTG characters
  if ($dna =~ /[^acgt]/i){
    $dna =~ s/[acgt]//g;
    push(@error2, "Gene error - $gene: DNA sequence contains the following non-ATCG characters: $dna\n"); 
    print "Gene error - $gene: DNA sequence contains the following non-ATCG characters: $dna\n" if $verbose;
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
other acedb) database.  The script will only analyse objects in the 'Predicted_gene' class.
The script emails the top 20 problems each day (sorted by severity).


=over 4

=item MANDATORY arguments: -database <database_path>

This argument must be a path to a valid acedb database, e.g.
check_predicted_genes.pl -database /wormsrv2/camace

=back

=over 4

=item OPTIONAL arguments: -log <logfile>, -verbose

If the file specified by -log already exists, the script will append the output to that file.
Otherwise it will attempt to write to a new file by that name.  If no log file is specified,
the script will generate a log file in /wormsrv2/logs.

If -verbose is specified, output will be written to screen as well as to the log file

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

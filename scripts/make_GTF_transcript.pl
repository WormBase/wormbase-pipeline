#!/software/bin/perl -w
#
# make_GTF_transcript.pl
#
# Ceates a GTF file for reading into Cufflinks of the transcripts of all expressed genes, pseudogene, non-coding transcripts, transposons and their isoforms
#
#
# by Gary Williams
#
# Last updated by: $Author: pad $                      
# Last updated on: $Date: 2011-08-19 12:54:35 $        

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
use Cwd;
use Ace;
use Sequence_extract;
use Coords_converter;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($output, $dbdir);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "output=s"    => \$output,
	    "database=s"  => \$dbdir,
	    "species=s"   => \$species,
	    );

$test = 0;
$debug = "gw3";

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

$dbdir          = $wormbase->test ? $wormbase->database("current") : $wormbase->autoace unless $dbdir; # Database path
my $gffdir      = $wormbase->gff_splits;        # GFF directory
my @chromosomes = $wormbase->get_chromosome_names; # chromosomes

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


if (!defined($output)){
  die "No -output file\n";
}

# open an ACE connection to parse details for mapping to genome
my $tace            = $wormbase->tace;        # TACE PATH
print "Connecting to Ace\n";
my $db = Ace->connect (-path => $dbdir,
                       -program => $tace) || die "cannot connect to database at $dbdir\n";

#my $seq_obj = Sequence_extract->invoke($dbdir, undef, $wormbase);

my $coords  = Coords_converter->invoke($dbdir, undef, $wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################

open (OUTDAT, ">$output") || die "Can't open $output\n";

my $next_start = 1;
my $prev_start = 1;

# foreach gene's transcripts  (coding_transcript, non-coding-transcript, tRNA, ncRNA, etc)
my $transcripts = $db->fetch_many(-query => 'Find Transcript where (NOT Method = history) AND (Species = "'.$wormbase->full_name.'")');
while ( my $transcript = $transcripts->next ) {

  my $gene_id = Transcript_gene($transcript);
  my $transcript_id = $transcript->name;

  # get chromosomal position
  my ($clone, $start, $end) = &get_position($transcript);
  if (! defined $clone) {next;}
  $start--;
  $end--;

  # get the strand
  my $strand = '+';
  if ($start > $end) {
    my $tmp = $end;
    $end = $start;
    $start = $tmp;
    $strand = '-';
  }

  # get the clone's start position in the chromosome
  my ($chromosome, $offset) = $coords->CloneOffset($clone);
  $start += $offset;
  $end += $offset;

  # get exons
  my ($exon_starts, $exon_ends) = &get_exons($transcript);





# GTF format for Cufflinks
# 1	seqname	chrX	Chromosome or contig name
# 2	source	Cufflinks	The name of the program that generated this file (always 'Cufflinks')
# 3	feature	exon	The type of record (always either "transcript" or "exon".
# 4	start	77696957	The leftmost coordinate of this record (where 0 is the leftmost possible coordinate)
# 5	end	77712009	The rightmost coordinate of this record, inclusive.
# 6	score	77712009	The most abundant isoform for each gene is assigned a score of 1000. Minor isoforms are scored by the ratio (minor RPKM/major RPKM)
# 7	strand	+	Cufflinks' guess for which strand the isoform came from. Always one of "+", "-", "."
# 7	frame	.	Cufflinks does not predict where the start and stop codons (if any) are located within each transcript, so this field is not used.
# 8	attributes	...	See below.
# Each GTF record is decorated with the following attributes:
# Attribute	Example	Description
# gene_id	CUFF.1	Cufflinks gene id
# transcript_id	CUFF.1.1	Cufflinks transcript id




  # write transcript
  print OUTDAT "$chromosome\t";
  print OUTDAT "WormBase\t";
  print OUTDAT "transcript\t";
  print OUTDAT "$start\t";
  print OUTDAT "$end\t";
  print OUTDAT "0\t";
  print OUTDAT "$strand\t";
  print OUTDAT ".\t";
  print OUTDAT "gene_id \"$gene_id\";";
  print OUTDAT "transcript_id \"$transcript_id\";\n";


  if ($strand eq '-') {
    @{$exon_starts} = reverse @{$exon_starts};
    @{$exon_ends} = reverse @{$exon_ends};
  }

  # write exons
  for ( my $i = 0; $i < @{$exon_starts}; $i++) {
    my $exon_start;
    my $exon_end;
    if ($strand eq '+') {
      $exon_start = $start + $exon_starts->[$i] - 1;
      $exon_end = $start + $exon_ends->[$i] - 1;
    } else {
      $exon_end = $end - $exon_starts->[$i] + 1;
      $exon_start = $end - $exon_ends->[$i] + 1;      
    }
    print OUTDAT "$chromosome\t";
    print OUTDAT "WormBase\t";
    print OUTDAT "exon\t";
    print OUTDAT "$exon_start\t";
    print OUTDAT "$exon_end\t";
    print OUTDAT "0\t";
    print OUTDAT "$strand\t";
    print OUTDAT ".\t";
    print OUTDAT "gene_id \"$gene_id\";";
    print OUTDAT "transcript_id \"$transcript_id\";\n";
  }
}

my @pseudogenes = $db->fetch(-query => 'Find Pseudogene where species = "'.$wormbase->full_name.'" AND (NOT Method = history) AND (NOT Method = Transposon_Pseudogene)');


foreach my $pseudogene (@pseudogenes) {

  my $gene_id = Pseudogene_gene($pseudogene);
  my $transcript_id = $pseudogene->name;

  # get chromosomal position
  my ($clone, $start, $end) = &get_position($pseudogene);
  if (! defined $clone) {next;}
  $start--;
  $end--;

  # get the strand
  my $strand = '+';
  if ($start > $end) {
    my $tmp = $end;
    $end = $start;
    $start = $tmp;
    $strand = '-';
  }

  # get the clone's start position in the chromosome
  my ($chromosome, $offset) = $coords->CloneOffset($clone);
  #$offset--; # make the offset start at position 0
  $start += $offset;
  $end += $offset;

  # get exons
  my ($exon_starts, $exon_ends) = &get_exons($pseudogene);

  # write pseudogene
  print OUTDAT "$chromosome\t";
  print OUTDAT "WormBase\t";
  print OUTDAT "transcript\t";
  print OUTDAT "$start\t";
  print OUTDAT "$end\t";
  print OUTDAT "0\t";
  print OUTDAT "$strand\t";
  print OUTDAT ".\t";
  print OUTDAT "gene_id \"$gene_id\";";
  print OUTDAT "transcript_id \"$transcript_id\";\n";


  if ($strand eq '-') {
    @{$exon_starts} = reverse @{$exon_starts};
    @{$exon_ends} = reverse @{$exon_ends};
  }

  # write exons
  for ( my $i = 0; $i < @{$exon_starts}; $i++) {
    my $exon_start;
    my $exon_end;
    if ($strand eq '+') {
      $exon_start = $start + $exon_starts->[$i] - 1;
      $exon_end = $start + $exon_ends->[$i] - 1;
    } else {
      $exon_end = $end - $exon_starts->[$i] + 1;
      $exon_start = $end - $exon_ends->[$i] + 1;      
    }
    print OUTDAT "$chromosome\t";
    print OUTDAT "WormBase\t";
    print OUTDAT "exon\t";
    print OUTDAT "$exon_start\t";
    print OUTDAT "$exon_end\t";
    print OUTDAT "0\t";
    print OUTDAT "$strand\t";
    print OUTDAT ".\t";
    print OUTDAT "gene_id \"$gene_id\";";
    print OUTDAT "transcript_id \"$transcript_id\";\n";
  }
}

close(OUTDAT);

$log->mail();
print "Finished" if ($verbose);
exit(0);




###############################
# Prints help and disappears  #
###############################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


#################################################
sub Transcript_gene {
  my ($transcript) = @_;

  my $gene_name;

  my $parent;
  if (defined $transcript->Corresponding_CDS) {
    $parent = $transcript->Corresponding_CDS;
    $gene_name = $parent->Gene;
    return $gene_name;

  } elsif (defined $transcript->Gene) {
    $gene_name = $transcript->Gene;
    return $gene_name;
  }  else {
    die "What is $transcript\n";
  }
  
}


#################################################
sub Pseudogene_gene {
  my ($pseudogene) = @_;

  return $pseudogene->Gene;

}


#################################################
sub Transcript_description {

  my ($transcript) = @_;

  my ($description, $cgc_name, $gene_name, $gene, $type, $method);
  my $full_string = "";

  if($transcript->Method(1)){
    $method = $transcript->Method(1);
  }

  if($transcript->Transcript){
    $type = $transcript->Transcript;
    ($description = $transcript->Transcript(2)) if ($transcript->Transcript(2)); # text field, not always present
  }
  else{
    $type = "";                 # non-coding transcript isoforms have no tag after 'Transcript'
  }

  # set empty text field if $description is empty to prevent -w warnings
  $description = "" if (!defined($description));

  if ($type eq '') { # non-coding transcript isoforms have no tag after 'Transcript'
    $full_string .= $wormbase->full_name('-short' => 1)." non-coding isoform ";
  }
  elsif ($method eq 'tRNAscan-SE-1.23') { # tRNAs
    $full_string .= $wormbase->full_name('-short' => 1)." tRNA ";
  }
  elsif ($type eq 'ncRNA') { # RNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." non-coding RNA ";
  }
  elsif ($type eq 'snRNA') { # snRNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." small nuclear RNA $description ";
  }
  elsif ($type eq 'snoRNA') { # snoRNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." small nucleolar RNA $description ";
  }
  elsif ($type eq 'miRNA') { # miRNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." microRNA ";
  }
  elsif ($type eq 'scRNA') { # scRNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." small cytoplasmic RNA ";
  }
  elsif ($type eq 'piRNA') { # piRNA genes
    $full_string .= $wormbase->full_name('-short' => 1)." piRNA ";
  }

  return $full_string;
}

###############################################################

sub Pseudogene_description {

  my ($pseudogene) = @_;

  my $full_string = "";

  my ($description, $cgc_name, $gene_name, $gene, $type, $method);

  $gene_name = $pseudogene->Gene;
  $gene = $db->fetch(Gene => $gene_name);

  if(defined($gene->CGC_name)){
    $cgc_name = $gene->CGC_name;
  }

  # get type of pseudogene
  $type = $pseudogene->Coding_pseudogene;

  if (($type) || ($cgc_name)) {
    if (($type) && ($cgc_name)) {
      $full_string .= $wormbase->full_name('-short' => 1)." $type pseudogene $cgc_name ";
    }
    elsif ($type) {
      $full_string .= $wormbase->full_name('-short' => 1)." $type pseudogene ";
    }
    elsif ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." pseudogene $cgc_name ";
    }
  }
  else {
    $full_string .= $wormbase->full_name('-short' => 1)." predicted pseudogene ";
  }

  return $full_string;
}



##########################################

#    my ($clone, $start, $end) = &get_position($pseudogene);
# get the clone that it is in
# look in the clone for the transcript or pseudogene
# get the start/end position

sub get_position {
  my ($obj) = @_;

  my $clone_obj = $obj->Sequence;
  if (!defined $clone_obj) {print "No Sequence tag for ",$obj->name,"\n"; return undef;}

  my @transcript_ids =  $clone_obj->at('Smap.S_child.Transcript[1]');
  my @transcript_starts = $clone_obj->at('Smap.S_child.Transcript[2]');
  my @transcript_ends = $clone_obj->at('Smap.S_child.Transcript[3]');

  for (my $i=0; $i<@transcript_ids; $i++) {
    if ($transcript_ids[$i]->name eq $obj->name) {return ($clone_obj->name, $transcript_starts[$i]->name, $transcript_ends[$i]->name)}
  }

  my @pseudogene_ids =  $clone_obj->at('Smap.S_child.Pseudogene[1]');
  my @pseudogene_starts = $clone_obj->at('Smap.S_child.Pseudogene[2]');
  my @pseudogene_ends = $clone_obj->at('Smap.S_child.Pseudogene[3]');

  for (my $i=0; $i<@pseudogene_ids; $i++) {
    if ($pseudogene_ids[$i]->name eq $obj->name) {return ($clone_obj->name, $pseudogene_starts[$i]->name, $pseudogene_ends[$i]->name)}
  }

  die "can't find the position of ",$obj->name,"\n";
}

##########################################
# get exons
#    my ($exons_start, $exons_end) = &get_exons($pseudogene);
# positions start at 1
sub get_exons {
  my ($obj) = @_;

  my @exons_start = $obj->Source_exons;
  my @exons_end = $obj->Source_exons(2);


  return (\@exons_start, \@exons_end);
}

##########################################

__END__

=pod

=head2 NAME - find_intergenic.pl

=head1 USAGE

=over 4

=item find_intergenic.pl [-options]

=back

find_intergenic.pl writes out intergenic sequences in Fasta format. By
default it will write out all sequences between gene transcripts that
do not overlap other genes on either strand. It can optionally be
restricted to write out only up to a specified length of region around
genes. It can be restricted to only write the region around the 5' or
around the 3' end of genes. If only the proximal region around the
genes are being printed out, then the sequence will be output in the
same sense as that gene. By default, if the whole intergenic region is
being printed, only the forward sense sequence is output. (What should
the sense of a sequence spanning the region between a forward and a
reverse sense gene be?).

find_intergenic.pl mandatory arguments:


=over 4

=item -output, Specifies the output file name

=back

find_intergenic.pl optional arguments:

=over 4

=item -proximity integer, Restricts the output to proximal regions.

=over 4

If the intergenic region is larger than the integer parameter, then
the output will be restricted to the length of that integer, if the
intergenic region is smaller than the integer parameter, then the
smaller value will be used. The default value is 0 which turns the
proximity restriction off. If the proximal intergenic regions of
adjacent genes overlap, then that overlapping section of the genome
will be output twice, once for each gene.

=back

=item -side "both" | "5" | "3", Output one side or other of the gene.

=over 4

Restricts the output to either just
the 5' side of the gene, or the 3' side, or (the default) writes both
sides of the gene. This only has effect if the '-proximity' option is
set to a non-zero size.

=back

=item -operons "include" | "only" | "no", Include operons or not.

=over 4

If the parameter is set to "include" (the default), then operons have
no effect on the output - they are ignored. If it is set to "only",
then only those genes that are within an operon will be considered, so
only intergenic regions of genes within operons will be output. If it
is set to "no", then intergenic regions within operons are ignored and
not output.

=back

=item -debug, Verbose/Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=back

=cut

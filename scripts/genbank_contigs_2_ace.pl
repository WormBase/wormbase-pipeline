#!/software/bin/perl -w

# got the Brugia data from
# http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=154235816
# from this page, get the WGS contig data in Genbank format
# and the WGS_SCAFLD supercontig data in Genbank format
# (Don't select the latter with sequences, you will loose the CONTIG line)
# put into the files 'contigs.genbank' and 'supercontigs.genbank'
#
# the mitochondrial genome sequencde by NITE is available from
# http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome&cmd=Retrieve&dopt=Overview&list_uids=16713


use strict;
use Bio::Seq;
use Bio::SeqIO;
use IO::String;
use Carp;
 
my $test = 0;			# only read small test files

my $species_prefix = "Bmal_";	# species-specific prefix for contig and supercontig names

my $contigs_file = 'contigs.genbank';
my $supercontigs_file = 'supercontigs.genbank';
if ($test) {
  $contigs_file = 'c.test';
  $supercontigs_file = 'sc.test';
}

my %acc_start;			# start positions of the contigs
my %acc_stop;			# stop positions of the contigs
my %acc_strand;			# strand of the contig
my %acc_len;			# length of the contig

my %acc_2_id;			# contig accession to id conversion
$acc_2_id{gap} = 'gap';	# initialise the accesion to contig conversion hash with gap to gap


my $agp_file = "supercontigs.agp";
open (AGP, ">$agp_file") || croak "Can't open $agp_file\n";

my $ace_file = "output.ace";
open (ACE, ">$ace_file") || croak "Can't open $ace_file\n";

print "parsing supercontig CONTIG data\n";
&parse_supercontigs_file();

print "parsing and writing contigs\n";
&parse_contigs();

print "parsing and writing supercontigs\n";
&write_supercontigs();


close(ACE);
close(AGP);
exit(0);


#################################
#
# Subroutines
#
#################################

#################################
# parse supercontig data for the start/stop/strand information in the
# CONTIG lines

sub parse_supercontigs_file {

  # get the sequences from the supercontig genbank file
  my $seqio_supercontigs = Bio::SeqIO->new(-file => $supercontigs_file,
					   -format => "genbank" );

  # read the supercontig data
  while(my $seq_obj = $seqio_supercontigs->next_seq) {

    # convert the ID from the genbank to the Brugia supercontig name
    my $supercontig_name = &get_supercontig_name($seq_obj);

# get the contigs in this supercontig
#tagname : CONTIG
# Value: join(AAQA01007080.1:1..1460,gap(7342),AAQA01008861.1:1..1290)
# or
#tagname : CONTIG
# Value: join(AAQA01003469.1:1..3114)


    my $contig_line = "";
    my $anno = $seq_obj->annotation;
    my @annotations = $anno->get_Annotations('CONTIG');
    for my $value ( @annotations ) {
      $contig_line .= $value->as_text;
    }
    $contig_line =~ s/Value: //g;

    # parse the CONTIG line

    # make the dummy Genbank entry, with ID 'gap' to help highlight the gaps
    my $str="LOCUS       gap
FEATURES             Location/Qualifiers
     contig          $contig_line
ORIGIN
";

    my $io = IO::String->new(\$str);
    my $seqio_dummy = Bio::SeqIO->new(-format => 'genbank',
				      -fh     => $io );

    my @overlap_output;		# output for ACE overlap details

    while(my $seq_dummy = $seqio_dummy->next_seq) {
      my $start = 0;		# start position of contig in supercontig
      my $end = 0;		# end position of contig in supercontig

      for my $feature ($seq_dummy->top_SeqFeatures) {
	my $ft_type = $feature->primary_tag;

	for my $location ( $feature->location->sub_Location ) {
	  my $acc = $location->seq_id();
	  # remove the version number from the end of the contig's accession number
	  $acc =~ s/\.\d+$//;
	  next if ($acc eq 'gap');
	  if (exists $acc_start{$acc}) {print "ERROR: found $acc used in $supercontig_name for second time\n";};		# check all contigs are only used once
	  # store the results
	  $acc_strand{$acc} = $location->strand;
	  $acc_start{$acc} = $location->start;
	  $acc_stop{$acc} = $location->end;
	}
      }
    }
  }
}

#################################
# creates the supercontig ID name from the information in the DESCRIPTION line

sub get_supercontig_name {
  my ($seq_obj) = @_;

  my $desc = $seq_obj->desc;
  my ($supercontig_number) = ($desc =~ /scf_(\d+)/);
  # Brugia also defines the scaffold name like this when the contigs are dgenerate (low quality or low coverage):
  if (! defined $supercontig_number) {
    ($supercontig_number) = ($desc =~ /degen_scaff_(\d+)/);
    # make the supercontig_number unique and add indication that this
    # is a low quality sequence
    $supercontig_number = 'Degenerate' . $supercontig_number; 
  }
  my $supercontig_name = $species_prefix . "supercontig" . $supercontig_number; # explicitly start the name with supercontig

  return $supercontig_name;
}

#################################
# creates the contig ID name from the information in the DESCRIPTION line

sub get_contig_name {
  my ($seq_obj) = @_;

  my $desc = $seq_obj->desc;
  my ($contig_number) = ($desc =~ /ctg_(\d+)/);
  my $contig_name = $species_prefix . "contig" . $contig_number; # explicitly start the name with 'contig'
  
  return $contig_name;
}

#################################
# parse contig data

sub parse_contigs {

  # get the sequences from the contig genbank file
  my $seqion_contigs = Bio::SeqIO->new(-file => $contigs_file,
				       -format => "genbank" );


  # read the contig data
  while(my $seq_obj = $seqion_contigs->next_seq) {

    # convert the ID from the genbank to the Brugia contig name
    my $id = $seq_obj->display_id;
    my $sv = $seq_obj->seq_version;
    my $acc = $seq_obj->accession;
    my $contig_name = &get_contig_name($seq_obj);
    # get just the sub-sequence we need to make the supercontig without problematical overlaps which might not overlap!
    my $start = $acc_start{$acc};
    my $stop  = $acc_stop{$acc};
    my $dna = $seq_obj->subseq($start, $stop); 
    my $dna_len = length $dna;
    # sanity check
    if ($dna_len != $stop - $start + 1) {print "ERROR contig $id is $dna_len bp long but we need $start..$stop for the supercontig\n";}
    $acc_len{$acc} = $dna_len;
    my $species = $seq_obj->species->binomial;

    $acc_2_id{$acc} = $contig_name;

    # write the contig
    print ACE "\n";
    print ACE "Sequence  :  \"$contig_name\"\n";
    print ACE "DNA          \"$contig_name\" $dna_len\n";
    print ACE "Database     \"EMBL\" \"NDB_ID\" \"$id\"\n";
    print ACE "Database     \"EMBL\" \"NDB_AC\" \"$acc\"\n";
    print ACE "Database     \"EMBL\" \"NDB_SV\" \"${acc}.${sv}\"\n";
    print ACE "Species      \"$species\"\n";
    print ACE "Clone        \"$contig_name\"\n";
    print ACE "Genomic_canonical\n";
    print ACE "Method       \"Genomic_canonical\"\n";
    
    $dna = reformat($dna);
    print ACE "\n";
    print ACE "DNA : $contig_name\n";
    print ACE "$dna\n";
  }
}


##########################################

sub write_supercontigs {

  # get the sequences from the supercontig genbank file
  my $seqio_supercontigs = Bio::SeqIO->new(-file => $supercontigs_file,
					   -format => "genbank" );

  # read the supercontig data
  while(my $seq_obj = $seqio_supercontigs->next_seq) {

    # convert the ID from the genbank to the Brugia supercontig name
    my $supercontig_name = &get_supercontig_name($seq_obj);

    my $species = $seq_obj->species->binomial;
    
    # get the contigs in this supercontig
    my $contig = "";
    my $anno = $seq_obj->annotation;
    my @annotations = $anno->get_Annotations('CONTIG');
    for my $value ( @annotations ) {
      $contig .= $value->as_text;
    }
    $contig =~ s/Value: //g;
    #print "CONTIG: $contig\n";
    
    # write the supercontig
    print ACE "\n";
    print ACE "Sequence : \"$supercontig_name\"\n";
    print ACE "Species  \"$species\"\n";
    print ACE "Method   \"Link\"\n";
    if ($supercontig_name =~ /Degenerate/) {
      print ACE "Remark   \"This is assembled from degenerate (low quality or low coverage) contigs.\"\n";
    }
    # parse the CONTIG line and write the Supercontig's Subsequence data and AGP file
    # and the Contigs' Source data
    &parse_contig($supercontig_name, $contig, %acc_2_id); # this completes the supercontig Sequence ACE output block and adds others

    ####################################################################
    #
    # write the genes, transcripts, CDS, operons, tRNA, Transposon_CDS, 
    # Proteins etc. from the FEATURES part of the GenBank entry
    #
    ####################################################################
    for my $feature ($seq_obj->top_SeqFeatures) {
      my $ft_type = $feature->primary_tag;
      
      if ($ft_type eq 'source') {next}
      if ($ft_type eq 'gap') {next}
      print ACE "\n// $ft_type\n";

      my $parent_tag;		# the tag in the supercontig that the feature is added to
      my $id;			# the ID of the feature
      my $gene_id;		# the ID of the gene it belongs to
      my $codon_start = 0;	# the default value for a CDS codon_start tag

      # gene
      if ($ft_type eq "gene") {
	$parent_tag = "Gene_child";
	$gene_id = &get_gene_id($supercontig_name, $feature);
	$id = $gene_id;
	print ACE "\n";
	print ACE "Gene : \"$gene_id\"\n";
	print ACE "Sequence \"$supercontig_name\"\n";
	print ACE "Sequence_name \"$gene_id\"\n";
	print ACE "Public_name \"$gene_id\"\n";
	print ACE "Species \"$species\"\n";
	print ACE "Live\n";
	print ACE "Method Gene\n";
	print ACE "Remark \"Defined by TIGR.\"\n";
      }

      # transcript
      if ($ft_type eq "mRNA") {
	$parent_tag = "Transcript";
	$gene_id = &get_gene_id($supercontig_name, $feature);
	$id = &get_id($supercontig_name, $feature);
	print ACE "\n";
	print ACE "Transcript : \"$id\"\n";
	print ACE "Sequence \"$supercontig_name\"\n";
	print ACE "Corresponding_CDS $id\n";
	print ACE "Species \"$species\"\n";
	print ACE "Method Coding_transcript\n";
	print ACE "Remark \"Defined by TIGR.\"\n";
	&write_exons($feature);
      }

      # tRNA (a type of Transcript)
      if ($ft_type eq 'tRNA') {
	$parent_tag = "Transcript";
	$gene_id = &get_gene_id($supercontig_name, $feature);
	$id = $gene_id;
	my $product = &get_tag('product', $feature);
	my $note = &get_tag('note', $feature);
	my ($aa) = ($product =~ /tRNA-(\S+)/);
	if (! defined $aa) {print "problem with aa in $supercontig_name tRNA\n";}
	my ($anticodon) = ($note =~ /anticodon\s+(\S+)/);
	if (! defined $anticodon) {print "problem with anticodon in $supercontig_name tRNA\n";}
	print ACE "\n";
	print ACE "Transcript : \"$id\"\n";
	print ACE "Sequence \"$supercontig_name\"\n";
	print ACE "Gene $id\n";
	print ACE "Brief_identification $product\n";
	print ACE "Transcript tRNA Type $aa\n";
	print ACE "Transcript tRNA Anticodon $anticodon\n";
	print ACE "Species \"$species\"\n";
	print ACE "Method tRNAscan-SE-1.23\n";
	print ACE "Remark \"Defined by TIGR.\"\n";
	&write_exons($feature);
      }

      # CDS
      if ($ft_type eq 'CDS') {
	$parent_tag = "CDS_child";
	$gene_id = &get_gene_id($supercontig_name, $feature);
	$id = &get_id($supercontig_name, $feature);
	my $product = &get_tag('product', $feature);
	$codon_start = &get_tag('codon_start', $feature);
	$codon_start--;		# change to 0, 1 or 2
	my $method;
	if ($product =~ /transposon/i) {
	  # Transposons don't have gene_id's, so remove this transposon's gene
	  print ACE "\n";
	  print ACE "-D Gene : \"$gene_id\"\n";
	  $method = "Transposon_CDS";
	} else {
	  $method = "JIGSAW";
	}
	print ACE "\n";
	print ACE "CDS : \"$id\"\n";
	print ACE "Sequence \"$supercontig_name\"\n";
	print ACE "Corresponding_transcript $id\n";
	print ACE "CDS\n";
	print ACE "Gene \"$gene_id\"\n";
	print ACE "Species \"$species\"\n";
	print ACE "Method $method\n";
	# +++ need to add the DB_info stuff for the protein
	print ACE "Remark \"Defined by TIGR.\"\n";
	&write_exons($feature, 'cds', $codon_start);
      }
      
      # repeat_region
      if ($ft_type eq 'repeat_region') {
	
      }
      

      # write the overall start and end positions in the parent sequence
      my $overall_start = $feature->location->start;
      my $overall_end = $feature->location->end;
      my $overall_strand = $feature->strand;
      if ($overall_strand == -1) { # swap start and end if reverse
	my $temp = $overall_start;
	$overall_start = $overall_end;
	$overall_end = $temp;
	$overall_start -= $codon_start; # adjust CDS codon start
      } else {
	$overall_start += $codon_start; # adjust CDS codon start
      }
      if ($ft_type ne 'repeat_region') {
	print ACE "\n";
	print ACE "Sequence : \"$supercontig_name\"\n";
	print ACE "$parent_tag $id $overall_start $overall_end\n";
      }
    }
    
    print ACE "\n";
    
  }				# iterate over supercontig entries

}


##########################################
# reformat a DNA string into lines of 60 characters

sub reformat {
    my $in_string = shift;
    my $out_string = "";

    my $string_len = length ($in_string);
    my $lines = int ($string_len / 60) ;

    for (my $i = 0; $i <= $lines; $i++) {
        $out_string = $out_string . substr($in_string,($i*60),60) . "\n";
    }
    return ($out_string);
}

##########################################
# parse the CONTIG line by creating a dummy GenBank-format string
# containing the CONTIG line as a Feature and reading it into bioperl
# using the bioperl parsing code to load it into a sequence's feature
# object.

sub parse_contig {
  my ($supercontig_name, $contig_line, %acc_2_id) = @_;

  my @contigs;			# list of the contigs in this supercontig

  # make the dummy Genbank entry, with ID 'gap' to help highlight the gaps
  my $str="LOCUS       gap
FEATURES             Location/Qualifiers
     contig          $contig_line
ORIGIN
";

  my $io = IO::String->new(\$str);
  my $seqio_dummy = Bio::SeqIO->new(-format => 'genbank',
				    -fh     => $io );

  while(my $seq_dummy = $seqio_dummy->next_seq) {
    my $pos1 = 0;		# start position of contig in supercontig
    my $pos2 = 0;		# end position of contig in supercontig
    my $number = 1;		# ordinal number of the contig or gap in this supercontig

    for my $feature ($seq_dummy->top_SeqFeatures) {
      my $ft_type = $feature->primary_tag;
      #print "found a contig line feature type: $ft_type\n";

      for my $location ( $feature->location->sub_Location ) {
	my $acc = $location->seq_id();
	# remove the version number from the end of the contig's accession number
	$acc =~ s/\.\d+$//;
	my $contig = $acc_2_id{$acc}; # get the contig name from the accession number
	push @contigs, $contig;
	my $strand = $location->strand;
	my $start = $location->start; 
	my $end = $location->end;     
	if ($contig eq 'gap') {
	  $start = 1;
	  $strand = 1;		# gaps are on '-' strand for some reason. force them to be '+'
	}
	#print "$contig $start $end $strand\n";

# write the AGP file
	my $len = $end - $start + 1;
	$pos1 = $pos2 + 1;
	$pos2 = $pos1 + $len - 1;
	my $type = ($contig eq 'gap') ? 'N' : 'W';  # 'N' indicates a gap in the AGP file, 'W' is contig
	$strand = ($strand eq '1') ? '+' : '-';
	if ($contig eq 'gap') {
	  print AGP "$supercontig_name\t$pos1\t$pos2\t$number\t$type\t$len\tfragment\tyes\n";

	} else {
	  print AGP "$supercontig_name\t$pos1\t$pos2\t$number\t$type\t$contig\t$start\t$len\t$strand\n";

	  # write the supercontig Subsequence data for this contig (not done for gaps)
	  if ($strand eq '+') { 
	    print ACE "Subsequence \"$contig\" $pos1 $pos2\n";
	  } else {
	    print ACE "Subsequence \"$contig\" $pos2 $pos1\n";
	  }
	}
	$number++;
      }
    }
  }

  #  write the Source data for all the contigs in the Supercontig
  foreach my $contig (@contigs) {
    if ($contig ne 'gap') {
      print ACE "\n";
      print ACE "Sequence : \"$contig\"\n";
      print ACE "Source $supercontig_name\n";
    }
  }
}

##########################################
# output the source exon positions
# $cds is only defined if this is a CDS
# $codon_start is only defined if there is a codon_start tag

sub write_exons {
  my ($feature, $cds, $codon_start) = @_;

  if (! defined $codon_start) {$codon_start = 0;}

  # is it a split location, or a simple one?
  my $start_exon;
  my $end_exon;
  my $strand = $feature->strand;

  if ( $feature->location->isa('Bio::Location::SplitLocationI')) {  
    my $origin;			# position of first exon's start
    for my $location ( $feature->location->sub_Location ) {
      my $start_pos_type = $location->start_pos_type();
      my $end_pos_type = $location->end_pos_type();
      if ($start_pos_type eq 'BEFORE' && defined $cds) {
	print ACE "Start_not_found\n";
      }
      if ($end_pos_type eq 'AFTER' && defined $cds) {
	print ACE "End_not_found\n";
      }
      my $start = $location->start;
      my $end = $location->end;

      if (! defined $origin) {
	if ($strand == 1) {
	  $start += $codon_start;
	  $origin = $start;
	} else {
	  $origin = $feature->end; # the end of the gene
	  $origin -= $codon_start;
	}
      }
      
      if ($strand == 1) {
	$start_exon = $start - $origin + 1;
	$end_exon = $start_exon + $end - $start; # start_exon + len
	print ACE "Source_exons $start_exon $end_exon\n";
      } else {
	if ($end > $origin) {$end = $origin;} # adjust the end if we have changed the origin for codon_start
	$end_exon = $origin - $end + 1;
	$start_exon = $end_exon + ($end - $start); # end_exon + len
	print ACE "Source_exons $end_exon $start_exon\n";
      }
    }
  } else {
    my $start = $feature->start;
    my $end = $feature->end;
    if ($strand == 1) {
      $start += $codon_start;
    } else {
      $end -= $codon_start;
    }
    $start_exon = 1;
    $end_exon = $end - $start + 1; # len
    print ACE "Source_exons $start_exon $end_exon\n";
  }

}  

##########################################
# get the value of the /locus_tag
sub get_gene_id {
  my ($supercontig_name, $feature) = @_;
  my $gene_id;
  if ($feature->has_tag('locus_tag')) {
    for my $val ($feature->get_tag_values('locus_tag')){
      $gene_id = $val;
    }
  }
  if (!defined $gene_id) {carp "can't find in $supercontig_name\n";}
  return $gene_id;
}
##########################################
# get ID from /note
sub get_id {
  my ($supercontig_name, $feature) = @_;
  my $id;
  if ($feature->has_tag('note')) {
    for my $val ($feature->get_tag_values('note')){
      if ($val =~ /transcript\s+(\S+)/) {
	$id = $1;
	$id =~ s/;//;		# remove any ';' at the end
      }
    }
  }
  if (!defined $id) {carp "can't find /note in $supercontig_name\n";}
  return $id;
}    
##########################################
# get the value of a tag such as /product
sub get_tag {
  my ($tag, $feature) = @_;

  if ($feature->has_tag($tag)) {
    for my $val ($feature->get_tag_values($tag)){
      return $val;
    }
  }
}    
##########################################


__END__



# overlap_right = name of clone to right of this one, position in this
#                 clone that the next clone's finished sequence starts at
#
# clone_right_end = name of clone to the left of this one, position in
#                   this clone that the previous clone ends at 
#                   - but this is unfinished sequence so we don't want to use this value to find overlaps
#
#                   and/or, name of this clone, length of this clone
#
# clone_left_end = name of clone to right of this one, position in
#                  this clone that the next clone starts at
#                  - but this is unfinished sequence so we don't want to use this value to find overlaps
#
#                   and/or, name of this clone, 1
#
# overlap_left = name of clone to the left of this one
#                - look in this clone to get its length and overlap_right value
#		   so that we can see the amount of overlap with this clone
#
# ------ = finished sequence
# ...... = unfinished sequence
#                 
#                 
#                 
#                 
#------------------------|..........|                                                                         previous clone
#                 
#          |.........|--------------V--------------------------V-------------------V-------|.......|          this clone
#                                   clone                      clone               overlap
#                                   right                      left                right
#                                   end                        end
#                                                              |...................|-----------------------   next clone
#                                                             

#!/usr/local/perl -w

use lib "/nfs/farm/Worms/ensembl/modules";
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::Runnable::Blat;
#use Bio::EnsEMBL::Analysis;
#use Bio::EnsEMBL::Root;

my $results = shift;  

my $dummy1 = Bio::Seq->new( -display_id => 'dummy1',
			     -seq => "acgt");

my $dummy2 = Bio::Seq->new( -display_id => 'dummy2',
			     -seq => "attatatatattata");

my @sequences = ($dummy1, $dummy2);
my $blat_dummy = Bio::EnsEMBL::Pipeline::Runnable::Blat->new(
							     -database    => "dummy",
							     -query_seqs  => \@sequences,
							     -blat        => "/nfs/farm/Worms/bin/blat"
							    );

my @features_within_features;
open (BLAT ,"<$results") or die "$results file\n";

while (<BLAT>) {

  ############################################################
  #  PSL lines represent alignments and are typically taken from files generated 
  # by BLAT or psLayout. See the BLAT documentation for more details. 
  #
  # 1.matches - Number of bases that match that aren't repeats 
  # 2.misMatches - Number of bases that don't match 
  # 3.repMatches - Number of bases that match but are part of repeats 
  # 4.nCount - Number of 'N' bases 
  # 5.qNumInsert - Number of inserts in query 
  # 6.qBaseInsert - Number of bases inserted in query 
  # 7.tNumInsert - Number of inserts in target 
  # 8.tBaseInsert - Number of bases inserted in target 
  # 9.strand - '+' or '-' for query strand. In mouse, second '+'or '-' is for genomic strand 
  #10.qName - Query sequence name 
  #11.qSize - Query sequence size 
  #12.qStart - Alignment start position in query 
  #13.qEnd - Alignment end position in query 
  #14.tName - Target sequence name 
  #15.tSize - Target sequence size 
  #16.tStart - Alignment start position in target 
  #17.tEnd - Alignment end position in target 
  #18.blockCount - Number of blocks in the alignment 
  #19.blockSizes - Comma-separated list of sizes of each block 
  #20.qStarts - Comma-separated list of starting positions of each block in query 
  #21.tStarts - Comma-separated list of starting positions of each block in target 
  ############################################################

  # first split on spaces:
  chomp;

  my (
      $matches,      $mismatches,    $rep_matches, $n_count, $q_num_insert, $q_base_insert,
      $t_num_insert, $t_base_insert, $strand,      $q_name,  $q_length,     $q_start,
      $q_end,        $t_name,        $t_length,    $t_start, $t_end,        $block_count,
      $block_sizes,  $q_starts,      $t_starts
     )
    = split;

  my $superfeature = Bio::EnsEMBL::SeqFeature->new();

  # ignore any preceeding text
  unless ( $matches and $matches =~/^\d+$/ ){
    next;
  }

  #print STDERR $_."\n";

  # create as many features as blocks there are in each output line
  my (%feat1, %feat2);
  $feat1{name} = $t_name;
  $feat2{name} = $q_name;

  $feat2{strand} = 1;
  $feat1{strand} = $strand;

  # all the block sizes add up to $matches + $mismatches + $rep_matches

  # percentage identity =  ( matches not in repeats + matches in repeats ) / ( alignment length )
  #print STDERR "calculating percent_id and score:\n";
  #print STDERR "matches: $matches, rep_matches: $rep_matches, mismatches: $mismatches, q_length: $q_length\n";
  #print STDERR "percent_id = 100x".($matches + $rep_matches)."/".( $matches + $mismatches + $rep_matches )."\n";
  my $percent_id = sprintf "%.2f", ( 100 * ($matches + $rep_matches)/( $matches + $mismatches + $rep_matches ) );

  # or is it ...?
  ## percentage identity =  ( matches not in repeats + matches in repeats ) / query length
  #my $percent_id = sprintf "%.2d", (100 * ($matches + $rep_matches)/$q_length );

  # we put basically score = coverage = ( $matches + $mismatches + $rep_matches ) / $q_length
  #print STDERR "score = 100x".($matches + $mismatches + $rep_matches)."/".( $q_length )."\n";

  unless ( $q_length ){
    $self->warn("length of query is zero, something is wrong!");
    next;
  }
  my $score   = sprintf "%.2f", ( 100 * ( $matches + $mismatches + $rep_matches ) / $q_length );

  # size of each block of alignment (inclusive)
  my @block_sizes     = split ",",$block_sizes;

  # start position of each block (you must add 1 as psl output is off by one in the start coordinate)
  my @q_start_positions = split ",",$q_starts;
  my @t_start_positions = split ",",$t_starts;

  $superfeature->seqname($q_name);
  $superfeature->score( $score );
  $superfeature->percent_id( $percent_id );

  # each line of output represents one possible entire aligment of the query (feat1) and the target(feat2)
  for (my $i=0; $i<$block_count; $i++ ) {

    #### working out the coordinates: #########################
    #
    #                s        e
    #                ==========   EST
    #   <----------------------------------------------------| (reversed) genomic of length L
    #
    #   we would store this as a hit in the reverse strand, with coordinates:
    #
    #   |---------------------------------------------------->
    #                                   s'       e'
    #                                   ==========   EST
    #   where e' = L  - s  
    #         s' = e' - ( e - s + 1 ) + 1
    #
    #   Also, hstrand will be always +1
    ############################################################

    my ($query_start,$query_end);

    if ( $strand eq '+' ) {
      $query_start = $q_start_positions[$i] + 1;
      $query_end   = $query_start + $block_sizes[$i] - 1;
    } else {
      $query_end   = $q_length  - $q_start_positions[$i];
      $query_start = $query_end - $block_sizes[$i] + 1;
    }

    #$feat2 {start} = $q_start_positions[$i] + 1;
    #$feat2 {end}   = $feat2{start} + $block_sizes[$i] - 1;
    $feat2 {start} = $query_start;
    $feat2 {end}   = $query_end;
    if ( $query_end <  $query_start ) {
      $self->warn("dodgy feature coordinates: end = $query_end, start = $query_start. Reversing...");
      $feat2 {end}   = $query_start;
      $feat2 {start} = $query_end;
    }

    $feat1 {start} = $t_start_positions[$i] + 1;
    $feat1 {end}   = $feat1{start} + $block_sizes[$i] - 1;

    # we put all the features with the same score and percent_id
    $feat2 {score}   = $score;
    $feat1 {score}   = $feat2 {score};
    $feat2 {percent} = $percent_id;
    $feat1 {percent} = $feat2 {percent};

    # other stuff:
    $feat1 {db}         = undef;
    $feat1 {db_version} = undef;
    $feat1 {program}    = 'blat';
    $feat1 {p_version}  = '1';
    $feat1 {source}     = 'blat';
    $feat1 {primary}    = 'similarity';
    $feat2 {source}     = 'blat';
    $feat2 {primary}    = 'similarity';

    my $feature_pair = $blat_dummy->create_FeaturePair(\%feat1,\% feat2);
    $superfeature->add_sub_SeqFeature( $feature_pair,'EXPAND');
  }
  push(@features_within_features, $superfeature);
}
close BLAT;

print "\n";
print "Features created:\n";
foreach my $superf ( @features_within_features ){
  foreach my $subf ( $superf->sub_SeqFeature ){
    print $subf->gffstring."\n";
  }
}


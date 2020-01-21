#!/software/bin/perl -w
#
# transcript_builder.pl
# 
# by Gary Williams
#
# Script to make ?Transcript objects
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-12-23 13:14:17 $
use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Data::Dumper;
use Coords_converter;
use Wormbase;
use Modules::SequenceObj;
use Modules::Transcript;
use Modules::CDS;
use Modules::Strand_transformer;
use Modules::Overlap;
use File::Path;
use Storable;
use Ace;

my ($debug, $store, $help, $verbose, $really_verbose, @test_est,
    $database, $test, @chromosomes, $chunk_total, $chunk_id,
    $gff_dir, $transcript_dir, $ace_fname, $problem_fname,
    @test_cds, $wormbase, $species, $detailed_debug);

my $gap = 15;			# $gap is the gap allowed in an EST alignment before it is considered a "real" intron
my $COVERAGE_THRESHOLD = 95.0;  # the alignment score threshold below which any cDNAs that overlap two genes will not be considered.

# to track failings of system calls
my $errors = 0;

GetOptions ( "debug:s"      => \$debug,
	     "help"             => \$help,
	     "verbose"          => \$verbose,
	     "really_verbose"   => \$really_verbose,
	     "est=s@"           => \@test_est, # only use the specified EST sequence - for debugging
	     "gap:s"            => \$gap,
	     "database:s"       => \$database,
	     "test"             => \$test,
	     "chromosome:s"     => \@chromosomes,
             "chunktotal:s"     => \$chunk_total,
             "chunkid:s"        => \$chunk_id,
	     "cds=s@"           => \@test_cds, # only use the specified CDS object - for debugging
	     "store:s"          => \$store,
	     "species:s"        => \$species,
             "gffdir:s"         => \$gff_dir,
             "transcriptdir:s"  => \$transcript_dir,
             "acefname:s"       => \$ace_fname,
             "problemfname:s"   => \$problem_fname,
             "detaileddebug"    => \$detailed_debug,
	   ) ;

#@test_cds = ("F55A4.13", "F55A4.2"); # testing bog-standard single-exon gene with nearby multi-exon gene with overlapping ESTs between them
#@test_cds = ("W04G3.6a", "W04G3.6b", "W04G3.6c", "W04G3.6d"); # 'a' isoform has multiple transcripts reverted back to the CDS structure

#@test_cds = ("B0025.1c", "B0025.1a"); 
#@test_cds = ("T03D8.1a", "T03D8.1c"); # different START codons in the same CDS exon specified by publication
#@test_cds = ("F31D4.1"); # bog-standard CDS with lots of matching ESTs - no TSLs
#@test_cds = ("F13H10.2b"); # TSL Feature - but the current coding transcript based on EST doesn't get up to the TSL Feature
#@test_cds = ("ZK593.5c"); # TSL Feature - but the START is in the next exons after an intron after the Features
#@test_cds = ("B0024.9a", "B0024.9b"); # two isoforms sharing the same Transcript structure
#@test_cds = ("C14A4.4a", "C14A4.4b"); # two isoforms sharing the same Transcript structure
#@test_cds = ("ZK970.2a", "ZK970.2b", "ZK970.2d"); # two isoforms (a/d) sharing the same Transcript structure - isoform 'a' has an incorrect Feature so should end up with same Transcript as 'd' (correct this before next Build!)
#@test_cds = ("C09G5.8a", "C09G5.8b"); # two isoforms 'a' is full-length, 'b' is defined by a SL1 but currently the 'b' coding transcript is extended up by three exons that look bogus.
#@test_cds = ("R03D7.6a", "R03D7.6b", "R03D7.6c"); # three isoforms 'a' & 'b' are full-length and bogusly short, 'c' is defined by a SL1 but currently the 'c' coding transcript is extended up by two exons that look bogus.
#@test_cds = ("C08G5.4a", "C08G5.4b", "C08G5.4c"); # bog-standard normal set of isoforms


# want to check if a CDS has a TSL asserted as evidence and when cDNAs are added to the structire, truncate the cDNA structure so the transcript is never longer than the TSL site.
# add asserted Isoform Feature_evidence to the $cds objects as $cds->SL data
# when new transcript objects are made from the $cds objects, add the  $cds->SL data to the transcript object
# when a transcript object is extended, check to see if there is a SL feature and only add exons after this - truncate exons overlapping the SL site


#$debug = "gw3";
##$detailed_debug = 1;
#$test = 1;
if (defined $test) {
  print "WORKING IN TEST_BUILD\n";
}


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
                           );
}

print "ACE connect\n" if ($debug);
my $db;
if (not defined $database or $database eq "autoace") {
  $db = $wormbase->autoace;
} else {
  $db = $database
}
my $ace = Ace->connect(-path => $db) || die(Ace->error);

my $tace = $wormbase->tace;

my %cds2gene = $wormbase->FetchData('cds2wbgene_id');
$species = $wormbase->species if not defined $species;

# other variables and paths.
@chromosomes = split(/,/,join(',',@chromosomes));


# Log Info
if ($detailed_debug) {
  # this uses a class variable in SequenceObj - not sure of Storable impact at mo' - this will still work.
  SequenceObj->debug(1);
}

my $log = Log_files->make_build_log($wormbase);

&check_opts;


# awful hack to override the location of GFF_SPLITS, needed in test mode
if (defined $gff_dir) {
  if (not $test) {
    $log->log_and_die("You can only override the GFF directory in test mode");
  } elsif (not -d $gff_dir) {
    $log->log_and_die("Given GFF directory does not exist");
  }
  $wormbase->{gff_splits} = $gff_dir;
}

# get coords obj to return clone and coords from chromosomal coords
print "Coords converter\n" if ($debug);
my $coords = Coords_converter->invoke($db, undef, $wormbase);

# Load in Feature_data : cDNA associations from COMMON_DATA
print "Load Feature_data\n" if ($debug);
my %feature_data;
&load_features( \%feature_data );
my %Features; # all TSL/Poly-A Features whether or not they are associated with a cDNA - # (182772,  182773, +, SL1) keyed by Feature_id

# process chromosome at a time
if (not @chromosomes) {
  if (defined $chunk_total and defined $chunk_id) {
    @chromosomes = $wormbase->get_chunked_chroms(-prefix => 1, 
                                                 -chunk_total => $chunk_total,
                                                 -chunk_id    => $chunk_id);
  } else {
    @chromosomes = $wormbase->get_chromosome_names('-prefix' => 1) unless @chromosomes;
  }
}

my $contigs = 1 if ($wormbase->assembly_type eq 'contig');

$ace_fname = sprintf("transcripts.%s.ace", $chromosomes[0]) if not defined $ace_fname;
$problem_fname = sprintf("problems.%s.txt", $chromosomes[0]) if not defined $problem_fname;
$transcript_dir = $wormbase->transcripts if not defined $transcript_dir;

my $out_file = sprintf("%s/%s", $transcript_dir, $ace_fname);
$log->write_to("Going to write output to $out_file\n") if ($verbose);
open (my $out_fh,">$out_file") or $log->log_and_die("cant open $out_file\n");

my $prob_file = sprintf("%s/%s", $transcript_dir, $problem_fname);
$log->write_to("Going to write problems to $prob_file\n") if ($verbose);
open (my $prob_fh,">$prob_file") or $log->log_and_die("cant open $prob_file\n");

print "Get Read_coverage\n" if ($debug);
my %Read_coverage;
#if (! $debug) { # this takes a long time to load, so ignore this stuff during debugging
  %Read_coverage = get_read_coverage();
#} else {
#  print "Not read in because debugging and it takes a long time\n";
#}

print "Get ignored EST data\n" if ($debug);
my %ignored_EST_data = &get_ignored_EST_data();

print "Get Feature evidence\n" if ($debug);
my %feature_evidence;;
#if (! $debug) { # this takes a long time to load, so ignore this stuff during debugging
  %feature_evidence = get_feature_evidence();
#} else {
#  print "Not read in because debugging and it takes a long time\n";
#}


foreach my $chrom ( @chromosomes ) {

  print "Chromosome $chrom\n" if ($debug);

  my ($link_start,$link_end);
  my %genes_exons;
  my %genes_span;
  my %cDNA;
  my %cDNA_span;
  my %cDNA_index;
  my %features;
  my $transformer;
  my @cdna_objs;
  my @cds_objs;
  my $index = 0;

  #
  # parse GFF file to get CDS and exon info
  #
  print "parse GFF file to get CDS and exon info\n" if ($debug);

  my $GFF = $wormbase->open_GFF_file($chrom, 'curated',$log);
  while (<$GFF>) {
    my @data = split;
    next if( $data[1] eq "history" );
    #  GENE STRUCTURE
    if ( $data[1] eq "curated" ) {
      $data[9] =~ s/\"//g;
      next if @test_cds and not grep { $data[9] eq $_ } @test_cds;
      if ( $data[2] eq "CDS" ) {
    	# GENE SPAN
    	$genes_span{$data[9]} = [($data[3], $data[4], $data[6])];
      } elsif ($data[2] eq "exon" ) {
    	# EXON 
    	$genes_exons{$data[9]}{$data[3]} = $data[4];
      }
    }
  }
  close $GFF;
  
  # read BLAT data
  print "read BLAT data\n" if ($debug);
  my @BLAT_methods = qw( BLAT_EST_BEST BLAT_mRNA_BEST BLAT_OST_BEST BLAT_RST_BEST BLAT_Trinity_BEST BLAT_IsoSeq_BEST BLAT_Nanopore_BEST);
  foreach my $method (@BLAT_methods) {
    # need to check that the GFF file for this method exists; not all species
    # have all methods
    my ($gff_file) = glob($wormbase->gff_splits . "/*${method}.gff");
    next if not defined $gff_file;

    my $GFF = $wormbase->open_GFF_file($chrom, $method, $log);
    while ( <$GFF> ) {
      next if (/^\#/); 		 # miss header
      next unless (/BEST/);

      my @data = split;

      $data[9] =~ s/\"//g;
      $data[9] =~ s/Sequence:// ;
      next if @test_est and not grep { $_ eq $data[9] } @test_est;
      
      next if ($ignored_EST_data{$data[9]}); # ignore anything with a Ignore tag set

      $cDNA{$data[9]}{$data[3]} = $data[4];

      # keep min max span of cDNA
      if ( !(defined($cDNA_span{$data[9]}[0])) or ($cDNA_span{$data[9]}[0] > $data[3]) ) {
    	$cDNA_span{$data[9]}[0] = $data[3];
    	$cDNA_span{$data[9]}[2] = $data[6]; #store strand of cDNA
      } 
      if ( !(defined($cDNA_span{$data[9]}[1])) or ($cDNA_span{$data[9]}[1] < $data[4]) ) {
    	$cDNA_span{$data[9]}[1] = $data[4];
      }
      
      $cDNA_span{$data[9]}[5] = $data[5]; # coverage score of the alignment
    }
    close $GFF;
  }
  
  # Chromomsome info
  print "read Chromosome data\n" if ($debug);

  $GFF = $wormbase->open_GFF_file($chrom, 'Link', $log);
  #create Strand_transformer for '-' strand coord reversal
  CHROM:while( <$GFF> ){
    my @data = split;
    my $chr = $wormbase->chromosome_prefix;
    if ( ($data[1] eq "Link" and $data[9] =~ /$chr/) or
    	 ($data[1] eq "Genomic_canonical" and $data[9] =~ /$chr/) ){
      $transformer = Strand_transformer->new($data[3],$data[4]);
      last CHROM;
    }
  }
  close $GFF;
  
  # add feature_data to cDNA (and store all Feature data whether or not it is associated with a cDNA)
  # CHROMOSOME_I  SL1  SL1_acceptor_site   182772  182773 .  -  .  Feature "WBsf016344"

  print "add feature_data to cDNA\n" if ($debug);
  my @feature_types = qw(SL1 SL2 polyA_site polyA_signal_sequence);
  foreach my $Type (@feature_types){
    # need to check that the GFF file for this method exists; not all species
    # have all methods
    my ($gff_file) = glob($wormbase->gff_splits . "/*${Type}.gff");
    next if not defined $gff_file;

    my $GFF = $wormbase->open_GFF_file($chrom, $Type, $log);
    while( <$GFF> ){
      my @data = split;
      if ( $data[9] and $data[9] =~ /(WBsf\d+)/) { # Feature "WBsf003597"
    	my $feat_id = $1;
	$Features{$feat_id} = [$data[3], $data[4], $data[6], $data[1]]; # 182772  182773 + SL1
    	my $dnas = $feature_data{$feat_id};
    	if ( $dnas ) {
    	  foreach my $dna ( @{$dnas} ) {
    	    # print "$dna\t$data[9]  --- $data[6] ---  ",$cDNA_span{"$dna"}[2],"\n";
    	    next unless ( $cDNA_span{"$dna"}[2] and $data[6] eq $cDNA_span{"$dna"}[2] ); # ensure same strand
    	    $cDNA_span{"$dna"}[3]{"$data[1]"} = [ $data[3], $data[4], $1 ]; # 182772  182773 WBsf01634 # add Feature to cDNA
    	  }
    	}
      }
    }
    close $GFF;
  }


  # need to sort the cds's into ordered arrays + and - strand genes are in distinct coord space so they need to be kept apart
  print "sort the cds's into ordered arrays\n" if ($debug);
  my %fwd_cds;
  my %rev_cds;
  foreach ( keys %genes_span ) {
    if ( $genes_span{$_}->[2] eq "+" ) {
      $fwd_cds{$_} = $genes_span{$_}->[0];
    } else {
      $rev_cds{$_} = $genes_span{$_}->[0];
    }
  }

  close $GFF;
  
  print "load EST data\n" if ($debug);
  &load_EST_data(\%cDNA_span, $chrom);  
  # &checkData(\$gff,\$%cDNA_span, \%genes_span); # this just checks that there is some BLAT and gene data in the GFF file
  &eradicateSingleBaseDiff(\%cDNA);

  # add SL Feature to $cds if it has an Isoform Feature asserted
  #create transcript obj for each CDS
  # fwd strand cds will be in block first then rev strand
  print "add SL Feature to cds\n" if ($debug);
  foreach (sort { $fwd_cds{$a} <=> $fwd_cds{$b} } keys  %fwd_cds ) {
    my $SL_ids = $feature_evidence{$_}; # add Feature_evidence to CDS
    my @SLs;
    if (defined $SL_ids) {
      foreach my $SL_id (@{$SL_ids}) {
	if ($Features{$SL_id}->[3] && $Features{$SL_id}->[3] =~ /^SL/) {
	  push @SLs, [$Features{$SL_id}->[0], $Features{$SL_id}->[1], $SL_id, $Features{$SL_id}->[3]]; # [start, end, WBsfID, type];
	}
      }
    } 
    my $cds = CDS->new( $_, $genes_exons{$_}, $genes_span{$_}->[2], $chrom, $transformer, \@SLs );
    push( @cds_objs, $cds);
    $cds->array_index("$index");
    $index++;
  }
  foreach ( sort { $rev_cds{$b} <=> $rev_cds{$a} } keys  %rev_cds ) {
    my $SL_ids = $feature_evidence{$_}; # add Feature_evidence to CDS
    my @SLs;
    if (defined $SL_ids) {
      foreach my $SL_id (@{$SL_ids}) {
	if ($Features{$SL_id}->[3] && $Features{$SL_id}->[3] =~ /^SL/) {
	  push @SLs, [$Features{$SL_id}->[0], $Features{$SL_id}->[1], $SL_id, $Features{$SL_id}->[3]]; # [start, end, WBsfID, type];
	}
      }
    } 
    my $cds = CDS->new( $_, $genes_exons{$_}, $genes_span{$_}->[2], $chrom, $transformer, \@SLs );
    push( @cds_objs, $cds);
    $cds->array_index("$index");
    $index++;
  }


  print "Make EST objects\n" if ($debug);
  my $count0 = 0;
  foreach my $cdna_id ( keys %cDNA ) {

    my $cdna = SequenceObj->new( $cdna_id, $cDNA{$cdna_id}, $cDNA_span{$cdna_id}->[2] );

    if ( $cDNA_span{$cdna_id}->[3] ) {
      foreach my $feat ( keys %{$cDNA_span{$cdna_id}->[3]} ) {
	$cdna->$feat( $cDNA_span{$cdna_id}->[3]->{"$feat"} ); # add feature to cDNA 
      }
    }
    # add paired read info
    if ( $cDNA_span{$cdna_id}->[4] ) {
      $cdna->paired_read( $cDNA_span{$cdna_id}->[4] );
    }

    # add coverage score
    if ( $cDNA_span{$cdna_id}->[5] ) {
      $cdna->coverage( $cDNA_span{$cdna_id}->[5] );
    }

    $cdna->transform_strand($transformer,"transform") if ( $cdna->strand eq "-" );

    #check for and remove ESTs with internal SL's 
    if (&sanity_check_features( $cdna )) {
      push @cdna_objs, $cdna;
    } else {
      $count0++;
    }
  }

  ######
  # sort the cDNAs such that the ones with features are dealt with first
  # This is necessary to ensure deterministic behaviour
  ######
  print "Sort the ESTs\n" if ($debug);
  @cdna_objs = sort {
    my $a_score = 0;
    $a_score += 2 if defined $a->SL;
    $a_score += 1 if defined $a->polyA_site;

    my $b_score = 0;
    $b_score += 2 if defined $b->SL;
    $b_score += 1 if defined $b->polyA_site;

    return $b_score - $a_score;
    
  } @cdna_objs;
  

  ######
  # Index for later rapid retrieval
  ######

  print "Index EST objects\n" if ($debug);
  for(my $i=0; $i < @cdna_objs; $i++) {
    my $cdna_obj = $cdna_objs[$i];
 
    my $cdna_id = $cdna_obj->name;

    $cdna_obj->array_index($i);
    $cDNA_index{$cdna_id} = $i;

  }


  # these are no longer needed so free memory !
  %genes_exons = ();
  %genes_span= ();
  %cDNA = ();
  %cDNA_span = ();

  ##########################################################
  # DATA LOADED - START EXTENDING THE TRANSCRIPTS          #
  ##########################################################

  # First round
  #
  # remove any cDNA that overlaps two CDSs and which has a score of
  # less than $COVERAGE_THRESHOLD for the alignment coverage score
  #
  # At present this does a very simple check to see if the EST
  # overlaps with two or more genes, but it should be improved to
  # check whether the EST exons overlap with the CDS exons because at
  # the moment genes in the introns of other genes in the same sense
  # have their weak EST rejected.

  my $round = 'First (poor quality) round:';
  my $count1 = 0;

  $log->write_to("\nTB : $round\n") if ($verbose);
  print "\nTB: first round\n\n";

  # now get the cDNAs that match introns of CDS structures this allows
  # us to change the program to iterate over CDSs first, instead of
  # iterating over all cDNAs and running out of memory. We can then
  # just add the matched cDNAs to the CDS structures and process one
  # CDS at a time instead of holding all CDS matches to all cDNAs in
  # memory at once (> 50 Gb!)
  print "\nGet CDS introns\n";

  my %CDS_introns; # key = {start}{end} value = [CDS_name]
  foreach my $cds ( @cds_objs ) {
    foreach my $intron_start (keys %{$cds->{introns}}) {
      my $intron_end = $cds->{introns}{$intron_start};
      push @{$CDS_introns{$intron_start}{$intron_end}}, $cds->name;
    }
  }

  print "\nGet matches of cDNA introns to CDS introns\n";

  my %cDNA_CDS; # key = {cDNA_name}{CDS_name} value = no. of times match is observed
  my %CDS_cDNA; # key = {CDS_name}{cDNA_name} value = no. of times match is observed
  foreach my $cdna (@cdna_objs) {
    foreach my $intron_start (keys %{$cdna->{introns}}) {
      if (exists $CDS_introns{$intron_start}) {
	my $intron_end = $cdna->{introns}{$intron_start};
	if (exists $CDS_introns{$intron_start}{$intron_end}) {
	  foreach my $CDS_name (@{$CDS_introns{$intron_start}{$intron_end}}) {
	    $cDNA_CDS{$cdna->name}{$CDS_name}++;
	    $CDS_cDNA{$CDS_name}{$cdna->name}++;
	  }
	}
      }
    }
  }

  print "Ignore any cDNA that matches the introns of two or more genes";
  foreach my $cdna ( @cdna_objs) {
    next if ( @test_est and not grep { $cdna->name eq $_ } @test_est); #debug line

    if (exists $cDNA_CDS{$cdna->name}) {
      my @matches = (keys %{$cDNA_CDS{$cdna->name}});
      
      my @matching_genes = ($species eq 'elegans') 
	? list_of_matched_genes_by_seqname($wormbase->seq_name_regex, \@matches)
	  : list_of_matched_genes(\%cds2gene, \@matches);
      
      if (scalar(@matching_genes) > 1 and $cdna->coverage < $COVERAGE_THRESHOLD) { 
	if ($verbose) {
	  $log->write_to("TB: $round : cDNA " . 
			 $cdna->name . 
			 " overlaps two or more genes and has an alignment score of less than $COVERAGE_THRESHOLD so it will not be used in transcript-building:");
	  map { $log->write_to("TB : $round :\t$_") } @matching_genes;
	  $log->write_to("\n");
	}
	
	print($prob_fh "$round cDNA " .
	      $cdna->name .
	      " overlaps two or more genes and has an alignment score of less than $COVERAGE_THRESHOLD so it will not be used in transcript-building: @matching_genes\n");
	$count1++;
	$cdna->mapped(1);  # mark the cDNA as used with a dummy CDS reference
      }
    }
  }



  ##########################################################
  # Second round.
  #
  # here we go through all of the cDNAs looking for matches of their
  # introns with the introns of the CDS structures. We store all
  # instances of a match in the cDNA structure so that we can go
  # through them all later and check first of all whether there are
  # any cDNAs with introns matching two or more genes - these are
  # candidates for merging genes or chimeric ESTs and the cDNA is not
  # used to build the transcripts as they will almost certainly
  # produce erroneous structures. Where the cDNA matches just one
  # gene, the cDNA may match several CDSs equally well. All of the
  # CDSs that match with an equal number of consecutive introns are
  # then given the option of adding the cDNA to their transcripts.

  # changed program to do one CDS at a time, to save memory
  # sort the CDSs by their names 
  my $prev_gene_name = "";
  my @cds_objs_in_gene = ();
  foreach my $next_cds (sort {$a->name cmp $b->name} @cds_objs) { 
    
    $round = "Second (intron) round (".$next_cds->name."):";
    
    $log->write_to("\nTB : $round\n") if ($verbose);
    print "\nTB: $round round\n\n";
    
    # want to check if the cDNA has introns that match one and only one gene
    # for each cDNA that has matching introns with the CDS
    foreach my $cdna_name (keys %{$CDS_cDNA{$next_cds->name}}) {
      my $cdna = $cdna_objs[$cDNA_index{$cdna_name}];
      
      next if ( defined($cdna->mapped) );
      next if ( @test_est and not grep { $cdna->name eq $_ } @test_est); #debug line
      
      # want to see which fresh set of CDSs this matches
      $cdna->reset_probably_matching_cds;
      
      if ($next_cds->map_introns_cDNA($cdna) ) { 
	# note each CDS and gene that this cDNA matches, together with the number of contiguous CDS introns matched
	$log->write_to("TB : $round : Registered intron match between " . $next_cds->name . " and " . $cdna->name ."\n") if ($verbose);
      } else {
	$log->write_to("TB : $round : No intron match between " . $next_cds->name . " and " . $cdna->name . "\n") if $verbose;
      }
      
      # now that this cDNA has information on which CDS's introns it
      # matches, we can add any cDNA that matches just one gene at the
      # introns level to the transcripts.  There may be more than one
      # CDS in a gene that matches the same number of CDS introns with
      # no mismatches
      
      my @matching_genes = ($species eq 'elegans') 
        ? $cdna->list_of_matched_genes_by_seqname($wormbase->seq_name_regex)
	  : $cdna->list_of_matched_genes(\%cds2gene);
      
      if (scalar(@matching_genes) == 1) { 
	my @best_cds = &get_best_CDS_matches($cdna); # get those CDSs for this cDNA that have the most introns matching
	foreach my $cds (@best_cds) {
	  if ( $cds->map_cDNA($cdna) ) { # add it to the transcript structure
	    $cdna->mapped($cds);  # mark the cDNA as used
	    $log->write_to("TB : $round : Used intron match of " . $cdna->name . " in a transcript of " . $cds->name ."\n") if ($verbose);
	  }
	}
      } elsif (scalar(@matching_genes) > 1) { 
	# we want to report those ESTs that have introns that match two or more CDSs as this may indicate required gene mergers.
	if ($verbose) {
	  $log->write_to("TB : $round : cDNA " .
			 $cdna->name . 
			 " matches introns in two or more genes and will not be used in transcript-building:");
	  foreach my $gene (@matching_genes) {
	    $log->write_to("TB : $round :\t" . $gene);
	  }
	  $log->write_to("\n");
	}
	print($prob_fh "$round cDNA " .
	      $cdna->name .
	      " matches introns in two or more genes and will not be used in transcript-building: @matching_genes\n");
	$cdna->mapped(1);  # mark the cDNA as used with a dummy CDS reference
      }
    }
    
    # Here we use the cDNAs that were not used in the intron round
    
    ##########################################################
    $round = "Third (extend transcripts) round:";
    $log->write_to("\nTB : $round\n") if ($verbose);
    print "\nTB: $round round\n\n";
    
    
    foreach my $cdna_name (keys %{$CDS_cDNA{$next_cds->name}}) {
      my $cdna = $cdna_objs[$cDNA_index{$cdna_name}];
      
      next if ( @test_est and not grep { $cdna->name eq $_ } @test_est); #debug line
      next if ( defined($cdna->mapped) );
      
      # want to see which fresh set of CDSs this matches
      $cdna->reset_probably_matching_cds;
      
      # here we are now looking for overlapped transcripts and want to
      # avoid using cDNAs that overlap two genes with no intron evidence
      # as to which one it should be added to.
      
      if ($cdna->overlap($next_cds)) {
	$cdna->probably_matching_cds($next_cds, 1);
      }
      
      # now check how many genes the cDNA overlaps - we only want those that overlap one
      my @matching_genes = ($species eq 'elegans') 
        ? $cdna->list_of_matched_genes_by_seqname($wormbase->seq_name_regex)
	  : $cdna->list_of_matched_genes(\%cds2gene);
      
      if (scalar(@matching_genes) == 1) { # just one matching gene
	foreach my $cds_match (@{$cdna->probably_matching_cds}) {
	  my $cds = $cds_match->[0];
	  if ( $cds->map_cDNA($cdna) ) { # add it to the transcript structure
	    $cdna->mapped($cds );  # mark the cDNA as used
	    $log->write_to("TB : $round : Have used first round addition of " . $cdna->name . " in the transcript " . $cds->name ."\n") if ($verbose);
	  }
	}
      } elsif (scalar(@matching_genes) > 1) { 
	$log->write_to("TB : $round : cDNA ",$cdna->name," overlaps two or more genes and will not be used in transcript-building:") if ($verbose);
	foreach my $gene (@matching_genes) {
          $log->write_to("TB : $round :\t" . $gene) if ($verbose);
	}
	$log->write_to("\n") if ($verbose);
	print($prob_fh "$round cDNA " .
	      $cdna->name .
	      " overlaps two or more genes and will not be used in transcript-building: @matching_genes\n");
	$cdna->mapped(1);  # mark the cDNA as used with a dummy CDS reference
      }
    }
    
    
    ##########################################################
    # Fourth round - use read-pair information to extend with cDNAs that do not overlap
    
    $round = "Fourth (read pairs) round:";
    $log->write_to("TB : $round\n") if ($verbose);
    print "\nTB: $round round\n\n";
    
  PAIR:   foreach my $cdna_name (keys %{$CDS_cDNA{$next_cds->name}}) {
      my $cdna = $cdna_objs[$cDNA_index{$cdna_name}];
      
      next if ( @test_est and not grep { $cdna->name eq $_ } @test_est); #debug line
      next if $cdna->mapped;
      
      # get name of paired read
      next unless (my $mapped_pair_name = $cdna->paired_read );
      
      # retrieve object from array
      next unless ($cDNA_index{$mapped_pair_name} and (my $partner = $cdna_objs[ $cDNA_index{$mapped_pair_name} ] ) );
      
      # get cds that paired read maps to 
      if (my $cds = $partner->mapped) {
	if ($cds != 1) { # don't want to start using the dummy 'cds = 1' flags we used earlier to mark cDNAs we don't wish to use in a transcript
	  my $index = $cds->array_index;
	  
	  # find next downstream CDS - must be on same strand
	  my $downstream_CDS;
	DOWN: while (! defined $downstream_CDS ) {
	    $index++;
	    if ( $downstream_CDS = $cds_objs[ $index ] ) {
	      
	      unless ( $downstream_CDS ) {
		last;
		$log->write_to("TB : $round : last gene in array\n") if ($verbose);
	      }
	      # dont count isoforms
	      my $down_name = $downstream_CDS->name;
	      my $name = $cds->name;
	      $name =~ s/[a-z]//;
	      $down_name =~ s/[a-z]//;
	      if ( $name eq $down_name ) {
		undef $downstream_CDS;
		next;
	      }
	      # @cds_objs is structured so that + strand genes are in a block at start, then - strand
	      last DOWN if( $downstream_CDS->strand ne $cds->strand );
	      
	      # check unmapped cdna ( $cdna ) lies within 1kb of CDS that paired read maps to ( $cds ) and before $downstream_CDS
	      
	      #print "$round trying ",$cds->name, " downstream is ", $downstream_CDS->name," with ",$cdna->name,"\n" if ($verbose);
	      if ( ($cdna->start > $cds->gene_end) and ($cdna->start - $cds->gene_end < 1000) and ($cdna->end < $downstream_CDS->gene_start) ) {
		$log->write_to("TB : $round : adding 3' cDNA " . $cdna->name . " to " . $cds->name . "\n") if ($verbose);
		$cds->add_3_UTR($cdna);
		$cdna->mapped($cds);
		last;
	      }
	    } else {
	      last DOWN;
	    }
	  }
	}
      }
    }
    
    ##########################################################=
    #  $round = "Fifth (extend transcripts again) round:";
    #  $log->write_to("$round\n") if ($verbose);
    #
    #  foreach my $CDNA ( @cdna_objs) {
    #    next if ( defined $est and $CDNA->name ne "$est"); #debug line
    #    next if ( defined($CDNA->mapped) );
    #    #sanity check features on cDNA ie SLs are at start
    #    next if ( &sanity_check_features( $CDNA ) == 0 );
    #    foreach my $cds ( @cds_objs ) {
    #      if ( $cds->map_cDNA($CDNA) ) {
    #	$CDNA->mapped($cds);
    #	print "$round ",$CDNA->name," overlaps ",$cds->name,"\n" if ($verbose);
    #      }
    #    }
    #  }
    
    
    my $seq_name_regexp = $wormbase->seq_name_regex;
    my ($gene_name) = ($next_cds->name =~ /($seq_name_regexp)/);
    if ($prev_gene_name eq $gene_name) {
      push @cds_objs_in_gene, $next_cds;
    } else {
      
      # find transcript isoforms which have weak evidence compared to their siblings and remove them
      purge_weakly_supported_transcript_structures(\@cds_objs_in_gene, \%Read_coverage);
      
      # find duplicate transcripts in a gene and make then non-identical so that the ENA will accept them
      purge_duplicates(\@cds_objs_in_gene, \%cds2gene);
      
      # output
      foreach my $cds (@cds_objs_in_gene) {
	$cds->report($out_fh, $coords, $wormbase->full_name, \%cds2gene);
      }

      # delete objects
      foreach my $cds (@cds_objs_in_gene) {
	$cds->delete;
      }
      
      @cds_objs_in_gene = ($next_cds);
      $prev_gene_name = $gene_name;
    }
    
    
  } # foreach $cds (@cds_objs)
  
  # find transcript isoforms which have weak evidence compared to their siblings and remove them
  purge_weakly_supported_transcript_structures(\@cds_objs_in_gene, \%Read_coverage);
  
  # process last gene
  # find duplicate transcripts in a gene and make then non-identical so that the ENA will accept them
  purge_duplicates(\@cds_objs_in_gene, \%cds2gene);
  
  # output
  foreach my $cds (@cds_objs_in_gene) {
    $cds->report($out_fh, $coords, $wormbase->full_name, \%cds2gene);
  }
  
  $log->write_to("$count0 cDNAs rejected in round 0 (inconsistent attached features)\n");
  $log->write_to("$count1 cDNAs rejected in round 1 (low quality and overlaps two or more genes)\n");
  #  $log->write_to("$count2 cDNAs rejected in round 2 (introns matched in two or more genes)\n");
  #  $log->write_to("$count3 cDNAs rejected in round 3 (overlaps two or more genes)\n");

} # foreach (@chromosomes)

print $out_fh "\n\n// Finished.\n";

close($prob_fh);
close($out_fh);

# cleanup things
$coords=undef;
$ace->close;

$log->mail();
exit(0);


######################################################################################################
#
#
#                           T  H  E        S  U  B  R  O  U  T  I  N  E  S
#
#
#######################################################################################################


# find duplicate transcripts within a gene locus and make them non-identical so that they are accepted by the ENA
sub purge_duplicates {
  my ($cds_objs, $cds2gene) = @_;

  my %genes; # hash of array of CDS objects, keyed by geneID

  foreach my $cds (@{$cds_objs}) {
    my $gene = $cds2gene->{$cds->name};
    push @{$genes{$gene}}, $cds; # add the current CDS to the CDSs to check against, because we had an instance of two Transcripts in the same isoform with the same structure! (T24B8.3c.1 and T24B8.3c.2)
    my $reverted = 0; # flag for already reverted a transcript in this CDS to the CDS structure
    # go through %genes checking every cds and every transcript for duplicates of the transcripts in $cds
    foreach my $transcript ($cds->transcripts) {
      foreach my $sibling_cds (@{$genes{$gene}}) { # look at other CDSs in this gene that we have already processed
	foreach my $sibling_transcript ($sibling_cds->transcripts) {
	  if ($sibling_transcript->name eq $transcript->name) {next} # don't check against itself
	  if (exists $sibling_transcript->{'ignore'}) {next}
	  if ($transcript->duplicate($sibling_transcript)) {
	    print $prob_fh "\nDuplicate Transcript: ",$transcript->name," ",$sibling_transcript->name,"\n";
	    $log->write_to("\nDuplicate Transcript: ".$transcript->name," ".$sibling_transcript->name."\n");
	    if ($reverted) { # don't want to revert to the CDS structure if we have already done this for another transcript in this CDS
	      print $prob_fh "\nWARNING: Two transcripts were reverted to the same CDS structure: ",$reverted," NOT CREATING: ",$transcript->name,"\n";
	      $log->write_to("\nWARNING: Two transcripts were reverted to the same CDS structure: ".$reverted." NOT CREATING: ".$transcript->name."\n");
#	      $log->error;
	      $transcript->{'ignore'} = 1; # ignore, don't report this.
	      last;
	    } else {
	      # force the structure back to being the original CDS structure
	      # set start, end, sorted_exons
	      $transcript->start($cds->start);
	      $transcript->end($cds->end);
	      $transcript->{'sorted_exons'} = $cds->sorted_exons;
	      $reverted = $transcript->name;
	      # now retest the new structure for duplicates
	      foreach my $sibling_cds2 (@{$genes{$gene}}) { # look at other CDSs in this gene that we have already processed
		foreach my $sibling_transcript2 ($sibling_cds2->transcripts) {
		  if (exists $sibling_transcript2->{'ignore'}) {next}
	          if ($sibling_transcript2->name eq $transcript->name) {next} # don't check against itself
		  if ($transcript->duplicate($sibling_transcript2)) {
		    print $prob_fh "\nWARNING: Duplicate Transcript found after reversion to CDS structure: ",$sibling_transcript2->name," NOT CREATING: ",$transcript->name,"\n";
		    $log->write_to("\nWARNING: Duplicate Transcript found after reversion to CDS structure: ".$sibling_transcript2->name." NOT CREATING: ".$transcript->name."\n");
#		    $log->error;
		    $transcript->{'ignore'} = 1; # ignore, don't report this.
		    last;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

}


# find transcript isoforms which have weak evidence compared to their siblings and remove them
# 'lowest_exon_count' is the lowest count (weakest evidence for an exon) of supporting cDNAs of all exons in the transcript
# as calculated by CDS::increment_exon_counts()
sub  purge_weakly_supported_transcript_structures {
  my ($cds_objs, $Read_coverage) = @_;

  foreach my $cds (@{$cds_objs}) {
    my $best_weak_count = 0;

    print "In purge_weakly_supported_transcript_structures\n" if ($debug);
    print "\n\tCDS: ".$cds->name."\n" if ($debug);
    # get the strongest evidence count of all the lowest evidence counts of all transcripts
    foreach my $transcript ($cds->transcripts) {
      print "\n\tTranscript: ".$transcript->name."\n" if ($debug);
      if (exists $transcript->{'ignore'}) {
	print "\tTranscript ".$transcript->name." has an Ignore set already\n" if ($debug);
	next}
      if (!exists $transcript->{'evidence'}) {
	print "\tTranscript ".$transcript->name." has no exon evidence values\n" if ($debug);
	next} # this happens if the structure has been reverted back to the CDS structure
      
      foreach my $cdna (@{$transcript->{'matching_cdna'}}) {  # debug info - not needed for the code to work
	my $name = $cdna->name;
	print "\tDebug: adding $name to exon evidence counts\n";
      }

      my $lowest_count = 1000000;
      foreach my $exon (@{$transcript->sorted_exons}) {
	my $exon_count = 0;
	foreach my $exon_start (keys %{$transcript->{'evidence'}}) {
	  my $cdna_name = $transcript->{'evidence'}{$exon_start};
	  if (exists $Read_coverage->{$cdna_name}) {
	    $exon_count += $Read_coverage->{$cdna_name}; # add the Nanopore Read_coverage to the count of evidence supporting this exon
	  } else {
	    $exon_count++; # just increment the count of EST/mRNA/trinity transcripts supporting this exon
	  }
	}
	if ($exon_count < $lowest_count) {$lowest_count = $exon_count}
	print " and now lowest_count = $lowest_count\n";
      }
      $transcript->{'lowest_exon_count'} = $lowest_count;
      if ($transcript->{'lowest_exon_count'} > $best_weak_count) {$best_weak_count = $transcript->{'lowest_exon_count'}}
    }
    print "\tbest_weak_count=$best_weak_count\n" if ($debug);
    my $best_weak_count_30_percent = $best_weak_count / 3;
    my $count = 0;
    foreach my $transcript (sort {$b->{'lowest_exon_count'} <=> $a->{'lowest_exon_count'}} $cds->transcripts) { # sort by reverse lowest_exon_count
      if (exists $transcript->{'ignore'}) {next}
      if (!exists $transcript->{'lowest_exon_count'}) {next} # this happens if the structure has been reverted back to the CDS structure
# C15H9.6a has lots of isoforms from the introns in the 3' UTR, but also has single good isoforms from 5' whispy bits we want to keep
# try removing the ' $transcript->{'lowest_exon_count'} == 1 && $best_weak_count > 1' to see what effect this has
# what happens in the cases of?: daf-16, C24A3.2a, C25H3.9a, C34F11.6, Y65B4A.6a/b, K08E7.2, F21H12.5
      print "\tChecking transcript ".$transcript->name."\n" if ($debug);
      print "\tTranscript lowest_exon_count = $transcript->{'lowest_exon_count'}\n" if ($debug);
      if ($transcript->{'lowest_exon_count'} < $best_weak_count_30_percent) {
	print $prob_fh "\nWeak evidence for Transcript: ",$transcript->name," ",$transcript->{'lowest_exon_count'}," / $best_weak_count\n";
	$log->write_to("\nWeak evidence for Transcript: ".$transcript->name." ".$transcript->{'lowest_exon_count'}." / $best_weak_count\n");
	$transcript->{'ignore'} = 1; # ignore, don't report this.
	print "\t*** Ignore this transcript because its evidence is too weak ".$transcript->name." (< $best_weak_count / 3) *********\n" if ($debug);
      } elsif ($transcript->{'lowest_exon_count'} == 1 && $best_weak_count > 1) {
	print $prob_fh "\nSingle EST evidence for Transcript: ",$transcript->name," ",$transcript->{'lowest_exon_count'}," / $best_weak_count\n";
	$log->write_to("\nSingle EST evidence for Transcript: ".$transcript->name." ".$transcript->{'lowest_exon_count'}." / $best_weak_count\n");
	$transcript->{'ignore'} = 1; # ignore, don't report this.
	print "\t*** Single EST evidence for Transcript: ".$transcript->name." *********\n" if ($debug);
      } else {
	$count++; # only want the two best ones
	if ($count > 2) {
	  print $prob_fh "\nOnly want the two best Transcripts, so ignore: ",$transcript->name,"\n";
	  $log->write_to("\nOnly want the two best Transcripts, so ignore: ".$transcript->name."\n");
	  $transcript->{'ignore'} = 1; # ignore, don't report this.
	  print "\t*** Only want the two best Transcripts, so ignore: ".$transcript->name."\n" if ($debug);
	}
      }
    }
  }
}




sub eradicateSingleBaseDiff {
  my $cDNAh = shift;
  $log->write_to( "\nMerging bitty alignment fragments caused by single bp mismatches\n\n") if ($verbose);
  foreach my $cdna_id (keys %{$cDNAh} ) {
    my $last_key;
    my $altered = 0;

    $log->write_to("\n############### $cdna_id #############\n") 
        if $really_verbose;
    foreach my $exon_start (sort { $cDNAh->{$cdna_id}->{$a} <=> $cDNAh->{$cdna_id}->{$b} } keys %{$cDNAh->{$cdna_id}}) {
      $log->write_to(" $exon_start -> $cDNAh->{$cdna_id}->{$exon_start}\n") 
          if $really_verbose;

      my $new_last_key = $exon_start;
      if ( $last_key ) {
        if ( $cDNAh->{$cdna_id}->{$last_key} >= $exon_start - $gap ) { #allows seq error gaps up to $gap bp
          $cDNAh->{$cdna_id}->{$last_key} = $cDNAh->{$cdna_id}->{$exon_start};
          delete $cDNAh->{$cdna_id}->{$exon_start};
          $new_last_key = $last_key;
          $altered = 1;
        }
      }
      $last_key = $new_last_key;
    }
    if ( $altered ) {
      $log->write_to(" Single base diffs removed from $cdna_id:\n") if $really_verbose;
      foreach my $exon_start (sort keys  %{$cDNAh->{$cdna_id}}) {
        $log->write_to(" $exon_start -> $cDNAh->{$cdna_id}->{$exon_start}\n") if $really_verbose;
      }
    }
  }
}

sub check_opts {
  # sanity check options
  if ( $help ) {
    system("perldoc $0");
    exit(0);
  }
}

sub checkData
  {
    my $file = shift;
    my $cDNA_span = shift;
    my $genes_span = shift;
    die "There's no BLAT data in the gff file $$file\n" if scalar keys %{$cDNA_span} == 0;
    die "There are no genes in the gff file $$file\n" if scalar keys %{$genes_span} == 0;
  }



###################################################################################

sub load_EST_data {
  my $cDNA_span = shift;
  my $chrom = shift;
  my %est_orient;

  $wormbase->FetchData("estorientation",\%est_orient) unless (5 < scalar keys %est_orient);
  
  foreach my $EST ( keys %est_orient ) {
    if ( exists $$cDNA_span{$EST} && defined $$cDNA_span{$EST}->[2]) {
      my $GFF_strand = $$cDNA_span{$EST}->[2];
      my $read_dir = $est_orient{$EST};
      CASE:{
        ($GFF_strand eq "+" and $read_dir eq "5") && do {
          $cDNA_span->{$EST}->[2] = "+";
          last CASE;
        };
        ($GFF_strand eq "+" and $read_dir eq "3") && do {
          $cDNA_span->{$EST}->[2] = "-";
          last CASE;
        };
        ($GFF_strand eq "-" and $read_dir eq "5") && do {
          $cDNA_span->{$EST}->[2] = "-";
          last CASE;
        };
        ($GFF_strand eq "-" and $read_dir eq "3") && do {
          $cDNA_span->{$EST}->[2] = "+";
          last CASE;
        };
      }
    }
  }

  # load paired read info
  $log->write_to("Loading EST paired read info\n") if ($verbose);
  my $pairs = $wormbase->common_data."/EST_pairs.txt";
  
  if ( -e $pairs ) {
    open ( PAIRS, "<$pairs") or $log->log_and_die("cant open $pairs :\t$!\n");
    while ( <PAIRS> ) {
      chomp;
      s/\"//g;#"
          s/Sequence://g;
      next if( ( $_ =~ /acedb/) or ($_ =~ /\/\//) );
      my @data = split;
      $cDNA_span->{$data[0]}->[4] = $data[1];
    }
    close PAIRS;
  } else {
    my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
    
    open (TACE, "echo '$cmd' | $tace $db |") or die "cant open tace to $db using $tace\n";
    open ( PAIRS, ">$pairs") or die "cant open $pairs :\t$!\n";
    while ( <TACE> ) { 
      chomp;
      s/\"//g;#"
          my @data = split;
      next unless ($data[0] && $data[1]);
      next if $data[0]=~/acedb/;
      $$cDNA_span{$data[0]}->[4] = $data[1];
      print PAIRS "$data[0]\t$data[1]\n";
    }
    close PAIRS;
  }
}


sub get_ignored_EST_data {

  my %ignored_EST_data;
  my $cmd = "select cdna from cdna in class cDNA_sequence where exists_tag cdna->Ignore";

  open (TACE, "echo '$cmd' | $tace $db |") or die "cant open tace to $db using $tace\n";
  while ( <TACE> ) { 
    chomp;
    s/\"//g;#"
    my @data = split;
    next unless ($data[0]);
    next if $data[0]=~/acedb/;
    $ignored_EST_data{$data[0]} = 1;
  }
  return %ignored_EST_data;
}

sub load_features {
  my $features = shift;
  my %tmp;
  $wormbase->FetchData("est2feature",\%tmp);
  foreach my $seq ( keys %tmp ) {
    my @feature = @{$tmp{$seq}};
    foreach my $feat ( @feature ) {
      push(@{$$features{$feat}},$seq);
    }
  }
}

sub sanity_check_features {
  my ($cdna) = @_;
  
  if ( my $sl = $cdna->SL ) {
    if( abs($sl->[0] - $cdna->start) != 1 ) {
      $log->write_to("${\$cdna->name} failed sanity check because attached SL feature (" . $sl->[2] . ") not at expected position " . $sl->[0] . ":".$sl->[1] ." " . $cdna->start . ":" . $cdna->end . "\n") if $verbose;
      return 0;
    }
  }
  
  if ( my $polyA = $cdna->polyA_site ) {
    if($polyA->[0] != $cdna->end){
      $log->write_to("${\$cdna->name} failed sanity check because attached polyA feature (" . $polyA->[2] . ") not at expected position \n") if $verbose;
      return 0;
    }
  }
  
  return 1;
}

# get the list of CDS objects from the $cdna->probably_matching_cds
# and return those with the highest number of matching introns
sub get_best_CDS_matches {
  my ($cdna) = @_;

  my @cds_matches = @{$cdna->probably_matching_cds};
  @cds_matches = sort {$b->[1] <=> $a->[1]} @cds_matches; # reverse sort by number of matching introns
  my $max_introns = $cds_matches[0][1]; # get the highest number of matching introns
  my @result;
  foreach my $next_cds (@cds_matches) {
    if ($next_cds->[1] == $max_introns) {
      push @result, $next_cds->[0];
    }
  }
  return @result;
}

# get CDS where the Isoform has one or more Properties.Isoform.Feature_evidence supporting evidence tags
sub get_feature_evidence {
  my %Feature_evidence;

#my @cds = $ace->fetch( -query => "Find CDS Where Isoform AND Method = \"curated\"");

  my $cdsIt = $ace->fetch_many(-query => "Find CDS Where Isoform AND Method = \"curated\"");
  
  while (my $cds = $cdsIt->next()) {
    # The 'Feature_evidence' values are in an Evidence hash attached to the 'Isoform' tag. You have to creep up on it slowly.
    if (defined $cds->at('Properties.Isoform.Feature_evidence')) {
      foreach my $ft ($cds->at('Properties.Isoform.Feature_evidence')->col(1)) {
	push @{$Feature_evidence{$cds->name}}, $ft->name;
      }
    }
  }
  return %Feature_evidence;
}


# where the Nanopore transcripts are equal in sequence, they have been collapsed to one Sequence with a Read_coverage holding the count
# get the Read_coverage value from each of the Nanopore transcripts
sub get_read_coverage {
  my %Read_coverage;
  my $NanoIt = $ace->fetch_many(-query => "Find Sequence Where Method = \"Nanopore\" AND Read_coverage > 1");
  
  while (my $nano = $NanoIt->next()) {
    # The 'Read_coverage' values should all be populated, but check just in case
    if (defined $nano->at('Origin.Read_coverage')) {
      my $coverage = $nano->Read_coverage->name;
      if ($coverage > 1) {
	$Read_coverage{$nano->name} = $coverage;
      }
    }
  }
  return %Read_coverage;
}



=head2 list_of_matched_genes

    Title   :   list_of_matched_genes
    Usage   :   intron_matched_genes(\%cds2gene, \@matching_cds_names)
    Function:   returns the list of genes derived from the names of the CDSs in $cdna->probably_matching_cds
    Returns :   list of gene names
    Args    :   
               

=cut

sub list_of_matched_genes {
  my ($cds2gene, $matches) = @_;

  my %genes;
  #print "Checking probably_matching_cds for ",$self->name,"\n";
  if (! @{$matches}) {return ()}
  foreach my $cds_name (@{$matches} ) {

    if (exists $cds2gene->{$cds_name}) {
      $genes{$cds2gene->{$cds_name}} = 1;
    }
  }

  return keys %genes;
}



sub list_of_matched_genes_by_seqname {
  my ($seq_name_regexp, $matches) = @_;

  my %genes;

  return () if not @{$matches};

  foreach my $cds_name (@{$matches}) {
    my ($gene) = ($cds_name =~ /($seq_name_regexp)/);
  
    if (defined $gene) {
      $genes{$gene} = 1;
    }
  }

  return keys %genes;
}




__END__

=pod

=head2 NAME - transcript_builder.pl

=head1 USAGE

=over 4

=item transcript_builder.pl  [-options]

=back

This script "does exactly what it says on the tin". ie it builds transcript objects for each gene in the database

To do this it ;

1) Determines matching_cDNA status for each gene. Goes through gff files and examines each cDNA to see if it matches any gene that it overlaps.

2) For each gene that has matching cDNAs it then confirms that every exon of the gene that lies within the region covered by the cDNA is covered correctly.  So if a gene has an extra exon that is in the intron of a cDNA, that cDNA will NOT be linked to that gene.


=back

=head2 transcript_builder.pl arguments:

=over 4

=item * verbose and really_verbose  -  levels of terminal output
  
=item * est:s     - just do for single est 

=item * gap:s      - when building up cDNA exon structures from gff file there are often single / multiple base pair gaps in the alignment. This sets the gap size that is allowable [ defaults to 5 ]

=item * gff_dir:s  - pass in the location of chromosome_*.gff files that have been generated for the database you are generating Transcripts for.

=item * gff:s         - pass in a gff file to use

=item * database:s      - either use autoace if used in build process or give the full database path. Basically retrieves paired read info for ESTs from that database.

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

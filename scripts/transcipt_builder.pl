#!/usr/local/bin/perl5.6.1 -w

use strict;
use lib glob("~ar2/wormbase/scripts/");
use Ace;
use Getopt::Long;

my ($debug, $verbose, $really_verbose, $est, $count, $report, $gap, $transcript);
$gap = 1; # $gap is the gap allowed in an EST alignment before it is considered a "real" intron

GetOptions ( "debug" => \$debug,
	     "verbose" => \$verbose,
	     "really_verbose" => \$really_verbose,
	     "est:s" => \$est,
	     "count" => \$count,
	     "report" => \$report,
	     "gap:s"  => \$gap,
	     "transcript" => \$transcript
	   ) ;

my $database = glob("~wormpub/DATABASES/TEST_DBs/ANT_matchingace");
my $gff = glob("~wormpub/DATABASES/TEST_DBs/ANT_matchingace/dump/T23G11.gff");

open (GFF,"<$gff") or die "gff\n";
# C43G2	  curated	  exon         10841   10892   .       +       .       Sequence "C43G2.4"
# C43G2   curated         CDS          10841   10892   .       +       0       Sequence "C43G2.4"
# C43G2   curated         Sequence     10841   13282   .       +       .       Sequence "C43G2.4"
# C43G2   BLAT_EST_BEST   similarity   8877    9029    100     +       .       Target "Sequence:yk1125e01.5" 37 189

my %genes_exons;
my %genes_span;
my %cDNA;
my %cDNA_span;
my $superlink;

# parse GFF file to get Sequence, exon and cDNA info
while (<GFF>) {
  my @data = split;
  $superlink = $data[0];

  #  GENE STRUCTURE
  if( $data[1] eq "curated" ) {
    $data[9] =~ s/\"//g;
    if( $data[2] eq "Sequence" ) {
      # GENE SPAN
      $genes_span{$data[9]} = [($data[3], $data[4], $data[6])];
    }
    elsif ($data[2] eq "exon" ) {
      # EXON 
      $genes_exons{$data[9]}{$data[3]} = $data[4];
    }
  }

  # BLAT DATA
  elsif( ($data[1] eq "BLAT_EST_BEST") or ($data[1] eq "BLAT_OST_BEST") or ($data[1] eq "BLAT_mRNA_BEST") ){
    $data[9] =~ s/\"//g;
    $data[9] =~ s/Sequence:// ;
    $cDNA{$data[9]}{$data[3]} = $data[4];
    # keep min max span of cDNA
    if( !(defined($cDNA_span{$data[9]}[0])) or ($cDNA_span{$data[9]}[0] > $data[3]) ) {
      $cDNA_span{$data[9]}[0] = $data[3];
    } 
    if( !(defined($cDNA_span{$data[9]}[1])) or ($cDNA_span{$data[9]}[1] < $data[4]) ) {
      $cDNA_span{$data[9]}[1] = $data[4];
    }
  }
}

close GFF;

&eradicateSingleBaseDiff;

my %gene2cdnas;
foreach my $CDNA (keys %cDNA_span) {
  $CDNA = $est if $est;
  my @genes = &findOverlappingGene(\$CDNA);
  print "\n=============================\n\n$CDNA overlaps @genes\n" if $verbose;
  foreach (@genes){
    print "____________________________\n\n\nchecking $_ against $CDNA\n";
    if(&checkExonMatch(\$_,\$CDNA) == 1 ){
      print "$CDNA matches $_\n";
      push( @{$gene2cdnas{$_}}, $CDNA);
    }
    else { print "$CDNA NOT match $_\n";}
  }
  last if $est;
}

open (ACE,">$database/transcripts.ace") or die "transcripts\n";
foreach (keys %gene2cdnas) {
  print "$_ matching cDNAs => @{$gene2cdnas{$_}}\n" if $report;
  print "$_ matches ",scalar(@{$gene2cdnas{$_}}),"\n" if $count;
  if( $transcript ) {
    next unless $genes_span{$_}->[2] eq "+"; #just do forward for now
    my %transcript;
    %transcript = %{$genes_exons{$_}};   # put the gene model in to the transcript

    print ACE "\nSequence : \"$_.trans\"\n";
    print ACE "Source T23G11\n";
    print ACE "method history\n";
    foreach my $cdna (@{$gene2cdnas{$_}}) { # iterate thru array of cDNAs ( 1 per gene )
      print ACE "matching_CDNA \"$cdna\"\n";
      foreach my $cExon(keys %{$cDNA{$cdna}} ){ #  then thru exons of the cDNA
	next if( defined $transcript{$cExon} and  $transcript{$cExon} == $cDNA{$cdna}->{$cExon} ); #skip if the cDNA exon is same as one in transcript
	
	# iterate thru existing transcript exons and check for overlaps
	foreach my $transExon (keys %transcript ) {
	  if( $cExon < $transExon and $cDNA{$cdna}->{$cExon} >= $transExon ) {
	
	    if( $cDNA{$cdna}->{$cExon} > $transcript{$transExon} ) {
	      $transcript{$cExon} = $cDNA{$cdna}->{$cExon};
	    }
	    else {
	      $transcript{$cExon} = $transcript{$transExon};
	    }
	    delete $transcript{$transExon};
	  }
	  elsif( $cExon > $transExon and $cExon < $transcript{$transExon} and $cDNA{$cdna}->{$cExon} >= $transExon ) {
	    $transcript{$transExon} = $cDNA{$cdna}->{$cExon};
	  }
	}
      }
    }
    
    my @exons = (sort { $transcript{$a} <=> $transcript{$b} } keys %transcript);
    foreach (@exons) {
      print ACE "source_exons ",$_ - $exons[0] + 1," ",$transcript{$_} - $exons[0]+1," \n";
    }

    print ACE "\nSequence : T23G11\n";
    print ACE "Subsequence $_.trans $exons[0] $transcript{$exons[-1]}\n";
  }
}
close ACE;
exit(0);

sub findOverlappingGene
  {
    my $cdna = shift;
    my @overlap_genes;
    foreach ( sort { $genes_span{$a}[0]<=>$genes_span{$b}[0]  } keys %genes_span ) {
      if ($cDNA_span{$$cdna}[0] > $genes_span{$_}[1] ) { next; } #
      elsif( ($cDNA_span{$$cdna}[0] < $genes_span{$_}[0]) and ($cDNA_span{$$cdna}[1] > $genes_span{$_}[0]) ) { # cDNA starts B4 gene and ends after gene start
	# overlaps 
	push( @overlap_genes, $_);
      }
      elsif( ($cDNA_span{$$cdna}[0] > $genes_span{$_}[0]) and ($cDNA_span{$$cdna}[1] < $genes_span{$_}[1]) ) { # cDNA contained within gene
	#overlaps
	push( @overlap_genes, $_);
      }
      elsif( ($cDNA_span{$$cdna}[0] < $genes_span{$_}[1]) and ($cDNA_span{$$cdna}[1] > $genes_span{$_}[1]) ) { # cDNA starts within and extends past gene end
	#overlaps
	push( @overlap_genes, $_);
      }
      elsif ($cDNA_span{$$cdna}[1] < $genes_span{$_}[0]) { last; }
    }
    return @overlap_genes;
  }


sub checkExonMatch
  {
    my $gene = shift;
    my $cdna = shift;
    my $cdna_exon;
    my $match = 1;
    if ($really_verbose) {
      foreach ( sort keys %{$genes_exons{$$gene}} ) {
	print "\t$_ $genes_exons{$$gene}->{$_}\n" if $really_verbose;
      }
    }

    #sort the cDNA and gene exon starts into ordered array
    my @cdna_exon_starts = sort { $cDNA{$$cdna}->{$a} <=> $cDNA{$$cdna}->{$b} } keys %{$cDNA{$$cdna} } ;
    my @gene_exon_starts = sort { $genes_exons{$$gene}->{$a} <=> $genes_exons{$$gene}->{$b} } keys %{$genes_exons{$$gene}};

    #check if cDNA exon fits with gene model
    foreach my $cExonStart (@cdna_exon_starts) {

      my $gExonS;
      # do cDNA and gene share exon start position
      if( $genes_exons{$$gene}->{$cExonStart} ) {
	if( $genes_exons{$$gene}->{$cExonStart} == $cDNA{$$cdna}->{$cExonStart} ) {
	  #exact match
	  print "\tExact Match\n";
	  next;
	}
	#is this final gene exon
	elsif( $cExonStart == $gene_exon_starts[-1] ) {
	  print "\tMatch - last gene exon\n";
	  next;
	}
	# or final cDNA exon
	elsif ( $cExonStart == $cdna_exon_starts[-1] ) {
	  # . . must terminate within gene exon
	  if ( $cDNA{$$cdna}->{$cExonStart} > $genes_exons{$$gene}->{$cExonStart} ) {
	    print "\tMISS - $$cdna $cExonStart => $cDNA{$$cdna}->{$cExonStart} extends over gene exon boundary\n";
	    $match = 0;
	  }
	  else {
	    print "\tMatch - last cDNA exon\n";
	    next;
	  }
	}
      }

      # do cDNA and gene share exon end position
      elsif ( ($gExonS = &geneExonThatEndsWith(\$$gene, $cDNA{$$cdna}->{$cExonStart}) ) and ($gExonS != 0) ) {
#	# shared exon end
	
	if( $gExonS == $gene_exon_starts[0] ) {	           #is this the 1st gene exon 
	  if ( $cExonStart == $cdna_exon_starts[0] ) {      # also cDNA start so always match
	    print "\tMatch - 1st exons overlap\n";
	    next;
	  }
	  elsif( $cExonStart < $gene_exon_starts[0] ) {      # cDNA exon overlap 1st gene exon
	    print "\tMatch - cDNA exon covers 1st exon\n";
	    next;
	  }
	  else {
	    print "\tMISS - cDNA exon splices in gene exon\n";
	    print "\t\t$$cdna $cExonStart => $cDNA{$$cdna}->{$cExonStart}\n";
	    print "\t\t$$gene $gExonS => $genes_exons{$$gene}->{$gExonS}\n";
	    $match = 0;
	  }
	}
	# exon matched is not 1st of gene
	elsif ( $cExonStart == $cdna_exon_starts[-1] ) {    # start of cDNA
	  print"\tMatch - 1st exon of cDNA starts in exon of gene\n";
	  next;
	}
      }

      # exon end or start not shared with gene model
      else{
	# cDNA exon wholly outside of gene model
	if($cExonStart > $genes_span{$$gene}->[1] or
	   $cDNA{$$cdna}->{$cExonStart} < $genes_span{$$gene}->[0] ) { 
	  print "\tMatch cDNA exon wholly outside of gene\n";
	  next;
	}
	# cDNA exon overlaps and starts in final gene exon and must be 1st of cDNA
	elsif(( $cExonStart == $cDNA_span{$$cdna}->[0] ) and #  1st exon of cDNA
	      ( $cExonStart > $gene_exon_starts[-1] and $cExonStart < $genes_exons{$$gene}->{$gene_exon_starts[-1]} ) and # cDNA exon starts in final exon of gene
	       ( $cDNA{$$cdna}->{$cExonStart} > $genes_exons{$$gene}->{$gene_exon_starts[-1]} )
	     )  {
	  print "cDNA exon start in final gene exon and protrudes past it\n";
	  next;
	}
	# cDNA exon overlaps 1st gene exon start and terminate therein
	elsif( ( $cExonStart == $cdna_exon_starts[-1]  ) and #  last exon of cDNA
	       ( $cExonStart < $genes_span{$$gene}->[0] ) and 
	       ( $cDNA{$$cdna}->{$cExonStart} > $genes_exons{$gene_exon_starts[0]} )
	     ) {
	  print "\tcDNA final exon overlaps first exon of gene and end therein\n";
	  next;
	}
	# single exon cDNA contained wholly within a single gene exon
	elsif( &cDNA_wholelyInExon(\$$gene, \$$cdna) == 1 ) {
	  print "\tsingle exon cDNA contained wholly within a single gene exon\n";
	  next;
	}
	else {
	  print "MISS - exon fell thru\n";
	  $match = 0;
	}
      }
    }

    return $match if $match == 0;

    # now confirm that each gene exon in the overlapping range is covered by the cDNA
    my @range = ($cDNA_span{$$cdna}->[0],$cDNA_span{$$cdna}->[1]);
    $range[0] = $genes_span{$$gene}->[0] if ($genes_span{$$gene}->[0] > $range[0]);
    $range[1] = $genes_span{$$gene}->[1] if ($genes_span{$$gene}->[1] < $range[1]);

  EXONS:
    foreach my $gene_exon (@gene_exon_starts) {

      next if ( ($genes_exons{$$gene}->{$gene_exon} < $range[0]) or ($gene_exon > $range[1]) ) ; # exon out of overlapping range
      my $exon_matched = 0;
      foreach (@cdna_exon_starts) {
	print "\tcomparing $_ : $cDNA{$$cdna}->{$_} with gene exon $gene_exon : $genes_exons{$$gene}->{$gene_exon}\n";
	if( ($_ <= $gene_exon or $_ == $range[0]) and ( $cDNA{$$cdna}->{$_} >= $genes_exons{$$gene}->{$gene_exon} or $cDNA{$$cdna}->{$_} == $range[1]) ) {
	  print "\t$$gene exon $gene_exon =>  $genes_exons{$$gene}->{$gene_exon} covered by $$cdna $_ => $cDNA{$$cdna}->{$_} \n"; 
	  $exon_matched = 1;
	  next EXONS;
	}
      }
      unless( $exon_matched == 1 ){
	$match = 0;
	last EXONS;  # a gene exon has not matched so no point carrying on
      }
    }

    return $match;
  }


sub cDNA_wholelyInExon
  {
    my $gene = shift;
    my $cdna = shift;
    my $cdna_extent0 = $cDNA_span{$$cdna}->[0];
    my $cdna_extent1 = $cDNA_span{$$cdna}->[1];
    foreach (keys %{$genes_exons{$$gene}} ) {
      return 1 if( $_ < $cdna_extent0 and $genes_exons{$$gene}->{$_} > $cdna_extent1);
    }
    return 0;
  }

sub geneExonThatEndsWith
  {
    my $gene = shift;
    my $exon_end = shift;
    foreach (keys %{$genes_exons{$$gene}} ) {
      return $_ if( $genes_exons{$$gene}->{$_} == $exon_end );
    }
    return 0;
  }

sub cDNA_ExonThatEndsWith
  {
    my $cdna = shift;
    my $exon_end = shift;
    foreach (keys %{$cDNA{$$cdna}} ) {
      return $_ if( $cDNA{$$cdna}->{$_} == $exon_end );
    }
    return 0;
  }

sub eradicateSingleBaseDiff
  {
    foreach my $cdna_hash (keys %cDNA ) {
      my $last_key;
      my $check;
      foreach my $exons (sort { $cDNA{$cdna_hash}->{$a} <=> $cDNA{$cdna_hash}->{$b} } keys %{$cDNA{$cdna_hash}}) {
	print "\n############### $cdna_hash #############\n" if $really_verbose;
	print "$exons -> $cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	my $new_last_key = $exons;
	if( $last_key ) {
	  if( $cDNA{$cdna_hash}->{$last_key} >= $exons - $gap ) { #allows seq error gaps up to 3 bp
	    $cDNA{$cdna_hash}->{$last_key} = $cDNA{$cdna_hash}->{$exons};
	    delete $cDNA{$cdna_hash}->{$exons};
	    $check = 1;
	    $new_last_key = $last_key;
	  }
	}
	$last_key = $new_last_key;
      }
      if ( $check ) {
	foreach my $exons (sort keys  %{$cDNA{$cdna_hash}}) {
	  print "single base diffs removed from $cdna_hash\n" if $really_verbose;
	  print "$exons -> $cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	}
      }
    }
  }

#!/usr/local/bin/perl5.6.1 -w

use strict;
use lib "/wormsrv2/scripts/"; 
use Getopt::Long;
use Data::Dumper;
use Coords_converter;
use Wormbase;

my $tace = &tace;
my $rundate = `date +%y%m%d`; chomp $rundate;

my $log = &make_log_file; 

#open(FH, ">/wormsrv2/logs/transcript_builder_$rundate") || die $!;
#my $log = (*FH);

my ($debug, $help, $verbose, $really_verbose, $est, $count, $report, $gap, $transcript, $gff, $show_matches, $database, $overlap_check, $load_matches, $load_transcripts, $build);

$gap = 5; # $gap is the gap allowed in an EST alignment before it is considered a "real" intron

GetOptions ( "debug" => \$debug,
	     "help" => \$help,
	     "verbose" => \$verbose,
	     "really_verbose" => \$really_verbose,
	     "est:s" => \$est,
	     "count" => \$count,
	     "report" => \$report,
	     "gap:s"  => \$gap,
	     "transcript" => \$transcript,
	     "gff:s"    => \$gff,
	     "show_matches"    => \$show_matches,
	     "database:s"   => \$database,
	     "overlap"      => \$overlap_check,
	     "load_transcripts" => \$load_transcripts,
	     "load_matches"    => \$load_matches,
	     "build"           => \$build
	   ) ;

&check_opts; # if -build set, this will set all relevant opts to works as if in build. Will NOT overwrite the others (eg -count)

$database = glob("~wormpub/DATABASES/TEST_DBs/transcripts") unless $database;


my %genes_exons;
my %genes_span;
my %cDNA;
my %cDNA_span;
my %transcript_span;
my $gff_file;
my @ordered_genes;

$gff_file = $gff if $gff;

my @chromosomes = qw(I II III IV V X);
foreach my $chrom ( @chromosomes ) {

  #Make sure these are empty
  %genes_exons = ();
  %genes_span= ();
  %cDNA = ();
  %cDNA_span = ();
  %transcript_span = ();
  @ordered_genes = ();
  
  unless ($gff) {
    $gff_file = "$database/CHROMOSOMES/CHROMOSOME_$chrom.gff";
  }

  open (GFF,"<$gff_file") or die "gff_file\n";
  print $log "reading gff file $gff_file\n";
  # C43G2	  curated	  exon         10841   10892   .       +       .       Sequence "C43G2.4"
  # C43G2   curated         CDS          10841   10892   .       +       0       Sequence "C43G2.4"
  # C43G2   curated         Sequence     10841   13282   .       +       .       Sequence "C43G2.4"
  # C43G2   BLAT_EST_BEST   similarity   8877    9029    100     +       .       Target "Sequence:yk1125e01.5" 37 189

  # parse GFF file to get Sequence, exon and cDNA info
  while (<GFF>) {
    my @data = split;
    
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

  &checkData(\$gff); # this just checks that there is some BLAT and gene data in the GFF file

  &eradicateSingleBaseDiff;

  # generate ordered array of genes to use as keys in sub findOverlappingGenes
  foreach ( sort { $genes_span{$a}[0]<=>$genes_span{$b}[0]  } keys %genes_span ) {
    push(@ordered_genes,$_);
  }

  my %gene2cdnas;
  foreach my $CDNA (keys %cDNA_span) {
    $CDNA = $est if $est;
    my @genes = &findOverlappingGene(\$CDNA);
    print "\n=============================\n\n$CDNA overlaps @genes\n" if $verbose;
    foreach (@genes){
      print "____________________________\n\n\nchecking $_ against $CDNA\n" if $verbose;
      if(&checkExonMatch(\$_,\$CDNA) == 1 ){
	print "$CDNA matches $_\n" if $verbose;
	push( @{$gene2cdnas{$_}}, $CDNA);
      }
      else { print "$CDNA NOT match $_\n" if $verbose;}
    }
    last if $est;
  }


  #This is the place to check for paired reads that dont overlap a gene model


  my $coords;
  # write out the transcript objects
  if( $transcript ) {
    system("mkdir $database/TRANSCRIPTS") unless -e "$database/TRANSCRIPTS";
    open (ACE,">$database/TRANSCRIPTS/transcripts_$chrom.ace") or die "transcripts\n";
    # get coords obj to return clone and coords from chromosomal coords
    $coords = Coords_converter->invoke($database);
  }
  
  foreach my $gene (keys %gene2cdnas) {
    print  "$gene matching cDNAs => @{$gene2cdnas{$gene}}\n" if $report;
    print "$_ matches ",scalar(@{$gene2cdnas{$gene}}),"\n" if $count;
    if( $transcript ) {
      # next unless $genes_span{$_}->[2] eq "-"; #just do forward for now

      # put the gene model in to the transcript
      my %transcript;
      %transcript = %{$genes_exons{$gene}};   
      my @transcript_span = @{ $genes_span{$gene} };

      #transcript object
      print ACE "\nTranscript : \"$gene\"\nSpecies \"Caenorhabditis elegans\"\n";
      foreach my $cdna (@{$gene2cdnas{$gene}}) { # iterate thru array of cDNAs ( 1 per gene )

	print ACE "Matching_cDNA \"$cdna\"\n";
	foreach my $cExon(keys %{$cDNA{$cdna}} ){ #  then thru exons of the cDNA

	  next if( defined $transcript{$cExon} and  $transcript{$cExon} == $cDNA{$cdna}->{$cExon} ); #skip if the cDNA exon is same as one in transcript
	
	  # if exon outside range of current transcript add it to the transcript UTR
	  if( $cExon > $transcript_span[1] or  $cDNA{$cdna}->{$cExon} < $transcript_span[0] ){
	    # add to transcript hash
           $transcript{$cExon} = $cDNA{$cdna}->{$cExon};

	   #expand transcript span to accomodate new exon
	   $transcript_span[0] =  $cExon if( $cExon < $transcript_span[0] );
	   $transcript_span[1] =  $cDNA{$cdna}->{$cExon} if( $cDNA{$cdna}->{$cExon} > $transcript_span[1] );
	  }

	  # iterate thru existing transcript exons and check for overlaps
	  foreach my $transExon (keys %transcript ) {
	    if( $cExon < $transExon and $cDNA{$cdna}->{$cExon} >= $transExon ) { # 

	      if( $cDNA{$cdna}->{$cExon} > $transcript{$transExon} ) {
		$transcript{$cExon} = $cDNA{$cdna}->{$cExon}; # ever possible ?
	      }
	      else {
		$transcript{$cExon} = $transcript{$transExon};
	      }
	      delete $transcript{$transExon};
	    }
	    elsif( $cExon >= $transExon and $cExon <= $transcript{$transExon} and $cDNA{$cdna}->{$cExon} > $transcript{$transExon} ) {
	      $transcript{$transExon} = $cDNA{$cdna}->{$cExon};
	    }
	  }
	}
      }

      my @exons = (sort { $transcript{$a} <=> $transcript{$b} } keys %transcript);
      foreach (@exons) {
	if( $genes_span{$gene}->[2] eq "+"){
	  print ACE "Source_exons ",$_ - $exons[0] + 1," ",$transcript{$_} - $exons[0]+1," \n"; # + strand
	}
	else {
	  print ACE "Source_exons ",$transcript{$exons[-1]} - $transcript{$_} + 1 ," ",$transcript{$exons[-1]} - $_ + 1 ," \n"; # - strand
	}
      }
      ($chrom) = $gff =~ /.*_(\w+)/ if $gff;
      my($source, $x, $y ) = $coords->LocateSpan("$chrom",$exons[0],$transcript{$exons[-1]});
      $transcript_span{"$gene.trans"} = [ ($exons[0],$transcript{$exons[-1]}) ];
      print ACE "Sequence \"$source\"\n";
      print ACE "Method \"transcript\"\n";

      # parent object
      print ACE "\nSequence : \"$source\"\n";
      if( $genes_span{$gene}->[2] eq "+"){
	print ACE "Transcript_child $gene $x $y\n"; # + strand
      }
      else {
	print ACE "Transcript_child $gene $y $x\n"; # - strand
      }   

      #link gene to transcript
      print ACE "\nSequence \"$gene\"\n";
      print ACE "Corresponding_transcript $gene\n" 
    }
  }
  close ACE if $transcript;

  if ($show_matches) { 
    system("mkdir $database/TRANSCRIPTS") unless -e "$database/TRANSCRIPTS";
    open(MATCHES,">$database/TRANSCRIPTS/chromosome${chrom}_matching_cDNA.dat");
    print $log "writing Matching_cDNA file $database/TRANSCRIPTS/chromosome${chrom}_matching_cDNA.dat\n";
    print MATCHES Data::Dumper->Dump([\%gene2cdnas]);
    close MATCHES;

    open (CDNAS, ">$database/chromosome${chrom}_matching_cDNA.ace") or die "output ace for chromosome $chrom\n";
    foreach my $gene ( keys %gene2cdnas ) {
      print CDNAS "\nSequence : \"$gene\"\n";
      foreach my $cdna ( @{$gene2cdnas{$gene}} ) {
	print CDNAS "Matching_cDNA \"$cdna\" Inferred_automatically \"transcript_builder.pl\"\n";
      }
    }
    close CDNAS;
  }


  &checkOverlappingTranscripts($chrom) if $overlap_check;

  last if $gff; # if only doing a specified gff file exit after this is complete
}

if( $load_transcripts ) {
  print $log "loading transcripts file to $database\n";
  system("cat $database/TRANSCRIPTS/ranscripts_*.ace > $database/TRANSCRIPTS/transcripts_all.ace");
  open(TACE, " | $tace -tsuser transcripts $database |") or warn "could open $database to load transcript file\n";
  print TACE "pparse $database/TRANSCRIPTS/transcripts_all.ace\nsave\nquit";
  close TACE;
}

if( $load_matches ) {
  print $log "loading matching_cDNA file to $database\n";
  system("cat $database/TRANSCRIPTS/chromosome*_matching_cDNA.ace > $database/TRANSCRIPTS/matching_cDNA_all.ace");
  open(TACE, " | $tace -tsuser matching_cDNA $database |") or warn "could open $database to load matching_cDNA file\n";
  print TACE "pparse $database/TRANSCRIPTS/matching_cDNA_all.ace\nsave\nquit";
  close TACE;
}
  
print $log "$0 finished at ",&runtime,"\n";
close $log;


exit(0);

sub findOverlappingGene
  {
    my $cdna = shift;
    my @overlap_genes;
    foreach ( @ordered_genes ) { # use this ordered array instead of sorting hash for every gene.
      print "testing overlap $_ $$cdna\n" if $verbose;
      if ($cDNA_span{$$cdna}[0] > $genes_span{$_}[1] ) { next; } #
      elsif( ($cDNA_span{$$cdna}[0] < $genes_span{$_}[0]) and ($cDNA_span{$$cdna}[1] > $genes_span{$_}[0]) ) { # cDNA starts B4 gene and ends after gene start
	# overlaps 
	push( @overlap_genes, $_);
      }
      elsif( ($cDNA_span{$$cdna}[0] >= $genes_span{$_}[0]) and ($cDNA_span{$$cdna}[1] <= $genes_span{$_}[1]) ) { # cDNA contained within gene
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
	  print "\tExact Match\n" if $verbose;
	  next;
	}
	#is this final gene exon
	elsif( $cExonStart == $gene_exon_starts[-1] ) {
	  print "\tMatch - last gene exon\n" if $verbose;
	  next;
	}
	# or final cDNA exon
	elsif ( $cExonStart == $cdna_exon_starts[-1] ) {
	  # . . must terminate within gene exon
	  if ( $cDNA{$$cdna}->{$cExonStart} > $genes_exons{$$gene}->{$cExonStart} ) {
	    print "\tMISS - $$cdna $cExonStart => $cDNA{$$cdna}->{$cExonStart} extends over gene exon boundary\n" if $verbose;
	    $match = 0;
	  }
	  else {
	    print "\tMatch - last cDNA exon\n" if $verbose;
	    next;
	  }
	}
      }

      # do cDNA and gene share exon end position
      elsif ( ($gExonS = &geneExonThatEndsWith(\$$gene, $cDNA{$$cdna}->{$cExonStart}) ) and ($gExonS != 0) ) {
#	# shared exon end
	
	if( $gExonS == $gene_exon_starts[0] ) {	           #is this the 1st gene exon 
	  if ( $cExonStart == $cdna_exon_starts[0] ) {      # also cDNA start so always match
	    print "\tMatch - 1st exons overlap\n" if $verbose;
	    next;
	  }
	  elsif( $cExonStart < $gene_exon_starts[0] ) {      # cDNA exon overlap 1st gene exon
	    print "\tMatch - cDNA exon covers 1st exon\n" if $verbose;
	    next;
	  }
	  else {
	    print "\tMISS - cDNA exon splices in gene exon\n" if $verbose;
	    print "\t\t$$cdna $cExonStart => $cDNA{$$cdna}->{$cExonStart}\n" if $verbose;
	    print "\t\t$$gene $gExonS => $genes_exons{$$gene}->{$gExonS}\n" if $verbose;
	    $match = 0;
	  }
	}
	# exon matched is not 1st of gene
	elsif ( ($cExonStart == $cdna_exon_starts[0] ) and # start of cDNA
		($cExonStart >$gExonS ) ) {                   # . . . is in gene exon
	  print"\tMatch - 1st exon of cDNA starts in exon of gene\n" if $verbose;
	  next;
	}
	else {
	  print "MISS - exon $$cdna : $cExonStart => $cDNA{$$cdna}->{$cExonStart} overlaps start of gene exon : $gExonS => $genes_exons{$$gene}->{$gExonS}\n" if $verbose;
	  $match = 0;
	}
      }

      # exon end or start not shared with gene model
      else{
	# cDNA exon wholly outside of gene model
	if($cExonStart > $genes_span{$$gene}->[1] or
	   $cDNA{$$cdna}->{$cExonStart} < $genes_span{$$gene}->[0] ) { 
	  print "\tMatch cDNA exon wholly outside of gene\n" if $verbose;
	  next;
	}
	# cDNA exon overlaps and starts in final gene exon and must be 1st of cDNA
	elsif(( $cExonStart == $cDNA_span{$$cdna}->[0] ) and #  1st exon of cDNA
	      ( $cExonStart > $gene_exon_starts[-1] and $cExonStart < $genes_exons{$$gene}->{$gene_exon_starts[-1]} ) and # cDNA exon starts in final exon of gene
	       ( $cDNA{$$cdna}->{$cExonStart} > $genes_exons{$$gene}->{$gene_exon_starts[-1]} )
	     )  {
	  print "cDNA exon start in final gene exon and protrudes past it\n" if $verbose;
	  next;
	}
	# cDNA exon overlaps 1st gene exon start and terminate therein
	elsif( ( $cExonStart == $cdna_exon_starts[-1]  ) and #  last exon of cDNA
	       ( $cExonStart < $genes_span{$$gene}->[0] ) and 
	       ( $cDNA{$$cdna}->{$cExonStart} < $genes_exons{$$gene}->{$gene_exon_starts[0]} )
	     ) {
	  print "\tcDNA final exon overlaps first exon of gene and end therein\n" if $verbose;
	  next;
	}
	# single exon cDNA contained wholly within a single gene exon
	elsif( &cDNA_wholelyInExon(\$$gene, \$$cdna) == 1 ) {
	  print "\tsingle exon cDNA contained wholly within a single gene exon\n" if $verbose;
	  next;
	}
	else {
	  print "MISS - exon fell thru\n" if $verbose;
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
	print "\tcomparing $_ : $cDNA{$$cdna}->{$_} with gene exon $gene_exon : $genes_exons{$$gene}->{$gene_exon}\n" if $verbose;
	if( ($_ <= $gene_exon or $_ == $range[0]) and ( $cDNA{$$cdna}->{$_} >= $genes_exons{$$gene}->{$gene_exon} or $cDNA{$$cdna}->{$_} == $range[1]) ) {
	  print "\t$$gene exon $gene_exon =>  $genes_exons{$$gene}->{$gene_exon} covered by $$cdna $_ => $cDNA{$$cdna}->{$_} \n" if $verbose; 
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
    print $log "removing small cDNA mismatches\n\n\n";
    foreach my $cdna_hash (keys %cDNA ) {
      my $last_key;
      my $check;
      foreach my $exons (sort { $cDNA{$cdna_hash}->{$a} <=> $cDNA{$cdna_hash}->{$b} } keys %{$cDNA{$cdna_hash}}) {
	print "\n############### $cdna_hash #############\n" if $really_verbose;
	print "$exons -> $cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	my $new_last_key = $exons;
	if( $last_key ) {
	  if( $cDNA{$cdna_hash}->{$last_key} >= $exons - $gap ) { #allows seq error gaps up to $gap bp
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




sub checkOverlappingTranscripts  {
  my $chrom = shift;
  print $log "checking overlapping transcripts for chromosome_$chrom\n";
  my @ordered_transcripts;
  foreach my $transcript ( sort { $transcript_span{$a}[0]<=>$transcript_span{$b}[0]  } keys %transcript_span) {
    push (@ordered_transcripts, $transcript);
  }

  my $trans_count = scalar (@ordered_transcripts);

  open (OLT,">$database/overlapping_transcripts_$chrom") or warn "cant open $database/overlapping_transcripts\n";
  TRANS:
    for( my $i = 0;$i < $trans_count; $i++) {
      for( my $j = 1; $j < 6; $j++ ) {
	next unless ($i-$j >= 0 and $transcript_span{ $ordered_transcripts[$i-$j] });
	if ( $transcript_span{ $ordered_transcripts[$i] }->[0] < $transcript_span{ $ordered_transcripts[$i-$j] }->[1] ) {
	  {
	    my $gene1 = substr ($ordered_transcripts[($i-$j)] ,0,-6);
	    my $gene2 = substr ($ordered_transcripts[$i] ,0,-6);
	    next TRANS if &checkIsos(\$gene1, \$gene2, *OLT) == 1;
	  }
	}
	next unless $transcript_span{ $ordered_transcripts[$i+$j] };
       if ( $transcript_span{ $ordered_transcripts[$i] }->[1] > $transcript_span{ $ordered_transcripts[$i+$j] }->[0] ) {
	 {
	   my $gene1 = substr ($ordered_transcripts[($i+$j)] ,0,-6);
	   my $gene2 = substr ($ordered_transcripts[$i] ,0,-6);
	   next TRANS if &checkIsos(\$gene1, \$gene2, *OLT) == 1;	 
	}
       }
      }
    }
  sub checkIsos
    {
      my $gene1 = shift;
      my $gene2 = shift;
      my $file  = shift;
      if( $genes_span{$$gene1}->[2] eq $genes_span{$$gene2}->[2] ){ # check same strand strand
	# dismiss isoforms
	my ($gene1_name) = $$gene1 =~ /(\w+\.\d+)\w?/;
	my ($gene2_name) = $$gene2 =~ /(\w+\.\d+)\w?/;
	
	unless( $gene1_name eq $gene2_name ) {
	  print $file "overlapping transcripts\t$$gene1\t$$gene2\n";
	  return 1;
	}
      }
      return 0;
    }
  close OLT;
}

sub check_opts {
  print $log "checking options . . . \n\n\n";
  # sanity check options
  if( $help ) {
    system("perldoc $0");
    exit(0);
  }

  unless ($transcript ) {
    if( $load_transcripts ) {
      die print "\n\n\nyou have chosen to load_transcripts without generating them\n Try adding -transcript if you want to build them first\n\n\n\n";
    }
    if( $overlap_check ) {
      die print "\n\n\nyou have chosen to check for overlapping transcripts without generating them\n Try adding -transcript if you want to build them first\n\n\n\n";
    }
    if( $build ) {
      $database = "/wormsrv2/autoace";
      $transcript = 1;
      $show_matches = 1;
      $load_transcripts = 1;
      $load_matches = 1;
      $overlap_check = 1;
    }
  }

  if( defined($load_matches) & !defined($show_matches ) ) {
    die print "\n\n\nyou have chosen to load_matches without generating them\n Try adding -show_matches if you want to generate them first\n\n\n\n";
  }
}

sub checkData
  {
    my $file = shift;
    die "There's no BLAT data in the gff file $$file\n" if scalar keys %cDNA_span == 0;
    die "There are no genes in the gff file $$file\n" if scalar keys %genes_span == 0;
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


The script outputs up to 3 files for each chromosome.

  chromosome*_transcripts.ace    -     transcript object ace file 

  chromosome*_matching_cDNA.ace  -     matching_cDNA assignments

  chromosome*_matching_cDNA.dat  -     matching_cDNA assignments as Data::Dumper hash   gene => [cdna, cdna]

=back

=head2 transcript_builder.pl arguments:

=over 4

=item * verbose and really_verbose  -  levels of terminal output
  
=item *  est:s     - just do for single est 

=item * count     - reports number of matching_cDNAs for each gene 

=item * report    - reports list of matching_cDNA for each gene

=item * gap:s      - when building up cDNA exon structures from gff file there are often single / multiple base pair gaps in the alignment. This sets the gap size that is allowable [ defaults to 5 ]

=item *  transcript   - write transcipt files

=item * gff:s         - pass in a gff file to use

=item * show_matches   - write the matching_cDNA ace files

=item * database:s      - database directory to use. Expects gff files to be in CHROMSOMES subdirectory.

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

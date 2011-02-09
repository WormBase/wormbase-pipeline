# functions used to map DNA_text to the genome, parse BLAT hits and identify overlapping genes


sub RunBlat {

    my $dna=shift;            # probe name that is stored under DNA_text tag in RNAi objects
    my $probe_name=shift;
    my $db=shift;             # acedb - used to in getOverlappingGenes together with GFF MySQL database

    my $first=1;

    $dna=~s/\n//g;
    $dna=~s/\s+//g;
	
    $temp_in_file=$temp_dir.int(rand(1000000));
    while (-e $temp_in_file) {
	$temp_in_file=$temp_dir.int(rand(1000000));
    }
    $temp_out_file=$temp_dir.int(rand(1000000));
    while (-e $temp_out_file) {
	$temp_out_file=$temp_dir.int(rand(1000000));
    }

    unless (-e $db_file) {  # $db_file is BLAT database (fasta)
	PrintError("Cannot run BLAT. Database file $db_file does not exist.");
	exit;
    }
    
    open OUT, ">$temp_in_file" || die $!;
    print OUT ">tmp_sequence\n";
    print OUT "$dna\n";
    close OUT;
    
    `blat $db_file $temp_in_file $temp_out_file -minIdentity=80 -maxIntron=10000 -noHead`;  # -minIdentity=80 -maxIntron=10000 are required to identify secondary targets

    open (IN, "<$temp_out_file") || die $!;
    my %blat_results=();
    while (<IN>) {
	chomp;
	next if /psLayout/;
	next if /match/;
	next if /------------------------------/;
	next unless $_;
	my @data=split('\t');
	push @{$blat_results{$data[9]}}, $_;
    }
    close IN;
    
    my $sqldb=Bio::DB::GFF->new(-dsn=>$db_name) || die "cannot open $db_name $!";   # $db_name is GFF MySQL database name - used to convert chromosomal coordinates to genomic canonical and in getOverlappingGenes


    my ($count,$ambigously_mapped_count,$ambiguous,$start,$stop,$strand,$mapped_count,$quality_match_count,$quality_partial_match_count,$quality_length_count,$quality_identity_count,$qstop,$qstart,$qsize,$unique_count);
    
    my ($seg, @features, $rel_name, $rel_start, $rel_stop, $gc, $tmp_name, $tmp_start, $tmp_stop);

    my ($primary_length, $primary_quality, $secondary_length, $secondary_quality, $max_hit_number, $block_size)  = (100, 95, 200, 80, 10, 1);  # parameters for primary/secondary targets
    my ($query_gap_cutoff, $target_gap_cutoff);
    
    foreach my $name (sort {$a cmp $b} keys %blat_results) {
	$count++;
	my $not_empty=0;
	my $hit_count=0;
	my @multiple_hits=();
	foreach (@{$blat_results{$name}}) {
	    my @blat=split('\t');
	    $hit_count++;
	    push @multiple_hits, \@blat;
	}
	
	my %best_hit_data_hash=%{findTargets(\@multiple_hits, $primary_length, $primary_quality, $secondary_length, $secondary_quality, $max_hit_number, $query_gap_cutoff, $target_gap_cutoff, $block_size)};  # find hits that satisfy primary/secondary criteria

	if ($probe_name=~/yk\d+/ || $probe_name=~/mv_/) {   # since true sequences of yk clones are not known, use just the best hit to avoid false hits from intronic sequences - a hack; same for Orfeome PCR products mv_
	    delete $best_hit_data_hash{'RNAi_secondary'};
	    splice @{$best_hit_data_hash{'RNAi_primary'}}, 1;
	}
		
	foreach my $method (keys %best_hit_data_hash) {
	    foreach (@{$best_hit_data_hash{$method}}) {
		    
	    my @best_hit_data=@{$_};
	    my $match_length=shift(@best_hit_data);
	    my $match_quality=shift(@best_hit_data);

	    $best_hit_data[13]=~/(CHROMOSOME_)(.*)/;
	    my $chrom=$2;
	    
	    $best_hit_data[15]++;  # blat coordinates are off-by-one
	    
	    ($start, $stop)=($best_hit_data[15], $best_hit_data[16]) if $best_hit_data[8] eq "+";
	    ($start, $stop)=($best_hit_data[16], $best_hit_data[15]) if $best_hit_data[8] eq "-";
	    ($qsize,$qstart,$qstop)=($best_hit_data[10], $best_hit_data[11], $best_hit_data[12]);
	    $strand=$best_hit_data[8];
	    
	    $mapped_count++;
	    
	    $seg=$sqldb->segment(-name=>$best_hit_data[13],-start=>$start,-end=>$stop, -absolute=>'1') || die "cannot fetch segment";  # start converting to genomic canonical using GFF database - pretty slow

	    my @f=$seg->contained_in(-type=>'region:Genomic_canonical');
	    if (!@f) {
		@f=grep {/super/i} $seg->contained_in(-type=>'region:Link');
	    }
	    if (!@f) {
		@f=$seg->contained_in(-type=>'region:Link');
	    }
	    my $rel_name=$f[0]->name;
	    $seg->refseq($rel_name);
	    my $rel_start=$seg->start;
	    my $rel_stop=$seg->stop;

            my @exon_starts=();
	    my @exon_ends=();
		
	    my @self_exon_starts=();   #within CDS itself
	    my @self_exon_ends=();
		
	    if ($best_hit_data[8] eq "+") {
		@exon_starts=split(',', $best_hit_data[20]);
		@self_exon_starts=split(',', $best_hit_data[19]);
		my @blocks=split(',', $best_hit_data[18]);
		my $diff=$start-$rel_start;
		for (my $e=0; $e<=$#exon_starts; $e++) {
		    $exon_starts[$e]++;
		    $exon_starts[$e]-=$diff;
		    $exon_ends[$e]=$exon_starts[$e]+$blocks[$e]-1;
		    
		    $self_exon_starts[$e]++;
		    $self_exon_ends[$e]=$self_exon_starts[$e]+$blocks[$e]-1;
		}
	    }
	    else {
		@exon_ends=split(',', $best_hit_data[20]);
		@self_exon_ends=split(',', $best_hit_data[19]);
		my @blocks=split(',', $best_hit_data[18]);
		my $diff=$stop-$rel_start;
		for (my $e=0; $e<=$#exon_ends; $e++) {
		    $exon_ends[$e]++;
		    $exon_ends[$e]-=$diff;
		    $exon_starts[$e]=$exon_ends[$e]+$blocks[$e]-1;
		    $self_exon_ends[$e]=$qsize-$self_exon_ends[$e];
		    $self_exon_starts[$e]=$self_exon_ends[$e]-($blocks[$e]-1);
		}
		    
		@exon_starts=sort {$a<=>$b} @exon_starts;
		@exon_ends=sort {$a<=>$b} @exon_ends;
		    
		@self_exon_starts=sort {$b<=>$a} @self_exon_starts;
		@self_exon_ends=sort {$b<=>$a} @self_exon_ends;
	    }

	   
	    
	    unless ($rel_name) {
		PrintError("$name not mapped correctly.");
		exit;
	    }
	    if ($rel_name) { 
		#do not merge adjacent blocks
		my $gc_length=$sqldb->segment(-name=>$rel_name)->length;

		print ACE "RNAi : \"$wbrnai\"\n";
		print ACE "Homol_homol\t\"$rel_name:RNAi\"\n";   # record homol_data object in RNAi
		if ($first) {                                    # record DNA_text - this can be done outside on RunBlat
		    if ($probe_name eq 'NULL') {
			print ACE "DNA_text\t\"".lc $dna."\"\t";
			print ACE "\"probe_$probe_count:$rel_name\"\n";
		    }
		    else {
			print ACE "DNA_text\t\"".lc $dna."\"\t";
			print ACE "\"$probe_name\"\n";
		    }
		    $first=0;
		}
		print ACE "\n";
     		
		print ACE "Homol_data : \"$rel_name:RNAi\"\n";   # add RNAi object data to Homol_data object
		print ACE "Sequence\t\"$rel_name\"\n";
		my %overlapping_genes;
		for (my $e=0; $e<=$#exon_starts; $e++) {
		    my $seg;
		    if ($best_hit_data[8] eq "+") {
			print ACE "RNAi_homol\t\"$wbrnai\"\t\"$method\"\t$match_quality\t$exon_starts[$e]\t$exon_ends[$e]\t$self_exon_starts[$e]\t$self_exon_ends[$e]\n";
			$seg=$sqldb->segment(-name=>$rel_name, -start=>$exon_starts[$e], -stop=>$exon_ends[$e]) || die "cannot fetch segment $rel_name:$exon_starts[$e]..$exon_ends[$e]:$!\n";
		    }
		    else {
			print ACE "RNAi_homol\t\"$wbrnai\"\t\"$method\"\t$match_quality\t$exon_ends[$e]\t$exon_starts[$e]\t$self_exon_ends[$e]\t$self_exon_starts[$e]\n";
			$seg=$sqldb->segment(-name=>$rel_name, -start=>$exon_starts[$e], -stop=>$exon_ends[$e]) || die "cannot fetch segment $rel_name:$exon_ends[$e]..$exon_starts[$e]:$!\n";
			
		    }
		    my $temp_ref=getOverlappingGenes($seg, $db);  # find overlapping gene - this is done at this point only for curator's reference and is commented out in the ace file
		    foreach (keys %{$temp_ref}) {
			$overlapping_genes{$_}=$$temp_ref{$_};
		    }
		}
		print ACE "\n";
		
		print ACE "Sequence : \"$rel_name\"\n";
		print ACE "Homol_data\t\"$rel_name:RNAi\"\t1\t$gc_length\n";
		print ACE "\n";

		if (%overlapping_genes) {
		    print ACE "\n";
		    print ACE "//Overlapping_genes\t", join("|", map {"$_($overlapping_genes{$_})"}sort {$a cmp $b} keys %overlapping_genes), "\n";
		    print ACE "\n";
		}
    
	    }
	}
    }
    }
  
    unlink $temp_in_file, $temp_out_file;

}



sub findTargets { # find hits that satisfy primary/secondary criteria
    my ($blat_ref, $primary_length, $primary_quality, $secondary_length, $secondary_quality, $max_hit_number, $query_gap_cutoff, $target_gap_cutoff, $block_size)=@_;

    my %targets=();
    
    my $hits;
    my @blat_data=();
    my $i=0;
    my $j=0;

    foreach (@$blat_ref) {
	$j=0;
	foreach (@$_) {
	    $blat_data[$i][$j]=$_;
	    $j++;
	}
	$i++;
    }
    $hits=$i;
    
    for ($i=0; $i<$hits; $i++) {
	my $target_type='';
	if ((defined($query_gap_cutoff) && $blat_data[$i][5] > $query_gap_cutoff) || (defined($target_gap_cutoff) && $blat_data[$i][7] > $target_gap_cutoff)) {
	    next;
	}

	my $match_length=$blat_data[$i][0]+$blat_data[$i][1];
	my $match_quality=100*$blat_data[$i][0]/$match_length;
	if ($match_quality >= $primary_quality && $match_length >= $primary_length) {
	    $target_type='RNAi_primary';
	}
	elsif ($match_quality >= $secondary_quality && $match_length >= $secondary_length) {
	    $target_type='RNAi_secondary';
	}
	else {
	    next;
	}
	
	push (@{$targets{$target_type}}, [$match_length, $match_quality, @{$blat_data[$i]}]); # store all hits by target type
    }

    my $result_count=0;
    foreach (keys %targets) {
	@{$targets{$_}}=sort {$$b[0] <=> $$a[0]} @{$targets{$_}}; # sort by match length - longest hit first
	$result_count+=scalar @{$targets{$_}};
    }
    
    if ($max_hit_number && $result_count > $max_hit_number) {     # return at most $max_hit_number of hits; remove secondary first, then shorted primary
	if ($targets{'RNAi_primary'} && scalar @{$targets{'RNAi_primary'}} >= $max_hit_number) {
	    delete $targets{'RNAi_secondary'};
	    splice @{$targets{'RNAi_primary'}}, $max_hit_number;
	}
	else {
	    splice @{$targets{'RNAi_secondary'}}, $max_hit_number-scalar @{$targets{'RNAi_primary'}};
	}
    }

    if ($result_count > 1 && defined($block_size)) {
	my %tmp_targets=();
	foreach my $target_type (keys %targets) {
	    for (my $i=0; $i <= $#{$targets{$target_type}}; $i++) {
		if ($target_type eq 'RNAi_primary' && $i == 0) {
		    push (@{$tmp_targets{$target_type}}, $targets{$target_type}[$i]);  # always return at least one (best) hit		    
		}
		else {
		    my @blat_data = @{$targets{$target_type}[$i]};
		    shift @blat_data;
		    shift @blat_data;
		    if ($blat_data[0] < $blat_data[10]) {  # match length is less than query length - go through blocks to make sure that at least one matches the length criteria (ensures that there is a continuous stretch of matched sequence that satisfies the lentgh criterion)
			my @blocks=split(/,/, $blat_data[18]);
			my $match;
			foreach (@blocks) {
			    if ($target_type eq 'RNAi_primary' && $_ >= $primary_length) {
				$match=1;
				last;
			    }
			    if ($target_type eq 'RNAi_secondary' && $_ >= $secondary_length) {
				$match=1;
				last;
			    }
			}
			if ($match) {
			    push (@{$tmp_targets{$target_type}}, $targets{$target_type}[$i]);
			}
		    }
		}
	    }
	}
	%targets=%tmp_targets;
    }
				
    return \%targets;
}




sub getOverlappingGenes { # find overlapping gene - this is done at this point only for curator's reference and is commented out in the ace file
    my $segment=shift;
    my $db=shift;
    my %genes;

    my %transcript_hash=();
    my %cds_hash=();
    my %pseudo_hash=();
    my %non_coding_transcript_hash=();
    my $non_coding_transcript_count=0;
    my $transcript_count=0;
    my $cds_count=0;
    my $pseudo_count=0;

    my @features=$segment->features(-types=>'exon');
    
    foreach (@features) {
	my $fname=$_->name;
	my $fsource=$_->source;
	if ($fsource eq "Coding_transcript") {
	    unless ($transcript_hash{$fname}) {
		$transcript_hash{$fname}=1;
		$transcript_count++;
	    }
	}
	if ($fsource eq "Non_coding_transcript") {
	    unless ($non_coding_transcript_hash{$fname}) {
		$non_coding_transcript_hash{$fname}=1;
		$non_coding_transcript_count++;
	    }
	}
	if ($fsource=~/curated/) {        #coding_exon:curated, exon:curated, intron:curated, CDS:curated
	    unless ($cds_hash{$fname}) {
		$cds_hash{$fname}=1;
		$cds_count++;
	    }
	}
	if ($fsource=~/Pseudogene/) {     #added for pseudogene mapping 6/8/04
	    unless ($pseudo_hash{$fname}) {
		$pseudo_hash{$fname}=1;
		$pseudo_count++;
	    }
	}
    }
    
    foreach (keys %transcript_hash) {
	my ($tr)=$db->fetch('Transcript', $_);
	my $cds=$tr->Corresponding_CDS if $tr;
	my $gene=$cds->Gene if $cds;
	$genes{$gene}=$gene->Public_name if $gene;
    }
    foreach (keys %non_coding_transcript_hash) {
	my ($tr)=$db->fetch('Transcript', $_);
	my $gene=$tr->Gene if $tr;
	$genes{$gene}=$gene->Public_name if $gene;
    }
    foreach (keys %pseudo_hash) {
	my ($ps)=$db->fetch('Pseudogene', $_);
	my $gene=$ps->Gene if $ps;
	$genes{$gene}=$gene->Public_name if $gene;
    }
    foreach (keys %cds_hash) {
	my ($cds)=$db->fetch('CDS', $_);
	my $gene=$cds->Gene if $cds;
	$genes{$gene}=$gene->Public_name if $gene;
    }

    if (%genes) {
	return \%genes;
    }
    else {
	return undef;
    }
}




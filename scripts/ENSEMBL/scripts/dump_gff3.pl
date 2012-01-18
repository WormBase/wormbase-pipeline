#!/usr/bin/perl
# based on the ZFish GFF3 VEGA dumper from Ian Sealy
# the Blast Filter is lifted from the WormBlast processing pipeline
#
# To lower the memory footprint and converts the rather large
# Ensembl objects into lighter hashes for later use


use strict;
use warnings;

use lib '../lib';
use lib '/software/worm/ensembl/ensembl/modules';
use lib $ENV{'CVS_DIR'}."/ENSEMBL/lib";

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use IO::File;
use DBI;
use Domain2Interpro;
use Getopt::Long;

my $debug =1;

# Vega database
my $vega_dbname = 'bmalayi2';
my $vega_dbhost = 'eagle';
my $vega_dbport =  3306;
my $vega_dbuser = 'wormadmin'; 
my $vega_dbpass = 'worms';
my @dump_slice ;
my $lsf;
my $dumpdir;
my $file;
my $slim;

GetOptions(
           'dbhost=s'     => \$vega_dbhost,
           'dbname=s'     => \$vega_dbname,
           'dbuser=s'     => \$vega_dbuser,
           'dbpass=s'     => \$vega_dbpass,
           'dbport=s'     => \$vega_dbport,
    	   'slice=s'      => \@dump_slice,
    	   'submit=i'     => \$lsf,
    	   'dumpdir=s'    => \$dumpdir,
    	   'file:s'       => \$file,
           'slim'         => \$slim,
          )or die ("Couldn't get options");


@dump_slice = split(/,/,join(',',@dump_slice));

my $vega_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-dbname  => $vega_dbname,
	-host    => $vega_dbhost,
	-port    => $vega_dbport,
	-user    => $vega_dbuser,
	-pass    => $vega_dbpass,
#	-species => 'vega',
);

my $mapper = Domain2Interpro->new();


#if single file passed print all to that rather than file per slice (chr or contig)
my  $pout;
if($file) {
    $pout = new IO::File (">$dumpdir/$file.gff3");
}

# Get all toplevel slices
my $sa = $vega_db->get_SliceAdaptor();
my $slices;
if ($dump_slice[0]){
    foreach (@dump_slice) {
	push @$slices, $sa->fetch_by_region('toplevel',$_);
    }
}
else {
   $slices = $sa->fetch_all('toplevel');
}
my $bsub_str = "bsub -R \"select[mem>3500] rusage[mem=3500]\" -o /dev/null -e /dev/null -J dump_$vega_dbname perl $0 -dbhost $vega_dbhost -dbname $vega_dbname -dbuser $vega_dbuser  -dbport $vega_dbport -dumpdir $dumpdir";
$bsub_str .= " -dbpass $vega_dbpass" if $vega_dbpass;
$bsub_str.=' -slice ';

my $slice_counter =0;
my $bsub_slices = "";
while( my $slice = shift @$slices) {
    my $slice_name = $slice->seq_region_name();
    my $slice_size = $slice->length;
    
    if ($lsf){
	$bsub_slices .= "$slice_name,";
	$slice_counter++;
	
	if($slice_counter == $lsf) {
	    print STDERR "$bsub_str $bsub_slices\n" if $debug;
	    print `$bsub_str $bsub_slices`;
	    $slice_counter =0;
	    $bsub_slices = "";
	}
    }else {
	
	print STDERR "Dumping $slice_name to $dumpdir/$slice_name.gff3 \n" if $debug;
	
	my $outf = $pout || new IO::File (">$dumpdir/$slice_name.gff3");
	
	print $outf "##gff-version 3\n";
	print $outf "##sequence-region $slice_name 1 $slice_size\n";
	
	# Get all the genes on this slice
	my $genes = $slice->get_all_Genes();
	while( my $gene=shift @$genes) {
	    my $gene_stable_id = 'gene:'.$gene->stable_id();
	    print STDERR " Dumping $gene_stable_id\n" if $debug;

	    my %gene_to_dump = (
				stable_id => $gene_stable_id,
				name      => $slice_name,
				start     => $gene->seq_region_start(),
				end       => $gene->seq_region_end(),
				strand    => $gene->strand(),
				note      => ($gene->status()||'PREDICTED' ). " " . $gene->biotype(),
				public_name => $gene->stable_id(),
				display     => $gene->stable_id(), # wrong but fixes db's without xref_mapping
				);

	    # get all transcripts of the gene
	    my $all_transcripts = $gene->get_all_Transcripts();
	    while( my $transcript = shift @{$all_transcripts}) {
		my $translation = $transcript->translation;

		my $transcript_stable_id = 'transcript:'.$transcript->stable_id();
		my $translation_id = 'cds:'.$transcript->translation->stable_id();

		print STDERR "  Dumping $transcript_stable_id\n" if $debug;
		print STDERR "  Dumping $translation_id\n" if $debug;

		push @{$gene_to_dump{'transcript'}}, {
		    stable_id => $transcript_stable_id,
		    translation_stable_id => $translation_id,
		    name      => $slice_name,
		    start     => $transcript->start(),
		    cds_start => $translation->genomic_start(),
		    cds_end   => $translation->genomic_end(),
		    end       => $transcript->end(),
		    strand    => $transcript->strand(),
		    note      => $transcript->biotype(),
		    info      => get_info($transcript),
		    public_name => $transcript->stable_id(),
		    display     => $transcript->stable_id(),
		};

		# get all exons and translatable_exons of the transcript
		my $all_exons=$transcript->get_all_Exons();
		my $all_t_exons = $transcript->get_all_translateable_Exons();

		while( my $exon = shift @{$all_exons}) {
		    my $exon_stable_id = 'exon:'.$exon->stable_id();
		    print STDERR "   Dumping $exon_stable_id\n" if $debug;
		    push @{${${$gene_to_dump{'transcript'}}[-1]}{'exon'}}, {
			stable_id => $exon_stable_id,
			name      => $slice_name,
			start     => $exon->seq_region_start(),
			end       => $exon->seq_region_end(),
			strand    => $exon->strand(),
		    }
		}
		
		$all_exons=undef;

		while (my $cds = shift @{$all_t_exons}) {
		    push @{${${$gene_to_dump{'transcript'}}[-1]}{'cds'}}, {
			stable_id => $translation_id,
			name      => $slice_name,
			start     => $cds->seq_region_start(),
			end       => $cds->seq_region_end(),
			strand    => $cds->strand(),
			phase     => (3-$cds->phase())%3, # phase/frame conversion to a sane system
		    }
		}
	    }
	    print $outf dump_gene(\%gene_to_dump);
	}
        
        next if $slim;
	# get all protein align features on the slice
        my %blastx_features;

        my $sth = $vega_db->prepare('SELECT seq_region_start,seq_region_end,seq_region_strand,evalue,hit_name,cigar_line,hit_start,hit_end,seq_region_strand,score,protein_align_feature_id,logic_name FROM protein_align_feature JOIN analysis USING(analysis_id) WHERE seq_region_id=?');
        my ($seq_region_start,$seq_region_end,$seq_region_strand,$e_value,$hit_name,$cigar,$hstart,$hend,$strand,$score,$dbID,$logic_name);
    	$sth->execute($slice->get_seq_region_id);

    	$sth->bind_col(1,\$seq_region_start);
	$sth->bind_col(2,\$seq_region_end);
	$sth->bind_col(3,\$seq_region_strand);
    	$sth->bind_col(4,\$e_value);
    	$sth->bind_col(5,\$hit_name);
    	$sth->bind_col(6,\$cigar);
	$sth->bind_col(7,\$hstart);
	$sth->bind_col(8,\$hend);
	$sth->bind_col(9,\$strand);
	$sth->bind_col(10,\$score);
	$sth->bind_col(11,\$dbID);
	$sth->bind_col(12,\$logic_name);

        while($sth->fetch){
	    my $cigar_line = flipCigarReference($cigar); # for Lincoln
	    my $stripped_feature = {
		hit_id      => $hit_name,
		target_id   => $slice->seq_region_name,
		target_start=> $hstart,
		target_stop => $hend,
		strand      => ($strand > 0?'+':'-'),
		hit_start   => $seq_region_start,
		hit_stop    => $seq_region_end,
		score       => $score,
		p_value     => $e_value,
		dbid        => $dbID,
		logic_name  => $logic_name,
		cigar       => ($strand > 0 ? $cigar_line : reverse_cigar($cigar_line)),
		feature_type=> 'protein_match',
	    };
	    push @{$blastx_features{$stripped_feature->{logic_name}}}, $stripped_feature;
	}

	while (my($k,$v)=each %blastx_features){
	    my @filtered_features=filter_features($v,$slice_size);
	    map {print $outf dump_feature($_)} @filtered_features;
	}

        # get all dna align features on the slice
	
	my $features=$slice->get_all_DnaAlignFeatures();
	print STDERR "dumping ${\scalar @$features} DNA Align Features\n" if $debug;
        foreach my $feature(@$features){
	    my $cigar_line = flipCigarReference($feature->cigar_string); # for Lincoln
	    my $stripped_feature = {
		hit_id      => $feature->hseqname,
		target_id   => $feature->slice->seq_region_name,
		target_start=> $feature->hstart,
		target_stop => $feature->hend,
		strand      => ($feature->strand > 0?'+':'-'),
		hit_start   => $feature->seq_region_start,
		hit_stop    => $feature->seq_region_end,
		score       => $feature->score,
		p_value     => $feature->p_value,
		dbid        => $feature->dbID,
		logic_name  => $feature->analysis->logic_name,
		cigar       => ($feature->strand > 0 ? $cigar_line : reverse_cigar($cigar_line)),
		feature_type=> 'nucleotide_match',
	    };
	    print $outf dump_feature($stripped_feature);
	}
	

	# get all repeat features on the slice
        my $repeats = $slice->get_all_RepeatFeatures;
        foreach my $feature (@$repeats){
	    my $stripped_feature = {
		target_id   => $feature->slice->seq_region_name,
		strand      => ($feature->strand > 0?'+':'-'),
		hit_start   => $feature->seq_region_start,
		hit_stop    => $feature->seq_region_end,
		score       => ($feature->score||'.'),
		dbid        => $feature->dbID,
		logic_name  => $feature->analysis->logic_name,
		feature_type=> 'repeat_region',
	    };
	    print $outf dump_feature($stripped_feature);
	}

	print $outf '#'x80;
	print $outf "\n";
	$outf->close unless $file;
    }
}
$pout->close if $file;

#if( $lsf) {
#    print STDERR "$bsub_str $bsub_slices\n" if $debug;
#    print `$bsub_str $bsub_slices`;
#}

# CIGAR to old GFF3 CIGAR format converter
sub cigar_to_almost_cigar {
    my $i=shift;
    $i=~s/(\d*)(\w)/$2$1 /g;
    return $i;
}

#convert refernece strand for Lincoln
sub flipCigarReference {
   my $i=shift;
   $i=~tr/ID/DI/;
   return $i;
}

# reverse cigar string
sub reverse_cigar {
    my $i=shift;
    my @pairs=$i=~/(\d*[MIDFR])/g;
    my $reversed_cigar = join('',reverse @pairs);
    return $i;
}

# print the feature using some funky template
sub dump_feature {
    my $i=shift;
    my %feature=%{$i};
    my $gff_line=
	"$feature{target_id}\t$feature{logic_name}\t$feature{feature_type}\t$feature{hit_start}\t$feature{hit_stop}\t".
	"$feature{score}\t$feature{strand}\t.".
	"\tID=$feature{logic_name}.$feature{dbid}".
	($feature{cigar}?";Name=$feature{hit_id};Target=$feature{hit_id} $feature{target_start} $feature{target_stop};Gap=$feature{cigar}\n":"\n");
    return $gff_line;
}


# build the info tag including protein features and interpro
sub get_info {
    my $transcript= shift;
    my $info='';
    # get all protein_features on the transcript
    my $features=$transcript->translation->get_all_ProteinFeatures();
    # get logic_name and hit_id
    
    my %plain_features;
    my $rest_features;
    map {
	if ($mapper->get_method2database($_->analysis->logic_name())) {
	    push @{$plain_features{$_->analysis->logic_name()}},
	        [$_->display_id(),$_->start(), $_->end(),$_->hstart(),$_->hend(),$_->score(),$_->p_value()]
	}
	else {
	    push @$rest_features,$_
	}
    } @$features;
    my @interpros=$mapper->get_mapping(\%plain_features);
    map {$info.=sprintf( "position:%d-%d method:%s accession:%s description:%s %%0A", $_->[1], $_->[2],
			 'InterPro', $_->[0] , $_->[7]) if $_->[1]} @interpros;

    while ( my $pfeature = shift @$rest_features ) {
	my $logic_name = $pfeature->analysis()->logic_name();
         $info.=sprintf( "position:%d-%d %s method:%s accession:%s %%0A", $pfeature->start(), $pfeature->end(), 
		 $pfeature->p_value(),$logic_name, $pfeature->id());
	}
	return $info;
}

# print the gene including transcripts and exons
sub dump_gene {
	my ($gene) = @_;
	
	my $output = '';
	
	# Dump gene
	$output .= "# Gene " . $gene->{'stable_id'} . "\n";
	$output .= gff_line(
		$gene->{'name'}, 'gene', $gene->{'start'}, $gene->{'end'},
		$gene->{'strand'}, $gene->{'stable_id'},undef, $gene->{'display'}, $gene->{'note'},undef,$gene->{'public_name'}
	);
	
	# Dump transcripts
	my $parent = $gene->{'stable_id'};
	my %exon_parent;
	foreach my $transcript (@{$gene->{'transcript'}}) {
		$output .= gff_line(
			$transcript->{'name'}, 'mRNA', $transcript->{'start'}, $transcript->{'end'},
			$transcript->{'strand'}, $transcript->{'stable_id'}, undef ,$transcript->{'display'} || undef, undef, 
			($transcript->{info}||undef),$transcript->{'public_name'},$parent);

		# Store the parent of this transcript's exons
		foreach my $exon (@{$transcript->{'exon'}}) {
			${$exon_parent{$exon->{'stable_id'}}}{$transcript->{'stable_id'}} = 1;
		}
	}
	
	# Dump exons
	foreach my $transcript (@{$gene->{'transcript'}}) {
		foreach my $exon (@{$transcript->{'exon'}}) {
			next if !$exon_parent{$exon->{'stable_id'}}; # If there are no parents then we've already dumped this exon
			my @parents = keys %{$exon_parent{$exon->{'stable_id'}}};
			delete $exon_parent{$exon->{'stable_id'}};
			$output .= gff_line(
				$exon->{'name'}, 'exon', $exon->{'start'}, $exon->{'end'},
				$exon->{'strand'}, $exon->{'stable_id'}, undef,undef, undef,undef,undef ,@parents);
		}
	}
	
	# Dump coding_exons
	foreach my $transcript (@{$gene->{'transcript'}}) {
		my $parent = $transcript->{'stable_id'};
		foreach my $cds (@{$transcript->{'cds'}}) {
			$output .= gff_line(
				$cds->{'name'}, 'CDS', $cds->{'start'}, $cds->{'end'},
				$cds->{'strand'}, $cds->{'stable_id'}, $cds->{'phase'},undef, 
				undef,undef ,undef, $transcript->{'stable_id'});
		}
	}
	
	$output .= "###\n";
	
	return $output;
}

# a template for a GFF line
sub gff_line {
	my ($seqid, $type, $start, $end, $strand, $stable_id, $phase,$name, $note, $info,$public,@parents) = @_;
	
	my $output = '';

	$phase='.' unless defined $phase;
        $strand= $strand>0? '+' : '-';

	$output .= "$seqid\tWormBase\t$type\t$start\t$end\t.\t$strand\t$phase\t";
	my @tags;
	push @tags, "ID=$stable_id" if defined $stable_id;
	if (@parents) {
		my $parent = join(',', @parents);
		push @tags, "Parent=$parent";
	}
	push @tags, "Name=$name" if defined $name;
	push @tags, "Note=$note" if defined $note;
	push @tags, "info=$info" if defined $info;
	push @tags, "public_name=$public" if defined $public;
	$output .= join(';', @tags);
	$output .= "\n";
	
	return $output;
}

# remove < 75% of evalue features from 100bp windows
sub filter_features {
	my ($features,$length)=@_;
        my %f_features;
	
	# 50bp bin size
	my $size=100;

	# should I bin them instead?
	my @bins;

	# put the feature in all bins between start and stop
	foreach my $f(@$features){
		my ($_start,$_end)=sort {$a <=> $b } ($f->{hit_start},$f->{hit_stop});
		for (my $n=int($_start/$size);$n <=int($_end/$size);$n++){
			push @{$bins[$n]},$f;
		}
	}

	# get the best 5 within 25% of the best hsp and add them to a hash
	foreach my $bin(@bins) {
		next unless $bin; # skip empty bins
		my $best=0;
		my $max_hsp=0;
		my @sorted_bin = sort {$a->{p_value} <=> $b->{p_value}} @$bin;
		$best=&p_value($sorted_bin[0]->{p_value});
		map { $f_features{$_->{dbid}}=$_ if (&p_value($_->{p_value}) > $best*0.75 && $max_hsp++ <5)} @$bin; # <= cutoff place
	}
	
	# flatten hash into an array
	my @_filtered=values %f_features;
	return @_filtered;
}

# p_value shorthand
sub p_value {
        my ($p)=@_;
        my $log = (eval($p) > 0 ? -(log(eval($p)))/log(10) : 999.9); # if e-value is zero set it to 999.9
        return $log;
}




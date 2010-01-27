#!/software/bin/perl -w
use strict;                              
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/ensembl/bioperl-live';
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Path;
use Bio::SeqIO;

my($alignments,$single,$merge, $transcripts, $output, $readsdir, $reads);

GetOptions ("alignments:s" => \$alignments,
	    "single"       => \$single,
	    "merge"        => \$merge,
	    "transcripts:s"=> \$transcripts,
	    "output:s"     => \$output,
	    "readsdir:s"   => \$readsdir,
	    "reads:s"      => \$reads
	    );

die "choose single | merge\n" unless ($single or $merge);
die "output required\n" if ($merge and !$output);

if ($single) {
    my($filename, $directories, $suffix) = fileparse($alignments);
    my $dat = "$directories/$filename.dat";
    &single($alignments, $dat);
}

my %store;
my %transcripts;
if($merge) {
    #merge each batch of alignments.
    my @files = glob ("$readsdir/chunks/exonerate_$reads.fa_*.fa/$reads.fa_*.fa_exonerate.dat");
    foreach (@files){
	&merge_hashes(\%store, $_);
    }
    #read in transcript sequences from file mapped against
    &get_sequences;

    my $mapped_reads;
    #work out how many reads mapped.
    foreach  my $gene (keys %transcripts) {
	$mapped_reads += $store{$gene}{'count'} if $store{$gene} ;
    }

    #report each gene
    open (OUT,">$output") or die "cant write to $output : $!\n";
    foreach my $gene (keys %transcripts) {
	if($store{$gene}) {
	    my $rpkm = sprintf("%.2f",$store{$gene}{'count'} / (((length $transcripts{$gene}))/1000) / ($mapped_reads/1000000));
	    my $coverage = sprintf("%.2f",($store{$gene}{'end'} - $store{$gene}{'start'}) / length $transcripts{$gene});
	    print OUT "$gene\treads=".$store{$gene}{'count'}."\tcoverage=".$store{$gene}{'start'}."-".$store{$gene}{'end'}." ($coverage)\tRPKM=$rpkm\n";
	}
	else {
	    print OUT "$gene\tno reads\n";
	}
    }
    close OUT;
}
sub merge_hashes {
    my ($store, $add) = @_;
    open( FH, "<$add" ) or die "can't open $add\t:$!";
    undef $/;
    my $VAR1;
    my $data = <FH>;
    eval $data;
    if($@) {
	warn "$add failed hash merge\n";
	return;
    }
    $/ = "\n";
    close FH;
    my $keycount = scalar keys %$VAR1;
    warn "retrieval through FetchData failed - dat file is empty\n" if $keycount == 0;

    foreach my $gene (keys %$VAR1) {
	foreach my $match ( @{$$VAR1{$gene}} ){
	    my($start, $end) = @{$match};
	    if($$store{$gene}) {
		if($start < $$store{$gene}{'start'}) {
		    $$store{$gene}{'start'} = $start;
		}
		if($end > $$store{$gene}{'end'}) {
		    $$store{$gene}{'end'} = $end;
		}
	    }else {
		$$store{$gene}{'start'} = $start;
		$$store{$gene}{'end'} = $end;
	    }
	    $$store{$gene}{'count'}++;
	}
    }
}


sub single {
    my $alignments = shift;
    my $dat = shift;

    my %store;
    open(VUL,"grep vulgar $alignments | ");
    while(<VUL>) {
	#vulgar: 1_147_363_275_2783276 0 35 + B0001.8.1|RefSeq 48 13 - 175 M 35 35
	my @data = split;
	my $read = $data[1];
	my $gene = $data[5];
	my ($start,$end)= &order_ends($data[6], $data[7]);
	
	$gene =~ s/\|.*//;

	push(@{$store{$gene}},[$start,$end]);
    }
    close VUL;    

    open (DUMP,">$dat");
    print DUMP Data::Dumper->Dump([\%store]);
    close DUMP;
    return 0;
}


sub order_ends {
    my ($start,$end) = @_;
    if($end < $start) {
	my $x = $start;
	$start = $end;
	$end = $x;
    }
    return($start,$end);
}


sub get_sequences {
    my $seqs = Bio::SeqIO->new('-file' => $transcripts, '-format' => 'fasta');
    while( my $seq = $seqs->next_seq){
	$transcripts{$seq->primary_id} = $seq->seq;
    }
}

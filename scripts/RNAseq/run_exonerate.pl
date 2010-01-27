#!/software/bin/perl -w
# Take RNAseq data and create clusters with
use strict;                                      
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Log_files;
use File::Basename;
use File::Path;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

my ($help,$index,$make_index,$genome,$reads,$readsdir,$outdir, $paired, $chunk);
my ($prepare_reads, $align, $process);

GetOptions ("help"         => \$help,
	    "genome:s"     => \$genome,
	    "reads:s" => \$reads,
	    "readsdir:s"   => \$readsdir,
	    "chunk:i"      => \$chunk,    #if set chunk reads in to this many per file.
	    "outdir:s"     => \$outdir,
	    "prep_reads"   => \$prepare_reads,
	    "align"        => \$align,
	    "process"      => \$process
	    ) or die;

my $log = Log_files->make_log("/tmp/exonerate.$$", 'ar2');


# output to subdir of sequences unless specified.
$outdir = "$readsdir" unless ($outdir);

my $lsf = LSF::JobManager->new();
#prepare reads from scratch (1 fastq file per lane)
if($prepare_reads) {
    #merge fastqs
    mkdir("$readsdir/chunks");
  #  system("zcat $readsdir/*txt.gz >> $readsdir/$reads.fastq"); 

    #chunk to 2million
    my $chunksdir = "$readsdir/chunks";
    &chunk_reads($readsdir, "$reads.fastq", 2000000, $chunksdir);
    $readsdir = $chunksdir;

    #convert to fasta
    &fastq_to_fasta($readsdir);
}


&align if $align;
&process($readsdir) if $process;

$log->mail;

sub process {
    my $readsdir = shift;
    #cat vulgar lines to one file.
    my $vulgar_file = "$readsdir/${reads}2cds.exonerate.vulgar";
    system("bsub \"grep -h $readsdir/chunks/exonerate*/*exonerate >> $vulgar_file\"");

    my $rpkm_script = $ENV{'CVS_DIR'}."/RNAseq/calc_RPKM.pl";
    my $err = "$readsdir/process.err";
    my $out = "$readsdir/process.out";

    #convert to hashes
    my @files = glob("$readsdir/chunks/exonerate*/*exonerate");
    foreach my $file (@files) {
	my @bopts = ("-J => $file", "-o => $out", "-e => $err");
	$lsf->submit( @bopts, "perl $rpkm_script -single -alignments $file -output");
    }
    $lsf->wait_all_children(history => 1);

    #merge hashes and calc RKPM
    my @bopts = ("-J RKPM", "-o $out", "-e $err");
    $lsf->submit( @bopts, "perl $rpkm_script -merge -readsdir $readsdir -reads $reads -output $readsdir/exonerate2cds.out -transcripts $genome");    
}

sub align {
#align chunks
    die "align to what!? -genome <target> \n" unless ($genome and -e $genome);
    my $command = "/software/vertres/bin/exonerate-2.2.0 -Q DNA -T DNA -s 174 ";

#single end reads 
    opendir(my $dh, $readsdir) || die "can't opendir $readsdir: $!";
    my @reads; 
    while (my $f = readdir($dh)){
	push (@reads,$f) if ($f =~ /\S+\.fa$/ and !(-d "$readsdir/$f"));
    }
    closedir $dh;

    foreach my $read (@reads) {
	my $output_dir = "$readsdir/exonerate_$read";
	mkdir $output_dir unless -e $output_dir;

	my @bsub_opts = ("-J $read", "-o => $output_dir/exonerate.out", "-e => $output_dir/exonerate.err");
	my $bsub = "$command -q $readsdir/$read -t $genome  > $output_dir/${read}_exonerate";
	$lsf->submit(@bsub_opts,$bsub); 
    }
    $lsf->wait_all_children(history => 1);
}
sub fastq_to_fasta {
    my $chunksdir = shift;
    my $script = $ENV{'CVS_DIR'}."/RNAseq/fq_allstd.pl fq2fa ";
    
    opendir(my $dh, $chunksdir) || die "can't opendir $chunksdir: $!";
    my @fastq;
    while (my $f = readdir($dh)){
	push (@fastq,$f) if ($f =~ /^\S+\.fastq$/ and !(-d "$chunksdir/$f"));
    }
    closedir $dh;
    foreach my $fq (@fastq){
	my $fasta = $fq;
	$fasta =~ s/fastq/fa/g;
	system ("bsub -J fq2fa_$fasta \"perl $script $chunksdir/$fq > $chunksdir/$fasta\"");
    }
}




sub chunk_reads {
    my ($readsdir, $reads, $chunk, $chunksdir) = @_;
    if($paired) {
    }
    else {
	my @bsub_opts = (-J => "chunk_$reads", -e => "$readsdir/chunk.err", -o => "$readsdir/chunk.out");
	$lsf->submit(@bsub_opts,"/software/bin/perl ".$ENV{'CVS_DIR'}."/shatter.pl -file $readsdir/$reads -bin $chunk -output $chunksdir/$reads -format fastq"); 
	$lsf->wait_all_children(history => 1);
    }
}
    

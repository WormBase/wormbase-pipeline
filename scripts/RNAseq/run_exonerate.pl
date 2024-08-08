#!/software/bin/perl -w
# Take RNAseq data and create clusters with
use strict;                                      
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Log_files;
use File::Basename;
use File::Path;
use Modules::WormSlurm
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
    WormSlurm::submit_job_and_wait("grep -h $readsdir/chunks/exonerate*/*exonerate >> $vulgar_file", 'production', '500m', '00:30:00', '/dev/null', '/dev/null');

    my $rpkm_script = $ENV{'CVS_DIR'}."/RNAseq/calc_RPKM.pl";
    my $err = "$readsdir/process.err";
    my $out = "$readsdir/process.out";

    #convert to hashes
    my %slurm_jobs;
    my @files = glob("$readsdir/chunks/exonerate*/*exonerate");
    foreach my $file (@files) {
	my $cmd = "perl $rpkm_script -single -alignments $file -output";
	my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '500m', '1:00:00', $out, $err, $file);
	$slurm_jobs{$job_id} = $cmd;
    }
    WormSlurm::wait_for_all(keys %slurm_jobs);

    #merge hashes and calc RKPM
    my @bopts = ("-J RKPM", "-o $out", "-e $err");
    my $cmd = "perl $rpkm_script -merge -readsdir $readsdir -reads $reads -output $readsdir/exonerate2cds.out -transcripts $genome";
    WormSlurm::submit_job_and_wait($cmd, 'production', '500m', '1:00:00', $out, $err);
}

sub align {
#align chunks
    die "align to what!? -genome <target> \n" unless ($genome and -e $genome);
    my $base_command = "/software/vertres/bin/exonerate-2.2.0 -Q DNA -T DNA -s 174 ";

#single end reads 
    opendir(my $dh, $readsdir) || die "can't opendir $readsdir: $!";
    my @reads; 
    while (my $f = readdir($dh)){
	push (@reads,$f) if ($f =~ /\S+\.fa$/ and !(-d "$readsdir/$f"));
    }
    closedir $dh;

    my %slurm_jobs;
    foreach my $read (@reads) {
	my $output_dir = "$readsdir/exonerate_$read";
	mkdir $output_dir unless -e $output_dir;

	my $cmd = "$base_command -q $readsdir/$read -t $genome  > $output_dir/${read}_exonerate";
	my $job_id = WormSlurm::submit_job_with_name($cmd , 'production', '500m', '1:00:00', "$output_dir/exonerate.out", "$output_dir/exonerate.err", $read);
	$slurm_jobs{$job_id} = $cmd;
    }
    WormSlurm::wait_for_jobs(keys %slurm_jobs);
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
    my $slurm_jobs;
    foreach my $fq (@fastq){
	my $fasta = $fq;
	$fasta =~ s/fastq/fa/g;
	my $cmd = "perl $script $chunksdir/$fq > $chunksdir/$fasta";
	my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '500m', '1:00:00', '/dev/null', '/dev/null', "fq2fa_$fasta");
	$slurm_jobs{$job_id} = $cmd;
    }
    WormSlurm::wait_for_jobs(keys %slurm_jobs);
}




sub chunk_reads {
    my ($readsdir, $reads, $chunk, $chunksdir) = @_;
    if($paired) {
    }
    else {
	WormSlurm::submit_job_and_wait("perl " . $ENV{'CVS_DIR'} . "/shatter.pl -file $readsdir/$reads -bin $chunk -output $chunksdir/$reads -format fastq",
				       'production', '500m', '1:00:00', "$readsdir/chunk.out", "$readsdir/chunk.err");
    }
}
    

#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2010-08-27 14:50:53 $

use lib '/software/worm/ensembl/bioperl-live';
use lib '/software/worm/lib';
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Bio::SearchIO;
use LSF::JobManager;
use Ace;
use Coords_converter;
use Bio::Tools::Blat;

################################
# command-line options         #
################################

my $help;          # Help perldoc
my $test;          # Test mode
my $debug;         # Debug mode, verbose output to user running script
my $load;          # use to automatically load file to autoace
my $ace;           # specify where to put the output ace file
my $store;         # specify a frozen configuration file
my $chrom;         # specify a chromosome
my $species;       
my $gffdir;        # use a specific GFF_SPLIT directory
my $dbdir;         # connect to a different acedb
my $dist;          # just creates and distributes jobs;
my $seqs;          # sequence file to align
my $genome;        # genome sequence file to align to (or file of files)
my $db;

GetOptions(
	   "debug=s"      => \$debug,
	   "test"         => \$test,
	   "species:s"    => \$species,
	   "help"         => \$help,
	   "load"         => \$load,
	   "acefile=s"    => \$ace,
	   'store=s'      => \$store,
	   'dbdir=s'      => \$dbdir,
	   "database:s"   => \$db,
	   "distribute"   => \$dist,
	   "seqs:s"       => \$seqs,
	   "genome:s"     => \$genome,
	   );

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
############################
my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test, -organsim => $species) }

my $log = Log_files->make_build_log($wormbase);
###########################################

if($dist) {
    use Ace;
    my $lsf = LSF::JobManager->new();
    my $ace = Ace->connect('-path' => $db, '-program' => $wormbase->tace ) or $log->log_and_die(Ace->error."\n");
    my $RNAis = $ace->fetch_many('-class' => "RNAi");
    my $count = 0;
    my $file = "/lustre/cbi4/work1/wormpub/RNAi/seqs_$count";
    my $seq;
    my $chunk = 5000;
    open($seq,">$file");
    while(my $rnai = $RNAis->next){
	if($rnai->DNA_text) {
	    print $seq "\n>",$rnai->name,"\n";
	    print $seq $rnai->DNA_text,"\n";
	}
	if((++$count % $chunk) == 0) {
	    ($seq, $file) = &submit_job($seq,$file, $lsf, $count);
	}
    }
    &submit_job($seq,$file, $lsf, $count);
}
else {
    $log->log_and_die("must have -seq and -genome") unless ($seqs and $genome and (-e $seqs) and (-e $genome) );
    my $blat = '/software/worm/bin/blat/blat ';
    my $out = "$seqs.blat";
    $wormbase->run_command("$blat $genome $seqs  $out",$log);

    my $coords = Coords_converter->invoke($db, 0, $wormbase);

    my $file;
    open($file, "<$out") or $log->log_and_die("cant open blat output $out : $!\n");
    my $parser = Bio::Tools::Blat->new('-fh' => $file);
    while(my $result = $parser->next_result) {
	my @blocks = $result->get_SeqFeatures;
	my $clone;
	my @data;
	# if all blocks lie in same clone great- otherwise need to convert to SL/supercontig
	my ($target,$res_srt, $res_end) = $coords->LocateSpan($blocks[0]->seq_id, $result->start, $result->end);
	my $convert = $coords->isa_clone($target);
	
	foreach my $block (@blocks){
	    my $tseq   = $block->seq_id;
	    my $tstart = $block->feature1->start;
	    my $tend   = $block->feature1->end;
	    my $qseq   = $block->hseq_id;
	    my $qstart = $block->feature2->start;
	    my $qend   = $block->feature2->end;
	    my $score  = $block->score;

	    ($tseq, $tstart, $tend) = $coords->LocateSpan($tseq,$tstart,$tend) if $convert; #convert from chrom -> clone / SL
	    push(@data,[$tseq,$qseq,$score, $tstart, $tend, $qstart, $qend]);
	}
	my $match = &match(\@data);
	if($match) {
	    print "\nHomol_data : \"".$target.":RNAi_homol\"\n";
	    foreach my $hsp (@data) {
		print "RNAi_homol ".$$hsp[1]," $match ",join("\t",@$hsp[2..6])."\n";
	    }  
	}
    }
}

$log->mail;
exit;



sub match {
    #return "RNAi_primary";
    my $data = shift;
    my $longest = 0;
    my $total =0;
    foreach my $hsp (@{$data}){
	my $length = abs($hsp->[4] - $hsp->[3]);
	$longest = $length if $length >  $longest;
	$total += $length;
    }
    
    #work out if Primary or Secondary
    if($longest >= 95  and $total > 100)  { # RNAi_primary match
	return "RNAi_primary";
    }
    elsif($longest >= 80 and $total > 200){
	return "RNAi_secondary";
    }
    else {
	print STDERR "no good\n";
	return undef;
    }
}


sub submit_job {
    my($seq, $file, $lsfm, $count) = @_;
    close $seq;
    $wormbase->remove_blank_lines( $file, $log);
    my $job_name = "worm_".$wormbase->species."_rnai";
    my @bsub_options = (-J => $job_name);
    my $cmd = $wormbase->build_cmd("RNAi2Genome.pl -seqs $file -genome $db/genome_seq");
    $lsfm->submit(@bsub_options, $$cmd);

    $file = "/lustre/cbi4/work1/wormpub/RNAi/seqs_$count";
    open($seq,">$file");
    return($seq, $file);
}

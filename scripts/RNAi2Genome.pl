#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2009-03-19 16:56:56 $

use lib '/software/worm/ensembl/bioperl-live';
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
    my $chunk = 8000;
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
    my $blat = '/software/worm/bin/blat/blat -out=blast  ';
    my $out = "$seqs.blat";
    #$wormbase->run_command("$blat $genome $seqs  $out",$log);

    my $coords = Coords_converter->invoke($db, 0, $wormbase);

    my $file;
    open($file, "<$out") or $log->log_and_die("cant open blat output $out : $!\n");
    my $parser = new Bio::SearchIO('-format' => 'blast', '-fh' => $file);
    while(my $result = $parser->next_result) {
	while(my $hit = $result->next_hit()) { # enter code here for hit processing
	    # this creates concatenated string of homology (|||   |||||||||) 
	    # use this to determine if there is a 200 bp run of identical bases.
	    my $homology_string = "";
	    map($homology_string .= $_->homology_string,$hit->hsps);
	    my $match = match($homology_string);
	    if(defined $match){
		my $score = 100 *( $homology_string =~ tr/\|/\|/ / length($homology_string));
		foreach my $hsp (sort {$a->start <=> $b->start} $hit->hsps) {
		    my($clone, $start, $stop) = $coords->LocateSpan($hit->name, $hsp->hit->start, $hsp->hit->end);
		    print "\nHomol_data : \"$clone:RNAi\"\nRNAi_homol ".$hsp->seq_id." $match $score $start $stop ".$hsp->start."\t".$hsp->end."\n";
#		    print join("\t",($match,$hit->name, $hsp->seq_id, $hsp->hit->start,$hsp->hit->end, $hsp->frac_identical, $hsp->evalue, $hsp->score))."\n";
		}     
	    }
	}
    }
}

$log->mail;
exit;



sub match {
    my $string = shift;    
    if($string =~ /\|{95}/ and length($string) > 100)  { # RNAi_primary match
	return "RNAi_primary";
    }
    elsif($string =~ /\|{80}/ and length($string) > 200){
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
    my $cmd = $wormbase->build_cmd("RNAi2Genome.pl -seqs $file -genome $db/genome_seq");
    $lsfm->submit($cmd);

    $file = "/lustre/cbi4/work1/wormpub/RNAi/seqs_$count";
    open($seq,">$file");
    return($seq, $file);
}

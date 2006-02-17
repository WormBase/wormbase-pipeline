#!/usr/local/bin/perl5.8.0 -w

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use File::Path;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $chrom);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "chromosome:s"=> \$chrom
	  );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);


# do a single chromosome if prompted else do the lot......

my @chromosome;

if ($chrom) {
    push (@chromosome,$chrom);
}
else {
    @chromosome = qw( I II III IV V X );
}

# other vars
my $workdir = $wormbase->wormpub."/analysis/UTR";
my $gff_dir = $wormbase->gff_splits;

my ($file1,$file2,$file3,$file4);

my ($gene,$gene_start,$gene_stop,$line);
my ($gene_strand,$exon5_stop,$exon3_start);
my @f;
my $CDS;
my $Transcript;
my %transcripts2process;
my $isoform;
my $dotno;

#
#




#################################
# Main Loop for each chromosome #
#################################


foreach my $chrom ( @chromosome ) {

    open (OUTPUT, ">$workdir/CHROMOSOME_${chrom}_UTR.gff") or $log->log_and_die("cant open $workdir/CHROMOSOME_${chrom}_UTR.gff : $!\n");

    #################################
    # coding_transcript split files #
    #################################

    $wormbase->run_command("grep exon $gff_dir/CHROMOSOME_${chrom}_Coding_transcript.gff > $workdir/CHROMOSOME_${chrom}_coding_transcript_exon.gff", $log);
    $wormbase->run_command("grep coding_exon $gff_dir/CHROMOSOME_${chrom}_curated.gff    > $workdir/CHROMOSOME_${chrom}_coding_exon.gff",            $log);

    #######################
    # curated split files #
    #######################

    $file1 = "$workdir/test_$chrom.coding_exon_exon.gff";
    $file2 = "$workdir/test_$chrom.coding_transcript_exon.gff";

    #############################################
    # Loop through each gene in this chromosome #
    #############################################

    open (GENES, "<$gff_dir/CHROMOSOME_${chrom}_curated.gff") or $log->log_and_die("cant open $gff_dir/CHROMOSOME_${chrom}_curated.gff : $!\n");
    while (<GENES>) {
	next if (/^\#/);
	($gene) = (/CDS \"(\S+)\"/);
	
	next unless ($gene);
	
	print "// Parsing gene $gene\n\n" if ($verbose);
	
	#####################################
        # make the small files to work with #
	#####################################

	$wormbase->run_command ("grep -iw $gene $workdir/CHROMOSOME_${chrom}_coding_exon.gff              > $file1",$log);
	$wormbase->run_command ("grep -iw $gene $workdir/CHROMOSOME_${chrom}_coding_transcript_exon.gff   > $file2",$log);
	
	###################
        # CDS coordinates #
	###################

	$line = 1;
	open (CDS, "<$file1") or $log->log_and_die("cant open $file1 $!\n");
	while (<CDS>) {
	    chomp;
	    @f = split /\t/;
	    if ($line == 1) {
		$gene_start = $f[3];
		$exon5_stop = $f[4];
	    }
	    $gene_stop   = $f[4];
	    $exon3_start = $f[3];
	    $gene_strand = $f[6];
	    $line++;
	}
	close CDS;
	
	print "GENE $gene [$gene_start -> $gene_stop] [FIRST EXON stop $exon5_stop : LAST EXON start $exon3_start]\n" if ($verbose);
	
	#################
        # New UTR exons #
	#################	

	open (NEW, "gff_overlap -unsorted -not $file1 $file2 | ")  or $log->log_and_die ("cant open gff_overlap $!\n");
	while (<NEW>) {
	    @f = split /\t/;

	    # Assign $Transcript and $CDS
	    ($Transcript) = $f[8] =~ (/\"(\S+)\"/);
	    $CDS = $Transcript;
	    $dotno = $CDS =~ tr /\./\./;
	    if ($dotno > 1) {
		chop $CDS;
		chop $CDS;
	    }
	    $transcripts2process{$Transcript} = 1;   # Ensure that we write coding_exons for this one
	    
	    if ($f[4] < $gene_start) {
		if ($f[6] eq "+") {
		    print OUTPUT "$f[0]\tCoding_transcript\tfive_prime_UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
		else {
		    print OUTPUT "$f[0]\tCoding_transcript\tthree_prime_UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
	    }
	    elsif ($f[3] > $gene_stop) {
		if ($f[6] eq "+") {
		    print OUTPUT "$f[0]\tCoding_transcript\tthree_prime_UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
		else {
		    print OUTPUT "$f[0]\tCoding_transcript\tfive_prime_UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
	    }
	}
	close NEW;
	
	##########################
	# Partially coding exons #
	##########################

	open (NEW, "gff_overlap -unsorted -minfrac1 1 -quiet $file1 $file2 | ") or $log->log_and_die ("cant open gff_overlap $!\n");
	while (<NEW>) {
	    @f = split /\t/;

	    # Assign $Transcript and $CDS
	    ($Transcript) = $f[8] =~ (/\"(\S+)\"/);
	    $CDS = $Transcript;
	    $dotno = $CDS =~ tr /\./\./;
	    if ($dotno > 1) {
		chop $CDS;
		chop $CDS;

	    }
	    $transcripts2process{$Transcript} = 1;   # Ensure that we write coding_exons for this one

	    
	    if ( ($f[4] == $exon5_stop) && ($f[3] != $gene_start)) {

		$gene_start--;
		$gene_stop++;

		if ($f[6] eq "+") {
		    print OUTPUT "$f[0]\tCoding_transcript\tfive_prime_UTR\t$f[3]\t$gene_start\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
		else {
		    print OUTPUT "$f[0]\tCoding_transcript\tthree_prime_UTR\t$f[3]\t$gene_start\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
				
		$gene_start++;
		$gene_stop--;

	    }
	    elsif ( ($f[3] == $exon3_start) && ($f[4] != $gene_stop) ) {
		
		$gene_start--;
		$gene_stop++;

		if ($f[6] eq "+") {
		    print OUTPUT "$f[0]\tCoding_transcript\tthree_prime_UTR\t$gene_stop\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
		else {
		    print OUTPUT "$f[0]\tCoding_transcript\tfive_prime_UTR\t$gene_stop\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$Transcript\"\n";
		}
		
		$gene_start++;
		$gene_stop--;

	    }
	}
	close NEW;
	

	#########################################
	# Coding_exon lines for the transcripts #
	#########################################

	# CHROMOSOME_I    curated           coding_exon     11641   11689   .       +       0       CDS "Y74C9A.2"
	#
	# becomes 
	#
	# CHROMOSOME_I    Coding_transcript coding_exon     11641   11689   .       +       0       Transcript "Y74C9A.2" ; CDS "Y74C9A.2"

	foreach $isoform (keys %transcripts2process) {
	    
	    next if ($isoform eq "");
	    
	    open (CODINGEXONS, "<$file1") or $log->log_and_die("cant open $file1 $!\n");
	    while (<CODINGEXONS>) {
		chomp;
		@f = split /\t/;

		# Assign $CDS
		($CDS) = $f[8] =~ (/\"(\S+)\"/);

		print OUTPUT"$f[0]\tCoding_transcript\tcoding_exon\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tTranscript \"$isoform\" ; CDS \"$CDS\"\n";
	    }
	}
	close CODINGEXONS;
	


	############
	# clean up #
	############
	
	undef $gene;
	%transcripts2process = ();

	unlink $file1;
	unlink $file2;

	
    #_ end of gene line
    }
#_ end of chromosome   
}

exit(0);

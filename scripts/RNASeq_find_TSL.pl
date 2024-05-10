#!/usr/bin/env perl
#
# RNASeq_find_TSL.pl
# 
# by Gary Williams                        
#
# Script to find TSL sites from the unaligned short read data
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-06-17 09:41:26 $      

# TSL sequences from 
# PLOS Genetics
# Nov 2006
# Operon Conservation and the Evolution of trans-Splicing in the Phylum Nematoda
# David B. Guiliano, Mark L. Blaxter
# http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.0020198
# Figure 4


# To make the ace file to load into geneace, find the next free Feature number (e.g. WBsf123456) 
# concatenate all of the SRA/*/TSL/TSL_new_features.ace files 
# acezip.pl -file concatenated_file.ace 
# to make each feature definintion at a position unique in the file, then
# cat concatenated_file.ace | perl -ne 'BEGIN{$FEATURE = "WBsf919697" ; $count=0; %hash=()}{if (s/Feature\s+:\s+"(\S+)"/Feature : "$FEATURE"/){$hash{$1}=$FEATURE; $FEATURE++ ; $count++};if ($_=~/Feature_object\s+"(\S+)"/) {$o=$1;$r=$hash{$o};$o=~s/\+/\\\+/g;s/Feature_object\s+"$o"/Feature_object "$r"/}; print "$_"} END{print STDERR "Number of Features: $count\n"}' > ! file_to_load.ace



use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Coords_converter;
use Feature_mapper;
use Modules::PWM;
use Sequence_extract;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $experiment_accession);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "expt:s"     => \$experiment_accession, # specify a single experiment to do, default is to do all in this species
	   );


$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


######################################
# variables and command-line options # 
######################################

my $status;
my $data;
my $coords;
my $mapper;
my $pwm;
my $seq_obj;

#my $database = $wormbase->database('current');
my $database = $wormbase->autoace;

$coords = Coords_converter->invoke($database, 0, $wormbase);
$mapper = Feature_mapper->new($database, undef, $wormbase);
$pwm = PWM->new;
$seq_obj = Sequence_extract->invoke($database, 0, $wormbase);
my %chromosomes;


#OUT of sequence when trying to find flanking regions in OVOC_OM2_83, 564, 565
#    my $is_zero_length = 1;
#    my $min_flank_length = 30;
#    my $no_unique_check = 0;
#    my $short_end_flanks_allowed = 1;
#    my ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature('OVOC_OM2_83', 564, 565, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed); 


my $RNASeq = RNASeq->new($wormbase, $log, 0, 1); # not -new_genome, -check

$log->write_to("Get experiments from config\n");
my $data = $RNASeq->get_transcribed_long_experiments();

my @experiment_accessions = keys %{$data};
if (defined $experiment_accession) {@experiment_accessions = ($experiment_accession)} # just do one experiment if it is specified


# get the TSL sequences for this species and reduce them to just the last 8 bases
my %TSL = $wormbase->TSL;
my @TSL;
foreach my $TSL_type (values %TSL) {
  $TSL_type = substr($TSL_type, -8);
  push @TSL, $TSL_type;
}

foreach my $experiment_accession (@experiment_accessions) {
  print "Experiment: $experiment_accession\n";

  my $experiment = $data->{$experiment_accession};
    
  my $experiment_accession=$experiment->{'experiment_accession'};
  my $library_strategy=$experiment->{'library_strategy'};
  my $library_selection=$experiment->{'library_selection'};
  my $library_layout=$experiment->{'library_layout'};
  my $library_source=$experiment->{'library_source'};
  my $analysis=$experiment->{'analysis'};
  if (!defined $analysis) {next}
  $log->write_to("Running alignment job with experiment=$experiment_accession $library_source $library_strategy $library_selection $library_layout\n");
  
  my $RNASeqGenomeDir = $RNASeq->{RNASeqGenomeDir};
  my $RNASeqSRADir = $RNASeq->{RNASeqSRADir};
  my $TSL_Dir =  "$RNASeqSRADir/$experiment_accession/TSL";
  my $TSL_Dir_tmp =  "$RNASeqSRADir/$experiment_accession/TSL-tmp";
  my $genomeDir = $RNASeqGenomeDir."/STAR";
  my $srr_done_file = "$RNASeqSRADir/$experiment_accession/SRR/fastq.done";
  my $solid_done_file = "$RNASeqSRADir/$experiment_accession/SRR/solid.done";
  my $TSL_done_file = "$TSL_Dir/TSL.ace.done";
  my $tsl_fastq = "$TSL_Dir/sl.fq";
  
  if (-e $TSL_done_file) {
    $log->write_to("$experiment_accession has been analysed for TSL sites already - not repeating this\n");
    
  } else { # no TSL_done_file
    if (!-e $srr_done_file) { 
      $RNASeq->get_SRA_files($experiment_accession);
    }
    
    mkdir "$TSL_Dir_tmp", 0777;
    mkdir "$TSL_Dir", 0777;
    
    chdir "$RNASeqSRADir/$experiment_accession"; 

    if (!-e "$tsl_fastq") {
      my $outfile="TSL1";
      foreach my $file (glob "SRR/*.fastq") {
	$log->write_to("Running get_TSL_RNASeq_reads.pl on $file\n");
	$wormbase->run_script("get_TSL_RNASeq_reads.pl -infile $file -outfile $TSL_Dir_tmp/$outfile",$log);
	$outfile++;
      }
      
      # now concatenate all the TSL files to make one
      $wormbase->run_command("cat TSL-tmp/* > $tsl_fastq", $log);
      $wormbase->run_command("rm -rf TSL-tmp", $log);
    }

    chdir $TSL_Dir; 

    if (!-e "$TSL_Dir/Aligned.out.bam") {
      # useful descriptions of parameters are found in the program's STAR/parametersDefault file
      
      #/pathToStarDir/STAR --genomeDir /path/to/GenomeDir --readFilesIn /path/
      #to/read1 [/path/to/read2] --runThreadN <n> --<inputParameterName> <input
      #parameter value(s)>
      
      #my $options = "";
      # outReadsUnmapped Fastx   : output in separate fasta/fastq files, Unmapped.out.mate1/2
      # alignIntronMin 25
      # alignIntronMax maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins (500Kb)
      # alignMatesGapMax maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins (500Kb)
      # outFilterMultimapNmax 10 int: read alignments will be output only if the read maps fewer than this value, otherwise no alignments will be output (1=Only uniquely mapping reads)
      # outFilterMismatchNoverLmax 0.05 float: alignment will be output only if its ratio of mismatches to mapped length is less than this value
      # chimSegmentMin int>0: minimum length of chimeric segment length, if ==0, no chimeric output
      # chimJunctionOverhangMin int>0: minimum overhang for a chimeric junction
      # outSAMstrandField None string: Cufflinks-like strand field flag None : not used intronMotif : strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.
      
      # uses 3.2 Gb memory for single-end elegans - 4 threads takes 15 mins (tophat takes 2 hours)
      # uses 4.3 Gb memory for paired-end remanei - 4 threads takes 90 mins - 3670 contigs
      
      my $options = "--alignIntronMin 25 --outReadsUnmapped Fastx --alignIntronMax 15000 --alignMatesGapMax 15000 --outFilterMismatchNoverLmax 0.02";
      
      my $threads = 4;
      my $STARcmd = $RNASeq->{Software}."/star/STAR --genomeDir $genomeDir --readFilesIn $tsl_fastq $options --runThreadN $threads";
      $log->write_to("Running STAR command:\n$STARcmd\n\n");
      $status = $wormbase->run_command($STARcmd, $log);
      if ($status != 0) {$log->log_and_die("Didn't run STAR to do the alignment successfully for $experiment_accession\n");}
      
      $status = $RNASeq->sam_to_bam('Aligned.out.sam', 'Aligned.out.bam');
      if ($status != 0) {$log->log_and_die("sam_to_bam() failed for $experiment_accession\n");}
      $status = $RNASeq->sort_bam_file('Aligned.out.bam');
      if ($status != 0) {$log->log_and_die("sort_bam_file() failed for $experiment_accession\n");}
      if (-e glob("Aligned.out.bam.sorted*")) {$log->log_and_die("sort_bam_file() failed for $experiment_accession - a BAM sort file is still existant\n");}
      $status = $RNASeq->index_bam_file('Aligned.out.bam');
      if ($status != 0) {$log->log_and_die("index_bam_file() failed for $experiment_accession\n");}
    }

# This is a normal hit to an SL1 site
# samtools view Aligned.out.bam|ty
# ERR225732.4719707.+.SL1 0       OVOC_OM1a       80997   255     94M     *       0       0       AGTAAATTAAACTTTGGTTGCAATTTCTTCAAACGGTATATAAAAATTCAAGTTTTACTTTCAAAATTTCTTGTGGATTATTTTTGATAAATCG  IIIKIHKMEK>GMMLKLIKLJKLLNKLHJGMIMGJJEKJJLLLKLKJNKJNJGNCLMJKLGLKIJJJIKMKCIMLLHKKLKLLNJFIHEBBAA8  NH:i:1  HI:i:1  AS:i:92 nM:i:0

# typical example of alignment across an intron
#ERR225732.1301342.+.SL1 0       OVOC_OM2        10249326        255     81M220N10M      *       0       0       AGGAACGACGACAGCACGGCAACACATCAGTATGTTGTTGCAGCTGGTAATCAATCTACGTCAGCCTCGTCTTTCCCTCAGTATCAACCAT     IILEEMEBIEDGKCJL::DDK8IAJJFCJKGLIJCEIAKJJDHHKHKGANLLLGD6HFELBIFHFLFHGH=KMIIKKFHCFGHFGDCBB@3     NH:i:1  HI:i:1  AS:i:89 nM:i:0

# alignment in reverse sense - the pos field gives the pos of the 3' end of the short read, so we must add on the cigar length to find the 5' end
# so why is it 
# ERR225732.15761390.+.SL1        16      OVOC_OM2        310068  255     92M     *       0       0       AATTACATGGCATAGCAGGCATCGCATTTTGAACAATTTCTTGCCCCTCTCGTAGTCACCTTCCGGAACTGACATAACGATATGTAAGTTCT    2AABBEHHHNGGL,JKKDIHKIKKJAKJKJMDHELKMMLI@KKNMKJNJDKLDLJJL@JLJHKIMMGJLKLLLKLKPMKMMLGKKJIKLKII    NH:i:1  HI:i:1  AS:i:88 nM:i:1


    # now parse the results looking for TSL sites
    my $TSLfile = "Aligned.out.bam";
    my %results; # hash of TSL sites found by RNAseq alignment
    my %example; # hash of example read IDs for debugging
    $log->write_to("get the TSL sites found by alignment\n");

    my $samtools =  $RNASeq->{Software}."/samtools/samtools view $TSLfile";
    open (HITS, "$samtools |") || $log->log_and_die("can't run samtools\n");
    while (my $line = <HITS>) {
      my $sense = '+';
      my ($id, $flags, $chrom, $pos, $cigar) = ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\S+\s+(\S+)/);
      if ($id =~ /\-.SL\d$/) {next} # the reverse sense reads appear to all be spurious, so ignore them
      my ($tsl) = ($id =~ /(SL\d)/);
      my $cigar_length=0;
      if ($flags & 0x10) { # find and deal with reverse alignments
	$sense = '-';
	$cigar =~ s/^\d+s//; # strip digits-S from the start of the cigar as these are not counted as part of the positioning 
	my @cigar_lengths = split /\D+/, $cigar;
	map { $cigar_length += $_ } @cigar_lengths; # sum the elements of the array i.e. lengths of the things in the cigar string;
	$pos = $pos + $cigar_length - 2; # go past the 'A' of the 'AG' we added on to ensure we hit a splice site
      } else {
	$pos = $pos + 1; # go past the 'A' of the 'AG' we added on to ensure we hit a splice site
      }

      if (!defined $chrom) {print "undefined chrom in ID $id\n";}
      if (!defined $pos) {print "undefined pos in ID $id\n";}
      if (!defined $sense) {print "undefined sense in ID $id\n";}
      if (!defined $tsl) {print "undefined tsl in ID $id\n";}

      $chrom = $wormbase->chromosome_prefix . $chrom;

      $results{$chrom}{$pos}{$sense}{$tsl}++;
      $example{$chrom}{$pos}{$sense}{$tsl} = "$id $sense $cigar";
      #print "store: ${chrom} ${pos} ${sense} ${tsl} value: $results{$chrom}{$pos}{$sense}{$tsl} $seq\n";
    }
    close(HITS);
    
    # now check that the TSL sites match known TSL Feature objects and add evidence to the Features
    $log->write_to("check that the TSL sites match known TSL Feature objects\n");
    my %found; # hash of TSL Feature object IDs with evidence from RNASeq - holds no. of supporting reads
    my %not_found; # hash of TSL Feature object IDs with no evidence from this RNASeq experiment
    
    print "Using database $database to look up Features\n";
    my $table_def = &write_feature_def;
    my $table_query = $wormbase->table_maker_query($database, $table_def);
    while(<$table_query>) {
      chomp;
      s/\"//g;  #remove "
      next if (/acedb/ or /\/\//);
      my @data = split("\t",$_);
      my ($feature, $clone, $start, $end, $method) = @data;
      if (!defined $feature) {next}
      if (!defined $clone) {$log->write_to("Feature $feature does not have a position mapped\n"); next}
      if (!defined $start) {$log->write_to("Feature $feature does not have a position mapped\n"); next}
      my ($TSL_chrom, $TSL_coord);
      ($TSL_chrom, $TSL_coord) = $coords->Coords_2chrom_coords( $clone, $start );
      if ($start < $end) { # forward sense
	print "looking for: ${TSL_chrom} ${TSL_coord} + ${method}\n";
	
	if (exists $results{$TSL_chrom}{$TSL_coord}{'+'}{$method}) {
	  $found{$feature} = $results{$TSL_chrom}{$TSL_coord}{'+'}{$method};
	  delete $results{$TSL_chrom}{$TSL_coord}{'+'}{$method};
	  print "found $feature +\n";
	} else {
	  $not_found{$feature} = 1;
	  print "NOT found $feature +\n";
	}
      } else { # reverse sense
	print "looking for: ${TSL_chrom} ${TSL_coord} - ${method}\n";
	if (exists $results{$TSL_chrom}{$TSL_coord}{'-'}{$method}) {
	  $found{$feature} = $results{$TSL_chrom}{$TSL_coord}{'-'}{$method};
	  delete $results{$TSL_chrom}{$TSL_coord}{'-'}{$method};
	  print "found $feature -\n";
	} else {
	  $not_found{$feature} = 1;
	  print "NOT found $feature -\n";
	}
      }    
    }
    
    
    # write out evidence for matched existing Feature objects
    print "Writing evidence for existing Features\n";
    my $aceout = "TSL_evidence.ace";
    open (ACE, ">$aceout") || $log->log_and_die("Can't open ace file $aceout\n");
    foreach my $feature (keys %found) {
      print ACE "\n\nFeature : \"$feature\"\n";
      my $reads = $found{$feature};
      print ACE "Defined_by_analysis $analysis $reads\n";
      print ACE "Remark \"Defined by RNASeq data from $analysis with $reads reads\"\n";
    }
    close(ACE);
    

    # the only data left in %results now are the matches to the genome
    # where there is no defined TSL Feature object, so these are all new
    # TSL sites that we could define a new Feature for.
    
    my $total_matches = keys %results; # so we can work out the no. of alignments per million (or billion?) matches
    my $full_species = $wormbase->full_name;
    # find all RNASeq evidence for TSL sites that do not have Features
    my $newaceout = "TSL_new_features.ace";
    
    print "Start output\n";
    open (NEWACE, ">$newaceout") || $log->log_and_die("Can't open ace file $newaceout\n");
    $log->write_to("find all RNASeq evidence for TSL sites that do not have Features\n");
    foreach my $TSL_chrom (keys %results) {
      foreach my $TSL_coord (keys %{$results{$TSL_chrom}}) {
	foreach my $sense (keys %{$results{$TSL_chrom}{$TSL_coord}}) {
	  
	  # get the flanking sequences and write the Feature_object coords for this new Feature
	  #print "$TSL_chrom, $TSL_coord, $sense ";
	  my ($feat_clone, $clone_start, $clone_end, $left_flank, $right_flank) = get_feature_flanking_sequences($TSL_chrom, $TSL_coord, $TSL_coord+1, $sense);
	  if (defined $feat_clone) {
	    if (substr($left_flank, -2) =~ /AG/i) { # ensure it looks like a splice site
	      foreach my $method (keys %{$results{$TSL_chrom}{$TSL_coord}{$sense}}) {

		# check to see if there is a Sl1/SL2 trans-splice sequence in the genome just upstream
		my $left_flank_end = substr($left_flank, -8);
		if (!grep /$left_flank_end$/, @TSL) {

		  # check to see if this is a splice site with a score > 1 
		  my $chrom_seq;
		  if (!exists $chromosomes{$TSL_chrom}) {
		    $chrom_seq = $seq_obj->Sub_sequence($TSL_chrom);
		    $chromosomes{$TSL_chrom} = $chrom_seq; # store the chromsomal sequence for re-use
		  } else {
		    $chrom_seq = $chromosomes{$TSL_chrom};
		  }
		  my $splice_score = $pwm->splice3($chrom_seq, $TSL_coord-1, $sense); # -1 to convert from acedb sequence pos to perl string coords	
		  $log->write_to("splice score = $splice_score sense=$sense\n");
		  if ($splice_score >= 1) {
		  
		    if (!defined $results{$TSL_chrom}{$TSL_coord}{$sense}{$method}) {next} # ignore the deleted entries
		    my $example = $example{$TSL_chrom}{$TSL_coord}{$sense}{$method};
		    my $ft_id = "WBsf#${TSL_chrom}#${TSL_coord}#${sense}#${method}"; # totally bogus but unique Feature ID, needs to be changed before adding to Build
		    print NEWACE "\n\nFeature : \"$ft_id\"\n";
		    print NEWACE "Sequence \"$feat_clone\"\n";
		    print NEWACE "Mapping_target \"$feat_clone\"\n";
		    print NEWACE "Flanking_sequences \"$left_flank\" \"$right_flank\"\n";
		    print NEWACE "Species \"$full_species\"\n";
		    print NEWACE "Description \"$method trans-splice leader acceptor site\"\n";
		    print NEWACE "SO_term SO:0000706\n";
		    print NEWACE "Method $method\n";
		    my $reads =  $results{$TSL_chrom}{$TSL_coord}{$sense}{$method}; # reads
		    print NEWACE "Defined_by_analysis $analysis $reads\n";
		    print NEWACE "Remark \"Defined by RNASeq data (example read: $example) from $analysis with $reads reads\"\n";
		    print NEWACE "\n\n";
		    print NEWACE "Sequence : \"$feat_clone\"\n";
		    print NEWACE "Feature_object \"$ft_id\" $clone_start $clone_end\n\n";
		  } else {
		    $log->write_to("Ignoring site chromosome: $TSL_chrom, position: $TSL_coord sense: $sense because the splice score is too low ($splice_score)\n");
		    
		  }
		} else {
		  $log->write_to("Ignoring site chromosome: $TSL_chrom, position: $TSL_coord sense: $sense because it is just downstream of the last 8 bases of the standard SL1 and so it probably a spurious hit and not trans-spliced\n");
		}
	      }
	    } else {
	      #$log->write_to("Ignoring site chromosome: $TSL_chrom, position: $TSL_coord sense: $sense because it doesn't look like a splice site (".substr($left_flank, -2).")\n");
	    }
	  } else {
	    $log->write_to("Couldn't find flanking sequence for chromosome: $TSL_chrom, position: $TSL_coord sense: $sense\n");
	  }
	}
      }
    }
    close(NEWACE);
    print "End output\n";
    
    # show all Feature objects that do not have evidence from this RNASeq data-set
    #  $log->write_to("show all Feature objects that do not have evidence from this RNASeq data-set\n");
    #  foreach my $feature (keys %not_found) {
    #    $log->write_to("Not found: $feature\n");
    #  }
    
    $log->write_to("TSL stuff done\n");
    
    $wormbase->run_command("touch ${aceout}.done" ,$log); # set flag to indicate we have finished this
    
    $status = $wormbase->run_command("touch $TSL_done_file", $log); # set flag to indicate we have finished this


  } # no TSL_done_file
} # foreach my $experiment_accession

$log->mail();
print "Finished.\n" if ($verbose);

exit(0);


############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_feature_def {
  my $def = "/tmp/Features_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END;

Sortcolumn 1

Colonne 1
//Subtitle Feature
Width 20
Optional
Visible
Class
Class Feature
From 1
Condition ((Method = "SL1") OR (Method = "SL2")) AND (Species = "$species")

Colonne 2
//Subtitle Clone
Width 12
Mandatory
Visible
Class
Class Sequence
From 1
Tag Sequence

Colonne 3
//Subtitle Feature2
Width 12
Optional
Hidden
Class
Class Feature
From 2
Tag Feature_object
Condition IS \\%1

Colonne 4
//Subtitle Start_pos
Width 12
Optional
Visible
Integer
Right_of 3
Tag HERE

Colonne 5
//Subtitle End_pos
Width 12
Optional
Visible
Integer
Right_of 4
Tag HERE

Colonne 6
//Subtitle Type
Width 12
Optional
Visible
Class
Class Method
From 1
Tag Method

END

  print TMP $txt;
  close TMP;
  return $def;
}

############################################################################
# get the clone coords and the flanking sequences


sub get_feature_flanking_sequences {
  my ($chromosome, $start, $end, $sense) = @_;

# get the chromosome location of the Feature
# map the Feature onto the chromosome
# get the flanking sequence start & end and find the minimal sequence object that will allow these positions to map to it (LocateSpan) and return this

  # returned data
  my ($left_flank, $right_flank);

  my @clone_coords = $coords->LocateSpan($chromosome, $start, $end ); 
  my ($target_clone, $clone_start, $clone_end) = @clone_coords;
    
  # use reverse sense coords
  if ($sense eq '-') {
    ($clone_end, $clone_start) = ($clone_start, $clone_end);
    $clone_end--;
    $clone_start--;
  }

  # get the left and right flanking sequences for this position
  # fourth argument says this is a zero-length feature
  # seventh argument allows non-unique short sequences at the contig sequence ends
  my $is_zero_length = 1;
  my $min_flank_length = 30;
  my $no_unique_check = 0;
  my $short_end_flanks_allowed = 0; # don't allow truncated flanks at the end of the clone
  ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($target_clone, $clone_start, $clone_end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed);

  if (defined $left_flank && defined $right_flank) {
    return ($target_clone, $clone_start, $clone_end, $left_flank, $right_flank)
  }
   


  # if the search for unique flanks in the clone fails, search in the whole chromosome
  # this can be quite slow, which is why we look for possible unique flanks in the clone first

  # get the left and right flanking sequences for this position
  # fourth argument says this is a zero-length feature
  # seventh argument allows non-unique short sequences at the contig sequence ends
  $is_zero_length = 1;
  $min_flank_length = 30;
  $no_unique_check = 0;
  $short_end_flanks_allowed = 1; # allow truncated flanks at the end of the chromosoem/contig
  ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($chromosome, $start, $end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed); 
  
  if (! defined $left_flank) {
    $log->write_to( "ERROR when trying to find flanking regions at $chromosome, $start, $end\n");
    return (undef, undef, undef);
  }

  return ($chromosome, $start, $end, $left_flank, $right_flank)
    
}

#!/usr/local/bin/perl

use strict;
use Carp;
$| = 1;

##################
# variables etc. #
##################

my $dir      = "/wormsrv1/newcamace2";
my @files    = glob("$dir/dump*.ace");
#my @files    = "$dir/lala.ace"; # for debug
my $tace     = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace /wormsrv1/camace";
my (%feature,%CDS,%hist,%pseudo,%RNA,%NDB,@lines,$name);

# what we don't want anymore
my @notwantedclass   = qw(KeySet LongText Locus Transgene 2_point_data Allele Anatomy Author Genetic_map gMap GO_term Grid Journal Lineage Cell Other_locus Pos_neg_data Rearrangement Gene_class Strain Expr_pattern);
my @notwantedtagsCDS = qw(Locus_genomic_seq Annotated Allele misc_feature Keyword Paired_read EMBL_feature Subsequence Status Map Matching_Genomic Spliced_cDNA Feature Title DB_annotation Homol In_alignment Clone);
my @notwantedtagsRNA = (@notwantedtagsCDS,'CDS_predicted_by'); 


########
# MAIN #
########


#######################################
# get a list of genes to be converted #
#######################################

&make_lists();

###################
# read dump files # 
###################

my $sep     = $/; # set record separator to slurp in chapters
$/          = "";

foreach my $file (@files) {

	my $newfilename = $file.".new";
	open (F,$file);
	open (OUT,">$newfilename");

	my ($skipclass,$class);
	
	OUTER:while (<F>) {
					
		###################
		# general changes #
		###################

		s/From_Laboratory/From_laboratory/g;
		s/Source_Exons/Source_exons/g;
		s/Main_Locus/Main_gene/g;

		#################################
		# check whether class is wanted #
		#################################
		
		if ($skipclass == 1) {
			next OUTER unless (/\/\/\s+Class/);
		} 
		if (/\/\/\s+Class\s+(\S+)/) {
			$class = $1;
			$skipclass = 0;
			foreach my $try (@notwantedclass) {
				if ($class =~ /^$try$/i) {
					$skipclass = 1;
				}
				else {
				}
			}
			if ($skipclass == 1) { 
				next OUTER;
			}
			else {
				next OUTER;
			}
		}		

		##############################
		# classes with minor changes #
		##############################
		
		if ($class =~ /^Accession_number$/i) {
			s/Sequence/Primary_for_sequence/g;
			print OUT $_;
			next OUTER;
		}
		if ($class =~ /Transposon/i) {
			s/Canonical_parent/Sequence/g;
			print OUT $_;
			next OUTER;
		}

		##############################################
		# major changes: convert objects and entries #
		##############################################

		elsif ($class =~ /Sequence/i) {
			
			@lines = split /\n/;

			if (/Sequence\s+:\s+\"(\S+)\"/) {

				$name = $1;
				# check whether Sequence is CDS or history object
				if (exists $CDS{$name} || exists $hist{$name} || exists $pseudo{$name}) {
					&convert_object('CDS');
					next OUTER;
				}	

				# check whether Sequence is RNA = Transcript object
				elsif (exists $RNA{$name}) {
					&convert_object('RNA');
					next OUTER;
				}	

				# ognore NDB_CDS objects
				elsif (exists $NDB{$name}) {
					next OUTER;
				}
				
				# all other Sequence objects
				else {
					foreach my $line (@lines) {
						if ($line =~ /^Matching_Genomic/i) {
							my ($entry) = ($line =~ /Matching_Genomic\s+\"(\S+)\"/);
							if (exists $CDS{$entry} || exists $hist{$entry} || exists $pseudo{$entry}) {
								$line =~ s/Matching_Genomic/Matching_CDS/g;
								print OUT "$line\n";
							}
							elsif (exists $RNA{$entry}) {
								$line =~ s/Matching_Genomic/Matching_transcript/g;
								print OUT "$line\n";
							}
							else {
								print STDERR "Could not relocate line $line of $name in file $file\n";
							}
						}
						# unwanted tags including EMBL_features
						elsif ($line !~ /^Allele|^Brief_identification|^Locus_other_seq|^Locus_genomic_seq|^EMBL_feature|^Matching_Genomic|^repeat_region|^CAAT_signal|^GC_signal|^TATA_signal|^allele_seq|^conflict|^mat_peptide|^misc_binding|^misc_feature|^misc_signal|^misc_recomb|^modified_base|^old_sequence|^polyA_signal|^polyA_site|^prim_binding|^prim_transcript|^promoter|^repeat_region|^repeat_unit|^satellite|^sig_peptide|^variation|^enhancer|^protein_bind|^stem_loop|^primer_bind|^transit_peptide|^misc_structure|^precursor_RNA|^LTR/) { 
							print OUT "$line\n";
						}
					}
					print OUT "\n";
					next OUTER;
				}		
			}
		}

		####################
		# whatever is left #
		####################
		
		else {
			print OUT $_;
		}
	}
}
	
##############################	
# connect CDS etc to parents #
##############################

# get subseq coordinates
$/ = $sep;
my %data;

my $command1 =<<EOF;
Table-maker -p "/wormsrv1/camace/wquery/subseq_coords.def"
quit
EOF

open (TACEI, "echo '$command1' | $tace | ");
while (<TACEI>) {
	chomp;
	next if ($_ eq "");
	next if (/\/\//);
	s/acedb\>\s//g;
	s/\"//g;
	if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
		my ($cosmid,$name,$start,$end) = ($1,$2,$3,$4); 
		$data{$name} = [$cosmid,$start,$end];
	}
}

foreach my $gene (sort keys %feature) {

	################################
	# generate parent SMap entries #
	################################
	 
	if ($feature{$gene}->{'source'} eq $data{$gene}->[0]) {
		print OUT "\nSequence \"$feature{$gene}->{'source'}\"\n";
	}
	else {
		print STDERR "Problem with $gene: ".$feature{$gene}->{'source'}." ne ".$data{$gene}->[0]."\n";
		exit(0);
	}	
	if ($data{$gene}->[1] && $data{$gene}->[2]) {
		if ((exists $CDS{$gene}) || (exists $hist{$gene}) || (exists $NDB{$gene}) || (exists $pseudo{$gene})) {
			print OUT "S_Child CDS_child \"$gene\"  $data{$gene}->[1] $data{$gene}->[2]\n";
		}
		elsif (exists $RNA{$gene}) {
			print OUT "S_Child Transcript_child \"$gene\"  $data{$gene}->[1] $data{$gene}->[2]\n";
		}
		else {
			print STDERR "Problem: $gene doesn't fit into any class\n";
		}
	}
	else {
		print STDERR "Problem: no coordinates for $gene\n";
		exit(0);
	}

	#########################
	# generate delete lines #
	#########################
	
	print OUT "\n-D Sequence \"$gene\"\n";
	print OUT "\nSequence \"$feature{$gene}->{'source'}\"\n";
	print OUT "-D Subsequence \"$gene\"\n";
}
$/ = "";











################################################################################
#                                subroutines                                   #
################################################################################

sub make_lists {

# genes (CDSs)
my $command2 = <<EOF;
find Predicted_gene
list 
quit
EOF

	open (TACEII, "echo '$command2' | $tace | ");
	while (<TACEII>) {		
		chomp;
		next if ($_ eq "");
		next if (/\/\//);
		s/acedb\>\s//g;
		s/\"//g;
		if (/(\S+)/) {
			$CDS{$1} = $1;		
		}	
	}

# history objects
my $command3 = <<EOF;
find Sequence 
query Method = "history"
list
quit
EOF

	open (TACEIII, "echo '$command3' | $tace | ");
	while (<TACEIII>) {
		chomp;
		next if ($_ eq "");
		next if (/\/\//);
		s/acedb\>\s//g;
		s/\"//g;
		if (/(\S+)/) {
			$hist{$1} = $1;		
		}	
	}

# pseudogenes
my $command4 = <<EOF;
find Sequence
query Method = "Pseudogene" 
list
quit
EOF

	open (TACEIV, "echo '$command4' | $tace | ");
	while (<TACEIV>) {
		chomp;
		next if ($_ eq "");
		next if (/\/\//);
		s/acedb\>\s//g;
		s/\"//g;
		if (/(\S+)/) {
			$pseudo{$1} = $1;		
		}	
	}

# RNA genes (transcripts)
my $command5 = <<EOF;
find Sequence
query RNA
list
quit
EOF

	open (TACEV, "echo '$command5' | $tace | ");
	while (<TACEV>) {
		chomp;
		next if ($_ eq "");
		next if (/\/\//);
		s/acedb\>\s//g;
		s/\"//g;
		if (/(\S+)/) {
			$RNA{$1} = $1;		
		}	
	}

# NDB_CDS objects
my $command6 = <<EOF;
find Sequence
query Method = "NDB_CDS"
list
quit
EOF

	open (TACEV, "echo '$command6' | $tace | ");
	while (<TACEV>) {
		chomp;
		next if ($_ eq "");
		next if (/\/\//);
		s/acedb\>\s//g;
		s/\"//g;
		if (/(\S+)/) {
			$NDB{$1} = $1;		
		}	
	}
}

################################################################################

sub convert_object {

	my $type  = shift;
	
	foreach my $line (@lines) {
		chomp($line);
		
		# convert "Sequence"
		if ($line =~ /^Sequence\s+:/) {
			if ($type eq 'CDS') {
				$line =~ s/Sequence/CDS/g;
			}
			elsif ($type eq 'RNA') {
				$line =~ s/Sequence/Transcript/g;
			}		
			print OUT "$line\n";	
		}
		
		# convert "Source"
		elsif ($line =~ /Source\s+.*?\"(\S+)\"/) {
			$feature{$name}->{'source'} = $1;
			$line =~ s/Source/SMap S_parent Sequence/g;
			print OUT "$line\n";
		}
		
		# convert RNA -> Transcript
		elsif (($type eq 'RNA') && ($line =~ /^RNA/)) {
			$line =~ s/RNA/Transcript/g;
			print OUT "$line\n";
		}
		
		# get rid of unwanted tags
		else {
			my $notwantedline = 0;
			if ($type eq 'CDS') {
				foreach my $tag (@notwantedtagsCDS) {
					$notwantedline = 1 if ($line =~ /^$tag/);	
				}
			}
			if ($type eq 'RNA') {
				foreach my $tag (@notwantedtagsRNA) {
					$notwantedline = 1 if ($line =~ /^$tag/);	
				}
			}
			print OUT "$line\n" unless ($notwantedline == 1);
		}
	}
	print OUT "\n";
}

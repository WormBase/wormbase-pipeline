#!/usr/local/bin/perl5.6.0 -w
#
# find_utrs.pl
#
# usage find_utrs.pl -d <database> -r <output dir>
#
# Ashwin Hajarnavis ah3@sanger.ac.uk  August 2002
#
# Last updated by: $Author: dl1 $                 
# Last updated on: $Date: 2002-08-16 10:30:36 $   

use strict;
use Getopt::Std;
use Data::Dumper;
use vars qw($opt_d $opt_r);

$|=0;

 ##############################
 # command-line options       #
 ##############################

getopts ('d:r:');

my $usage = $0;
$usage .= " -d [WS\#\#] -r [path to output directory](for Ace file, gff file, error file)";

unless ($opt_d) {
    die "$usage";
}

my $gff_dir            = "/wormsrv2/autoace/GFF_SPLITS/$opt_d";
my $best_results_files = $opt_r;

 ##############################
 # Paths etc                  #
 ##############################

my $dbdir ="/wormsrv2/$opt_d/";
if ($opt_d eq "current_DB") {
    my $WS = get_wormbase_version($dbdir);
    $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/$WS";
}

my $tace      = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";
my $query_def = "/nfs/team71/phd/ah3/UTR_DATA/Find_UTRs/cDNA_CDS_EST5_EST3_METHOD.def";
my $genes     = "genes.gff";
my $est_file  = "BLAT_EST_BEST.gff";
my $mrna_file = "BLAT_mRNA_BEST.gff";

log_file($opt_r, $dbdir) if $opt_r;

my %link_coordinate;
my %link_chrom;
my @chromosome = ('CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X');
my $chrom2parse;

# Link data
open (LINK, "</wormsrv2/autoace/BLAT/chromosome.ace") || die "whoops $!";
while (<LINK>) {
    if (/^Sequence \: \"CHROMOSOME_(\S+)\"/) {
	$chrom2parse = $1;
	next;
    }

    if (/Subsequence\s+\"(\S+)\"\s+(\d+)/) {
	$link_coordinate{$1} = $2;
	$link_chrom{$1}      = $chrom2parse;
    }
}
close LINK;

# check to see if EST hash data exists
# make it via tablemaker queries if absent

my %ESTacc2ESTname;
my %EST_name;
my %EST_dir;

unless (-e "/wormsrv2/autoace/BLAT/EST.dat") {
    (%EST_name,%EST_dir) = &make_EST_hash;
}
# else read it into memory
else {
    print "Read hash from file EST.dat\n";
    open (FH, "</wormsrv2/autoace/BLAT/EST.dat") or die "EST.dat : $!\n";
    undef $/;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}

print "Reverse key/value for EST_name hash .....";
%ESTacc2ESTname = reverse %EST_name;
print "OK.\n";

#
# Create hashes for gene positions and transcript mappings

my ($gene_est_map, $est_data) = parse_table($dbdir, $query_def);

# write gene gff positions
my $gene_data = read_gff($gff_dir, $genes);

# write EST gff positions
$est_data = read_gff($gff_dir, $est_file, $est_data);

# write mRNA gff positions
$est_data = read_gff($gff_dir, $mrna_file, $est_data);

 ##############################
 # Write output               #
 ##############################

if ($opt_r)  {
    open(OUTFILE, ">$opt_r/ALL.UTRs.tmp")|| die "Cannot open outfile $opt_r\n";
    open (ERRORFILE, ">$opt_r/data_errors.txt")|| die "Cannot dump errors\n";
    open(ACEFILE, ">$opt_r/UTRs.ace")||die "Cannot open ace file\n";
}

find_transcript($gene_est_map, $gene_data, $est_data);

if ($opt_r) {
    close (OUTFILE);
    gff_sort();
}
 

exit(0);

#####################################################################################################
#                                                                                                   #
# Subroutines:                                                                                      #
#                                                                                                   #
#####################################################################################################



sub find_transcript {
    my ($gene_est_map,$gene_data,$est_data) = @_;
    my $a = "a";
    my ($gene_id, $gene_chrom, $gene_start, $gene_end, $gene_strand, $three_start, $three_end, $five_start, $five_end);
    
    my $ests;
    while (($gene_id, $ests) = each (%{$gene_est_map})) {
	unless (defined($gene_data->{$gene_id})) {
	    print ERRORFILE "$gene_id is not in GFF\n" if $opt_r;  ### gene_id appears in tablemaker query but not in GFF
	    next;
	}
	
	if ($gene_id =~ /(.*)([a-z])$/) {
	    if ($2 ne "a") {
		unless (defined($gene_data->{$1.$a})) {
		    print ERRORFILE "$gene_id exists but $1"."a does not\n" if $opt_r;  ### for cases where there appears to be a .b variant, but no .a variant in GFF.
		    next;
		}
	    }
	}
	
	$gene_chrom  = $gene_data->{$gene_id}->[0];
	$gene_start  = $gene_data->{$gene_id}->[1];
	$gene_end    = $gene_data->{$gene_id}->[2];
	$gene_strand = $gene_data->{$gene_id}->[3];
	
	
	my ($transcript_start, $transcript_end, $min_cdna, $max_cdna) = find_extent($ests, $est_data, $gene_chrom, $gene_start, $gene_end, $gene_id, $gene_strand);
	my $comments;
	if (($transcript_start == 1.e12) || ($transcript_end == -1.e12)) {
	    print ERRORFILE "cDNAs associated with $gene_id not found in BEST.GFF\n" if $opt_r; ### if none of the cDNAs associated with a CDS are BEST matches.
	    next;
	}
		
	if(($transcript_start > $gene_start)&&($transcript_end < $gene_end)) {
	    $comments = "cds $gene_id spans transcript";
	}
	else {
	    $comments = "$gene_id";
	}

	if ($gene_strand eq "+") {
	    $five_start  = $transcript_start;
	    $five_end    = $gene_start - 1;
	    $three_start = $gene_end +1;
	    $three_end   = $transcript_end;
	    if ($five_start < $gene_start) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t5_UTR\t$five_start\t$five_end\t.\t$gene_strand\t.\t$comments\n" if $opt_r;
		print_ace("5_UTR", $gene_id, $min_cdna, $gene_chrom, $five_start, $five_end) if $opt_r;
	    }
	    if ($three_end > $gene_end) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t3_UTR\t$three_start\t$three_end\t.\t$gene_strand\t.\t$comments\n" if $opt_r;
		print_ace("3_UTR", $gene_id, $max_cdna, $gene_chrom, $three_start, $three_end) if $opt_r;
	    }
	}   
	
	if ($gene_strand eq "-") {
	    $five_start  = $gene_end+1;
	    $five_end    = $transcript_end;
	    $three_start = $transcript_start;
	    $three_end   = $gene_start -1;
	    if ($five_end > $gene_end) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t5_UTR\t$five_start\t$five_end\t.\t$gene_strand\t.\t$comments\n" if $opt_r;
		print_ace("5_UTR", $gene_id, $max_cdna, $gene_chrom, $five_end, $five_start) if $opt_r;
	    }
	    if ($three_start < $gene_start) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t3_UTR\t$three_start\t$three_end\t.\t$gene_strand\t.\t$comments\n" if $opt_r;
		print_ace("3_UTR", $gene_id, $min_cdna, $gene_chrom, $three_end, $three_start) if $opt_r;
	    }
	}   
    }
}

sub print_ace {
   my ($type,$gene_id,$evidence,$gene_chrom,$five_end,$three_end) = @_;

   print ACEFILE "\nUTR :\t\"$type:$gene_id\"\n";
   print ACEFILE "Matching_CDS\t$gene_id\n";
   if ($ESTacc2ESTname{$evidence} ne "") {
       print ACEFILE "$type\tAccession_evidence\t\"$ESTacc2ESTname{$evidence}\"\n";
   }
   else {
       print ACEFILE "$type\n";
   }
   print ACEFILE "Species\t\"Caenorhabditis elegans\"\n";
   print ACEFILE "Method\t\"UTR\"\n\n";

   my $link; my $rel_five_end; my $rel_three_end;
   foreach $link (keys %link_chrom) {
       next unless ($link_chrom{$link} eq $gene_chrom); 
#       print "Checking $link [$link_chrom{$link} $link_coordinate{$link}]\n";
       next if ($five_end  < $link_coordinate{$link});  # start is before link begins
       $rel_five_end  = $five_end  - $link_coordinate{$link} + 1;
       $rel_three_end = $three_end - $link_coordinate{$link} + 1;
       print ACEFILE "Sequence :\t\"$link\"\n";
       print ACEFILE "UTR\t\"$type:$gene_id\"\t$rel_five_end\t$rel_three_end\n\n";
       print ACEFILE "//Found match - UTR maps to $link at $rel_five_end - $rel_three_end [abs: $five_end $three_end]\n\n";
#       print "Found match - UTR maps to $link at $rel_five_end - $rel_three_end\n\n";
       last;
   }
} 
      

sub find_extent {
    my ($ests,$est_data,$gene_chrom,$gene_start,$gene_end,$gene_id,$gene_strand) = @_;
    my $current_start = 1.e12;
    my $current_end = -1.e12;
    my $min_cdna;
    my $max_cdna;
    my $i;
    foreach (@{$ests}) {
	unless (defined($est_data->{$_}->[0])) {
	    print ERRORFILE"$_ is not found in BLAT_EST_BEST or BLAT_mRNA_BEST\n" if $opt_r;
	    next;
	}
	unless(defined($est_data->{$_}->[3])) {
	    print ERRORFILE "$_ has not been allocated a method\n" if ($opt_r);
	    next;
	}
	unless ($gene_chrom eq $est_data->{$_}->[0]) {
	    print ERRORFILE "$gene_id on chromosome $gene_chrom, best $_ match is on  $est_data->{$_}->[0]\n" if $opt_r;
	    next;
	}
	 
	if ($est_data->{$_}->[3] eq "EST_3") {
	    if($est_data->{$_}->[2] eq  $gene_strand) {
		print ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the same strand as $gene_id ($gene_strand)\n" if ($opt_r);
		next;
	    }
	}
	elsif ($est_data->{$_}->[3] eq "EST_5") {
	    if($est_data->{$_}->[2] ne  $gene_strand) {
		print  ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the other strand to $gene_id ($gene_strand)\n" if ($opt_r);
		next;
	    }
	}
	elsif ($est_data->{$_}->[3] eq "NDB") {
	    if($est_data->{$_}->[2] ne $gene_strand) {
		print ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the other strand to $gene_id ($gene_strand)\n" if ($opt_r);
		next;
	    }
	} 
	
	for($i=0; $i<scalar(@{$est_data->{$_}->[1]}); $i++) {
	#    print "$_\t $est_data->{$_}->[0] $est_data->{$_}->[1]->[$i]->[0]\t $est_data->{$_}->[1]->[$i]->[1]\t $i\n";
	    if ((($gene_start - $est_data->{$_}->[1]->[$i]->[0]) > 10000) || ($est_data->{$_}->[1]->[$i]->[1] - $gene_end) >10000) {
		print ERRORFILE "cDNA $_\t$est_data->{$_}->[0]\t$est_data->{$_}->[1]->[$i]->[0]\t$est_data->{$_}->[1]->[$i]->[1]\tis matched to $gene_id\t$gene_start\t$gene_end\n" if $opt_r;
		next;
	    }
	    
	    if ($current_start > $est_data->{$_}->[1]->[$i]->[0]) {
		$current_start = $est_data->{$_}->[1]->[$i]->[0];
		$min_cdna = $_;
	    }
	    if ($current_end < $est_data->{$_}->[1]->[$i]->[1]) {
		$current_end = $est_data->{$_}->[1]->[$i]->[1];
		$max_cdna = $_;
	    }
	}
    }
    return ($current_start, $current_end, $min_cdna, $max_cdna);
}

sub read_gff {
    my $gff_dir = shift;
    my $gff_type = shift;
    my @chromosomes = qw(I II III IV V X);
    my $rh;
    if (@_) {
	print "second pass\n";
	$rh = shift;
    }
    else {
	print "first pass\n";
    }
	
    print "Reading $gff_type.\n";
    foreach(@chromosomes) {
	open (GFFFILE, "<$gff_dir/CHROMOSOME_$_\.$gff_type") || die "Cannot open $gff_type file $gff_dir/CHROMOSOME_$_\.$gff_type\n";
	while (<GFFFILE>) {
	    chomp;
	    next if /^\#/;
	    my ($chrom, $curated, $sequence, $start, $end, $score, $strand, $frame, $comments) = split /\t/;
	    
	    if ($curated eq "curated") {
		$comments =~ s/.*"(.*)".*/$1/;
	
		$chrom =~ s/CHROMOSOME_//;
		
		$rh->{$comments}->[0] = $chrom;
		$rh->{$comments}->[1] = $start;
		$rh->{$comments}->[2] = $end;
		$rh->{$comments}->[3] = $strand;
	    }
	    elsif (($curated eq "BLAT_EST_BEST")||($curated eq "BLAT_mRNA_BEST")) {
		$comments =~ s/.*"(.*)".*/$1/;
		$comments =~ s/.*://;
		$chrom =~ s/CHROMOSOME_//;
		$rh->{$comments}->[0] = $chrom;
		$rh->{$comments}->[2] = $strand;
		
		if (defined(@{$rh->{$comments}->[1]})) {
		    my @coords = ($start, $end);
		    push (@{$rh->{$comments}->[1]}, \@coords);
		}
		else {
		    $rh->{$comments}->[1]->[0]->[0] =$start;
		    $rh->{$comments}->[1]->[0]->[1] =$end;
		}
	    }
	    else {
		#	print "is this a trna $curated $comments\n";
		next;
	    }
	}
	close (GFFFILE);
    }
    print "Finished reading $gff_type.\n";
    return $rh;
}



    

sub parse_table {
    $ENV{'ACEDB'} = shift;
    my $query_def = shift;
    my %cds_est = ();
    my $cds_est = \%cds_est;
    my $est_data;
    my $q = 0;
    my $command = "Table-maker -p $query_def\nquit\n";
    open (TACE, "echo '$command' | $tace | ") || die "Cannot query acedb. $command  $tace\n";
    while (<TACE>) {
	chomp;
	s/acedb\>//g;
	next if ($_ eq "");
	next if (/\/\//);
	next if (/^\*/); ###removes lines starting with * 
	next if (/@/);   ###removes lines with email addresses (@)
	next if !(/\d/); ###ACHTUNG - this removes all lines with no numbers in them.
	s/\"//g;
	(/^(\S+)\s/);
	my ($est, $cds, $five_tag, $three_tag, $method) = split /\t/;
	
	if (defined(@{$cds_est->{$cds}})) {
	    push (@{$cds_est->{$cds}}, $est);
	}
	else {
	    $cds_est->{$cds}->[0] = $est;
	}
	if ($method eq "NDB") {
	    $est_data->{$est}->[3] = $method;
	}
	elsif ($method eq "EST_elegans") {
	    $est_data->{$est}->[3] = $five_tag if ($five_tag); 
	    $est_data->{$est}->[3] = $three_tag if ($three_tag); 
	    
	    unless(($five_tag)xor($three_tag)) {
		print ERRORFILE  "Check 3 and 5 tags on $est\n";
		next;
	    }
	}
    }
    print "Tablemaker query complete.\n";
    return (\%cds_est, $est_data);
}

sub printtest {
    my $cds_est = shift;
    my ($cds, $est, $i);
    while (($cds, $est) = each(%{$cds_est})) {
	print "$cds";
	if (scalar(@{$cds_est->{$cds}})>1) {
	    for ( $i=0; $i<scalar(@{$cds_est->{$cds}}); $i++) {
		print  "\t$cds_est->{$cds}->[$i]\n";
	    }
	}
	else {
	    print "\t$cds_est->{$cds}->[0]\n";
	}
    }
}
    

sub log_file {
    my ($file,$dbdir) = @_;
    my $WS = get_wormbase_version($dbdir);
    open (README, ">>$file/resultslog.txt")|| die "Cannot open README.txt\n";
    print README "\n\n\n##################################################\n";
    print README "$file/CHROMOSOME_n_.transcripts.gff  CREATED BY $0 \n\n";
    print README "Connecting to WS$WS\nperl $0 -d $opt_d -r $opt_r\n\nON ";
    print README  `date`;
    print README "BY ";
    print README `whoami`;
}



sub gff_sort {
    
    my (@a, @n, @s, @e, @st, @f, $f, $i, $fh);
    
    open(OUTFILE, "<$opt_r/ALL.UTRs.tmp")|| die "Cannot open outfile $opt_r\n";

    open(ONEFILE,   ">$opt_r/CHROMOSOME_I.UTRs.gff")   || die "Cannot open final file in  $opt_r\n";
    open(TWOFILE,   ">$opt_r/CHROMOSOME_II.UTRs.gff")  || die "Cannot open final file in  $opt_r\n";
    open(THREEFILE, ">$opt_r/CHROMOSOME_III.UTRs.gff") || die "Cannot open final file in  $opt_r\n";
    open(FOURFILE,  ">$opt_r/CHROMOSOME_IV.UTRs.gff")  || die "Cannot open final file in  $opt_r\n";
    open(FIVEFILE,  ">$opt_r/CHROMOSOME_V.UTRs.gff")   || die "Cannot open final file in  $opt_r\n";
    open(XFILE,     ">$opt_r/CHROMOSOME_X.UTRs.gff")   || die "Cannot open final file in  $opt_r\n";
    
    while (<OUTFILE>) {
	s/\#.*//;
	next unless /\S/;
	@f = split /\t/;
	push @a,  $_;
	push @n,  $f[0];
	push @s,  $f[3];
	push @e,  $f[4];
	push @st, $f[6];
    }
    
    foreach $i (sort { $n[$a] cmp $n[$b] or  $st[$a] cmp $st[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0..$#a) {
	print ONEFILE   "$a[$i]" if ($a[$i] =~ /CHROMOSOME_I\t/);
	print TWOFILE   "$a[$i]" if ($a[$i] =~ /CHROMOSOME_II\t/);
	print THREEFILE "$a[$i]" if ($a[$i] =~ /CHROMOSOME_III\t/);
	print FOURFILE  "$a[$i]" if ($a[$i] =~ /CHROMOSOME_IV\t/);
	print FIVEFILE  "$a[$i]" if ($a[$i] =~ /CHROMOSOME_V\t/);
	print XFILE     "$a[$i]" if ($a[$i] =~ /CHROMOSOME_X\t/);
    }     
}

sub get_wormbase_version {
  my $dbdir = shift;
  my $WS_version = `grep "NAME WS" $dbdir/wspec/database.wrm`;
  chomp($WS_version);
  $WS_version =~ s/.*WS//;    
  return($WS_version);
}






__END__;


// Spread sheet definition for the ACeDB software 
// User: ah3
// Date: 2002-08-13_14:36:19

// %n (%%n in the graphic) are parameter to be given on the command line in tace
// or by default by the Parameters given in this file
// \%n (%n in the graphic) are substituted by the value of column n at run time
// Line starting with // are ignored, starting with # are comments

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Condition (Method = "EST_elegans") OR (Method = "NDB")
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Sequence 
From 1 
Tag Matching_Genomic  
 
Colonne 3 
Width 12 
Optional 
Visible 
Show_Tag 
From 1 
Tag EST_5  
 
Colonne 4 
Width 12 
Optional 
Visible 
Show_Tag 
From 1 
Tag EST_3  
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Method 
From 1 
Tag Method 
 
 

// End of these definitions


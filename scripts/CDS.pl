#!/usr/local/bin/perl
#
# This program gets all CDS sequences from an EMBL file
# and writes them in FastA format.
# Author: Antonio Miranda
# Nov/99
# Version 1


###################################################################
# opening the file and parsing its lines
###################################################################

my $inputfile = shift;

open (INPUTFILE, "/usr/local/pubseq/bin/pfetch -F $inputfile |") || die "cannot open $inputfile: $!";
while (<INPUTFILE>) {

    chomp;

    if (/^\s+[acgtunACGTUN]/) {
	s/\d//g;
	s/\s//g;
	$seq .= $_;
    }
	  
    if (/^SV\s+(\S+)\.(\d+)/) {
	$acc = $1, $sv = $2;
    }
  
    if (/\/protein_id=\"(\S+.+)\"/) {
	$CDS_line  .= "\t$1\n";
	push (@protein_id,$1);
    }

    if (/\//) {
	$report = 0; #turn off;
	next;
    }

    if ((/^FT\s+CDS\s+/) || ($report == 1)) {
	$report = 1; #turn on;
	$CDS_line .= substr($_,21);
#	print "// Found a gene \"$CDS_line\"\n";
	next;
    }
}

###########################################################################
# Creation of separate arrays for one exon for'd,
# one exon compl, multi-exon for'd and multi exon compl.
###########################################################################

@init_list = split(/\n/,$CDS_line);

foreach (@init_list) {

    if ((/join/) && (/complem/)) {
#	print "#1 Processing $_ from entry\n";
	push(@multiexoncomp, $_);
    } 
    elsif (/join/) {
#	print "#2 Processing $_ from entry\n";
	push(@multiexonforw, $_);
    }
    elsif (/complem/) {
#	print "#3 Processing $_ from entry\n";
	push(@oneexoncomp, $_);
    }
    else {
#	print "#4 Processing $_ from entry\n";
	push(@oneexonforw, $_);
    }
}

########################################################################
# Creation of hashes containing the following data: CDS begin/end and EMBL_id

foreach (@oneexonforw) {
    @seg1 = split(/\t/,$_);
    @seg2 = (@seg2, @seg1);
}

foreach (@oneexoncomp) {
    @seg3 = split(/\t/,$_);
    @seg4 = (@seg4, @seg3);
}

foreach (@multiexonforw) {
    @seg5 = split(/\t/,$_);
    @seg6 = (@seg6, @seg5);
}

foreach (@multiexoncomp) {
    @seg7 = split(/\t/,$_);
    @seg8 = (@seg8, @seg7);
}

%oneef   = @seg2;
%oneec   = @seg4;
%multief = @seg6;
%multiec = @seg8;
 
###########################################################
# Cleaning the data: removal of words, characters, etc.
###############################################################

foreach (keys (%oneef)) {
    s/join//;
    s/\(//;
    s/\<//;
    s/\>//;
    s/\)//;
    push(@list1,$_);
#    print "#1 $_\n";
    @EMBL_list1 = values(%oneef);
}

foreach (keys (%oneec)) {
    s/join//;
    s/\(//;
    s/\)//;
    s/\<//;
    s/\>//;
    s/complement//;
    push(@list2,$_);
#    print "#3 $_\n";
    @EMBL_list2 = values(%oneec);
}

foreach (keys (%multief)) {
    s/join//;
    s/\(//;
    s/\)//;
    s/\<//;
    s/\>//;
    push(@list3,$_);
#    print "#2 $protein_id $_\n";
    @EMBL_list3 = values(%multief);
}

foreach (keys (%multiec)) {
    s/join//;
    s/\(//g;
    s/\)//g;
    s/\<//;
    s/\>//;
    s/complement//g;
    push(@list4,$_);
#    print "#4 $_\n";
    @EMBL_list4 = values(%multiec);
}

##############################################################
# Processing begin/end data, getting the substrings,
# constructing the whole CDS again (for multi exon
# genes)
###############################################################

foreach (@list1) { 
    while (/(\d+)\.\.(\d+)/g) {
	$x = $1, $y = $2;
	$span = $y - $x + 1;
	$oneexon = substr($seq, $x - 1, $span);
	chomp($oneexon);

	$CDS = $oneexon;
	&print_lines_sixty;
    }
}

foreach (@list2) { 
    while (/(\d+)\.\.(\d+)/g) {
	$x = $1, $y = $2;
	$span = $y - $x + 1;
	$oneexonc = substr($seq, $x - 1, $span);
	chomp($oneexonc);

	$oneexonc =~ s/(\w)/$1\./g;
	@oneexonc = split(/\./g, $oneexonc);
	@oneexonc = reverse(@oneexonc);
	$oneexonc = join("",@oneexonc);
	$oneexonc =~ tr/acgtu/tgcaa/;

	$CDS = $oneexonc;
	&print_lines_sixty;
    }
}

foreach (@list3) {
    while (/(\d+)\.\.(\d+)/g) {
	$x = $1, $y = $2;
	$span = $y - $x + 1;
	$exon = substr($seq, $x - 1, $span);
	chomp($exon);
	push(@exons,$exon);
    }
    $gene = join("",@exons);
    $CDS = $gene;
    &print_lines_sixty;
    @exons = ();
}

foreach (@list4) {
    while (/(\d+)\.\.(\d+)/g) {
	$x = $1, $y = $2;
	$span = $y - $x + 1;
	$exonc = substr($seq, $x - 1, $span);
	chomp($exonc);
	$exonc =~ s/(\w)/$1\./g;
	@exonc = split(/\./g, $exonc);
	@exonc = reverse(@exonc);
	$exonc = join("",@exonc);
	$exonc =~ tr/acgtu/tgcaa/;
	push(@exonsc,$exonc);
    }
    $genec = join("",@exonsc);
    $CDS = $genec;
    &print_lines_sixty;
    @exonsc = ();
}

###############################################################
# print the sequence in lines of sixty characters
###############################################################

sub print_lines_sixty {
    
    $CDS_length = length($CDS);
    $lines = $CDS_length / 60;
    if (($CDS_length % 60) != 0) {
	$lines = $lines + 1;
    }
    $x = 0;
    $span = 60;

    print ">$protein_id[$out_CDS] acc=$acc ver=$sv\n";
    for ($a = 1; $a <= $lines; $a++) {
	$subs = substr($CDS, $x, $span);
	print "$subs\n";
	$x = $x + 60;
    }

    $out_CDS++;
}

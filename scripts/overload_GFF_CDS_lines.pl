#!/usr/local/bin/perl5.8.0
#
# overload_GFF_CDS_lines.pl
# 
# by Dan Lawson
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2005-04-15 12:41:02 $
#

#
#    1. Brief_identification
#    2. Locus names when available
#    3. WormPep IDs
#    4. Confirmed_by status
#    5. WBGene IDs
#    6. Genbank accesions (for region:Genomic_canonical)
#
#	        [1]                                                                      [3]                    [2]             [4]                   [5]
# CDS "JC8.10a" ; Note "Inositol polyphosphate phosphatase, catalytic domain homologues" ; WormPep "WP:CE28239" ; Note "unc-26" ; Confirmed_by "cDNA" ; Gene "WBGene00006763"
#

use Getopt::Long;

my $chrom;

GetOptions (
	    "chrom:s"     => \$chrom
	    );

$verbose = 1;

my $release = shift;

# check that the supplied $release is numeric
unless ($release =~ /\d+/) {
    exit(1);x
}

print "// Working with wormpep release $release\n" if ($verbose);


# get data from wormpep release

open (WORMPEP, "</wormsrv2/WORMPEP/wormpep${release}/wormpep${release}");
while (<WORMPEP>) {
    if (/^>(\S+) (\S+) (WBGene\d+) (\S+.+)\s+status\:(\S+)/) {
	$CDS           = $1;
	$wormpep{$CDS} = $2;
	$geneID{$CDS}  = $3;
	$status{$CDS}  = $5;
	
	$line = $4;
	
	if ($line =~ /locus:(\S+)\s+(\S+.+)/) {
	    $locus{$CDS}   = $1;
	    $briefID{$CDS} = $2;
	}
	elsif ($line =~ /locus:(\S+)$/) {
	    $locus{$CDS}   = $1;
	}
	else {
	    $briefID{$CDS} = $line;
	}
    }
    elsif (/^>(\S+) (\S+) (WBGene\d+) status\:(\S+)/) {
	$CDS           = $1;
	$wormpep{$CDS} = $2;
	$geneID{$CDS}  = $3;
	$status{$CDS}  = $4;
    }
    
}
close WORMPEP;

# test output from wormpep
#foreach $i (sort keys %wormpep) {
#    print "CDS: $i\t$wormpep{$i}\t$geneID{$i}\t$locus{$i}\t$status{$i}\t\'$briefID{$i}\'\n";
#}

# parse GFF lines

our $splitdir = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
our $gffdir   = "/wormsrv2/autoace/CHROMOSOMES";

if (defined($chrom)){
    unless (grep { $chrom eq $_ } ('I','II','III','IV','V','X','MtDNA')) {
	die "ERROR: $chrom is an incorrect chromosome number, please use I, II, III etc.\n";
    }
    @files = (
	      "CHROMOSOME_${chrom}"
	      );
}
else {
    @files = (
              'CHROMOSOME_I',
              'CHROMOSOME_II',
              'CHROMOSOME_III',
              'CHROMOSOME_IV',
              'CHROMOSOME_V',
              'CHROMOSOME_X'
              );
}

our @gff_files = sort @files;
undef @files;

foreach my $file (@gff_files) {
    
    # Check for existance of GFFfile
    unless (-e "$gffdir/$file.gff") {
        &usage("No GFF file");
    }
 
    open (OUT, ">$gffdir/${file}.CSHL.gff");

    open (GFF, "<$gffdir/${file}.gff")  || die "Cannot open $file.gff\n";
    while (<GFF>) {
	chomp;
	
	#skip header lines of file
	unless  (/^\S+\s+curated\s+CDS/) {
	    print OUT "$_\n";
	    next;
	}
	
	($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name) = split /\t/;
	
	($i) = $name =~ (/CDS \"(\S+)\"/);
	
	print OUT "$chromosome\t$source\t$feature\t$start\t$stop\t$score\t$strand\t$other\tCDS \"$i\" ;";
	print OUT " Note \"$briefID{$i}\" ;"        if ($briefID{$i} ne "");
	print OUT " WormPep \"WP:$wormpep{$i}\" ; " if ($wormpep{$i} ne "");
	print OUT " Note \"$locus{$i}\" ; "         if ($locus{$i} ne "");
	print OUT " Status \"$status{$i}\" ; "      if ($status{$i} ne "");
	print OUT " Gene \"$geneID{$i}\" ; "        if ($geneID{$i} ne "");
	print OUT "\n";
    }
    close GFF; #_ end of input GFF file

    close OUT; #_ end of output GFF file
}


# hasta luego
exit(0);


sub usage {
 my $error = shift;
 
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  elsif ($error eq "No GFF file") {
      # No GFF file to work from
      print "One (or more) GFF files are absent from $gffdir\n\n";
      exit(0);
  }
}

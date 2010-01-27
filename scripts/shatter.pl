#!/usr/local/bin/perl
#
# shatter fasta file into separate sequences
#
# Usage : shatter <filename> <No. of entries per output file> <name for output files>
#

use strict;
use Getopt::Long;

my( $file,$bin_size,$output_name,$format );

GetOptions( "file:s" => \$file,
	    "bin:i"  => \$bin_size,
	    "output:s" => \$output_name,
	    "format:s" => \$format);


die "$0 -file <file> -bin <int> -output <dir/id> -format <fasta|fastq>\n\n" unless ($file and $bin_size and $output_name and $format);
my $count = 0;
my $output_count = 1;
my $rec_line_count;
my $line;


&fastq if ($format eq 'fastq');
&fasta if ($format eq 'fasta');


sub fastq {
    print "trying to shatter $file in to $bin_size called $output_name\n";
    open (OUTPUT, ">${output_name}_${output_count}.$format") || die "can't open file ${output_name}_${output_count}.$format :$! \n" ;

    open (FILE, "<$file") or die "cant open $file :$!\n";
    while (<FILE>) {
	$line = $_;
	
	s/^\n$//; #cleanup emplty lines
	print OUTPUT;
	$rec_line_count++;
	$count++ if (($rec_line_count % 4 )==0);
	#print STDERR "$count\n";
	if ($count == $bin_size) {
	    close OUTPUT;
	    
	    $output_count++;
	    
	    open (OUTPUT, ">${output_name}_${output_count}.$format") || die "can't open file ${output_name}_${output_count}.$format :$!\n" ;
	    
	    $count = 0;
	}
    }
    
    close FILE;
    close OUTPUT;
}


sub fasta {
    print "trying to shatter $file in to $bin_size called $output_name\n";

    open (OUTPUT, ">${output_name}_${output_count}") || die "can't open file ${output_name}_${output_count} :$! \n" ;
    
    open (FILE, "<$file") or die "cant open $file :$!\n";
    while (<FILE>) {
	$line = $_;
	
	s/^\n$//; #cleanup emplty lines

	if (/^>(\S+)/) {

	    if ($count == $bin_size) {
		close OUTPUT;

		$output_count++;

		open (OUTPUT, ">${output_name}_${output_count}") || die "can't open file ${output_name}_${output_count} :$!\n" ;

		$count = 0;
	    }
	    $count++;
	}
	print OUTPUT;
    }
    close FILE;

    close OUTPUT;
}

#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

use strict;
use Getopt::Std;
use GDBM_File;
use vars qw($opt_s $opt_t);

getopts ("st");
my $release = shift;
die "please enter release version for new dataset\n" unless $release;

my $usage = "cat swissprot/trembl .fasta file | swiss_trembl2slim.pl -release_no\n";
$usage .= "-s for swissprot\n";
$usage .= "-t for trembl\n";


my %exclude;
$exclude{'Caenorhabditis elegans'}       = 1;
$exclude{'Drosophila melanogaster'}      = 1;
$exclude{'Saccharomyces cerevisiae'}     = 1;
$exclude{'Homo sapiens'}                 = 1;
$exclude{'Human immunodeficiency virus'} = 1;
$exclude{'Caenorhabditis briggsae'}      = 1;

our $output; # file to write
my $output_dir = "/acari/work2a/wormpipe/swall_data";
my $input_dir = "/acari/work2a/wormpipe/swall_data";

my %HASH;

if ($opt_s && $opt_t) {
    die "$usage";
}
elsif ($opt_s) {
    unless (-s "$input_dir/swissprot2org") {
        die "$input_dir/swiss2org not found or empty";
    }
    tie %HASH,'GDBM_File', "$input_dir/swissprot2org",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
    $output = "$output_dir/slimswissprot";
}
elsif ($opt_t) {
    unless (-s "$input_dir/trembl2org") {
        die "$input_dir/trembl2org not found or empty";
    }
    tie %HASH,'GDBM_File', "$input_dir/trembl2org",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
    $output = "$output_dir/slimtrembl";
}
else {
    die "$usage";
}

read_fasta (\*STDIN);
untie %HASH;

sub read_fasta {
    local (*FILE) = @_;
    my ($id, $acc, $seq);
    open (OUT,">$output") or die "cant write to $output\n";
    while (<FILE>) {
        chomp;
        if (/^>(\S+)\.\d+\s+(\S+)/) {
            my $new_acc = $1;
	    my $new_id = $2;
            if ($acc) {
                $seq =~ tr/a-z/A-Z/;
                my $org;
	 	if (exists($HASH{$acc})){
		  $org = $HASH{$acc};
		  if( (exists $exclude{$org}) || ($org =~ /Human immunodeficiency virus/) ){
                    	$id = $new_id; $acc = $new_acc; $seq = "" ;
                    	next;
		    }
                    my $count = 0;
                    $seq = reverse $seq;
                    print OUT ">$acc";
                    while (my $base = chop $seq) {
                    	if ($count++ % 50 == 0) {
                        	print OUT "\n";
                    	}
                    	print OUT $base;
                    }
                    print OUT "\n";
            	}
            }
	    $id = $new_id ;$acc = $new_acc; $seq = "" ;
	} 
        elsif (eof) {
            if ($acc) {
                $seq .= $_ ;
                $seq =~ tr/a-z/A-Z/;
                my $org;
	 	if (exists($HASH{$acc})){
	            $org = $HASH{$acc};
                    if (exists $exclude{$org}) {
                    	next;
		    }
                    my $count = 0;
                    $seq = reverse $seq;
                    print OUT ">$id $acc ($org)";
                    while (my $base = chop $seq) {
                    	if ($count++ % 50 == 0) {
                        	print OUT "\n";
                    	}
                    	print OUT $base;
                    }
                    print OUT "\n";
            	}
            }
        }
        else {
            $seq .= $_ ;
        }
    }
}

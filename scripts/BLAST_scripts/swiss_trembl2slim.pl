#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

use strict;
use Getopt::Std;
use DB_File;
use vars qw($opt_s $opt_t);

getopts ("st");
my $release = shift;
die "please enter release version for new dataset\n" unless $release;

my $usage = "cat swissprot/trembl .fasta file | swiss_trembl2slim.pl -release_no\n";
$usage .= "-s for swissprot\n";
$usage .= "-t for trembl\n";


my %exclude = (
        'Caenorhabditis elegans'       => 1,
        'Drosophila melanogaster'      => 1,
        'Saccharomyces cerevisiae'     => 1,
        'Homo sapiens'                 => 1,
        'Human immunodeficiency virus' => 1,
        'Caenorhabditis briggsae'      => 1,
        'Caenorhabditis remanei'       => 1,
        'Pristionchus pacificus'       => 1,
	'Caenorhabditis japonica'      => 1,
	'Caenorhabditis brenneri'      => 1,
	'Pristionchus pacificus'       => 1,
	'Brugia malayi'                => 1,
        'Onchocerca volvulus'          => 1,
	) ;

our $output; # file to write
my $output_dir = $ENV{'PIPELINE'} . "/swall_data";
my $input_dir = $ENV{'PIPELINE'} . "/swall_data";


my %HASH;

if ($opt_s && $opt_t) {
    die "$usage";
}
elsif ($opt_s) {
    unless (-s "$input_dir/swissprot2org") {
        die "$input_dir/swiss2org not found or empty";
    }
    tie (%HASH, 'DB_File', "$input_dir/swissprot2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or die "cannot open $input_dir/swissprot2org DBM file\n";
    $output = "$output_dir/slimswissprot";
}
elsif ($opt_t) {
    unless (-s "$input_dir/trembl2org") {
        die "$input_dir/trembl2org not found or empty";
    }
    tie (%HASH, 'DB_File', "$input_dir/trembl2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or die "cannot open $input_dir/trembl2org DBM file\n";
    $output = "$output_dir/slimtrembl";
}
else {
    die "$usage";
}

read_fasta ($ENV{'PIPELINE'} . '/blastdb/Supported/uniprot'); # downloaded by update_blastDBs.pl
untie %HASH;

sub read_fasta {
    my $uniprot = shift;
    my ($id, $acc, $seq);
    open (OUT,">$output") or die "cant write to $output\n";
    open (FILE,"<$uniprot") or die "cant read $uniprot\n";
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
    close FILE;
    close OUT;
}

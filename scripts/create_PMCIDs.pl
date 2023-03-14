#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Const::Fast;
use LWP::Simple;
use IO::String;
use Path::Class;
use Log_files;
use Wormbase;

const my $WB_XREF_URL => 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref';
const my $EUROPE_PMC_FTP_URL => 'https://europepmc.org/pub/databases/pmc/DOI/PMID_PMCID_DOI.csv.gz';
const my $EUROPE_PMC_SAVE_FILE => $ENV{'BUILD_HOME'} . '/BUILD_DATA/MISC_DYNAMIC/PMID_PMCID_DOI.csv.gz';

my ($out_file, $in_file, $type, $debug, $test);
GetOptions(
    "out|o=s"   => \$out_file,
    "in|i=s"    => \$in_file,
    "type|t=s"  => \$type,
    "debug|d:s" => \$debug,
    "test"      => \$test
    ) or die ("USAGE: perl create_PMCIDs.pl --in <in_file> --out <out_file> --type <ace|list>\n");

my $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
    );

# establish log file.
my $log = Log_files->make_build_log($wormbase);

unless ($type eq 'ace' || $type eq 'list') {
    $log->log_and_die("Expecting type argument of either 'ace' or 'list'\n");
}


my %pmid2wb;
if ($type eq 'list') {
    # Get map of WB to PubMed IDs from Caltech
    my $wb2pmid = get_caltech_xref_map();

    # Create reversed map for entries in input file
    my $fh = file($in_file)->openr;
    while (my $wb_id = $fh->getline()) {
	chomp $wb_id;
	if (exists $wb2pmid->{$wb_id}) {
	    $pmid2wb{$wb2pmid->{$wb_id}} = $wb_id;
	}
    }
} else {
    # Parse ace file to get map of PubMed to WB IDs
    my $fh = file($in_file)->openr;
    my $current_wb_id;
    my $wbcount = 0;
    my $unique_count;
    while (my $line = $fh->getline()) {
	if ($line =~ /^Paper : "(WBPaper\d+)"/) {
	    $current_wb_id = $1;
	} elsif ($line =~ /^DB_info/) {
	    my @col = split /\s+/, $line;
	    if ($col[9] eq "\"PMID\"") {
		my $pmid = $col[12];
		$pmid =~ s/"//g;
		$unique_count++ unless exists $pmid2wb{$pmid};
		$pmid2wb{$pmid} = $current_wb_id;
		$wbcount++;
	    }
	}
    }
}

# Download PMCID ID mapping file from FTP site and save in build data folder
my $rc = getstore($EUROPE_PMC_FTP_URL, $EUROPE_PMC_SAVE_FILE);
if (is_error($rc)) {
  $log->log_and_die("getstore of <$EUROPE_PMC_FTP_URL> failed with $rc");
}

# Create outfile
my $out_fh = file($out_file)->openw or $log->log_and_die("Cannot write to $out_file\n");
$out_fh->print('Database : "Europe_PMC"' . "\n" .
	       'Name "Europe PMC"' . "\n" .
	       'URL "https://europepmc.org"' . "\n" .
	       'URL_constructor "https://europepmc.org/search?query=%s"' . "\n\n" .
	       'Database_field : "PMCID"' . "\n\n");

open(MAP, "gunzip -c $EUROPE_PMC_SAVE_FILE |") or $log->log_and_die("Cannot open gzip stream for $EUROPE_PMC_SAVE_FILE |");
my $first_line = 1;
while (<MAP>) {
    chomp;
    my @col = split(",", $_);
    if ($first_line) {
	die ("Unexpected format for file $EUROPE_PMC_SAVE_FILE\n") if ($col[1] ne 'PMCID' || $col[0] ne 'PMID');
	$first_line = 0;
	next;
    }
    my $pmid = $col[0];
    my $pmcid = $col[1];
    next unless exists $pmid2wb{$pmid};
    $out_fh->print('Paper : "' . $pmid2wb{$pmid} . "\"\n" . 'Database "Europe_PMC" "PMCID" "' . $pmcid . "\"\n\n") if length $pmcid > 0;
}
close (MAP);

$log->mail;
exit(0);

sub get_caltech_xref_map {
    my %xref_map;

    my $handle = IO::String->new(get($WB_XREF_URL));
    while (defined (my $line = <$handle>)) {
	chomp $line;
	$line =~ s/<BR>$//;
	my ($wbpaper_id, $xref_id) = split(" ", $line);
	if ($xref_id =~ /pmid(\d+)/) {
	    $xref_map{$wbpaper_id} = $1;
	}
    }
    close $handle;

    return \%xref_map;
}

#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

use strict;
use DBI;
use Getopt::Std;
use vars qw($opt_h $opt_u $opt_p $opt_d $opt_a $opt_s $opt_g);

getopts ("h:u:p:d:a:sg:");

my $usage .= "waba2gff.pl\n";
$usage .= "-h [db host]\n";
$usage .= "-u [db user]\n";
$usage .= "-d [db name]\n";
$usage .= "-p [db pass]\n";
$usage .= "-a [analysisId]\n";
$usage .= "-s dump states, and not whole waba matches\n";
$usage .= "-g [gff file needed to remap the briggsae sequences to their parents, e.g. supercontigs]\n";

unless ($opt_h && $opt_u && $opt_d && $opt_a) {
    die "$usage";
}

my $analysis = $opt_a;

# connect to the mysql database
my $dbh = DBI -> connect("DBI:mysql:$opt_d:$opt_h", $opt_u, $opt_p, {RaiseError => 1})
    || die "cannot connect to db, $DBI::errstr\n";

###################################
# read the gff file if required
my %gff;
if ($opt_g) {
    open (GFF, "$opt_g") || die "cannot open gff file $opt_g\n"; 
    while (<GFF>) {
	chomp;
	my @a = split /\t/;
	$a[8] =~ /\"(\S+)\"/;
	my $child = $1;
	$gff{$child} = [$a[0], $a[3], $a[4]];
    }
    close GFF;
}

###################################
# get the mapping of internal_id to id
my %internal2id;
my $sth_c = $dbh->prepare ( q{ SELECT id, internal_id, length
                                 FROM contig
			     } );
$sth_c->execute;
my $ref = $sth_c->fetchall_arrayref;
foreach my $aref (@$ref) {
    my ($id, $internal_id, $length) = @$aref;
    $internal2id{$internal_id} = $id;
}

###################################
# get the features, and write the gff file
my $sth;
if ($opt_s) {
    $sth = $dbh->prepare ( q{ SELECT contig, seq_start, seq_end,
                                     hid, hstart, hend,
                                     score, perc_id, strand,
                                     state, cigar
                                FROM waba_feature
                               WHERE analysis = ?
                            ORDER BY hid,hstart,hend
			    } );
}
else {
    $sth = $dbh->prepare ( q{ SELECT contig, seq_start, seq_end,
                                     hid, hstart, hend,
                                     score, perc_id, strand
                                FROM waba_fset
                               WHERE analysis = ?
                            ORDER BY hid,hstart,hend
                            } );
}

$sth->execute ($analysis);
while (my ($internal_id, $start, $end, $hid, $hstart, $hend, $score, $pid, $strand, $state, $cigar) = $sth->fetchrow_array) {
#while (<>) {
#    chomp;
#    my ($internal_id, $start, $end, $hid, $hstart, $hend, $score, $pid, $strand, $state, $cigar) = split /\s+/;
    $strand = ($strand > 0) ? "+" : "-";
    # get id
    my $id;
    if (exists $internal2id{$internal_id}) {
        $id = $internal2id{$internal_id};
    }
    else {
        print STDERR "no mapping of internal_id to id for $internal_id\n";
    }
    # process id's from ensembl names
    if ($id =~ /^(.+)\.[^\.]+\.[^\.]+\.[^\.]+/) {
        $id = $1;
    }
    my $new_id;
    if ($opt_g) {
	unless (exists ($gff{$id})) {
	    die "no parent info for $id\n";
	}
	$start += $gff{$id}->[1]-1;
	$end += $gff{$id}->[1]-1;
	$new_id = $gff{$id}->[0];
	my @cigars = split (/:/, $cigar);
	for (my $i = 0; $i < scalar(@cigars); $i++) {
	    my @coors = split (/,/, $cigars[$i]);
	    $coors[0] += $gff{$id}->[1]-1;
	    $cigars[$i] = join (",", @coors);
	}
	$cigar = join (":", @cigars);
    }
    else {
      $new_id = $id;
    }
    unless ($opt_s) {
        $state = "all";
    }
    print "$hid\twaba_$state\tsimilarity\t$hstart\t$hend\t$score\t$strand\t.\tTarget \"Sequence:$new_id\" $start $end $cigar\n";
}

$sth_c->finish;
$sth->finish;
$dbh->disconnect;







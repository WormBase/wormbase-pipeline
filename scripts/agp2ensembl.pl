#!/usr/local/bin/perl -w
#
# agp2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# modified for reading in .agp files for worm ensembl

=pod

=head1 NAME

agp2ensembl.pl

=head1 SYNOPSIS

agp2ensembl.pl -agp chromX.agp <db options> -write

=head1 DESCRIPTION

Load clone into EnsEMBL database from raw sequence. Takes an agp file
and checks to see which clones need to be loaded. Calls fetch on the
clone.version and loads clone as a single contig.

=head1 OPTIONS

    -dbhost  DB host
    -dbuser  DB user
    -dbname  DB name
    -dbpass  DB pass
    -phase   clone phase
    -agp     agp file
    -write   write clone
    -v       print info about clones and contigs

=head1 CONTACT

Simon Potter: scp@sanger.ac.uk
Kerstin Jekosch : ar2@sanger.ac.uk

=head1 BUGS

Insert list of bugs here!


=cut

use lib '/nfs/disk100/humpub/modules/PerlModules';

use strict;
use Getopt::Long;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::PerlDB::Clone;
use Bio::EnsEMBL::PerlDB::Contig;
use Hum::NetFetch qw( wwwfetch );
use Hum::EMBL;

my($id, $acc, $ver, $phase, $contigs);
my($agp, $single, $seqio, $seq);
my($dbname, $dbhost, $dbuser, $dbpass);
my($help, $info, $write, $replace, $verbose);


$Getopt::Long::autoabbrev = 0;   # personal preference :)
$dbuser = 'wormadmin';           # default

my $ok = &GetOptions(
    "phase=s"   => \$phase,
    "agp=s"     => \$agp,
    "dbname=s"  => \$dbname,
    "dbhost=s"  => \$dbhost,
    "dbuser=s"  => \$dbuser,
    "dbpass=s"  => \$dbpass,
    "help"      => \$help,
    "info"      => \$info,
    "write"     => \$write,
    "v"         => \$verbose
);

if ($help || not $ok) {
    &usage;
    exit 2;
} elsif ($info) {
    exec("perldoc $0");
}

unless ($dbname && $dbuser && $dbhost) {
    print STDERR "Must specify all DB parameters\n";
    exit 1;
}

unless ($agp) {
    print STDERR "Must specify apg file\n";
    exit 1;
}

if ($phase < 0 && $phase > 4) {

    print STDERR "Phase should be 1, 2, 3 or 4\n";
    exit 1;
}
$phase = -1 unless defined $phase;


# open connection to EnsEMBL DB
my $dbobj;
if ($dbpass) {
    $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-pass'   => $dbpass,
        '-dbname' => $dbname

    ) or die "Can't connect to DB $dbname on $dbhost as $dbuser"; # Do we need password as well?
}
else {
    $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-dbname' => $dbname
    ) or die "Can't connect to DB $dbname on $dbhost as $dbuser"; # Do we need password as well?
}


open AGP, "< $agp" or die "Can't open AGP file $agp";
while (<AGP>) {

    my @fields = split;
    my $sv = $fields[5]; # col 6 is SV - do we need anything else?
    ($acc, $ver) = $sv =~ /(\w+)\.(\d+)/;
    unless ($acc && $ver) {
        print STDERR "Invalid $sv\n";
        next;
    }

    if (&is_in_db($dbobj, $sv)) {
        print "Found $sv; skipping\n";
        next;
    }
    $seq = fetch_seq($acc, $ver);
    unless ($seq) {
        print STDERR "Can't fetch $sv\n";
        next;

    }

    my $clone = new Bio::EnsEMBL::PerlDB::Clone;
    my $contig = new Bio::EnsEMBL::PerlDB::Contig;
    my $length = $seq->length;

    $clone->id          ($acc);
    $clone->htg_phase   ($phase);
    $clone->embl_id     ($acc);
    $clone->version     (1);
    $clone->embl_version($ver);
    $contig->id         ("$acc.$ver.1.$length");
    $contig->embl_offset(1);
    $contig->length     ($length);
    $contig->seq        ($seq);
    $contig->version    (1);
    $contig->embl_order (1);

    print "Clone ", $clone->id, "\n";
    if ($verbose) {
        print "\tembl_id     ", $clone->embl_id, "\n";
        print "\tversion     ", $clone->version, "\n";
        print "\temblversion ", $clone->embl_version, "\n";
        print "\thtg_phase   ", $clone->htg_phase, "\n";
    }
    print "Contig ", $contig->id, "\n";
    if ($verbose) {
        print "\toffset: ", $contig->embl_offset, "\n";
        print "\tlength: ", $contig->length, "\n";
        print "\tend:    ", ($contig->embl_offset + $contig->length - 1), "\n";
        print "\tversion:", $contig->version, "\n";
        print "\torder:  ", $contig->embl_order, "\n";
        print "\tlength: ", $contig->length, "\n";
    }
    print "\n";

    $clone->add_Contig($contig);

    if ($write) {
        eval {
            $dbobj->write_Clone($clone);
        };
        if ($@) {
            print STDERR "Error writing clone $sv\n"; 
        }
    }
}



sub usage {
    print <<EOF
$0 [options]
Options:
  -dbname
  -dbhost
  -dbuser
  -dbpass
  -agp      agp file
  -phase    clone phase
  -write    write clone
  -v        print info about clones and contigs
EOF
}



sub fetch_seq {
    my ($acc, $sv) = @_;

    my $embl_str = wwwfetch($acc)
        or die "Can't fetch '$acc'";
    my $embl = Hum::EMBL->parse(\$embl_str);
    
    my $embl_sv = $embl->SV->version;
    die "EMBL SV '$embl_sv' doesn't match SV '$sv' for clone '$acc'" unless $sv == $embl_sv;
    
    my $seq_str = $embl->Sequence->seq;
    my $seq = Bio::Seq->new(
        -id     => $acc,
        -seq    => $seq_str,
        );
    return $seq;
}    


sub is_in_db {
    my ($dbobj, $sv) = @_;
    my $clone;

    eval {
        $clone = $dbobj->get_Clone($sv);
    };
    if (defined $clone) {
        return 1;
    }
    else {
        return 0;
    }
}

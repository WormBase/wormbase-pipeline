#!/usr/local/ensembl/bin/perl -w
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
Kerstin Jekosch : kj2@sanger.ac.uk

=head1 BUGS

Insert list of bugs here!


=cut

use lib '/nfs/farm/Worms/Ensembl/ensembl-pipeline/modules';
use lib '/nfs/farm/Worms/Ensembl/ensembl/modules';
use lib '/nfs/disk100/humpub/modules/PerlModules';

use strict;
use Getopt::Long;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;
use Hum::NetFetch qw( wwwfetch );
use Hum::EMBL;

my($id, $acc, $ver, $phase, $contigs);
my($agp, $single, $seqio, $seq, $fasta, $strict);
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
    "v"         => \$verbose,
    "fasta=s"   => \$fasta,
    "strict"    => \$strict
);

if ($help || not $ok) {
    &usage;
    exit 2;
} elsif ($info) {
    exec("perldoc $0");
}

unless ($dbname && $dbuser && $dbhost) {
    print STDERR "Must specify all DB parameters\n";
    &usage;
    exit 1;
}

unless ($agp) {
    print STDERR "Must specify apg file\n";
    exit 1;
}

if (defined $phase && ($phase < 0 || $phase > 4)) {

    print STDERR "Phase should be 1, 2, 3 or 4\n";
    exit 1;
}
$phase = -1 unless defined $phase;


# open connection to EnsEMBL DB
my $dbobj;
if ($dbpass) {
    $dbobj = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-pass'   => $dbpass,
        '-dbname' => $dbname

    ) or die "Can't connect to DB $dbname on $dbhost as $dbuser"; # Do we need password as well?
}
else {
    $dbobj = Bio::EnsEMBL::DBSQL::Pipeline::DBAdaptor->new(
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-dbname' => $dbname
    ) or die "Can't connect to DB $dbname on $dbhost as $dbuser"; # Do we need password as well?
}

my $clone_adaptor = $dbobj->get_CloneAdaptor();
my $sic = $dbobj->get_StateInfoContainer();
my $analysis_adaptor = $dbobj->get_AnalysisAdaptor();
my $submitted_analysis = $analysis_adaptor->fetch_by_dbID(1); #1 is dummy analysis to mark addition

my %seqs;
if ($fasta) {
    open (FH , "$fasta") || die "cannot open file $fasta";
    %seqs = read_fasta (\*FH);
}


open AGP, "< $agp" or die "Can't open AGP file $agp";
while (<AGP>) {
    chomp;
    my @fields = split;
    my $sv = $fields[5]; # col 6 is SV - do we need anything else?
    ($acc, $ver) = $sv =~ /(\w+)\.(\d+)/;
    unless ($acc && $ver) {
        if ($strict) {
            print "Invalid $sv: $_\n";
            next;
	}
        else {
            $acc = $sv;
            $ver = 1;
	}
    }

    if (&is_in_db($clone_adaptor, $sv)) {
        print "Found $sv; skipping\n";
        next;
    }
    elsif ( &update_existing_clone($clone_adaptor, $sv) ) {
      print "Found old version of $sv; updated\n";
      next;
    }

    if ($fasta) {
        my $seq_str = $seqs{$acc};
        $seq = Bio::Seq->new(
            -id     => $acc,
            -seq    => $seq_str,
            );
        unless ($seq) {
            print "Can't fetch $sv\n";
            next;
	}
    }
    else {
        $seq = fetch_seq($acc, $ver);
        unless ($seq) {
            print "Error fetching $sv\n";
            next;
	}
    }

    my $clone = new Bio::EnsEMBL::Clone;
    my $contig = new Bio::EnsEMBL::RawContig;
    my $length = $seq->length;

    $clone->id       ($acc);
    $clone->htg_phase   ($phase);
    $clone->embl_id     ($acc);
    $clone->version     (1);
    $clone->embl_version($ver);
    my $time = time;
    $clone->created($time);

    $contig->name         ("$acc.$ver.1.$length");
    $contig->embl_offset(1);
    $contig->length     ($length);
    $contig->seq        ($seq->seq);
  #  $contig->version    (1);
  #  $contig->embl_order (1);

    print "Clone ", $clone->id, "\n";
    if ($verbose) {
        print "\tembl_id     ", $clone->embl_id, "\n";
        print "\tversion     ", $clone->version, "\n";
        print "\temblversion ", $clone->embl_version, "\n";
        print "\thtg_phase   ", $clone->htg_phase, "\n";
    }
    print "Contig ", $contig->id, "\n";
    if ($verbose) {
     #   print "\toffset: ", $contig->embl_offset, "\n";
        print "\tlength: ", $contig->length, "\n";
     #   print "\tend:    ", ($contig->embl_offset + $contig->length - 1), "\n";
     #   print "\tversion:", $contig->version, "\n";
     #   print "\torder:  ", $contig->embl_order, "\n";
    }
    print "\n";

    $clone->add_Contig($contig);
    $clone->modified($time);

    if ($write) {
        eval {
            $clone_adaptor->store($clone);
        };
        if ($@) {
            print "Error writing clone $sv\n"; 
        }
	else {
	  $sic->store_input_id_analysis($contig->name,$submitted_analysis) ;
	  if( $@ ) {
	    print "Error update input_id_analysis table for ",$contig->name,"\n";
	  }
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
  -fasta    fasta file of sequences, no wwwfetch
  -strict   enforces accession.version syntax
EOF
}



sub fetch_seq {
    my ($acc, $sv) = @_;

    my $embl_str = wwwfetch($acc)
        or die "Can't fetch '$acc'";

    my $embl = Hum::EMBL->parse(\$embl_str);
    
    my $embl_sv = $embl->SV->version;
    unless ($sv == $embl_sv) {
        print "EMBL SV '$embl_sv' doesn't match SV '$sv' for clone '$acc', skip it\n";
        return 0;
    }
    my $seq_str = $embl->Sequence->seq;
    my $seq = Bio::Seq->new(
        -id     => $acc,
        -seq    => $seq_str,
        );
    return $seq;
}

# If clone already exists in database we can update rather than add the new version
sub update_existing_clone
  {
    my $clone_adaptor = shift;
    my $sv = shift;
    $sv =~ /(\w+)\.(\d+)/;

    my $acc = $1;
    my $ver = $2;

    my $clone = $clone_adaptor->_generic_sql_fetch(qq{ WHERE name = '$acc' });# I know this is bad :)
    if ($clone ) {
      my $old_ver = $clone->embl_version();

      if ($old_ver and ($ver > $old_ver) ) {
	# update clone in database

	#update the clone version and embl_version
	my $db_ver = $clone->version();
	$db_ver++;
	my $time = time;

	&make_SQL_query($dbobj,"UPDATE clone SET version = $db_ver WHERE name = \"$acc\"");
	&make_SQL_query($dbobj,"UPDATE clone SET embl_version = $ver WHERE name = \"$acc\"");
	&make_SQL_query($dbobj,"UPDATE clone SET modified = FROM_UNIXTIME($time) WHERE name = \"$acc\"");

	#update contig
	my $contigs = $clone->get_all_Contigs();
	#each clone has only one contig
	foreach my $contig ( @{$contigs} ) { 
	
	  # update the DNA 
	  my $dna = $seqs{$acc};
	  my $dna_id = $contig->dna_id();
	  &make_SQL_query($dbobj,"UPDATE dna SET sequence = \"$dna\" WHERE dna_id = $dna_id");

	  # update the contig
	  my $length = length($dna);
	  my $name = "$acc."."$ver."."1."."$length";
	  my $contig_id = $contig->dbID;

	  &make_SQL_query($dbobj,"UPDATE contig SET name = \"$name\" WHERE contig_id = $contig_id");
	  $contig->length($length);
	  &make_SQL_query($dbobj,"UPDATE contig SET length = $length WHERE contig_id = $contig_id");
	  return 1;
	}
      }
      else {
	print "clone $acc already has an equal or newer verion (database ver = $old_ver ) in the database\n";
	return 1;
      }
    }
    else {
      return 0;
    }
  }

sub make_SQL_query
  {
    my ($dbobj,$sql) = @_;
    my $sth = $dbobj->prepare($sql);
    $sth->execute;

    print "done\n";
  }



sub is_in_db {
    my ($dbobj, $sv) = @_;
    my $clone;

    eval {
      /(\w+)\.(\d+)/;
      my $acc = $1;
      my $ver = $2;
      $clone = $dbobj->fetch_by_accession_version($acc, $ver);
    };
    if (defined $clone) {
        return 1;
    }
    else {
        return 0;
    }
}

sub read_fasta {
    local (*FILE) = @_;
    my ($id , $seq , %name2seq);
    while (<FILE>) {
        chomp;
        if (/^>(\S+)/) {
            my $new_id = $1;
            if ($id) {
                $seq =~ tr/A-Z/a-z/;
                $name2seq{$id} = $seq;
            }
            $id = $new_id ; $seq = "" ;
        } 
        elsif (eof) {
            if ($id) {
                $seq .= $_ ;
                $seq =~ tr/A-Z/a-z/;
                $name2seq{$id} = $seq;
            }
        }
        else {
            $seq .= $_ ;
        }
    }
    return %name2seq;
}

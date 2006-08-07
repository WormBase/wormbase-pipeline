#!/usr/bin/perl -w
#
# perl compara_test.pl <CDS name>
#

use lib '/nfs/acari/wormpipe/ensembl/bioperl-live/';
use lib '/nfs/acari/wormpipe/ensembl/ensembl-compara/modules/';
use lib '/nfs/acari/wormpipe/ensembl/ensembl/modules/';

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::AlignIO;

# setup database connection
my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => 'ecs1f',        # change that
    -user   => 'wormro',       # and that
    -dbname => 'worm_compara'
);

my $member_adaptor = $compara_db->get_MemberAdaptor();

# fetch a Member
my $member = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE', $ARGV[0] );

$member->print_member;

# then you get the homologies where the member is involved
my $homology_adaptor = $compara_db->get_HomologyAdaptor();
my $homologies       = $homology_adaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );

# the bioperl alignment formatter
my $alignIO = Bio::AlignIO->newFh(
    -interleaved => 0,
    -fh          => \*STDOUT,
    -format      => "clustalw",
    -idlength    => 20
);

my %seqs;    # to store peptide sequences

foreach my $homology ( @{$homologies} ) {

    printf( "homology(%d) %s\n", $homology->dbID, $homology->description );

    # some statistics from codeml
    print join " ", ( map { $homology->{"_$_"} } qw(dn ds n s lnl threshold_on_ds) );
    printf " dS/dN %.4f \n", $homology->{_ds} / $homology->{_dn};

    # print the alignment
    print $alignIO $homology->get_SimpleAlign;

    # short version of the member pair
    foreach my $this ( @{ $homology->get_all_PeptideAlignFeature } ) { $this->display_short() }

    # get the peptide sequences
    foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
        my ( $me, $at ) = @{$ma};
        foreach $pepm ( @{ $me->get_all_peptide_Members() } ) { $seqs{ $pepm->stable_id } = $pepm->sequence }
    }
}

# print the fasta file of all members
while ( my ( $k, $v ) = each %seqs ) { print ">$k\n$v\n" }

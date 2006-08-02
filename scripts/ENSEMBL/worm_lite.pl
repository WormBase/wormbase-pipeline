#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  worm_lite.pl
#
#        USAGE:  ./worm_lite.pl
#
#  DESCRIPTION:
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  07/07/06 15:37:31 BST
#     REVISION:  ---
#===============================================================================

#####################################################################
# needs some makefile love to pull together all GFFs and Fastas
####################################################################

use strict;
use YAML;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use WormBase;
use DBI qw(:sql_types);

my ( $debug, $species, $setup, $dna, $genes );

GetOptions(
    'species=s'  => \$species,
    'setup'      => \$setup,
    'load_dna'   => \$dna,
    'load_genes' => \$genes,
    'debug'      => \$debug,
);

my $config = ( YAML::LoadFile("ensembl_lite.conf") )->{$species};
my $cvsDIR = '/nfs/acari/wormpipe/ensembl/';

our $gff_types="curated coding_exon";

&setupdb( $config->{database}, $config->{version} ) if $setup;
&load_dna($config)   if $dna;
&load_genes($config) if $genes;

# create database from schema
sub setupdb {
    my ( $db, $version ) = @_;
    my $status = 0;
    print ">>creating new database $db->{dbname}\n" ;
    my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password}";
    system("$mysql -e \"DROP DATABASE IF EXISTS $db->{dbname};\"") && die;
    system("$mysql -e \"create database $db->{dbname};\"") && die;
    system( "$mysql $db->{dbname} < " . $cvsDIR . "ensembl-pipeline/sql/table.sql" ) && die;
    system( "$mysql $db->{dbname} < " . $cvsDIR . "ensembl/sql/table.sql" ) && die;
    system("$mysql -e 'INSERT INTO coord_system VALUES (1,\"scaffold\",\"$version\",1,\"default_version,top_level\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO coord_system VALUES (2,\"superlink\",\"$version\",2,\"default_version,sequence_level\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.version\",\"$version\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO analysis (created,logic_name,module) VALUES ( NOW(),\"wormbase\",\"wormbase\");INSERT INTO analysis_description VALUES(1,\"imported from WormBase\",\"WormGene\");' $db->{dbname}") && die;
    system("$mysql $db->{dbname} <attrib_type.sql") && die;
    system("/nfs/team71/worm/mh6/bin/perl /nfs/acari/wormpipe/ensembl/ensembl-pipeline/scripts/load_taxonomy.pl -name \"Caenorhabditis $species\" -taxondbhost ecs2 -taxondbport 3365 -taxondbname ncbi_taxonomy -lcdbhost ecs1f -lcdbport 3306 -lcdbname worm_ensembl_$species -lcdbuser wormadmin -lcdbpass worms") && die;

    if ($status) { die("Error while building the database.") }
}

# load genome sequences
sub load_dna {
    my ($config) = @_;
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $config->{database}->{host},
        -user   => $config->{database}->{user},
        -dbname => $config->{database}->{dbname},
        -pass   => $config->{database}->{password},
        -port   => $config->{database}->{port},
    );
    my $cs = $db->get_CoordSystemAdaptor()->fetch_by_dbID(1);    # chromosome
    my $ss = $db->get_CoordSystemAdaptor()->fetch_by_dbID(2);    # superlink
    my $sa = $db->get_SliceAdaptor();

    foreach my $file ( glob $config->{fasta} ) {
        my $seqs = Bio::SeqIO->new( '-format' => 'Fasta', '-file' => $file );
        next if $file =~ /masked|CSHL/;
        print "processing $file\n" if $debug;
        while ( my $seq = $seqs->next_seq() ) {
            print $seq->display_id(), "\t", $seq->length(), "\n" if $debug;

            my $cslice=&store_slice( $db, $seq->display_name(), 1, $seq->length, 1, $cs );

            # do something else if size > 16MB
	    if ($seq->length > 16000000){
		my $part1=substr($seq->seq(),0,($seq->length()/2));
		my $part2=substr($seq->seq(),($seq->length()/2));
		my $bslice1=&store_slice( $db, $seq->display_name().'_1', 1, (length $part1), 1, $ss, $part1);
		my $bslice2=&store_slice( $db, $seq->display_name().'_2', 1, (length $part2), 1, $ss, $part2);

		&insert_agp_line($cslice->get_seq_region_id($cslice), (length $part1)+1,$seq->length, $bslice2->get_seq_region_id($bslice2), 1,(length $part2) , 1, $db);
		&insert_agp_line($cslice->get_seq_region_id($cslice), 1, (length $part1), $bslice1->get_seq_region_id($bslice1), 1,(length $part1) , 1, $db);
	    }
	    else {
		    my $bslice=&store_slice( $db, $seq->display_name(), 1, $seq->length, 1, $ss, $seq->seq());
		    &insert_agp_line($cslice->get_seq_region_id($cslice), 1, $seq->length, $bslice->get_seq_region_id($bslice), 1,$seq->length , 1, $db);
	    }

        }
    }

    $db->get_MetaContainer()->store_key_value( 'assembly.mapping', $cs->name . ":" . $cs->version . "|" . $ss->name );
    $db->dbc->do('INSERT INTO seq_region_attrib (seq_region_id,attrib_type_id,value) SELECT seq_region_id,11,5 FROM seq_region WHERE Name LIKE "%Mt%"');
    $db->dbc->do('INSERT INTO seq_region_attrib (seq_region_id,attrib_type_id,value) SELECT seq_region_id,6,6 FROM seq_region WHERE coord_system_id=1');
}

sub load_genes {
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $config->{database}->{host},
        -user   => $config->{database}->{user},
        -dbname => $config->{database}->{dbname},
        -pass   => $config->{database}->{password},
        -port   => $config->{database}->{port},
    );
    my $analysis = $db->get_AnalysisAdaptor()->fetch_by_logic_name('wormbase');

    foreach my $file ( glob $config->{gff} ) {
        next if $file =~ /masked|CSHL|BLAT_BAC_END|briggsae/;
        $file=~/.*\/(.*)\.gff/;
	print "parsing $1 from $file\n";
        my $slice = $db->get_SliceAdaptor->fetch_by_region('Scaffold',$1);
        my $genes = &parse_gff( $file, $slice, $analysis );
        &write_genes( $genes, $db );
    }
    $db->dbc->do('UPDATE gene SET biotype="protein_coding"');
}
#  /(\w+\.\d+)[a-z A-Z]*/
package WormBase;

# redefine subroutine to use different tags
sub process_file {
    my ($fh) = @_;
    my ( %genes, $transcript, %five_prime, %three_prime );

  LOOP: while (<$fh>) {
        chomp;

        my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene ) = split;
        my $element = $_;

        next LOOP if ( /^#/ || $chr =~ /sequence-region/ || ( !$status && !$type ) );
        my $line = $status . " " . $type;
        $gene =~ s/\"//g;
        if ( ( $line eq 'Coding_transcript five_prime_UTR' ) or ( $line eq 'Coding_transcript three_prime_UTR' ) ) {
            $transcript = $gene;

            #remove transcript-specific part: Y105E8B.1a.2
            $gene =~ s/(\.\w+)\.\d+$/$1/;
            my $position = $type;
            if ( $position =~ /^five/ ) {
                $five_prime{$gene} = {} if ( !$five_prime{$gene} );
                $five_prime{$gene}{$transcript} = [] if ( !$five_prime{$gene}{$transcript} );
                push( @{ $five_prime{$gene}{$transcript} }, $element );
            }
            elsif ( $position =~ /^three/ ) {
                $three_prime{$gene} = {} if ( !$three_prime{$gene} );
                $three_prime{$gene}{$transcript} = [] if ( !$three_prime{$gene}{$transcript} );
                push( @{ $three_prime{$gene}{$transcript} }, $element );
            }
            next LOOP;
        }
        elsif ( $line ne $gff_types ) { next LOOP }    # <= here goes the change needs tp become $line eq "$bla $blub"

        if ( !$genes{$gene} ) { $genes{$gene} = []; push( @{ $genes{$gene} }, $element ) }
        else { push( @{ $genes{$gene} }, $element ) }
    }
    print STDERR "Have " . keys(%genes) . " genes (CDS), " . keys(%five_prime) . " have 5' UTR and " . keys(%three_prime) . " have 3' UTR information\n";
    return \%genes, \%five_prime, \%three_prime;
}

1;

#!/usr/local/ensembl/bin/perl -w
#===============================================================================
#
#         FILE:  worm_lite.pl
#
#        USAGE:  ./worm_lite.pl
#
#  DESCRIPTION:
#
#       AUTHOR:   (Michael Han), <mh6@sanger.ac.uk>
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
use Storable;

use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose warning);
verbose('OFF');
use FindBin;
use lib "$FindBin::Bin/../lib";
use WormBase;
use DBI qw(:sql_types);

my ( $debug, $species, $setup, $dna, $genes, $test,$store, $yfile, $agp );

GetOptions(
  'species=s'  => \$species,
  'setup'      => \$setup,
  'load_dna'   => \$dna,
  'load_genes' => \$genes,
  'load_agp'   => \$agp,
  'debug'      => \$debug,
  'test'       => \$test,
  'store=s'    => \$store,
  'yfile=s'    => \$yfile,

) || die("bad commandline parameter\n");

if ($store){
	my $storable=Storable::retrieve($store);
	$species= ref $storable;
	$species=~tr/[A-Z]/[a-z]/;
}

die "You must supply a valid YAML config file\n" if not defined $yfile or not -e $yfile;

my $config = ( YAML::LoadFile($yfile) )->{$species};
if ($test) {
  $config = ( YAML::LoadFile($yfile) )->{"${species}_test"};
}
my $cvsDIR = $test
  ? ( YAML::LoadFile($yfile) )->{test}->{cvsdir}
  : ( YAML::LoadFile($yfile) )->{generics}->{cvsdir};

$WormBase::Species = $species;
our $gff_types = ($config->{gff_types} || "curated coding_exon");


&setupdb( $config->{database}, $config->{version} ) if $setup;
&load_dna($config)   if $dna;
&load_agp($config) if $agp;
&load_genes($config) if $genes;


##################################
# create database from schema
#
# hardcoded paths:
#      /nfs/acari/wormpipe/ensembl/ensembl-pipeline/scripts/DataConversion/wormbase/attrib_type.sql
#      /nfs/acari/wormpipe/ensembl/ensembl-pipeline/scripts/load_taxonomy.pl
#
# taxondb: ia64f -taxondbport 3365 -taxondbname ncbi_taxonomy / ens-livemirror 
#  if the taxondb goes down bully Abel

sub setupdb {
    my ( $db, $version ) = @_;
    my $status = 0;
    print ">>creating new database $db->{dbname} on $db->{host}\n";
    my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password}";
    system("$mysql -e \"DROP DATABASE IF EXISTS $db->{dbname};\"") && die;
    system("$mysql -e \"create database $db->{dbname};\"")         && die;
    print "loading table.sql from ensembl\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "ensembl/sql/table.sql" ) && die;
    print "loading table.sql from ensembl-pipeline\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "ensembl-pipeline/sql/table.sql" ) && die;

    system("$mysql -e 'INSERT INTO coord_system VALUES (1,1,\"chromosome\",\"$version\",1,\"default_version,top_level\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO coord_system VALUES (2,1,\"superlink\",\"$version\",2,\"default_version,sequence_level\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO coord_system VALUES (3,1,\"clone\",\"$version\",3,\"default_version\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.version\",\"$version\");' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.start_date\",NOW());' $db->{dbname}") && die;
    system("$mysql -e 'INSERT INTO analysis (created,logic_name,module) VALUES ( NOW(),\"wormbase\",\"wormbase\");' $db->{dbname}")          && die;
    system("$mysql -e 'INSERT INTO analysis_description (analysis_id,description,display_label) VALUES (1,\"imported from WormBase\",\"WormGene\");' $db->{dbname}") && die;
    system("$mysql $db->{dbname} <$cvsDIR/ensembl-pipeline/scripts/DataConversion/wormbase/master_attrib_type.sql") && die;
    system("$mysql $db->{dbname} <$cvsDIR/ensembl-pipeline/scripts/DataConversion/wormbase/attrib_type.sql") && die;
    system("perl $cvsDIR/ensembl-pipeline/scripts/load_taxonomy.pl -name \"$config->{species}\" -taxondbhost ens-livemirror -taxondbport 3306 -taxondbname ncbi_taxonomy -lcdbhost $db->{host} -lcdbport $db->{port} -lcdbname $db->{dbname} -lcdbuser $db->{user} -lcdbpass $db->{password}"
      )
      && die("cannot run taxondb update");

    # in case the taxondb gives up uncomment below
    #    system("$mysql $db->{dbname} </nfs/acari/wormpipe/ensembl/ensembl-pipeline/scripts/DataConversion/wormbase/taxonomy.sql") && die;
    #    system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"species.classification\",\"$species\");' $db->{dbname}") && die;
    #    system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"species.taxonomy_id\",\"$config->{taxon_id}\");' $db->{dbname}") && die;

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
    my $cl = $db->get_CoordSystemAdaptor()->fetch_by_dbID(3);    # clone
    my $sa = $db->get_SliceAdaptor();

    foreach my $file ( glob $config->{fasta} ) {
        my $seqs = Bio::SeqIO->new( '-format' => 'Fasta', '-file' => $file );
        next if $file =~ /masked|CSHL/;
        print "processing $file\n" if $debug;
        while ( my $seq = $seqs->next_seq() ) {
            print $seq->display_id(), "\t", $seq->length(), "\n" if $debug;
            $seq->{'primary_seq'}->{'seq'}=~s/[^acgtnACGTN]/n/g; # removes ambiguity codes
            my $cslice = &store_slice( $db, $seq->display_name(), 1, $seq->length, 1, $cs );

            # do something else if size > 16MB
            if ( $seq->length > 16000000 ) {
                my $part1 = substr( $seq->seq(), 0, ( $seq->length() / 2 ) );
                my $part2 = substr( $seq->seq(), ( $seq->length() / 2 ) );
                my $bslice1 = &store_slice( $db, $seq->display_name() . '_1', 1, ( length $part1 ), 1, $ss, $part1 );
                my $bslice2 = &store_slice( $db, $seq->display_name() . '_2', 1, ( length $part2 ), 1, $ss, $part2 );

                &insert_agp_line( $cslice->get_seq_region_id($cslice), ( length $part1 ) + 1,
                    $seq->length, $bslice2->get_seq_region_id($bslice2), 1, ( length $part2 ), 1, $db);
                &insert_agp_line( $cslice->get_seq_region_id($cslice), 1, ( length $part1 ),
                    $bslice1->get_seq_region_id($bslice1), 1, ( length $part1 ), 1, $db);
            }
            else {
                my $bslice = &store_slice( $db, $seq->display_name(), 1, $seq->length, 1, $ss, $seq->seq() );
                &insert_agp_line( $cslice->get_seq_region_id($cslice), 1, $seq->length, $bslice->get_seq_region_id($bslice), 1, $seq->length, 1, $db );
            }

        }
    }

    $db->get_MetaContainer()->store_key_value( 'assembly.mapping', $cs->name . ":" . $cs->version . "|" . $ss->name .':'.$cs->version );
    $db->get_MetaContainer()->store_key_value( 'assembly.mapping', $cs->name . ":" . $cs->version . "|" . $cl->name .':'.$cs->version );
    $db->get_MetaContainer()->store_key_value( 'assembly.mapping', $cl->name . ":" . $cs->version . "|" . $cs->name .':'.$cs->version."|".$ss->name .':'.$cs->version );

    $db->dbc->do('INSERT INTO seq_region_attrib (seq_region_id,attrib_type_id,value) SELECT seq_region_id,11,5 FROM seq_region WHERE Name LIKE "%Mt%"');
    $db->dbc->do('INSERT INTO seq_region_attrib (seq_region_id,attrib_type_id,value) SELECT seq_region_id,6,6 FROM seq_region WHERE coord_system_id=1');
    
#    $db->close;
    undef $db;
    undef $cs;
    undef $ss;
    undef $cl;
    undef $sa;

    # agp fun
    load_agp($config) if ($config->{agp}); 

}

sub load_agp {
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
    my $cl = $db->get_CoordSystemAdaptor()->fetch_by_dbID(3);    # clone
    my $sa = $db->get_SliceAdaptor();

    if ($config->{agp}) {
        foreach my $file ( glob $config->{agp} ) {
	    my $bzhook= $file=~/.bz2$/ ? "bzcat $file|" : $file;
            my $fh = new IO::File $bzhook or die "couldn't open " . $file . " $!";

	    #### need a slice of a chromosome to get the id :o
	    $file=~/(\w+)\.agp/;
	    print "processing chromosome $1\n";
	    my $chromosome = $db->get_SliceAdaptor()->fetch_by_region('toplevel',$1); # stupid assumption about the name of the file being the chromosome name
            my %chr_hash = %{ &agp_parse( $fh, $chromosome->adaptor->get_seq_region_id($chromosome), $config->{version} ) };

            foreach my $id ( keys(%chr_hash) ) {
		my $length = $chr_hash{$id}->{'contig_end'}-$chr_hash{$id}->{'contig_start'}+1;
		print "adding $1 $id - $chr_hash{$id}->{contig_start} - $chr_hash{$id}->{contig_end} length:$length\n" if $debug;
		##############################
                my $contig = &store_slice( $db, $id, 1, $length, 1, $cl );
                $contig = $contig->adaptor->get_seq_region_id($contig);
                my $agp    = $chr_hash{$id};
                &insert_agp_line($agp->{'chromosome_id'},$agp->{'chr_start'},$agp->{'chr_end'},$contig,$agp->{'contig_start'},$agp->{'contig_end'},$agp->{'contig_ori'},$db);
            }
    	}
    }

}

sub load_genes {
  my ($config) = @_;

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $config->{database}->{host},
    -user   => $config->{database}->{user},
    -dbname => $config->{database}->{dbname},
    -pass   => $config->{database}->{password},
    -port   => $config->{database}->{port},
      );
  my $analysis = $db->get_AnalysisAdaptor()->fetch_by_logic_name('wormbase');
    
  my (%slice_hash, @path_globs, @gff_files, $genes); 

  foreach my $slice (@{$db->get_SliceAdaptor->fetch_all('toplevel')}) {
    $slice_hash{$slice->seq_region_name} = $slice;
    if ($species eq 'elegans') {
      my $other_name;
      if ($slice->seq_region_name !~ /^CHROMOSOME/) {
        $other_name = "CHROMOSOME_" . $slice->seq_region_name; 
      } else {
        $other_name = $slice->seq_region_name;
        $other_name =~ s/^CHROMOSOME_//;
      }
      $slice_hash{$other_name} = $slice;
    }
  }
  
  @path_globs = split(/,/, $config->{gff});

  foreach my $fglob (@path_globs) {
    push @gff_files, glob($fglob);
  }
  
  open(my $gff_fh, "cat @gff_files |") or die "Could not create GFF stream\n";
  $genes = &parse_gff_fh( $gff_fh, \%slice_hash, $analysis);
  &write_genes( $genes, $db );
  
  $db->dbc->do('UPDATE gene SET biotype="protein_coding"');
  $db->dbc->do('INSERT INTO meta (meta_key,meta_value) VALUES ("genebuild.start_date",NOW())');
}

package WormBase;

# redefine subroutine to use different tags
sub process_file {
    my ($fh) = @_;
    my ( %genes, $transcript, %five_prime, %three_prime, %parent_seqs );

  LOOP: while (<$fh>) {
        chomp;

        my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene ) = split;
        my $element = $_;

        next LOOP if ( /^#/ || $chr =~ /sequence-region/ || ( !$status && !$type ) );
        my $line = $status . " " . $type;
        $gene =~ s/\"//g if $gene;
        if ( ( $line eq 'Coding_transcript five_prime_UTR' ) or ( $line eq 'Coding_transcript three_prime_UTR' ) ) {
            $transcript = $gene;

            #remove transcript-specific part: Y105E8B.1a.2
            $gene =~ s/(\.\w+)\.\d+$/$1/ unless $species eq 'brugia'; # for that if i will go to hell :-(
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

        if (not exists $genes{$gene}) {
          $genes{$gene} = [];
          $parent_seqs{$gene} = $chr;
        }
        push( @{ $genes{$gene} }, $element );
    }
    print STDERR "Have " . keys(%genes) . " genes (CDS), " . keys(%five_prime) . " have 5' UTR and " . keys(%three_prime) . " have 3' UTR information\n";
    return \%genes, \%five_prime, \%three_prime, \%parent_seqs;
}

1;

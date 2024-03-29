#!/usr/bin/env perl
# Populate these meta keys with interesting enough results:
# sample.gene_param    Gene ID (ex: maker-uti_cns_0003289-snap-gene-0.4)
# sample.gene_text  Gene ID or name (ex: fkbp6)
# sample.location_param   Location around gene (ex: uti_cns_0003289:9000-13000)
# sample.location_text  Location around gene (ex: uti_cns_0003289:9000-13000)
# sample.search_text  Random word from the gene description (ex: isomerase)
# sample.transcript_param Transcript ID (ex: maker-uti_cns_0003289-snap-gene-0.4-mRNA-1)
# sample.transcript_text  Transcript ID (ex: maker-uti_cns_0003289-snap-gene-0.4-mRNA-1)
use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use List::MoreUtils qw(uniq);
use Carp;


use Getopt::Long;

my (
  $dbname, $compara_dbname,
  $host, $compara_host,
  $user, $compara_user, 
  $port, $compara_port,
  $pass,
  $gene_id,
);



&GetOptions(
  'dbname=s' => \$dbname,
  'user=s'   => \$user,
  'host=s'   => \$host,
  'port=s'   => \$port,
  'pass=s'   => \$pass,
  'compara_dbname=s' => \$compara_dbname,
  'compara_user=s'   => \$compara_user,
  'compara_host=s'   => \$compara_host,
  'compara_port=s'   => \$compara_port,
  'gene_id=s' => \$gene_id,
);

# get a compara dba
my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => sprintf( 'mysql://%s:%s@%s:%d/%s',
                                                                                $compara_user,
                                                                                "",
                                                                                $compara_host,
                                                                                $compara_port,
                                                                                $compara_dbname),
                                                               -species => 'Multi' ) unless ( ! defined $compara_dbname);
my $core_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbname => $dbname,
                                                   -host => $host,
                                                   -port => $port,
                                                   -user => $user, 
                                                   -pass => $pass );
my $gdb_adaptor = $compara_db->get_GenomeDBAdaptor();
my $genome_db = $gdb_adaptor->fetch_by_core_DBAdaptor($core_db);
my $gene = $gene_id ? lookup_gene($core_db, $gene_id) : choose_best_gene(get_candidate_genes($core_db, $compara_db));

store_gene_sample( $core_db, $gene);

sub lookup_gene {
  my ($core_db, $gene_id) = @_;
  return $core_db->get_GeneAdaptor->fetch_by_stable_id($gene_id);
}
sub get_candidate_genes {
  my ($core_db, $compara_db) = @_;
  my @all_genes;

  SEQ: for my $slice (sort { $b->length <=> $a->length } @{$core_db->get_SliceAdaptor->fetch_all('toplevel')}) {

    for my $g (sort { $a->start <=> $b->start } @{$slice->get_all_Genes_by_type('protein_coding')}) {
      my $stats = { 
        description     => [],
        domains         => [],
        model_orthologs => [],
      };

      push @{$stats->{description}}, $g->description if defined $g->description and $g->description !~ /Uncharacterized/i;
      push @{$stats->{description}}, "display_xref: ".$g->display_xref->display_id if defined $g->display_xref;
      
      my @domains;
      for my $feature (@{$g->canonical_transcript->translation->get_all_DomainFeatures}){
        push @domains, $feature->interpro_ac if $feature->interpro_ac;
      }
      @domains= uniq(@domains);
      $stats->{domains} = \@domains;
      my $gm = $compara_db->get_GeneMemberAdaptor->fetch_by_stable_id_GenomeDB($g->stable_id, $genome_db);
      for my $target_species ('homo_sapiens',
                                  'mus_musculus',
                                  'danio_rerio',
                                  'drosophila_melanogaster',
                                  'saccharomyces_cerevisiae',
                                  'caenorhabditis_elegans_prjna13758') {
        for my $homology (@{$compara_db->get_HomologyAdaptor->fetch_all_by_Member($gm, -TARGET_SPECIES => $target_species)}) {
          if ($homology->description() eq 'ortholog_one2one') {
            push @{$stats->{model_orthologs}}, $target_species;
          }
        }
      }
      
      push @all_genes, {
        gene => $g,
        stats => $stats,
        score => &score($stats)
      };
      last SEQ if @all_genes  > 1000;
    }
  }
  return @all_genes;
}

sub choose_best_gene {
  my ($best, $second, $third) = sort {
    $b->{score} <=> $a->{score}
  } @_;

  print "Chose gene " . &choice_to_string($best);
  return $best->{gene};
}




#######################################3

sub choice_to_string {
  my $o = shift;
  return $o->{gene}->stable_id . " score ".$o->{score}."\n" . Data::Dumper::Dumper($o->{stats})."\n";
}
sub score {
  my ($stats) = @_;
  my $desc = @{$stats->{description}};
  my $dom = @{$stats->{domains}};
  my $o = @{$stats->{model_orthologs}};
  # We want a description, three orthologs, and ten protein domains
  return (3 * $desc + 1) * ( 5 * (6-$o)*$o + 1 ) * ((20 - $dom) * $dom + 1 ) ; 

}
sub store_gene_sample {
    my ( $dba, $gene ) = @_;
    my $meta     = $dba->get_MetaContainer();
    my $sr_name  = $gene->seq_region_name();
    my $sr_start = $gene->seq_region_start();
    my $sr_end   = $gene->seq_region_end();
    
    # adjust bounds by 10%
    my $flank = int(( $sr_end - $sr_start + 1 )/10);
    $sr_start -= $flank;
    $sr_end += $flank;
    $sr_start = 1 if ( $sr_start < 0 );
    $sr_end = $gene->slice()->seq_region_length()
	if ( $sr_end > $gene->slice()->seq_region_length() );
    
    $meta->delete_key('sample.location_param');
    $meta->store_key_value( 'sample.location_param',
			    "$sr_name:${sr_start}-${sr_end}" );
    $meta->delete_key('sample.location_text');
    $meta->store_key_value( 'sample.location_text',
			    "$sr_name:${sr_start}-${sr_end}" );
    $meta->delete_key('sample.gene_param');
    $meta->store_key_value( 'sample.gene_param', $gene->stable_id() );
    $meta->delete_key('sample.gene_text');
    $meta->store_key_value( 'sample.gene_text',
			    ( $gene->external_name() || $gene->stable_id()
			    ) );
    my $transcript = $gene->canonical_transcript;
    $meta->delete_key('sample.transcript_param');
    $meta->store_key_value( 'sample.transcript_param',
			    $transcript->stable_id() );
    $meta->delete_key('sample.transcript_text');
    $meta->store_key_value( 'sample.transcript_text',
			    ( $transcript->external_name() ||
			      $transcript->stable_id() ) );
    $meta->delete_key('sample.search_text');
    $meta->store_key_value( 'sample.search_text', 'synthetase' );
    return;
} 

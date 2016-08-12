#!/usr/bin/env perl
  
use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Carp;


use Getopt::Long;

my (
  $dbname, $compara_dbname,
  $host, $compara_host,
  $user, $compara_user, 
  $port, $compara_port,
  $pass, $compara_pass,
  $test,
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
  'compara_pass=s'   => \$compara_pass,
  'test'             => \$test,
);

    
# get a compara dba

my $compara_url = sprintf( 'mysql://%s:%s@%s:%d/%s',
                           $compara_user, 
                           $compara_pass,
                           $compara_host,
                           $compara_port,
                           $compara_dbname);

my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => sprintf( 'mysql://%s:%s@%s:%d/%s',
                                                                                $compara_user, 
                                                                                $compara_pass,
                                                                                $compara_host,
                                                                                $compara_port,
                                                                                $compara_dbname),
                                                               -species => 'Multi' );
my $core_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbname => $dbname,
                                                   -host => $host,
                                                   -port => $port,
                                                   -user => $user, 
                                                   -pass => $pass );

# longest slice is bound to have multiple genes, right?

my @all_genes;

SEQ:foreach my $slice (sort { $b->length <=> $a->length } @{$core_db->get_SliceAdaptor->fetch_all('toplevel')}) {

  my @genes = sort { $a->start <=> $b->start } @{$slice->get_all_Genes_by_type('protein_coding')};
  
  print STDERR "There are ", scalar(@genes) . " gene\n";
  
  foreach my $g (@genes) {

    my $stats = { 
      gene            => $g , 
      description     => 0,
      domains         => 0,
      model_orthologs => 0,
    };

    $stats->{description}++ if defined $g->description and $g->description !~ /Uncharacterized/i;
    $stats->{description}++ if defined $g->display_xref;
        
    my @f = grep { $_->interpro_ac } @{$g->canonical_transcript->translation->get_all_ProteinFeatures};
    
    $stats->{model_domains} = 1 if @f;
    

    my $gm = $compara_db->get_GeneMemberAdaptor->fetch_by_stable_id($g->stable_id);
    foreach my $target_species ('homo_sapiens',
                                'mus_musculus',
                                'danio_rerio',
                                'drosophila_melanogaster',
                                'saccharomyces_cerevisiae',
                                'caenorhabditis_elegans_prjna13758') {
      foreach my $homology (@{$compara_db->get_HomologyAdaptor->fetch_all_by_Member($gm, -TARGET_SPECIES => $target_species)}) {
        if ($homology->description() eq 'ortholog_one2one') {
          $stats->{model_orthologs}++;
        }
      }
    }
    
    push @all_genes, $stats;

    last SEQ if @all_genes  > 1000;
  }


}

my ($best) = sort {
  $b->{description} <=> $a->{description} or $b->{model_orthologs} <=> $a->{model_orthologs}
} @all_genes;
      
print "CHOSE GENE ", $best->{gene}->stable_id, "\n";

&print_gene_sample( $best->{gene} );
if (not $test) {
  &store_gene_sample( $core_db, $best->{gene} );
}


sub print_gene_sample {
  my ($g) = @_;
  
  print "Selected gene :\n";
  printf "  acc  : %s\n", $g->stable_id;
  printf "  name : %s\n", $g->display_xref->display_id;
  printf "  desc : %s\n", $g->description;

}

#######################################3
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

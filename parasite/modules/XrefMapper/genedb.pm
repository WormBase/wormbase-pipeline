
package XrefMapper::genedb;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


# This module is activated by specifying "taxon=parasite" in the mapping input file
# It contains some common config for worms maintained by WormBase (i.e. having genes
# with WBGene ids etc)

sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
           ExonerateGappedBest5 => ['RefSeq_mRNA',
                                    'RefSeq_mRNA_predicted', 
                                    'RefSeq_ncRNA', 
                                    'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}


sub transcript_names_from_gene {
  return;
}

sub gene_display_xref_sources {
  my $self = shift;

  my ($list_ref, $ignore_ref) = $self->XrefMapper::DisplayXrefs::gene_display_xref_sources();


  my @list = qw(
                 Uniprot_gn
                 GeneDB
               );

  return [\@list, $ignore_ref];
}


sub transcript_display_xref_sources {
  my ($self) = @_;

  my ($list_ref, $ignore_ref) = $self->XrefMapper::DisplayXrefs::transcript_display_xref_sources();

  my @list = qw(
                Uniprot_gn
               );

  return [\@list, $ignore_ref];
}



sub gene_description_sources {

  return (
    "Uniprot/SWISSPROT",
    "GeneDB",
    "Uniprot/SPTREMBL",
      );

}



sub gene_description_filter_regexps {

  return (
           '^Uncharacterized protein\s*',
           '^Putative uncharacterized protein\s+\S.+',
           '^Putative uncharacterized protein\s*',
           '^Hypothetical protein\s*',
   );

}

1;

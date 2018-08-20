
package XrefMapper::parasite;

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


sub gene_description_sources {

  return ("RFAM",
          "RNAMMER",
          "TRNASCAN_SE",
	  "miRBase",
          "HGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
          "Uniprot/SPTREMBL",
	  "RefSeq_peptide",
      );

}


sub gene_description_filter_regexps {

  return (
           '^Uncharacterized protein\s*',
           '^Putative uncharacterized protein\s*',
           '^Hypothetical protein\s*',
   );

}

1;

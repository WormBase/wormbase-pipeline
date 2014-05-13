#!/software/bin/perl -w


use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $species = shift;

if (! defined $species) {
    die "Specify a species_name.\n";
}

Bio::EnsEMBL::Registry->load_all();

my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "$species", 'Core', 'gene' );

if (! defined $gene_adaptor) {
    die "can't get a gene adaptor for species, $species\n";
}

# Get all transcripts through all genes

my $max_intron_size = 0;

my $genes_aref = $gene_adaptor->fetch_all_by_biotype('protein_coding');

print STDERR "processing " . @$genes_aref . " genes\n";

foreach my $gene (@$genes_aref) {
    my $transcripts_aref = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts_aref) {
	my $introns_aref = $transcript->get_all_Introns();
	foreach my $intron (@$introns_aref) {
	    if ($intron->length() > $max_intron_size) {
		$max_intron_size = $intron->length();
	    }
	}
    }
}

print "Max Intron size: $max_intron_size\n";

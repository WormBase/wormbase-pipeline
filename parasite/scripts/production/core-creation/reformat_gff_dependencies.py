transcript_types = ["mRNA", "tRNA", "rRNA", "transcript","nontranslating_transcript",
                    "ncRNA", "snRNA", "tRNA_pseudogene",
                    "snoRNA", "scRNA", "piRNA", "lincRNA",
                    "asRNA", "miRNA_mature", "miRNA", "pseudogenic_transcript"]

coding_transcript_types = ["mRNA", "transcript"]

pseudogene_types = ["pseudogene", "Pseudogene"]

gene_types = ["gene"] + pseudogene_types

exon_types = ["exon"]

intron_types = ["intron"]

cds_types = ["CDS"]

allowed_types = list(set(transcript_types + coding_transcript_types + gene_types + pseudogene_types + exon_types + cds_types))

gff_column_names = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

wormbase_source = "WormBase_imported"
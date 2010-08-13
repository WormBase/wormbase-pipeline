#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  kill_transcripts.pl
#
#        USAGE:  ./kill_transcripts.pl <killlist>
#
#  DESCRIPTION: reads from a file transcript stable ids and kills them one by one
#               including exons, translations, protein features, ...
#               will also remove genes with zero transcripts
#
#===============================================================================

use strict;

while(<>){
	chomp;
	print "DELETE t.* FROM transcript t JOIN transcript_stable_id USING(transcript_id) WHERE stable_id = \"$_\";\n";
}

print <<HERE;
DELETE et.* FROM exon_transcript et LEFT JOIN transcript t USING(transcript_id) WHERE t.transcript_id IS NULL;
DELETE ex.* FROM exon ex LEFT JOIN exon_transcript e USING(exon_id) WHERE e.transcript_id IS NULL;
DELETE esi.* FROM exon_stable_id esi LEFT JOIN exon e USING(exon_id) WHERE e.exon_id IS NULL;
DELETE tsi.* FROM transcript_stable_id tsi LEFT JOIN transcript t USING(transcript_id) WHERE t.transcript_id IS NULL;
DELETE tl.* FROM translation tl LEFT JOIN transcript t using(transcript_id) WHERE t.transcript_id IS NULL;
DELETE tsi.* FROM translation_stable_id tsi LEFT JOIN translation t USING(translation_id) WHERE t.translation_id IS NULL;
DELETE pf.* FROM protein_feature pf LEFT JOIN translation t USING(translation_id) WHERE t.translation_id IS NULL;
DELETE gene.* FROM gene LEFT JOIN transcript t USING(gene_id) WHERE t.gene_id IS NULL;
DELETE gsi.* FROM gene_stable_id gsi LEFT JOIN gene g USING(gene_id) WHERE g.gene_id IS  NULL;
DELETE iia.* FROM input_id_analysis WHERE input_id_type='TRANSLATIONID' AND NOT EXISTS(SELECT * FROM translation t WHERE t.translation_id=iia.input_id);
HERE


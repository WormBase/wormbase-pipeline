#!/usr/bin/env ruby

data=Hash.new
ARGF.each{|line|
	(array_name,gene_name,cds_name,type,transcripts)=line.split
	next if transcripts.nil?
	trans=transcripts.split('|')
	data[array_name]||=Array.new
	data[array_name].push(trans)
	data[array_name].flatten!
}

data.each{|array_name,v|
	v.each{|hit|
		print <<HERE

INSERT IGNORE INTO object_xref (ensembl_object_type,ensembl_id,xref_id,linkage_annotation,analysis_id) SELECT 'ProbeSet', probe_set_id,xref_id,'[imported from WormBase]',5 FROM probe_set p, xref x  WHERE x.dbprimary_acc like "#{hit}%" AND p.name="#{array_name}";
INSERT IGNORE INTO object_xref (ensembl_object_type,ensembl_id,xref_id,xref_id,linkage_annotation,analysis_id) SELECT 'Probe', probe_id,xref_id,'[imported from WormBase]',5 FROM probe join probe_set p using (probe_set_id), xref x  WHERE x.dbprimary_acc like "#{hit}%" AND p.name="#{array_name}";
INSERT IGNORE INTO object_xref (ensembl_object_type,ensembl_id,xref_id,linkage_annotation,analysis_id) SELECT 'ProbeFeature', probe_feature_id,xref_id,'[imported from WormBase]',5 FROM probe_feature join probe using (probe_id) join probe_set p using (probe_set_id), xref x  WHERE x.dbprimary_acc like "#{hit}%" AND p.name="#{array_name}";
HERE
  }
}

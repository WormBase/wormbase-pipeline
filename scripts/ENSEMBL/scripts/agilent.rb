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

INSERT IGNORE INTO object_xref (ensembl_object_type,ensembl_id,xref_id) SELECT 'ProbeSet', probe_set_id,xref_id FROM probe_set p, xref x  WHERE x.dbprimary_acc like "#{hit}%" AND p.name="#{array_name}";
INSERT IGNORE INTO object_xref (ensembl_object_type,ensembl_id,xref_id) SELECT 'Probe', probe_id,xref_id FROM probe p, xref x  WHERE x.dbprimary_acc like "#{hit}%" AND p.name="#{array_name}";
HERE
	}
}

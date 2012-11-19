#!/usr/bin/env ruby
# that is a small one-off script to parse Ladeanna's C.briggsae orthologs

def read_common_data(file)
	eval File.new("/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/#{file}",'r')
	return $VAR1
end

def print_ace(gene1,gene2,organism2)
	puts "Gene : \"#{gene1}\""
	puts "Ortholog #{gene2} \"#{organism2}\" Paper_evidence WBPaper00038778",""
end

cds2gene = read_common_data("cds2wbgene_id.dat")

ARGF.each{|line|
	line.chomp!
	columns = line.split(',')
	briggsae_gene = cds2gene[columns[4]]
	elegans_gene = cds2gene[columns[10]]
	next if briggsae_gene.nil? ||elegans_gene.nil?
	print_ace(elegans_gene,briggsae_gene,'Caenorhabditis briggsae')
	print_ace(briggsae_gene,elegans_gene,'Caenorhabditis elegans')
}


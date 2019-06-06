#!/usr/bin/env ruby
#
# THIS IS JUST A DRAFT, DON'T USE IN ANGER
#

# YAML file can be created from a TSV, as example from Table Maker
# perl -anE 'BEGIN{use YAML}$h{$F[0]}{name}=$F[0];$h{$F[0]}{paper}{$F[1]} = $F[2]||"n/a" if $F[1]; END{say Dump(\%h)}'
# looks like:
# ---
# WBVar00000001:
#   name: WBVar00000001
#   paper:
#     WBPaper00000932: 2428682
#     WBPaper00014052: n/a
#
# usage: ruby make_agr_variations.rb <filter.yaml> <genomic fasta> <wormbase gff3>

require 'bio'
require 'date'
require 'json'
require 'yaml'

term2so = {
	'insertion_site' => 'SO:0000667', # changed insertion_site to insertion
	'deletion'       => 'SO:0000159',
	'point_mutation' => 'SO:1000008',
	'substitution'   => 'SO:1000008', # to hack our wrong SO-terms in
}

chrom2ncbi = {
	'I' => 'RefSeq:NC_003279.8',
	'II' => 'RefSeq:NC_003280.10',
	'III' => 'RefSeq:NC_003281.10',
	'IV' => 'RefSeq:NC_003282.8',
	'V' => 'RefSeq:NC_003283.11',
	'X' => 'RefSeq:NC_003284.9',
	'MtDNA' => 'RefSeq:NC_001328.1',
}

variations = Array.new
chromosomes= Hash.new

filter = YAML.load_file(ARGV.shift)

Bio::FlatFile.open(Bio::FastaFormat, ARGV.shift).each_entry{|e|
   chromosomes[e.entry_id.to_s] = e.seq
}

ARGF.each_line{|line|
  variation = Hash.new
  cols = line.split("\t")
  next unless term2so[cols[2]]

  variation[:alleleId] = cols[8][/variation=(WBVar\d+)/,1]

  next unless ['Allele','KO_consortium','NBP_knockout','Million_mutation'].include?(cols[1])
  next unless filter.has_key?(variation[:alleleId])

  variation[:start] = cols[3].to_i
  variation[:end] = cols[4].to_i
  variation[:chromosome] = cols[0]
  variation[:assembly] = 'WBcel235'
  variation[:type]=term2so[cols[2]]

  if cols[2].eql?('insertion_site')
	  next unless cols[8]=~/insertion=([^;]+)/
          variation[:paddedBase] = chromosomes[cols[0]][variation[:start]-2]
	  variation[:genomicReferenceSequence]='N/A'
	  s = Bio::Sequence::NA.new($1.to_s)
	  s.complement! if cols[6].eql?('-') # as all variations are supposed to be on the forward strand
	  variation[:genomicVariantSequence] = s.to_s.upcase
	  variation[:end]+=1 if variation[:end] == variation[:start] # to make it inline with the HGVS coordinates

  elsif cols[2].eql?('deletion')
          variation[:paddedBase] = chromosomes[cols[0]][variation[:start]-2]
	  variation[:genomicReferenceSequence] = chromosomes[cols[0]].subseq(variation[:start],variation[:end])
          variation[:genomicVariantSequence] = 'N/A'

  elsif ['point_mutation','substitution'].include?(cols[2]) # multi-bp substitutions need to become DELINS in the future
	  next unless cols[8]=~/substitution=/ # to skip the crispr/cas9 alleles
	  next unless variation[:start] == variation[:end] # only allow single bp "substitutions"
	  to = cols[8][/substitution=([^;]+)/,1].split('/')[1]
	  variation[:genomicReferenceSequence] = chromosomes[cols[0]].subseq(variation[:start],variation[:end])
	  s = Bio::Sequence::NA.new(to)
	  s.complement! if cols[6].eql?('-') # as all variations are supposed to be on the forward strand
	  variation[:genomicVariantSequence] = s.to_s.upcase
  end
  
  variation[:references] = filter[variation[:alleleId]]["paper"].map{|k,v| # creates the publication part 
	  v.eql?('n/a')? {:publicationId => "WB:#{k}",:crossReference => {:id => "WB:#{k}",:pages => ['reference']}}:{:publicationId => "PMID:#{v}",:crossReference => {:id => "WB:#{k}",:pages => ['reference']}}
  } if filter[variation[:alleleId]]["paper"]

  variation[:sequenceOfReferenceAccessionNumber]=chrom2ncbi[variation[:chromosome]]
  variation[:alleleId]=variation[:alleleId].prepend('WB:') # prefix it with WB:

  variations.push(variation)
}

meta = {
	:dataProvider => {
		:crossReference => {
			:id => 'WormBase',
			:pages => ['homepage'] 
		},
                :type => 'curated'
	},
	:dateProduced => DateTime.now.to_s,
	:release => ENV['WS_RELEASE']
}

bleh =  {:data => variations, :metaData => meta}

puts JSON.pretty_generate(bleh)

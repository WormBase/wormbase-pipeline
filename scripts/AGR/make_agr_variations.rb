#!/usr/bin/env ruby
#
# THIS IS JUST A DRAFT, DON'T USE IN ANGER
#

# usage: ruby make_agr_variations.rb <genomic fasta> <wormbase gff3>
#

require 'bio'
require 'date'
require 'json'

term2so = {
	'insertion_site' => 'SO:0000667', # changed insertion_site to insertion
	'deletion'       => 'SO:0000159',
	'point_mutation' => 'SO:1000008',
}

variations = Array.new
chromosomes=Hash.new

Bio::FlatFile.open(Bio::FastaFormat, ARGV.shift).each_entry{|e|
   chromosomes[e.entry_id.to_s] = e.seq
}

ARGF.each_line{|line|
  variation = Hash.new
  cols = line.split("\t")
  next unless term2so[cols[2]]

  variation[:alleleId] = cols[8][/variation=(WBVar\d+)/,1]
  variation[:start] = cols[3]
  variation[:end] = cols[4]
  variation[:chromosome] = cols[0]
  variation[:assembly] = 'WBCel235'
  variation[:sotype]=term2so[cols[2]]
  if cols[8]=~/insertion=([^;])/
          variation[:paddedBase] = chromosomes[cols[0]][cols[3].to_i-2]
	  variation[:genomicReferenceSequence]='N/A'
	  s = Bio::Sequence::NA.new($1.to_s)
	  s.complement! if cols[6].eql?('-') # as all variations are supposed to be on the forward strand
	  variation[:genomicVariantSequence] = s.to_s.upcase
  elsif cols[2].eql?('deletion')
          variation[:paddedBase] = chromosomes[cols[0]][cols[3].to_i-2]
	  variation[:genomicReferenceSequence] = chromosomes[cols[0]].subseq(cols[3].to_i,cols[4].to_i)
          variation[:genomicVariantSequence] = 'N/A'
  elsif cols[2].eql?('point_mutation')
	  next unless cols[8]=~/substitution=/ # to skip the crispr/cas9 alleles
	  to = cols[8][/substitution=([^;]+)/,1].split('/')[1]
	  variation[:genomicReferenceSequence] = chromosomes[cols[0]].subseq(cols[3].to_i,cols[4].to_i)
	  s = Bio::Sequence::NA.new(to)
	  s.complement! if cols[6].eql?('-') # as all variations are supposed to be on the forward strand
	  variation[:genomicVariantSequence] = s.to_s.upcase
  end
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
	:dateProduced => DateTime.now.to_s
}

bleh =  {:data => variations, :metaData => meta}

puts JSON.pretty_generate(bleh)

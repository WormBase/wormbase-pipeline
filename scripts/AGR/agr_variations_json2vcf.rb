#!/usr/bin/env ruby
# based on 
# 	* VCF v4.3 (http://samtools.github.io/hts-specs/VCFv4.3.pdf)
# 	* AGR JSON 1.0.0.8 (https://github.com/alliance-genome/agr_schemas/blob/release-1.0.0.8/)

require 'json'
require 'date'

puts '##fileformat=VCFv4.3'
puts "##fileDate=#{Date.today.strftime("%Y%m%d")}"
puts '##source=AllianceJSON'
puts '##INFO<ID=DP,Number=1,Type=Integer,Description="Total Depth">'
puts ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'].join("\t")

# payload
parsed = JSON.parse(ARGF.read,symbolize_names: true)
parsed[:data].sort_by{|a|[a[:chromosome],a[:start]]}.each{|v|

	refSeq = v[:genomicReferenceSequence].eql?('N/A') ? '':v[:genomicReferenceSequence]
	varSeq = v[:genomicVariantSequence].eql?('N/A') ? '':v[:genomicVariantSequence]
	pos = v[:start]

	if v.has_key?(:paddedBase)
		if pos == 1 # special case for end positions
	           refSeq = refSeq.append(v[:paddedBase])
        	   varSeq = varSeq.append(v[:paddedBase])
		else
	           refSeq = refSeq.prepend(v[:paddedBase])
        	   varSeq = varSeq.prepend(v[:paddedBase])
		   pos = pos -1 # POS should include the padding base
		end
	end

	# the DP=100 is to fill the INFO column, as else some parsers break
	puts [v[:chromosome],v[:start],v[:alleleId],refSeq,varSeq,'.','PASS','DP=100'].join("\t")
}

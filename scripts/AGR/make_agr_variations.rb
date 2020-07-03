#!/usr/bin/env ruby

require 'bio'
require 'date'
require 'json'
require 'optparse'

Options = Struct.new(:fasta,:gff,:db,:ws_version,:wquery)

class Parser
	def self.parse(options)
		args = Options.new("world")

		args.wquery = "agr_variations.def" # default

		opt_parser = OptionParser.new do |opts|

			opts.banner = "Usage: make_agr_variations.rb [options]"

			opts.on("-f","--file NAME","FASTA file") do |f|
				args.fasta = f
			end

			opts.on("-g","--gff NAME","GFF file") do |g|
				args.gff = g
			end

			opts.on("-d","--database NAME","ACeDB database directory") do |d|
				args.db = d
			end

			opts.on("-w","--ws_version NAME","WormBase release (a.e. WS277)") do |v|
				args.ws_version = v
			end

			opts.on("-q","--wquery [DEFFILE]","Table-Maker DEF file") do |t|
				args.wquery = t
			end

			opts.on("-h","--help","print help") do
				puts opts
				exit
			end
		end

		opt_parser.parse!(options)
		return args
	end
end

class TableMaker
	def initialize(giface,database)
		@giface = giface
		@database = database
	end

	# as example query="query find Variation *;Reference"
	def execute_wquery(query,wquery,outfile='/tmp/tablemaker.outfile')
		IO.popen("#{@giface} #{@database}","r+"){|pipe|
			unless query.nil? # will quey all variation
				pipe.puts query
				pipe.puts "Table-maker -active -o #{outfile} #{wquery}"
			else # query only the current keyset from the query
				pipe.puts "Table-maker -o #{outfile} #{wquery}"
			end
			pipe.puts "quit"
			Process.wait(pipe.pid)
		}

                return self.parse_output(outfile)
	end

	# that one is specific to that particular query, so probably not too generally useful
	def parse_output(outfile)
		results=Hash.new
		File.open(outfile).each{|line|
			line.chomp!
			line.gsub!('"','')
			c = line.split
			results[c[0]]=Hash.new # WBVarXXX
			results[c[0]]["name"] = c[0] # WBVarXXX
			if columns[1]
				results[c[0]]["paper"]=Hash.new # WBPaperXXX
				results[c[0]]["paper"][c[1]] = c[2] || 'n/a' # PubmedID
			end
		}
		File.unlink(outfile)
		return results
	end
end

term2so = {
	'insertion_site' => 'SO:0000667', # changed insertion_site to insertion
	'deletion'       => 'SO:0000159',
	'point_mutation' => 'SO:1000008',
	'substitution'   => 'SO:1000032', # to hack our wrong SO-terms in => DELIN
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

options = Parser.parse(ARGV)

variations = Array.new
chromosomes= Hash.new

# hardcoded giface ... which is probably not needed
tablemaker = TableMaker.new('/nfs/panda/ensemblgenomes/wormbase/software/packages/acedb/RHEL7/4.9.62/giface',options.db)
filter = tablemaker.execute_wquery("query find Variation *;Reference",options.wquery) 

Bio::FlatFile.open(Bio::FastaFormat, options.fasta).each_entry{|e|
   chromosomes[e.entry_id.to_s] = e.seq
}

# parse from the GFF the respective variation lines
File.open(options.gff).each{|line|
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

	  next unless variation[:start] == variation[:end] # only allow single bp delin "substitutions" due to AGR limitations

	  to = cols[8][/substitution=([^;]+)/,1].split('/')[1]
	  variation[:genomicReferenceSequence] = chromosomes[cols[0]].subseq(variation[:start],variation[:end])
	  s = Bio::Sequence::NA.new(to)
	  s.complement! if cols[6].eql?('-') # as all variations are supposed to be on the forward strand
	  variation[:genomicVariantSequence] = s.to_s.upcase

	  if to.length == 1 # WS276 hack for the missing term, so if from is 1bp and to is 1bp it is a point mutation
	          variation[:type]='SO:1000008' 	  
	  end
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
	:release => options.ws_version
}

bleh =  {:data => variations, :metaData => meta}

puts JSON.pretty_generate(bleh)

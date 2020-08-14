#!/usr/bin/env ruby

require 'bio'
require 'date'
require 'json'
require 'optparse'
require 'zlib'

Options = Struct.new(:fasta,:gff,:db,:ws_version,:wquery,:outfile,:all)


##
# = Parser Class
#
# to extend getopts and generate a autohelp

class Parser
	def self.parse(options)
		args = Options.new("world")

		args.wquery = "agr_variations.def" # default
                args.all = false

		opt_parser = OptionParser.new do |opts|

			opts.banner = "Usage: make_agr_variations.rb [options]"

			opts.on("-f","--fasta NAME","FASTA file") do |f|
				args.fasta = f
			end

			opts.on("-o","--outfile NAME","JSON output file") do |o|
				args.outfile = o
			end

			opts.on("-g","--gff NAME","gzipped GFF file") do |g|
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

                        opts.on("-a","--all","Retrieve all naturally occurring variants and induced mutations from Million Mutation Project") do |a|
                                args.all = a
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

##
# = TableMaker class
#
# connects to giface and runs tablemaker def files against a keyset

class TableMaker
	def initialize(giface,database)
		@giface = giface
		@database = database
	end

	##
	# as example query="query find Variation *;Reference"
	def execute_wquery(query,wquery,outfile='/tmp/tablemaker.outfile')
		IO.popen("#{@giface} #{@database}","r+"){|pipe|
			unless query.nil? # will query all variation
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

	##
	# that one is specific to that particular query, so probably not too generally useful
	def parse_output(outfile)
		results=Hash.new
		File.open(outfile).each{|line|
			line.chomp!
			line.gsub!('"','')
			c = line.split("\t")
                        if !results[c[0]]
			        results[c[0]]=Hash.new # WBVarXXX
			        results[c[0]]["name"] = c[0] # WBVarXXX
			if !c[1].empty?
				results[c[0]]["paper"]||=Hash.new # WBPaperXXX
				results[c[0]]["paper"][c[1]] = c[2] || 'n/a' # PubmedID
			end
			if c[3]
				results[c[0]][:strains]||=[] # WBStrains
				results[c[0]][:strains].push(c[3]) # adds WBStrainId
			end
		}
#		File.unlink(outfile)
		return results
	end
end

term2so = {
	'insertion_site' => 'SO:0000667', # changed insertion_site to insertion
	'deletion'       => 'SO:0000159',
	'point_mutation' => 'SO:1000008',
	'substitution'   => 'SO:1000032', # to hack our wrong SO-terms in => DELIN
        'SNP'            => 'SO:0000694',
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
# filter = tablemaker.execute_wquery("query find Variation *;Reference;SMap",options.wquery)
if options.all
  filter = tablemaker.execute_wquery('query find Variation WHERE Live AND (!Gene OR COUNT(Gene) < 3) AND SMap AND (Natural_variant OR Reference = WBPaper00042537)',options.wquery)
else
  filter = tablemaker.execute_wquery('query find Variation WHERE Live AND COUNT(Gene) == 1 AND Reference AND SMap AND (Phenotype OR Disease_info OR Interactor) AND NOT Natural_variant',options.wquery)
end

Bio::FlatFile.open(Bio::FastaFormat, options.fasta).each_entry{|e|
   chromosomes[e.entry_id.to_s] = e.seq
}

allele_included = Hash.new

# parse from the GFF the respective variation lines
Zlib::GzipReader.open(options.gff).each{|line|
  variation = Hash.new
  cols = line.split("\t")
  next unless term2so[cols[2]]

  variation[:alleleId] = cols[8][/variation=(WBVar\d+)/,1]
  next if allele_included.key?(variation[:alleleId]) # Deals with identical variations from different sources in GFF
  allele_included[variation[:alleleId]] = 1
  
  next unless options.all or ['Allele','KO_consortium','NBP_knockout','Million_mutation'].include?(cols[1])
  
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

  elsif ['point_mutation','substitution','SNP'].include?(cols[2]) # multi-bp substitutions need to become DELINS in the future
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

  if options.all # add strains only to the large JSON
	  variation[:strains]=filter[variation[:alleleId]][:strains]||['WBStrain00000001'] # defaults to N2
  end

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

File.write(options.outfile , JSON.pretty_generate(bleh))

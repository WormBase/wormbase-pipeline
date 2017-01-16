#!/usr/bin/env ruby
#== Synopsis
# eggNOG v4 parser to ACeDB
# needs a load of common_data files to map ids around
# * wormpep can be found in the usual place
# * the description, funccat and members at the eggNOG site
# * the flybase file to map uniprot to flybase gene names is at gadfly
#== Usage
# ruby eggnog.rb --memberfile euNOG.members.txt --commondata cds2wormpep243.txt --funcfile euNOG.funccat.txt --descriptionfile euNOG.description.txt --type euNOG --flyfile dmel-all-translation-r5.57.fasta

require 'optparse'

class Protein
	attr_accessor :name, :species, :dbxref
	def initialize(id)
		@name=id
		@dbxref=[]
	end
end

class Cluster
	attr_accessor :name, :members, :description, :type
	def initialize(id)
		@name=id
		@members=[]
		@type=[]
	end
end

# mapping between NCBI taxonomy ids and real names
tax2name = {
	'6238' => 'Caenorhabditis briggsae',
	'6239' => 'Caenorhabditis elegans',
	'6334' => 'Trichinella spiralis',
	'31234' => 'Caenorhabditis remanei',
	'54126' => 'Pristionchus pacificus',
	'135651' => 'Caenorhabditis brenneri',
	'281687' => 'Caenorhabditis japonica',
	'10090' => 'Mus musculus',
	'7955' => 'Danio rerio',
	'9606' => 'Homo sapiens',
	'7227' => 'Drosophila melanogaster',
	'4932' => 'Saccharomyces cerevisiae',
}

clusters=Hash.new

memfile=''
funcfile=''
descfile=''
type='nemNOG'
commonData='/homes/mh6/cds2wormpep243.txt'
flyfile='dmel-all-translation-r5.57.fasta'
opt=OptionParser.new
opt.on('-m','--memberfile INFILE'){|f|memfile=File.new(f,'r')}
opt.on('-f','--funcfile INFILE'){|f|funcfile=File.new(f,'r')}
opt.on('-d','--descfile INFILE'){|f|descfile=File.new(f,'r')}
opt.on('-c','--commondata INFILE'){|f|commonData=f}
opt.on('-b','--flyfile INFILE'){|f|flyfile=f}
opt.on('-t','--type TYPE'){|t|type=t}

opt.parse(ARGV)

# slurp CDS -> Protein
# a.) common_data
def slurp_commondata(f)
	var1=Hash.new
	File.new(f,'r').each_line{|l|
		l.chomp!
		l.gsub!('"','')
		cols= l.split(/\t/)
		var1[cols[0]]=cols[1]
	}
	return var1
end

cds2wormpep=slurp_commondata(commonData)
# b.) flybase
def slurp_flybase(f)
	flypep2flygene=Hash.new
	File.new(f,'r').each_line{|l|
		match = />(\w+).*FlyBase_Annotation_IDs:(\w+)-\w+/.match(l)
		next unless match
		flypep2flygene[match[1].to_s]=match[2].to_s
		# puts "#{match[1].to_s} #{match[2].to_s}"
	}
	return flypep2flygene
end
fp2fg=slurp_flybase(flyfile)

# slurp functions
funcfile.each_line{|l|
	cols = l.split(/\t/)
	clusters[cols[0]]||=Cluster.new(cols[0])
	clusters[cols[0]].type=cols[1].scan(/./)
}
funcfile.close

# slurp descriptions
descfile.each_line{|l|
	l.chomp!
	c = l.split(/\t/)
	clusters[c[0]]||=Cluster.new(c[0])
	clusters[c[0]].description=c[1] unless c[1].nil?
}
descfile.close

# slurp members
memfile.each_line{|line|
        next if line=~/^#/;
	cols = line.split(/\t/)
        hits=/(\d+)\.(\S+)/.match(cols[1])

	taxId=hits[1].to_s
	protId=hits[2].to_s
        next unless tax2name[taxId]
	protein = Protein.new(cds2wormpep[protId]||fp2fg[protId]||protId)

	protein.dbxref << "UniProt UniProtAcc #{protId}" 
	protein.dbxref << "EnsEMBL EnsEMBL protein_id #{protId}" unless cds2wormpep[protId]
	protein.dbxref << "EggNog cluster #{cols[0]}"

	# fly
	if fp2fg[protId]
	 protein.dbxref << "FlyBase FlyBase_ID #{protein.name}" if fp2fg[protId]
	 protein.name="FLYBASE:#{protein.name}" if fp2fg[protId]
	end

	# yeast
	if tax2name[taxId].eql?('Saccharomyces cerevisiae')
	  protein.dbxref << "SGD SGD_systematic_name #{protein.name}"
	  protein.name="SGD:#{protein.name}"
	end

	if ['Mus musculus','Danio rerio','Homo sapiens'].include?(tax2name[taxId])
	  protein.name="ENSEMBL:#{protein.name}"
	end

	raise "cannot find species for #{taxId} #{cols[1]}" unless tax2name[taxId]
	protein.species=tax2name[taxId]

        clusters[cols[0]]||=Cluster.new(cols[0]);
        clusters[cols[0]].members << protein;
}
memfile.close

# generic print
clusters.each{|k,v|
        next unless v.members.size > 0
	puts "Homology_group : \"#{k}\""
	puts "Title \"#{v.description}\"" if v.description
	puts "DB_info Database EggNog cluster #{k}"
	puts "Group_type eggNOG eggNOG_type #{type}"
	v.type.each{|t|
		puts "Group_type eggNOG eggNOG_code Code_#{t}"
	}
	v.members.each{|m|
		puts "Protein #{m.name}"
	}
	puts ""
	v.members.each{|m|
		puts "Protein : #{m.name}"
		puts "Species \"#{m.species}\""
		m.dbxref.each{|x|
			puts "Database #{x}"
		}
		puts ""
	}
}


#/usr/bin/env ruby
# small Inparanoid 6 XML parser to generate Tier I /II mappings to H.sapiens (through Ortholog_other)

require 'rubygems'
require 'yaml'
require 'hpricot'

def read_common_data(file)
  eval File.new("/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/#{file}",'r').read

  return $VAR1
end

@@seq2gene = read_common_data('worm_gene2geneID_name.dat')
@@gen2seq  = @@seq2gene.invert

def fix_name(name)
  puts "fixing #{name}" if $DEBUG
  if (@@gen2seq[name])
    return name
  elsif (@@seq2gene[name])
    return @@seq2gene[name]
  elsif (@@seq2gene[name+'.1'])
	  return @@seq2gene[name+'.1']
  else
    return nil
  end
end

def xml2ace(parent,childs,child_organism)
    parent.each{|p|
      pname=fix_name(p.attributes['geneid'])
      next if pname.nil?
      puts "Gene : \"#{pname}\""
      childs.each{|c|
        cname=c.attributes['protid'] # remanei hack
        next if cname.nil?
        puts "Ortholog_other ENSEMBL:#{cname} From_analysis Inparanoid_6"
        }
      puts ""
      childs.each{|c|
        cname=c.attributes['protid'] # remanei hack
        next if cname.nil?
        puts "Protein ENSEMBL:#{cname}"
	puts "Species \"#{child_organism}\""
	puts "DB_info EnsEMBL ENSEMBL_proteinID #{cname}"
	puts ''
      }
    }
end

doc = open(ARGV.shift) { |f| Hpricot(f) }
clusters=doc.search('cluster')

clusters.each{|c|
	bgenes = c.search('//gene[@species=ensHOMSA]')
#	egenes = c.search('//gene[@species=modCAEBR]')
 	egenes = c.search('//gene[@species=modCAEEL]')
 
#  xml2ace(bgenes,egenes,'Caenorhabditis elegans') # parent/child
  xml2ace(egenes,bgenes,'Homo sapiens')
}


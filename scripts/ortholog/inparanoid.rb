#/usr/bin/env ruby
# as they switched in version 6 to XML files,
# a quick XML parser combined with the common_data for creating Tier I / Tier II orthologs

require 'rubygems'
require 'yaml'
require 'hpricot'

def read_common_data(file)
#  eval File.new("/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/#{file}",'r').read
  eval File.new("/nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/#{file}",'r').read

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
        cname=fix_name(c.attributes['geneid']) # remanei hack
        next if cname.nil?
        puts "Ortholog #{cname} \"#{child_organism}\" From_analysis Inparanoid_6"
        }
      puts ""
    }
end

doc = open(ARGV.shift) { |f| Hpricot(f) }
clusters=doc.search('cluster')

clusters.each{|c|
	bgenes = c.search('//gene[@species=modCAERE]')
#	egenes = c.search('//gene[@species=modCAEBR]')
 	egenes = c.search('//gene[@species=modCAEEL]')
 
  xml2ace(bgenes,egenes,'Caenorhabditis elegans') # parent/child
  xml2ace(egenes,bgenes,'Caenorhabditis remanei')
}


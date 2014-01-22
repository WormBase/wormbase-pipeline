#/usr/bin/env ruby
# as they switched in version 6 to XML files,
# a quick XML parser combined with the common_data for creating Tier I / Tier II orthologs

require 'rubygems'
require 'yaml'
require 'hpricot'

# Hashes
@@gene2species = Hash.new
@@gene2id      = Hash.new

def read_common_data(file,species="autoace")
  eval File.new("/homes/mh6/project/BUILD/#{species}/COMMON_DATA/#{file}",'r').read
  return $VAR1
end

@@seq2gene = read_common_data('cgc_name2gene.dat')
@@seq2gene.merge!(read_common_data('worm_gene2geneID_name.dat'))
@@seq2gene.merge!(read_common_data('worm_gene2geneID_name.dat','remanei'))
@@seq2gene.merge!(read_common_data('worm_gene2geneID_name.dat','briggsae'))
@@seq2gene.merge!(read_common_data('worm_gene2geneID_name.dat','brenneri'))
@@seq2gene.merge!(read_common_data('worm_gene2geneID_name.dat','japonica'))

@@gen2seq  = @@seq2gene.invert

def fix_name(name)
  puts "fixing #{name}" if $DEBUG
  if (@@gen2seq[name])
    return name
  elsif (@@seq2gene[name])
    return @@seq2gene[name]
  elsif (@@seq2gene[name+'.1'])
	  return @@seq2gene[name+'.1']
  elsif (@@seq2gene[name+'a'])
	  return @@seq2gene[name+'a']
  elsif name=~/WBGene\d+/
    return name
  else
#    puts "//barfed on #{name}"
    return nil
  end
end

def xml2ace(a,b)
    a.each{|p|
      pname=fix_name(@@gene2id[p.attributes['id']])
      next if pname.nil?
      puts "Gene : \"#{pname}\""
      b.each{|c|
        cname=fix_name(@@gene2id[c.attributes['id']])
        next if cname.nil?
        puts "Ortholog #{cname} \"#{@@gene2species[c.attributes['id']]}\" From_analysis Inparanoid_8"
      }
      a.each{|paralog|
         next if paralog.attributes['id'].eql?p.attributes['id']
         cname=fix_name(@@gene2id[paralog.attributes['id']])
         next if cname.nil?
         puts "Paralog #{cname} \"#{@@gene2species[paralog.attributes['id']]}\" From_analysis Inparanoid_8"
      }
      puts ""
    }
end

doc = open(ARGV.shift) { |f| Hpricot.XML(f) }



species=doc.search('species')
species.each{|s|
     spec = s.attributes['name'].chomp(" ")
     genes = s.search('//gene')
     genes.each{|g|
      @@gene2species[g.attributes['id']]=spec
      geneId=g.attributes['geneId']
      if spec.eql?('Caenorhabditis briggsae')
           unless geneId=~/Cbr-/
            geneId="Cbr-#{geneId}" if geneId=~/-/
           end
      elsif spec.eql?('Caenorhabditis remanei')
           unless geneId=~/Cre-/
            geneId="Cre-#{geneId}" if geneId=~/-/
           end
      elsif spec.eql?('Caenorhabditis brenneri')
           unless geneId=~/Cbn-/
            geneId="Cbn-#{geneId}" if geneId=~/-/
           end
      elsif spec.eql?('Caenorhabditis japonica')
           unless geneId=~/Cjp-/
            geneId="Cjp-#{geneId}" if geneId=~/-/
           end
      elsif spec.eql?('Pristionchus pacificus')
           unless geneId=~/Ppa-/
            geneId="Ppa-#{geneId}" if geneId=~/-/
           end
      end

      geneId.sub!('CELE_','')
      geneId.sub!('CRE_','CRE')
      geneId.sub!('CAEBREN_','CBN')
      geneId.sub!(/\/\S+/,'')
      @@gene2id[g.attributes['id']]=geneId
     }
}

doc.search('orthologGroup').each{|c|
	genes = c.search('//geneRef')
       # split into C.elegans

       # mapper
       species2genes= Hash.new       
 
       genes.each{|g|
         species2genes[@@gene2species[g.attributes['id']]]||=[]
         species2genes[@@gene2species[g.attributes['id']]] << g
       }
       
       v=species2genes.values
       xml2ace(v[0],v[1])
       xml2ace(v[1],v[0])
}


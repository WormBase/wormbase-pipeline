#/usr/bin/env ruby
# to parse the Orthology Matrix files from the ETH
# will also map it to a.) current_db and b.) cb25 (so a mapping file is available)

require 'yaml'

def read_common_data(file)
  eval File.new("/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/#{file}",'r').read

  return $VAR1
end

def read_flatfile
  ensembl2seq=Hash.new
  File.new('cb25ensembl_ids.lst','r').each_line{|line|
    line.chomp!
    (ens,seq)=line.split
    ensembl2seq[ens]=seq
    }
    return ensembl2seq
end

@@seq2gene = read_common_data('worm_gene2geneID_name.dat')

puts YAML::dump(@@seq2gene) if $DEBUG

#@@ens2seq  = read_flatfile

def fix_name(name)
  if (@@seq2gene[name])
    return @@seq2gene[name]
  else
    return nil
  end
end

def h2ace(parent,child_organism)
    parent.each{|k,v|
      puts "Gene : \"#{k}\""
      v.each{|c|
        puts "Ortholog #{c} \"#{child_organism}\" From_analysis OMA"
        }
      puts ""
    }
end

egenes=Hash.new
bgenes=Hash.new
ARGF.each{|line|
  line.chomp!
  line.gsub!(/\s/,'')
  ids=line.split(',')
  egene=fix_name(ids[0])
  bgene=fix_name(ids[2])
  puts "#{ids[0]} => #{egene} : #{ids[2]} => #{bgene}" if $DEBUG
  next if (egene.nil? || bgene.nil?)
  egenes[egene]||=[]
  egenes[egene] << bgene 
  bgenes[bgene]||=[]
  bgenes[bgene] << egene
  }
h2ace(egenes,'Caenorhabditis remanei')
h2ace(bgenes,'Caenorhabditis elegans')

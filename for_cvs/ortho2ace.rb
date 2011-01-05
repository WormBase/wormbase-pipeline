#!/usr/bin/env ruby
#
# merge ranjana's list and the EnsMart export into a ace file for testing
require 'rubygems'
require 'oniguruma'

# get omim -> ensembl gene id from file
def get_ensmart_hash
  omim2ensp=Hash.new
  infile=File.new('mart_export2.txt','r')
  infile.each_line{|line|
    line.chomp!
    match_data=line.split(" ")
    omim2ensp[match_data[1]]=match_data[0]
    }
  return omim2ensp
end
  
# iterator for the email from ranjana
def process_omims
    omim2ensp=get_ensmart_hash
    ensp2omim=omim2ensp.invert
    infile=File.new('ranjana','r')
    reg = Oniguruma::ORegexp.new('(WBGene\d+).*OMIM:\s((\d+\s)+)')
    infile.each_line{|line|
      match_data=reg.match(line)
      #puts line,match_data.inspect
      next if match_data.nil?
      wbgene=match_data[1].to_s
      ensembl_proteins=match_data[2].chop.split(" ").map{|o|
        $stderr.puts "cannot find EnsEMBL id for WormBase #{wbgene} -> OMIM #{o}" if omim2ensp[o].nil?
        omim2ensp[o]
        }
      yield [wbgene,ensembl_proteins,ensp2omim]
    }
end

process_omims{|wbgene,ensembl_proteins,o|
  puts "Gene : \"#{wbgene}\""
  ensembl_proteins.each {|ep|
    puts "Ortholog_other ENSEMBL:#{ep} Inferred_automatically omim_orthologs" unless ep.nil?
    }
    puts ""
    ensembl_proteins.each {|ep|
      puts "Protein : ENSEMBL:#{ep}"
      puts 'Species "Homo sapiens"'
      puts "DB_info Database OMIM gene #{o[ep]}"
      puts "DB_info Database EnsEMBL ENSEMBL_proteinID #{ep}",''
      }    
}

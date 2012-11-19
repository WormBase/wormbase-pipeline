#!/usr/bin/env ruby
# this script takes Erich's OrthoMCL output,
# combines it with mappings from the current_DB/COMMON_DATA
# and creates ACE output for it

def _print_ace(a,b,species)
                ace=String.new
                a.each{|_a|
                        next if b.nil? || _a.nil? ||b.size<1
                        ace << "Gene : \"#{_a}\"\n"
                        b.each{|_b|
                                ace << "Ortholog #{_b} \"#{species}\" Inferred_automatically OrthoMCL\n"
                        }
                        ace << "\n"
                }
                return ace
end

# that is just a way to hold the clusters
class Orthomcl
        def initialize(id)
                @briggsae=Array.new
                @elegans=Array.new
		@remanei=Array.new
                @inp_id=id
        end



        def to_ace
                ace_line =  _print_ace(@elegans,@briggsae,'Caenorhabditis briggsae')
                ace_line << _print_ace(@briggsae,@elegans,'Caenorhabditis elegans')
                return ace_line
        end

        attr_accessor :briggsae , :elegans, :inp_id, :remanei
end

eval File.new('/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/cds2wbgene_id.dat','r').read
cds2gene=$VAR1

ARGF.each{|line|
  c=line.split("\t")
  genes=c[1].split

  # cluster_id ... actually never used except as object_id
  mcl_id=c[0].split('(')[0]
  ortho=Orthomcl.new(mcl_id)
  genes.each{|gene|

      # is like gene or sequence id (organism)
      /(\w+)\((\w+)\)/.match(gene)
      gene_id=$1.to_s
      next if gene_id.nil?
      organism_id=$2.to_s

      if organism_id.eql?('elegans')
        ortho.elegans << gene_id
      elsif organism_id.eql?('briggsae')
	next if cds2gene[gene_id].nil?
        ortho.briggsae << cds2gene[gene_id]
      end
    }
    puts ortho.to_ace
  }

#!/usr/bin/env ruby

# shorthand to print ace format
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

# to store OrthoMCL results
class Orthomcl

        attr_accessor :acc, :inp_id

        # creates a structure with a primary id and a hash of arrays to keep the data
        def initialize(id)
                @briggsae=Array.new
                @elegans=Array.new
                @remanei=Array.new
                @brenneri=Array.new
                @japonica=Array.new
                @inp_id=id
                @acc={'elegans' => @elegans,'brenneri'=>@brenneri,'remanei' => @remanei,'briggsae' => @briggsae,'japonica' => @japonica}
        end

        # MxN minus self loop to print the acefile
        def to_ace
                ace_line=''
                @acc.each{|k,v|
                    @acc.each{|k2,v2|
                      next if k.eql?(k2)
                      ace_line <<  _print_ace(v,v2,"Caenorhabditis #{k2}")
                    }
                }
                return ace_line
        end

end


# main part, to read the raw OrthoMCL file and parse the line
ARGF.each{|line|
  c=line.split("\t")
  genes=c[1].split
  mcl_id=c[0].split('(')[0]
  ortho=Orthomcl.new(mcl_id)
  genes.each{|gene|
      /(\w+)\((\w+)\)/.match(gene)
      gene_id=$1.to_s
      organism_id=$2.to_s
      ortho.acc[organism_id]<<gene_id
  }
  puts "//from cluster #{mcl_id}"
  puts ortho.to_ace
}
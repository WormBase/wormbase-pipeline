#!/usr/bin/env ruby
#== Synopsis
# dumps the longest translation and their Pfam motifs of a gene
#== Usage
# dump_treefam.rb -db DB_NAME -species "Full SPECIESNAME"
#
#== Author:
# Michael Han (mh6@sanger.ac.uk)
# Wellcome Trust Sanger Institute
# United Kingdom


require 'optparse'

require 'rubygems'
require 'rdoc/ri/ri_paths'
require 'rdoc/usage'
require_gem 'ensembl-api' # Jan Aert's EnsEMBL API
include Ensembl::Core

database='worm_ensembl_elegans'
species='Caenorhabdits elegans'

opt=OptionParser.new
opt.on("-d", "--database DATABASE",''){|d|database=d}
opt.on("-s","--species SPECIES",''){|s|species=s}
opt.parse(ARGV) rescue RDoc::usage('Usage')

# generate a species suffix from the two names
(spec_prefix,spec_suffix)=species.split
species_string=(spec_prefix[0..2].upcase) + (spec_suffix[0..1].upcase)

# connect to db -- hardcoded to ia64d, can always be changed if needed
CoreDBConnection.establish_connection(:adapter=>'mysql',:host => 'ia64d',:database =>database,:username=>'wormro',:password => '')

# get all genes ...
Gene.find_all().each{|g|
    $stderr.puts "processing gene #{g.stable_id}" if $DEBUG
    t=g.transcripts # get all transcripts for the gene g
    
    # sort them by size and pick the longest one
    t.sort!{|a,b|a.cds_seq.length <=> b.cds_seq.length} # else it sorts it alphabetically
    longest_transcript=t.pop

    # then take all Pfam features and create a nice text representation
    $stderr.puts "processing domains" if $DEBUG
    pfams=longest_transcript.translation.protein_features.map{|a| a if a.analysis.logic_name.eql?('Pfam')}.compact # removes non-Pfam features
    d_line=pfams.map{|p| "#{p.hit_id}:#{p.seq_start}..#{p.seq_end}"}.join('; ')

    # mapping coordinates
    strand=g.strand>0?'+':'-'
    m_line="#{g.seq_region.name}(#{strand}):#{longest_transcript.start}-#{longest_transcript.stop}"
    
    # and print it as Fasta file
    $stderr.puts "printing fasta" if $DEBUG
#   seq=Bio::Sequence::AA.new(longest_transcript.protein_seq.sub('*',''))
    seq=Bio::Sequence::NA.new(longest_transcript.cds_seq)

    puts seq.to_fasta("#{g.stable_id}\t#{g.stable_id}_#{species_string}\t#{longest_transcript.translation.stable_id}\t#{m_line}\t#{d_line}",60)
}

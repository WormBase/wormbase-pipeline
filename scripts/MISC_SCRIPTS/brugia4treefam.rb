#!/usr/bin/env ruby
# author: Michael Han (mh6@sanger.ac.uk)
# dumping script for TreeFam
# dumps the longest translation and their Pfam motifs of a gene


require 'rubygems'
require_gem 'ensembl-api' # Jan Aert's EnsEMBL API
include Ensembl::Core

# connect to db -- hardcoded to brugia and ia64b, can always be changed if needed
CoreDBConnection.establish_connection(:adapter=>'mysql',:host => 'ia64b',:database =>'worm_ensembl_brugia',:username=>'wormro',:password => '')

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

    # and print it as Fasta file
    $stderr.puts "printing fasta" if $DEBUG
    seq=Bio::Sequence::AA.new(longest_transcript.protein_seq.sub('*',''))
    puts seq.to_fasta("#{g.stable_id}\t#{g.stable_id}\t#{g.stable_id}\t#{d_line}",60)
    }

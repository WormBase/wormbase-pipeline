#!/usr/bin/env ruby
#== Synopsis
# checks the cdses for cases where a reversal of the CDS would solve in-frame stops
#== Usage
# check_reverse -db DB_NAME
#
#== Author:
# Michael Han (mh6@sanger.ac.uk)
# Wellcome Trust Sanger Institute
# United Kingdom
# 
#== Warnings:
# the default installed ensembl-api uses deprecated require_gem calls
# which were removd in the latest gem versions.
# it is possible to get a working version from Jan's github, but the gemspec doesn't work, so no gem magic

require 'optparse'

require 'rubygems'
require 'rdoc/ri/ri_paths'
require 'rdoc/usage'
# require_gem 'ensembl-api' # Jan Aert's EnsEMBL API
# gem 'ensembl-api'
require 'ensembl'
include Ensembl::Core

database='worm_ensembl_elegans'

opt=OptionParser.new
opt.on("-d", "--database DATABASE",''){|d|database=d}
opt.parse(ARGV) rescue RDoc::usage('Usage')

# connect to db -- hardcoded to ia64d, can always be changed if needed
CoreDBConnection.establish_connection(:adapter=>'mysql',:host => 'ia64d',:database =>database,:username=>'wormro',:password => '')

# get all genes ...
Gene.find(:all).each{|g|
   $stderr.puts "processing gene #{g.stable_id}" if $DEBUG
    t=g.transcripts # get all transcripts for the gene g
    
    t.each{|transcript|
    
      seq=Bio::Sequence::NA.new(transcript.cds_seq)
      fstops= seq.translate.scan('*').length
      rstops= seq.reverse_complement[1..-1].translate.scan('*').length
      if fstops>1 and rstops <=1
        puts "#{g.stable_id} needs to be reversed (#{fstops} forward stops / #{rstops} reverse (+1bp))"
      elsif fstops>1 and rstops >1
        puts "#{g.stable_id} is generically shite (#{fstops} forward stops / #{rstops} reverse (+1bp))"
      end
    }
}

#!/usr/bin/env ruby
#== Synopsis
# a script to create other_nematode ACeDB and Fasta files from an EMBL EST file
# including filtering of Nemagenetag and WashU entries
# uses a patched BioRuby (fixed in HEAD, broken in 1.0.0)
#
#== Usage
#   ruby process_other_ESTs.rb --embl|-e FILENAME -o|--outdir DIRECTORYNAME [-v|--verbose]
#
#== Author
# Michael Han (mh6@sanger.ac.uk)
# WellCome Trust Sanger Institute
# United Kingdom
#
#== Copyright
# as ruby

require 'optparse'
#require 'rdoc/ri/paths'
#require 'rdoc/usage'
require 'bio'

# to bring the bio rubygem version of embl.rb to HEAD
class EMBL

  # my patch for 1.0.0
  def id_line(key=nil)
    unless @data['ID']
      tmp = Hash.new
      idline = fetch('ID').split(/; +/)         
      tmp['ENTRY_NAME'], tmp['DATA_CLASS'] = idline.shift.split(/ +/)
      if idline.first =~ /^SV/
        tmp['SEQUENCE_VERSION'] = idline.shift.split(' ').last
        tmp['TOPOLOGY'] = idline.shift
        tmp['MOLECULE_TYPE'] = idline.shift
        tmp['DATA_CLASS'] = idline.shift
      else
        tmp['MOLECULE_TYPE'] = idline.shift
      end
      tmp['DIVISION'] = idline.shift
      tmp['SEQUENCE_LENGTH'] = idline.shift.strip.split(' ').first.to_i

      @data['ID'] = tmp
    end
    
    if key
      @data['ID'][key]
    else
      @data['ID']
    end
  end

  # and Mitsuteru's modification of my patch
  def sv
    if (v = field_fetch('SV').sub(/;/,'')) == ""
      [id_line['ENTRY_NAME'], id_line['SEQUENCE_VERSION']].join('.') 
    else
      v
    end  
  end
  def version
    (sv.split(".")[1] || id_line['SEQUENCE_VERSION']).to_i
  end
end


output_dir='/tmp/' #for debugging
infile = '/dev/null'
verbose=false

washu_ests_file = "/nfs/disk100/wormpub/BUILD_DATA/cDNA/washu/washu_contig_est_gbacc.full_list"
nembase_ests_file = "/nfs/disk100/wormpub/BUILD_DATA/cDNA/nembase/nembase_contig_est_gbacc.full_list" 

opt=OptionParser.new
opt.on('-e','--embl INFILE'){|f|infile=File.new(f,'r')}
opt.on('-o','--outdir OUTDIR'){|o|output_dir=o}
opt.on('-v','--verbose'){verbose=true}
opt.parse(ARGV) #rescue RDoc::usage('Usage')

$stderr.puts "using #{infile.path} , output to #{output_dir}" if verbose
$stderr.puts "Fetching other nematode EST sequences" if verbose

# Read in data about which Genbank/EMBL ESTs are used to construct the
# WashU NemaGene or Edinburgh NemBase nemotode EST contigs.
# We wish to exclude the WashU and NemBase set of ESTs from this Other Nematode
# set of ESTs
@@contig_est_gbacc = Hash.new # bugger. I need to check the global variables :-(

# get the EST ids of the ESTs to skip into a hash
def get_ids(filename)
 $stderr.puts "fetching ids from #{filename}"

 ests_file=File.new(filename,'r')
 ests_file.each_line {|line|
    line.chomp!
    next if /^;/.match(line)           # ignore comment lines
    @@contig_est_gbacc[line] = 1
 }
end

[washu_ests_file,nembase_ests_file].each{|id|get_ids(id)}


# open filehandles for output files
file="#{output_dir}/other_nematode_ESTs"
out_nem=File.new(file,'w')
out_ace=File.new("#{file}.ace",'w')
out_prob=File.new("#{file}.prob",'w')

# for the acefile
def pretty_seq(seq)
	s=seq.to_s
	s.gsub!(/(\w{60})/,"\\1\n")
	return s
end

# Grab all EST sequences (division EST) from the EMBL flatfile
sequences=Bio::FlatFile.auto(infile)
sequences.each {|entry|
    if (@@contig_est_gbacc.has_key?(entry.accession))
	$stderr.puts "skipping #{entry.accession}" if verbose
	next
    end
    
     
    $stderr.print "." if verbose
    entry.definition.gsub!('"','') # remove quotation marks for ACeDB
    out_nem.puts ">#{entry.accession}.#{entry.version} #{entry.entry_id} #{entry.definition}\n#{pretty_seq(entry.seq)}"
    out_ace.puts "Sequence : \"#{entry.accession}.#{entry.version}\""
    out_ace.puts "Database EMBL NDB_AC #{entry.accession}"
    out_ace.puts "Database EMBL NDB_ID #{entry.entry_id}"
    out_ace.puts "Database EMBL NDB_SV #{entry.version}"
    if defined? entry.os
       entry.os.each{|o|
       out_ace.puts "Species \"#{o['os']}\""
      }
    # else it is something else (not EMBL) so use default behaviour
    else
       out_ace.puts "Species \"#{entry.organism}\""
    end
    out_prob.puts("#{entry.accession} #{entry.entry_id} #{entry.definition} has #{entry.os.size} species entries") if (defined?(entry.os) && entry.os.size > 1)
    out_ace.puts "Title \"#{entry.definition}\"\nMethod EST_nematode"
    out_ace.puts "DNA \"#{entry.accession}\.#{entry.version}\"\n\n"
    out_ace.puts "DNA \"#{entry.accession}\.#{entry.version}\""
    out_ace.puts pretty_seq(entry.seq)
    out_ace.print "\n" 
}    
  
# close filehandles
out_nem.close
out_ace.close

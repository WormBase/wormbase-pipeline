#!/usr/local/bin/perl5.6.1 -w

# What it does: Download EMBL flatfile (complete entry) via getz
#               Ouput: EMBL flatfile with the tags: ID, AC, SV, KW, DE, OS
#                                    without DNA and protein sequences 

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2003-09-04 14:46:57 $ 

use strict;

use Getopt::Long;

my $version;
GetOptions ("version|v=s"  => \$version);

my $start = `date`;
my $download = "/wormsrv1/chaokung/EMBL/OUTPUT/embl";
my $output_dir = "/wormsrv1/chaokung/EMBL/OUTPUT";

print "Downloading current EMBL Release . . .\n";
`getz -e "((([embl-Organism:caenorhabditis*]&[embl-Division:inv]))>[embl-FtQualifier:gene]) > parent" > $download`;

print "\nDeleting protein and DNA sequences from EMBL flatfile...\n";
open(IN, $download) || die "Can't read in file!";

while(<IN>){
  chomp;

  my $output="S2_".$version."_EMBL_entry_gene_tag_no_seqs";
  open(OUT_embl, ">>$output_dir/$output") || die "Can't write to file!";
  if ($_ =~ /(^ID[\w\W.]+)|(^AC[\w\W.]+)|(^SV[\w\W.]+)|(^KW[\w\W.]+)|(^DE[\w\W.]+)|(^OS[\w\W.]+)|(^\/\/)/){
    print OUT_embl $_,"\n";
  } 
  if (($_ =~ /^FT[\w\W.]+/)&&($_ !~ /^FT\s+\/translation[\w\W.]+/)                                                               
      &&($_ !~ /^FT\s+[A-Z]+/)){
    print OUT_embl $_,"\n"; 
  }      
}
close IN;

system("rm -f $download");

my $end = `date`;

print "\n$0\nstarted at $start\n";
print "finished at $end\n";



__END__



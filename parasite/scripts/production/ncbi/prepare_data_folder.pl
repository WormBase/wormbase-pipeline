#!/usr/bin/env perl
use strict;
use warnings;
# Given a config in the from-NCBI format, create directory structure
use YAML;
use File::Path qw/make_path/;
use File::Basename;
use Net::FTP;

my $conf;
$conf = Load(do {local $/; <STDIN>} || "--- {}") unless -t STDIN;
die "Usage: PARASITE_DATA=... $0 <conf" unless $conf;

my $bp = $conf->{GB_BioProjects}{Bioproj}{BioprojectAccn};
my ($spe, $cies) = split (/\s+/, $conf->{SpeciesName});
my $isolate = $conf->{Biosource}{Isolate};

my $species = lc(join("_", $spe, $cies, $bp));
die $species unless $species =~ /^\w+_\w+_\w+\d+$/;

goto PRINT unless $ENV{PARASITE_DATA};

my $dir = join("/", $ENV{PARASITE_DATA}, $species);
make_path $dir;

my ($assembly_path) = glob("$dir/*genomic.fna");

if (not $assembly_path){
  my $ftp = Net::FTP->new("ftp.ncbi.nlm.nih.gov")
     or die "$@";
  (my $uri = $conf->{FtpPath_GenBank}) =~ s{ftp://ftp.ncbi.nlm.nih.gov}{};
  $ftp->login("anonymous",'-anonymous@') or die $ftp->message;
  $ftp->cwd($uri);
  my ($remote, @others) = grep {$_ !~ /from_genomic\.fna\.gz$/} grep {$_ =~ /genomic\.fna\.gz$/ } $ftp->ls or die $ftp->message;
  die "Multiple files on the FTP: ".  join("\t", $remote, @others) if @others;
  die $uri unless $remote;
  $assembly_path = "$dir/$remote";
  system("curl --silent ftp://ftp.ncbi.nlm.nih.gov/$uri/$remote > $assembly_path") and die "$@";
  system("gunzip $assembly_path") and die "Failed: gunzip $assembly_path";  
  $assembly_path =~ s/\.gz$//;
}
die unless -s $assembly_path;

my $assembly_destination_path = join("/", $dir, "$species.fa");
my $seq_region_synonyms_path = join("/", $dir, "$species.seq_region_synonyms.tsv");
if( not -s $assembly_destination_path or not -s $seq_region_synonyms_path){
   open(my $assembly_in_fh, "<", $assembly_path) or die "$!: $assembly_path";
   open(my $assembly_out_fh, ">", $assembly_destination_path) or die "$!: $assembly_destination_path";
   open(my $seq_region_synonyms_fh, ">", $seq_region_synonyms_path) or die "$!: $seq_region_synonyms_path";
   while(<$assembly_in_fh>){
      unless(/^>/){
         print $assembly_out_fh $_;
         next;
      }
      (my $l = $_ ) =~ s/^>//;
      chomp $l;
      my ($scaffold_name_ncbi) = $l =~ /(.*?) /;
      $l =~ s/^$scaffold_name_ncbi//;
      $l =~ s/, whole genome shotgun sequence$//;
      my ($scaffold_name_local) = $l =~ /.*? (\S+)$/;
      die "$scaffold_name_local | $_" if grep {$scaffold_name_local =~ /$_/i} ("$spe.$cies", $isolate);
      print $assembly_out_fh ">$scaffold_name_local\n";
      print $seq_region_synonyms_fh join("\t", "toplevel", $scaffold_name_local, $scaffold_name_ncbi, "INSDC")."\n"; 
   }
   close $_ for ($assembly_in_fh, $assembly_out_fh, $seq_region_synonyms_fh);
}

PRINT:

my %urls = (
  "University of Hull" => "https://www.hull.ac.uk",
  "WTSI" => "http://www.sanger.ac.uk",
  "University of Melbourne" => "http://www.gasserlab.org",
  "Oregon State University" => "https://oregonstate.edu",
);
print Dump({
  $species => {
     taxon_id => $conf->{SpeciesTaxid} // "?",
     assembly_version => $conf->{AssemblyName} // "?",
     gff_sources => "WormBase_imported",
     meta => {
        "assembly.accession" => $conf->{AssemblyAccession} // "?",
        "provider.name" => $conf->{SubmitterOrganization} // "?",
        "provider.url" => $urls{$conf->{SubmitterOrganization}} // "?",
        "species.biosample" => $conf->{BioSampleAccn} // "?",
        "species.nematode_clade" => $species =~ /meloidogyne|globodera|heterodera/ ? "IV" : $species =~ /pristionchus/ ? "V": "?",
     },
  }
});

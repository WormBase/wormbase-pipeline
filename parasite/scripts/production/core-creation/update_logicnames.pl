#! /usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $user, $port, $pass, $genome, $source, $logicname);

GetOptions(
  'host=s'      =>   \$host,
  'pass=s'      =>   \$pass,
  'port=i'      =>   \$port,
  'user=s'      =>   \$user,
  'genome=s'    =>   \$genome,
  'source=s'    =>   \$source,
  'logicname=s' =>   \$logicname
) || die ("check command line parameters\n") ;

my $usage = "$0 -host <host> \
        	-user <write user> \ 
		-port <port> \ 
		-pass <password> \
		-genome <genus_species_bioproject> \
		-source <source of genes/transcripts to be updated (as in gene/transcript table of core db)> 
		-logicname <new logic name>";

foreach my $param ($host, $user, $port, $pass, $genome, $source, $logicname){
  die "Parameters not defined:\n$usage" if not defined $param;
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  -host => $host,
  -user => $user,
  -port => $port,
  -pass => $pass
);

# get adaptors 
my $gene_adaptor = $registry->get_adaptor($genome, "core", "Gene");
my $transcript_adaptor = $registry->get_adaptor($genome, "core", "Transcript");
my $analysis_adaptor = $registry->get_adaptor($genome, "core", "Analysis");

# add the new analysis, if it doesn't already exist
my $analysis = Bio::EnsEMBL::Analysis->new(
  -logic_name => $logicname,
  -gff_source => $source,
  -module     => "WormBase"
   );

$analysis_adaptor->store($analysis);

# or update it
$analysis_adaptor->update($analysis);

# fetch all genes and transcript of our source of interest

my $genes = $gene_adaptor->fetch_all_by_source($source);
die "No genes found with source $source" if scalar(@$genes) == 0;
my $transcripts = $transcript_adaptor->fetch_all_by_source($source);
die "No transcripts found with source $source" if scalar(@$transcripts) == 0;

print "Found ".scalar(@$genes)." genes and ".scalar(@$transcripts)." transcripts to update\n";

# update analysis descs
foreach my $gene (@$genes){
  $gene->analysis($analysis);
  $gene_adaptor->update($gene);
}

foreach my $transcript(@$transcripts){
  $transcript->analysis($analysis);
  $transcript_adaptor->update($transcript);
}

print "Done.\n";

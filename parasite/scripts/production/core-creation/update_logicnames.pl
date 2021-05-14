#! /usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $user, $port, $pass, $genome, $source, $biotype, $logicname);

GetOptions(
  'host=s'      =>   \$host,
  'pass=s'      =>   \$pass,
  'port=i'      =>   \$port,
  'user=s'      =>   \$user,
  'genome=s'    =>   \$genome,
  'source=s'    =>   \$source,
  'biotype=s'   =>   \$biotype,
  'logicname=s' =>   \$logicname
) || die ("check command line parameters\n") ;

my $usage = "$0 -host <host> \
        	-user <write user> \ 
		-port <port> \ 
		-pass <password> \
		-genome <genus_species_bioproject> \
		-source <source of genes/transcripts to be updated (as in gene/transcript table of core db)> 
                -biotype <biotype of genes/transcripts to be updated (as in gene/transcript table of core db)> 
		-logicname <new logic name>";

foreach my $param ($host, $user, $port, $pass, $genome, $source, $biotype, $logicname){
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

# fetch all genes and transcripts of our source and biotype of interest
# update their analysis IDs

my $genes = $gene_adaptor->fetch_all_by_source($source);
my $updated_genes = 0;
foreach my $gene (@$genes){
  if ($gene->biotype() eq $biotype){
    $gene->analysis($analysis);
    $gene_adaptor->update($gene);
    $updated_genes++;
  }
}

print "Updated $updated_genes genes to new logic name $logicname\n";

my $transcripts = $transcript_adaptor->fetch_all_by_source($source);
my $updated_transcripts = 0;
foreach my $transcript(@$transcripts){
  if ($transcript->biotype() eq $biotype){
    $transcript->analysis($analysis);
    $transcript_adaptor->update($transcript);
    $updated_transcripts++;
  }
}

print "Updated $updated_transcripts transcripts to new logic name $logicname\n";

print "Done.\n";

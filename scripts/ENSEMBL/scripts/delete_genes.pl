#!/usr/local/ensembl/bin/perl -w
=head1 NAME

  delete_genes.pl

=head1 SYNOPSIS
 
  delete_genes.pl
  deletes genes from given database whose ids are passed in through a file

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbport    For RDBs, what port to connect to (port= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -stable_id ids are stable_ids instead of gene_id
    -infile    file with the ids (one per line)
    -help      summary of options


=head2 EXAMPLES

./delete_genes.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** -dbname rat_3_4_test -infile del.txt -stable_id

=cut

use strict; 
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;



my $dbhost = 'ia64d';
my $dbuser = 'wormadmin';
my $dbpass = 'XXX';
my $dbport = 3306;
my $dbname = 'worm_ensembl_cangaria';
my $infile;
my $stable_id;

&GetOptions( 
	    'dbname=s'     => \$dbname,
	    'dbuser=s'     => \$dbuser,
	    'dbpass=s'     => \$dbpass,
            'dbhost=s'     => \$dbhost,
            'dbport=s'     => \$dbport,
            'infile=s'     => \$infile,
	    'stable_id'    => \$stable_id,
           );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $dbhost,
					    -user => $dbuser,
					    -dbname => $dbname,
					    -pass  => $dbpass,
					    -port => $dbport,
					   );

open (INFILE, '<',$infile) or die "Cannot open $infile";

my @gene_ids;

while (<INFILE>){
  chomp $_;
  my $gene_id = $_;

  push @gene_ids, $gene_id;
}

my $gene_adaptor = $db->get_GeneAdaptor;

foreach my $gene_id (@gene_ids){

  my $gene;

  if ($stable_id){
    $gene = $gene_adaptor->fetch_by_stable_id($gene_id);
  }else{
    $gene = $gene_adaptor->fetch_by_dbID($gene_id);
  }

  print "GENE_ID: $gene_id\n";
  unless ($gene){
	  print "skipping: $gene_id\n";
	  next;
  }
  $gene_adaptor->remove($gene);
}


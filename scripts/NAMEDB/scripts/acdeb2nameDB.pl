#!/usr/local/bin/perl -w


use strict;
use Getopt::Long;

my ($init, $database);
GetOptions ( 
	    "init" => \$init,
	    'database:s' => \$database,
	   );

use constant USERNAME => 'wormpub';
use constant PASSWD   => 'wormpub';

use lib '../lib','./lib','./blib';
use lib $ENV{'CVS_DIR'};

use Wormbase;
use NameDB;

$database = glob("~wormpub/DATABASES/current_DB") unless $database;

my $db;
if( $init ) {
  $db = NameDB->connect('gene_id',USERNAME,PASSWD);
  $db->initialize(1) if $init;  # start with a completely clean slate
  # create domains to be stored
  $db->addDomain(Protein => 'P%010d');
  $db->addDomain(Gene => 'WBGene%08d');

  print $db->getTemplate('Gene'),"\n";

  $db->setDomain('Gene');
  my $result  = $db->addNameType('CGC',1,1);
  $result  = $db->addNameType('GenBank',1);
  $result  = $db->addNameType('EMBL',0,0);
  $result  = $db->addNameType('WormPep',0,1);
  $result  = $db->addNameType('CDS',1,1);

  print "Gene name types currently accepted\n";
  my $types = $db->getNameTypes;
  print join "\n",@$types,"\n";
  print join(",",$db->getNameTypeAttributes('CGC')),"\n";
}

my $ace = Ace->connect(-path => "$database") or die Ace->error(),"\n";

my $WBGenes = $ace->fetch_many('Gene' => "*" );

my %import_gene;
while ( my $gene = $WBGenes->next ) {
  #print $gene->name,"\n";
  my $id = $gene->name;

  next unless $gene->History;
  my @row = $gene->History->row();
  my( $version_change, $ver, $date, $who, $event, $event_detail, $comment) = map($_->name,@row);

  my $cgc_name  = $gene->CGC_name;
  my @CDS       = $gene->Corresponding_CDS;
  my $seq_names = $gene->Sequence_name;
  my @mol_names = $gene->Molecular_name;


  #primary_identifier needs object_id | object_public_id | domain_id | object_version | object_live
  #identifier_log | object_id | log_version | log_what | log_who | log_when       | log_related_object | log_name_type | log_name_value

  $import_gene{'object_id'} = $id;
  $import_gene{'who'}       = $who;
  $import_gene{'when'}      = $when;
  $import_gene{'name'}      = $seq_name;
  import
  

  $db->import('Gene', %import_gene )
}




$db->disconnect if $db;

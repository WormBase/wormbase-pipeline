#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

my ($init, $database, $target);
GetOptions ( 
	    "init" => \$init,
	    'database:s' => \$database,
	    'gene:s'     => \$target,
	   );

use constant USERNAME => 'wormpub';
use constant PASSWD   => 'wormpub';

use lib '../lib','./lib','./blib';
use lib $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'}."/NAMEDB/lib";

use Wormbase;
use NameDB;
use HTTP::Date;

$database = glob("~wormpub/DATABASES/current_DB") unless $database;

my $db = NameDB->connect('wbgene_id:mcs2a',USERNAME,PASSWD);
if( $init ) {
  $db->initialize(1) if $init;  # start with a completely clean slate
  # create domains to be stored
  $db->addDomain(Gene => 'WBGene%08d');

  print $db->getTemplate('Gene'),"\n";

  $db->setDomain('Gene');
  my $result  = $db->addNameType('CGC',1,1);
 # $result  = $db->addNameType('GenBank',1,0);
 # $result  = $db->addNameType('EMBL',0,0);
 # $result  = $db->addNameType('WormPep',0,1);
  $result  = $db->addNameType('CDS',0,1);
  $result  = $db->addNameType('Sequence',1,1);
 # $result  = $db->addNameType('Mol_name',1,1);
  $result  = $db->addNameType('Public_name',1,1);

  print "Gene name types currently accepted\n";
  my $types = $db->getNameTypes;
  print join "\n",@$types,"\n";
  #print join(",",$db->getNameTypeAttributes('CGC')),"\n";
}

my $ace = Ace->connect(-path => "$database") or die Ace->error(),"\n";

my $WBGenes;
if ($target) {
  my $q=<<END;
find Gene $target;
END

  $WBGenes = $ace->fetch_many( -query => "$q" );
}
else {
  $WBGenes = $ace->fetch_many('Gene' => "*" );
}

$db->setDomain('Gene');
print $db->getDomainId,"\n";

while ( my $gene = $WBGenes->next ) {

  my $id = $gene->name;
  next unless $gene->History;

  my $cgc_name  = $gene->CGC_name->name        if $gene->CGC_name;
  my $seq_name  = $gene->Sequence_name->name   if $gene->Sequence_name;
  my $live      = defined ($gene->at('Status.Live')) ? 1 : 0 ;
  my @mol_names = map($_->name,$gene->Molecular_name) if $gene->Molecular_name;
  my @CDS_names = map($_->name,$gene->Corresponding_CDS) if $gene->Corresponding_CDS;

  my @history = $gene->Version_change;
 
  foreach my $his_event ( @history ) {
    {
      my @row = $his_event->row();
      my( $ver, $date, $who, $event, $event_detail, $comment) = map($_->name,@row);

      #fix date format 2004-07-01_09:55:50  -> 200400701
      $date = &convert_date( $date );

      #primary_identifier needs object_id | object_public_id | domain_id | object_version | object_live
      #identifier_log | object_id | log_version | log_what | log_who | log_when       | log_related_object | log_name_type | log_name_value

      if ( $ver == 1 ) {
	#create new 
	my %import_gene;
	$import_gene{'id'}   = $id;
	$import_gene{'who'}     = $who;
	$import_gene{'when'}    = $date;
	push( @{$import_gene{'names'}->{'CGC'}}, $cgc_name ) if $cgc_name;
	push( @{$import_gene{'names'}->{'Public_name'}},( $cgc_name ? $cgc_name : $seq_name) );
	push( @{$import_gene{'names'}->{'Sequence'}}, $seq_name);
	#push( @{$import_gene{'names'}->{'Mol_name'}}, @mol_names);
	push( @{$import_gene{'names'}->{'CDS'}}, @CDS_names);
	$import_gene{'live'}    = $live;
	my $internal_id = $db->import_object(\%import_gene );
      } 
      else {
	#update existing
	$comment = "" unless $comment;
	print "$id, $ver, $date, $who, $event, $event_detail, $comment\n";
      }
    }
  }
}

$db->disconnect if $db;

sub convert_date
  {
    my $date = shift;
    my @fields = &HTTP::Date::parse_date( $date );
    $fields[1] =  sprintf("%1d$fields[1]");
    $fields[-1] = "";

    my $datestamp = join("",@fields);
    return $datestamp;
}



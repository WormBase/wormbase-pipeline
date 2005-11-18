#!/usr/local/bin/perl -w

use lib "lib";
use lib "/nfs/WWW/SANGER_docs/perl";
use SangerWeb;

use NameDB;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);

my $sw = SangerWeb->new( {
			  'title'  => "WormBase NameDB",
			  'banner' => "Results",
			  'inifile'=> "/nfs/WWW/SANGER_docs/Projects/C_elegans/header.ini"
			 });

my $USER = $sw->username();

$USER = 'wormpub' unless $USER;

print $sw->header();
print start_html;
print "logged in as $USER<br><br>";

print "@ARGV\n";

my $action;
my $type;
my $name;
my $gene_id;

$action  = param('action');
$type    = param('type');
$name    = param('new_name') || param('remove_name');
$gene_id = param('id');
$merge_id = param('id_2');
$lookup   = param('gene');

my $domain = 'Gene';
my $DB   = 'wbgene_id;mcs2a';
#my $USER = 'wormpub';
my $PASS = 'wormpub';

#verify if valid name

my $db = NameDB_handler->connect($DB,$USER,$PASS);
my $id;
$db->setDomain($domain);


# query current status
if ( $action eq 'query') {
  my $gene = $lookup;
  ($gene) = $db->idGetByTypedName($type=>$lookup) unless ("$type" eq "WBGene");
  &print_history($gene);
  &printAllNames($gene);
}

# Add a new gene
elsif ( $action eq 'new_gene') {
  &validate_name($name, $type);
  &check_pre_exists($name, $type);
  &isoform_exists($name, $type);
  &make_new_obj($name, $type);
}

# Add name to existing gene
elsif ( $action eq 'add_name'){
  &validate_name($name, $type);
  &validate_geneid($name, $type);
  &check_pre_exists($name, $type);
  &add_name($name, $type);
  print "Added $name to $gene_id<br>";
}

# Kill existing gene
elsif( $action eq "kill_gene"){
  &validate_geneid;
  print "$gene_id killed<br>" if $db->idKill($gene_id);
}

# split existing gene
elsif ($action eq "split_gene") {
  &validate_geneid;
  &validate_name;
  &check_pre_exists;
  my $id = $db->idSplit($gene_id);
  $type = "CDS";
  &add_name($id);

  print "Split $gene_id creating $id with CDS $name<br>"
}
elsif ($action eq "merge_genes") {
  &validate_geneid();
  &validate_geneid($merge_id);
  if ($db->idMerge($merge_id,$gene_id)) {
    print "merged $merge_id into $gene_id<br>";
  }
  else {
    print "merge failed<br>";
  }
}
elsif( $action eq "remove_name" ) {
  &validate_geneid();
  &validate_name;
  my ($exist_id) = $db->idGetByTypedName($type,$name);
  if( "$exist_id" ne "$gene_id" ) {
    &dienice("$name is not a name of $gene_id\n");
  }
  if( $db->delName($gene_id,$type,$name) ) {
    print "Removed $name as name for $gene_id<br><br>Remaining names are <br>";
  }
  else {
    print "name removal failed<br>Names for $gene_id";
  }
  #print remaining names
  &printAllNames;
}


print end_html;

print $sw->footer();

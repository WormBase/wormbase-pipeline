#!/usr/bin/env perl

use strict;
use warnings;
use XML::LibXML;
use XML::LibXML::Reader;
use Net::FTP;

use Getopt::Long;

use Wormbase;


my $PROVIDER_ID = 1386;
my $PROVIDER_NAME =  'WormBase';
my $PROVIDER_DESC = 'WormBase Publication pages';
my $PROVIDER_EMAIL = 'hinxton@wormbase.org';

my $LINK_BASE_URL = 'http://www.wormbase.org/resources/paper';
my $SCHEMA_URL = 'http://europepmc.org/docs/labslink.xsd';

my ($provider_xml,
    $links_xml,
    $database,
    $tace,
    $upload,
    $test,
    $debug,
    $store,
    );

&GetOptions("providerxml=s" => \$provider_xml,
            "linksxml=s"    => \$links_xml,
            "database=s"    => \$database,
            "tace=s"        => \$tace,
            "upload"        => \$upload, 
            "test"          => \$test,
            "debug=s"       => \$debug,
            "store=s"       => \$store,
    );

if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new(-debug => $debug, 
                            -test => $test,
      );
}

my $log = Log_files->make_build_log($wormbase);

$tace = $wormbase->tace of not defined $tace;
$database = $wormbase->autoace if not defined $database;

# Build the XML (it is very simple, so just substitute values into a string)
my $provider_details_xml_string = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<providers>\n";
$provider_details_xml_string .= "<provider>\n<id>$PROVIDER_ID</id>\n<resourceName>$PROVIDER_NAME</resourceName>\n";
$provider_details_xml_string .= "<description>$PROVIDER_DESC</description>\n<email>$PROVIDER_EMAIL</email>\n</provider>\n</providers>";

my $parser = XML::LibXML->new;
my $provider_xml_obj = $parser->parse_string($provider_details_xml_string);

my $papers = &paper_tm_query();

# Create a blank XML DOM object and add the root element
my $link_xml_obj = XML::LibXML::Document->new('1.0', 'UTF-8');
my $root_element = $link_xml_obj->createElement('links');
$link_xml_obj->addChild($root_element);

foreach my $paper (@$papers) {
  my ($wb_id, $id) = @$paper;

  my $url = "$LINK_BASE_URL/$wb_id";
  my $source = "MED";
  my $title = "WormBase Paper $wb_id";
  
  # Create an in-memory representation of the XML structure for each link
  my $link_element = $link_xml_obj->createElement('link');
  my $resource_element = $link_xml_obj->createElement('resource');
  my $url_element = $link_xml_obj->createElement('url');
  my $title_element = $link_xml_obj->createElement('title') if $title;
  my $record_element = $link_xml_obj->createElement('record');
  my $source_element = $link_xml_obj->createElement('source');
  my $id_element = $link_xml_obj->createElement('id');
  $link_element->addChild($resource_element);
  $resource_element->addChild($url_element);
  $resource_element->addChild($title_element);
  $link_element->addChild($record_element);
  $record_element->addChild($source_element);
  $record_element->addChild($id_element);
  
  # Populate the XML elements with values from the current line of the tab-delimited file
  $url_element->appendTextNode($url);
  $title_element->appendTextNode($title);
  $source_element->appendTextNode($source);
  $id_element->appendTextNode($id);
  
  # Associate the link with the provider details via the provider ID
  $link_element->setAttribute('providerId', $PROVIDER_ID);
  
  # Add the completed link element to the root node
  $root_element->addChild($link_element);	
}

my $schema = XML::LibXML::Schema->new(location => $SCHEMA_URL);

foreach my $pair ([$provider_xml_obj, $provider_xml],
                  [$link_xml_obj, $links_xml]) {
  my ($xml_obj, $xml_file) = @$pair;

  eval {
    $schema->validate($xml_obj) 
  };
  if ($@) {
    die "$xml_obj did not validate; exiting\n";
  }

  open(my $fh, ">$xml_file") or die "Cannot open $xml_file for writing\n";
  print $fh $xml_obj->toString(1);
  close($fh);
}
exit 0;

######################################################
sub paper_tm_def {
  
  my $tm_def = "/tmp/wormbase_lablinks.tm.$$.def";

  my $query = <<'EOF';
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Paper 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Database 
From 1 
Tag Database 
Condition "MEDLINE"
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Database_field 
Right_of 2 
Tag  HERE  
Condition "PMID"
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  
EOF

  open(my $fh, ">$tm_def") or die "Could not open $tm_def for writing\n";
  print $fh $query;
  close($fh);

  return $tm_def;
}


######################################################
sub paper_tm_query {

  my $def = &paper_tm_def();
  my $tm_cmd = "Table-maker -p \"$def\"\nquit\n";
 
  my @papers;
 
  open(my $tace_fh, "echo '$tm_cmd' | $tace $database |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"\S+\"\s+\"\S+\"\s+\"(\S+)\"/ and do {
      my ($wb_paper, $medline_id) = ($1, $2);

      push @papers, [$wb_paper, $medline_id];
    };
  }

  unlink $def;

  return \@papers;
}



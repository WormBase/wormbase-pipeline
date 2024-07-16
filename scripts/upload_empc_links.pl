#!/usr/bin/env perl

use strict;
use warnings;
use XML::LibXML;
use XML::LibXML::Reader;
use Net::FTPSSL;
use LWP::Simple;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Wormbase;


my $PROVIDER_ID = 1386;
my $PROVIDER_NAME =  'WormBase';
my $PROVIDER_DESC = 'WormBase is an international consortium of biologists and computer scientists dedicated to providing the research community with accurate, current, accessible information concerning the genetics, genomics and biology of Caenorhabditis elegans and other nematodes.';
my $PROVIDER_EMAIL = 'hinxton@wormbase.org';

my $LINK_BASE_URL = 'http://www.wormbase.org/resources/paper';
my $SCHEMA_URL = 'https://europepmc.org/docs/labslink.xsd';

my ($provider_xml,
    $links_xml,
    $database,
    $tace,
    $upload,
    $test,
    $debug,
    $store,
    $authfile,
    $wormbase,
    );

&GetOptions("providerxml=s" => \$provider_xml,
            "linksxml=s"    => \$links_xml,
            "database=s"    => \$database,
            "tace=s"        => \$tace,
            "upload"        => \$upload, 
            "test"          => \$test,
            "debug=s"       => \$debug,
            "store=s"       => \$store,
            "authfile=s"    => \$authfile,
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

$tace = $wormbase->tace if not defined $tace;
$database = $wormbase->autoace if not defined $database;
$provider_xml = $wormbase->acefiles . "/Wormbase_provider.xml" if not defined $provider_xml;
$links_xml = $wormbase->acefiles . "/Wormbase_links.xml" if not defined $links_xml;
$authfile = $wormbase->wormpub . "/ebi_resources/EPMCFTP.s" if not defined $authfile;


$log->write_to("Fetching data from WormBase...\n");
my $papers = &paper_tm_query();

$log->write_to("Building XML documents...\n");
print "Building XML documents\n";

# Build the XML (it is very simple, so just substitute values into a string)
my $provider_details_xml_string = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<providers>\n";
$provider_details_xml_string .= "<provider>\n<id>$PROVIDER_ID</id>\n<resourceName>$PROVIDER_NAME</resourceName>\n";
$provider_details_xml_string .= "<description>$PROVIDER_DESC</description>\n<email>$PROVIDER_EMAIL</email>\n</provider>\n</providers>";

print "Created XML header\n";

my $parser = XML::LibXML->new;
my $provider_xml_obj = $parser->parse_string($provider_details_xml_string);

# Create a blank XML DOM object and add the root element
my $link_xml_obj = XML::LibXML::Document->new('1.0', 'UTF-8');
my $root_element = $link_xml_obj->createElement('links');
$link_xml_obj->addChild($root_element);

print "Created blank XML object\n";

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
print "Finished parsing papers\n";

$log->write_to("Validating and writing  XML documents...\n");
my $tmp_schema_file = $wormbase->acefiles . "/epmc-labslink.xmd";
getstore($SCHEMA_URL, $tmp_schema_file);
my $schema = XML::LibXML::Schema->new(location => $tmp_schema_file);

foreach my $pair ([$provider_xml_obj, $provider_xml],
                  [$link_xml_obj, $links_xml]) {
  my ($xml_obj, $xml_file) = @$pair;

  eval {
    print "Validating $xml_file\n";
    $schema->validate($xml_obj) 
  };
  if ($@) {
    $log->log_and_die("$xml_obj did not validate; exiting\n");
  }

  open(my $fh, ">$xml_file") or $log->log_and_die("Cannot open $xml_file for writing\n");
  print $fh $xml_obj->toString(1);
  close($fh);
  print "Done XML file\n";
}

my $provider_xml_z = "${provider_xml}.gz";
my $links_xml_z    = "${links_xml}.gz";
print "gz filenames\n";
unlink $tmp_schema_file;
unlink $provider_xml_z if -e $provider_xml_z;
unlink $links_xml_z if -e $links_xml_z;
print "Check files\n";
$wormbase->run_command("gzip $provider_xml", $log) and 
    $log->log_and_die("Failed to gzip $provider_xml");
$wormbase->run_command("gzip $links_xml", $log) and
    $log->log_and_die("Failed to gzip $links_xml");
print "Done gzip\n";

if ($upload) {
  my ($ftp_host, $ftp_user, $ftp_pass, $ftp_dir);

  open (my $authin, $authfile) or $log->log_and_die("Could not open $authfile for reading\n");
  while(<$authin>) {
    /^HOST:(\S+)$/ and $ftp_host = $1;
    /^USER:(\S+)$/ and $ftp_user = $1;
    /^PASS:(\S+)$/ and $ftp_pass = $1;
    /^DIR:(\S+)$/  and $ftp_dir  = $1;
  }
 
 print "FTP details $ftp_host $ftp_user $ftp_pass $ftp_dir\n"; 
  $log->write_to("Connecting to FTP site...\n");

  my $ftp = Net::FTPSSL->new($ftp_host, ReuseSession => 1, Debug => 0) 
      or $log->log_and_die("Cannot connect to $ftp_host: $@");
  $ftp->login($ftp_user,"$ftp_pass")
      or $log->log_and_die ("Cannot login to $ftp_host using WormBase credentials\n". $ftp->message);
  $ftp->cwd($ftp_dir) 
      or $log->log_and_die ("Cannot change into to_ena dir for upload of files\n". $ftp->message);
  $ftp->binary();

  print "Starting FTP\n";
  foreach my $file ("$provider_xml_z", "$links_xml_z") {
    $log->write_to("Depositing $file on FTP site...\n");
    $ftp->put($file) or $log->log_and_die ("FTP-put failed for $file: ".$ftp->message."\n");
  }
  $ftp->quit;
}

$log->mail();
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

  open(my $fh, ">$tm_def") or $log->log_and_die("Could not open $tm_def for writing\n");
  print $fh "$query\n";
  close($fh);

  return $tm_def;
}

exit;

######################################################
sub paper_tm_query {

  my $def = &paper_tm_def();
  my $tm_cmd = "Table-maker -p \"$def\"\nquit\n";
  print "Doing tablemaker: $tm_cmd\n";
 
  my @papers;
 
  open(my $tace_fh, "echo '$tm_cmd' | $tace $database |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"\S+\"\s+\"\S+\"\s+\"(\S+)\"/ and do {
      my ($wb_paper, $medline_id) = ($1, $2);

      push @papers, [$wb_paper, $medline_id];
    };
  }

  unlink $def;

  print "Completed Tablemaker\n";

  return \@papers;
}



#!/software/bin/perl
##########################################################
#
# a.) there is a case to use the XML emitter from the KO-Consortium if the number of tags increase
# b.) and maybe just pushing it through XML-tidy instead of formatting it by hand
#
use Getopt::Long;
use Ace;
use Wormbase;
use Log_files;
use Storable;
use Digest::MD5 qw(md5_hex);
use IO::File;
use strict;
	
my ($database,$outfile,$store,$debug,$test);
GetOptions(
		'database=s' => \$database,
		'file=s'     => \$outfile,
                'store=s'    => \$store,
                'debug=s'    => \$debug,
                'test'       => \$test,
)||die(@!);

my $wormbase;
if ($store){
   $wormbase= retrieve($store) or die("Can't restore wormbase from $store\n");
}else{
   $wormbase= Wormbase->new(-debug => $debug,
                            -test  => $test );
}

my $log=Log_files->make_build_log($wormbase);

$outfile||= "${\$wormbase->autoace}/ReutersCitationIndex.xml.gz";
my $file = new IO::File "| gzip -9 -c > $outfile";

$database||=$wormbase->autoace;

my $db = Ace->connect(-path => $database)||die(Ace->error);

my @timelist=localtime();
my $year = $timelist[5]+1900;
my $time=sprintf('%02u/%02u/%u',$timelist[3],$timelist[4],$year);

$log->write_to("dumping XML file for WS${$wormbase->version} $time\n");

# header
print $file "<report name=\"CI-Report\" date_generated=\"$time\">\n";

my $gIt = $db->fetch_many(-query => 'Find Gene CGC_name AND Corresponding_CDS')||$log->log_and_die(Ace->error);

while (my $gene =$gIt->next){
  print_gene($file,$gene);
}

# footer
print $file "</report>\n";

# cleanup
$file->close;
$db->close;
$log->mail;

sub print_gene{
  my ($f,$g)=@_;

  my @c = split(/\s+/,"${\$g->History->right(2)}");
  my $creation_date = $c[2];
 
  print STDERR "processing $g\n" if $debug;
  my $data =
  "  <date_provided>$year</date_provided>\n".
  "  <date_created>$creation_date</date_created>\n".
  "  <repository>WormBase</repository>\n".
  "  <owner>WormBase</owner>\n".
  "  <title>".$g->CGC_name." (${\$g->Gene_class->Description})</title>\n";

  ### people block ######################
  # search through the CGC_name, the Concise_description and the Provisional_description
  my @people = $g->at('Identity.Name.CGC_name')->col(3);
  push (@people, $g->at('Identity.Name.Other_name')->col(3))
              if ($g->Other_name && $g->at('Identity.Name.Other_name')->col(3));
  push (@people, $g->at('Gene_Info.Structured_description.Concise_description')->col(3)) 
              if $g->at('Gene_Info.Structured_description.Concise_description');
  push (@people, $g->at('Gene_Info.Structured_description.Provisional_description')->col(3))
              if $g->at('Gene_Info.Structured_description.Provisional_description');

  # and then print them
  # if we get our hands on ORCID ids, they should go in there too
  foreach my $p (@people){
    if ($p && $p->class eq 'Person'){
     $data.="  <author>\n";
     $data.='   <full_name>'.$p->Full_Name."</full_name>\n";
     $data.='   <affiliation>'.get_address($p)."</affiliation>\n" if (get_address($p));
     $data.="  </author>\n";
    }
  }
  #######################################

  $data.="  <abstract><![CDATA[".$g->Concise_description."]]></abstract>\n" if $g->Concise_description;
  $data.="  <publisher>WormBase</publisher>\n".
  "  <organism>".$g->Species."</organism>\n".
  "  <data_type>Protein Coding</data_type>\n";

  # Paper block
  foreach my $p ($g->Reference){
      foreach $db($p->Database){
         my @row = $db->row;
         $data.="  <citation pubmed_id=\"$row[2]\">$p</citation>\n" if $row[0] eq 'MEDLINE';
         last;
      }
  }

  # that is a checksum to make it possible to look for changes
  print $f " <record record_id=\"$g\" source_url=\"http://wormbase.org/db/gene/gene?name=$g\" checksum=\"".
  &md5_hex($data)."\">\n$data </record>\n";
}

# that is to get the Institution affiliation from the Address block, which is in semi-random order
sub get_address{
  my ($person) = @_;
  return undef unless $person->Address;
  my $n=0;
  while ($n <= scalar $person->Address->tags){
     return $person->Address->down($n)->right if "${\$person->Address->down($n)}" eq 'Institution';
     $n++;
  }
}

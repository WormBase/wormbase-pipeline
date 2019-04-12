#!/usr/bin/env perl

#####################################################################################
# Even though the GO term definitions are imported to autoace from
# citace, there may be some GO terms that are referenced by InterPro
# etc. that have not yet been defined by Caltech. We therefore do a
# check for those GO_term object that do not have a fully populated
# set of tags and we populate them here.
#####################################################################################

use strict;
use Ace;
use Getopt::Long;
use Storable;
          
use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

my $CURRENT_GO_OBO =  "http://purl.obolibrary.org/obo/go.obo";
## Old url http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo";

my ($help, $debug, $test, $store, $acedbpath, $noload, $acefile, $wormbase);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
            "noload"     => \$noload,
            "acefile"    => \$acefile, 
	    )||die(@!);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $log  = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace;        # TACE PATH

$acefile = $wormbase->acefiles . "/go_defs.ace" if not defined $acefile;
$acedbpath = $wormbase->autoace if not defined $acedbpath;

my $cmd  = "Query Find GO_term Where !Type\nlist -a\nquit";

my (@go, %go);

open (my $tace_fh, "echo '$cmd' | $tace $acedbpath |");
while (<$tace_fh>) {
  chomp;
  next if (/acedb\>/);
  next if (/\/\//);
  if (/GO_term\s+:\s+\"(\S+)\"/) {
    push @go, $1;
    $go{$1} = 1;
  }
}
close($tace_fh);

if (@go == 0) {
  $log->write_to("No GO_terms are missing a Type tag.\nAll appear OK.\n");
  $log->mail;
  exit(0);
} 
  
$log->write_to("There are " . scalar @go . " GO_terms that are incomplete - fixing them\n");    
$log->write_to("GO_terms objects are:\n@go\n\n");

my $gocount = 0;

#[Term]
#id: GO:0000003
#name: reproduction
#namespace: biological_process
#alt_id: GO:0019952
#alt_id: GO:0050876
#def: "The production by an organism of new individuals that contain some portion of their genetic material inherited from that organism." [GOC:go_curators, GOC:isa_complete, ISBN:0198506732 "Oxford Dictionary of Biochemistry and Molecular Biology"]
#subset: goslim_generic
#subset: goslim_pir
#subset: goslim_plant
#subset: gosubset_prok
#synonym: "reproductive physiological process" EXACT []
#xref: Wikipedia:Reproduction
#is_a: GO:0008150 ! biological_process

# We only really want the 'namespace' field for the 'Type' tag, but
# grap a few other things as well

my $go_obo_file = "/tmp/GO_file_$$.".$wormbase->species;
$log->write_to("Attempting to wget the latest GO data\n");
system("wget -q -O $go_obo_file $CURRENT_GO_OBO") 
    and $log->log_and_die("Could not fetch latest version of GO file ($CURRENT_GO_OBO)\n");

open (my $go_fh,"<$go_obo_file") or $log->log_and_die("Cannot open $go_obo_file\n");
open (my $go_out_fh, ">$acefile") or $log->log_and_die("Cannot write to $acefile\n");

my ($namespace, $name, $def, @id);

while (<$go_fh>){
  chomp;
  if ($_ =~ /^id:\s+(\S+)/) {
    push @id, $1;
  }
  if ($_ =~ /^alt_id:\s+(\S+)/) {
    push @id, $1;
  }
  if ($_ =~ /name:\s+(.+)$/) {
    my $raw_name = $1;
    $raw_name =~ s/\s*$//;
    $name = $raw_name;
  }
  if ($_ =~ /namespace:\s+(\S+)/) {
    $namespace = $1;
  }
  if ($_ =~ /def:\s+\"(.+)\"$/ or $_ =~ /def:\s+\"(.+)\"\s+\[.+\]$/) {
    $def = $1;
  }
  if ($_ =~ /^\s*$/) {
    my @match = grep { $go{$_} } @id;
    
    if (scalar(@match) == 1) {
      my ($wbid) = @match;
      
      print $go_out_fh "\nGO_term : \"$wbid\"\n";
      print $go_out_fh "Definition \"$def\"\n";
      print $go_out_fh "Name \"$name\"\n";
      print $go_out_fh "Type \"".ucfirst($namespace)."\"\n";
      print $go_out_fh "\n";
      $gocount++;
    }
    @id = ();
  }
}

$log->write_to("added $gocount GO definitions\n");
close($go_out_fh);
    
# load file to autoace if -load specified
if (-s $acefile and not $noload) {
  $wormbase->load_to_database($wormbase->autoace, "$acefile", 'go_defs', $log);
}

# tidy up and exit
$wormbase->run_command("rm $go_obo_file",$log);

$log->mail();
exit(0);


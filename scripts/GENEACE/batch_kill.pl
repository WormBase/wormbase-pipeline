#!/software/bin/perl -w
use strict;
use lib '../blib/lib';
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib $ENV{'CVS_DIR'};
use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;
=pod

=head batch_kill.pl

=item Options:

  -file	     file containing genes to kill <Mandatory>

    FORMAT:
    WBGene00001234
    Remark "Removed"
    WBPerson1983

    WBGene00001235
    Remark "Killed"
    WBPerson1849 

    The blank line between entries is ESSENTIAL

  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans genes
  -load      loads the resulting .ace file into geneace.
  -test      use the test nameserver
  -ns        Kill's the gene in the nameserver as well as producing the .ace 
             file for geneace
  -user      username                 <Manditory if using -ns>
  -password  password                 <Manditory if using -ns>

e.g. perl batch_kill.pl -file deathrow.txt [simple example]

     perl batch_kill.pl -u fred -p secret -file deathrow.txt -ns -test -debug mt3

=cut

my ($USER,$PASS, $test, $file, $ns, $debug, $load);
GetOptions(
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'ns'         => \$ns,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
	  ) or die;

my $species = 'elegans';
my $log;
if (defined $USER) {$log = Log_files->make_log("NAMEDB:$file", $USER);}
elsif (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}
my $DB;
my $db;
my $ecount;
if ($test) {
    $DB = 'test_wbgene_id;mcs2a;3305';
  } else {
    $DB = 'wbgene_id;mcs2a;3305';
}
my $wormbase = Wormbase->new("-organism" =>$species);
my $database = "/nfs/disk100/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) killing genes in file [${file}]\n\n");
$log->write_to("TEST mode is ON!\n\n") if $test;

if ($ns) {
$db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/htdocs');
$db->setDomain('Gene');
}
my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");

my $output = $database."/NAMEDB_Files/batch_kill.ace";

##############################
# warn/notify on use of -load.
##############################
if (!defined$load) {$log->write_to("2) You have decided not to automatically load the output of this script\n\n");}
elsif (defined$load) { $log->write_to("2) Output has beed scheduled for auto-loading.\n\n");}

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$output") or $log->log_and_die("cant write output: $!\n");
my($gene,$person,$remark);
my $count;
while(<FILE>){
    chomp;
    unless (/\w/) {
	&kill_gene;
    }
    else { #gather info
	if   (/^(WBGene\d{8})/) { $gene = $1; } 
	elsif(/^(WBPerson\d+)/) { $person = $1; }
	elsif(/^Remark\s+\"(.*)\"/){$remark = $1}
	else { $log->error("malformed line : $_\n") }
    }
}
$log->write_to("3) $count genes in file to be killed\n\n");
$log->write_to("4) $count genes killed\n\n");
&kill_gene; # remember the last one!
&load_data if ($load);
$log->write_to("5) Check $output file and load into geneace.\n") unless ($load);
$log->mail();


sub kill_gene {
    if($gene and $person and $remark) {
	$count++;
	my $geneObj = $ace->fetch('Gene', $gene);
	if($gene) {
	    my $ver = $geneObj->Version->name;
	    $ver++;
	    print ACE "\nGene : $gene\nVersion $ver\nHistory Version_change $ver now $person Event Killed\nDead\nRemark \"$remark\" Curator_confirmed $person\n-D Sequence_name\n-D Method\n-D Map_info\n-D Allele\n";
	    #nameserver kill
	    $db->kill_gene($gene) if $ns;
	}
	else {
	    $log->error("no such gene $gene\n");
	}
    }
#    elsif (!defined($gene) && !defined($person)&& !defined($remark)) {
elsif (!defined($gene && $person && $remark)) {
  $ecount++;
      $log->error("Warning: additional blank ($ecount) line in input file has been ignored\n");
}
    else {
	$log->error("missing info on $gene : $person : $remark\n");
    }
    undef $gene; undef $person ;undef $remark;
}

sub load_data {
# load information to $database if -load is specified
$wormbase->load_to_database("$database", "$output", 'batch_kill.pl');
$log->write_to("5) Loaded $output into $database\n\n");
$wormbase->run_command("rm $output\n");
$log->write_to("6) Output file has been cleaned away like a good little fellow\n\n");
print "Finished!!!!\n";
}

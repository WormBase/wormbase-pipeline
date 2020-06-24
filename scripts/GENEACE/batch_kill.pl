#!/software/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;

=pod

=head batch_kill.pl

=item Options:

  -test      use the test nameserver
  -file	     file containing genes to kill <Mandatory>

    FORMAT:
    WBGene00001234
    Remark : "Removed"
    WBPerson1983

    WBGene00001235
    Remark : "Killed"
    WBPerson1849 

    The blank line between entries is ESSENTIAL
  
  -domain    set domain to whatever (Gene Variation Feature)
  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans genes
  -load      loads the resulting .ace file into geneace.
  -ns        Kill's the gene in the nameserver as well as producing the .ace 
             file for geneace

e.g. perl batch_kill.pl -file deathrow.txt [simple example]

     perl batch_kill.pl -file deathrow.txt -ns -debug pad

=cut

my ($test, $help, $debug, $verbose, $store, $wormbase);
my ($file, $ns, $load, $transposon, $domain);
if (!defined $domain) {$domain = 'Gene'}


my $species = 'elegans';
GetOptions(
	   "test"       => \$test,
	   'file:s'     => \$file,
	   'ns'         => \$ns,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
	   'species:s'  => \$species,
	   'domain:s'   => \$domain,
	  ) or die;

my $log;
if (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}

my $db;

my $ecount;

$wormbase = Wormbase->new(
			  -organism => $species, 
			  -debug => $debug,
			  -test => $test,
			 );

my $database = "/nfs/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) killing ${domain}s in file [${file}]\n\n");

##############################
# warn/notify on use of -load.
##############################
if (!defined $load) {$log->write_to("You have decided not to automatically load the output of this script\n\n");}
elsif (defined $load) { $log->write_to("Output has been scheduled for auto-loading.\n\n");}

if ($ns) {
  $log->write_to("Contacting NameServer.....\n");
  $db = NameDB_handler->new($wormbase, $test);
}

my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");

my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
my $outname = "batch_kill.ace";
my $output = "$outdir"."$outname";

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$output") or $log->log_and_die("cant write output: $!\n");
my($gene,$person,$remark,$tflag,);
my $count;
while(<FILE>){
  chomp;
  unless (/\w/) {
    &kill_gene;
    undef $tflag;
  } else { #gather info
    if (/(WBGene\d{8}|WBVar\d{8})|WBsf\d+|WBStrain\d{8}/) { 
      $gene = $1; 
    } 
    elsif(/(WBPerson\d+)/) { $person = $1; }
    elsif(/Remark\s+\:\s+\"(.*)\"/){
      $remark = $1;
      if (/Transposon_CDS/) {$tflag = "1";}
    }
    else { $log->error("malformed line : $_\n") }
  }
}
&kill_gene; # remember the last one!
$log->write_to("\n$count ${domain}s in file to be killed\n\n");
$log->write_to("$count ${domain}s killed\n\n");
&load_data if ($load);
$log->write_to("Check $output file and load into geneace.\n") unless ($load);
$log->mail();


sub kill_gene {
  if($gene and $person and $remark) {
    $count++;
    my $geneObj = $ace->fetch($domain, $gene);
    if($geneObj) {
      
      if ($domain eq 'Gene'){
	#nameserver kill
	$log->write_to("NS->kill $gene\n");
	$db->kill_gene($gene) if $ns;

	my $ver = $geneObj->Version->name;
	$ver++;
	# We no longer Kill Transposon_CDSs/Pseudogenes but suppress them.
	if ($tflag) {
	  print ACE "\nGene : $gene\nVersion $ver\nHistory Version_change $ver now $person Event Suppressed\nSuppressed\nRemark \"$remark\" Curator_confirmed $person\n-D Method\n-D Map_info\n-D Allele\n";
	} else {
	  print ACE "\nGene : $gene\nVersion $ver\nHistory Version_change $ver now $person Event Killed\nDead\nRemark \"$remark\" Curator_confirmed $person\n-D Sequence_name\n-D Method\n-D Map_info\n-D Other_name\n-D Allele\n";
	}
      } elsif($domain eq 'Variation'){
	#nameserver kill
	$log->write_to("NS->kill $gene\n");
	$db->kill_variations([$gene]) if $ns;

	print ACE "\nVariation : $gene\nStatus Dead Curator_confirmed $person\nRemark \"$remark\"\n";
      } elsif($domain eq 'Feature'){
	#nameserver kill
	$log->write_to("NS->kill $gene\n");
	$db->kill_features([$gene]) if $ns;

        print ACE "\nFeature : $gene\nDeprecated \"feature has been retired\" Curator_confirmed $person\nRemark \"$remark\"\n";
      }
    }
    else {
      $log->error("no such $domain $gene in $database\n");
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
  $wormbase->load_to_database("$database", "$output", 'batch_kill.pl', $log, undef, 1);
  $log->write_to("Loaded $output into $database\n\n");
  $wormbase->run_command("mv $output $backupsdir"."$outname". $wormbase->rundate. "\n");
  $log->write_to("Output file has been cleaned away like a good little fellow\n\n");
  print "Finished!!!!\n";
}

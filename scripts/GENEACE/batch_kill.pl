#!/usr/local/bin/perl -w
use strict;
use lib '../blib/lib';
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib $ENV{'CVS_DIR'};
use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;

=pod

=head batch_addname.pl

=item Options:

  -user      username
  -password  password
  -file	     file containing genes to kill
    FORMAT: 
    WBGene00001234
    Remark "Removed"
    WBPerson1983

    WBGene00001235
    Remark "Killed"
    WBPerson1849 

 The blank line between entries is ESSENTIAL

  -test      use the test nameserver

e.g. perl batch_kill.pl -u fred -p secret -file deathrow.txt -test

=cut

my ($USER,$PASS, $test, $file, $species, $ns);
GetOptions(
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'ns'         => \$ns
	  ) or die;

$species = 'elegans' unless $species;

my $log = Log_files->make_log("NAMEDB:$file", $USER);
my $DB;
if ($test) {
    $DB = 'test_wbgene_id;mcs2a;3305';
  } else {
    $DB = 'wbgene_id;mcs2a;3305';
}

$log->write_to("killing genes in $file\n\n");
$log->write_to("TEST mode is ON!\n\n") if $test;

my $db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/htdocs');
my $ace = Ace->connect('-path', '/nfs/disk100/wormpub/DATABASES/geneace') or $log->log_and_die("cant open geneace: $!\n");

$db->setDomain('Gene');

#open file and read

open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">batch_kill.ace") or $log->log_and_die("cant write output: $!\n");
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
&kill_gene; # remember the last one!

$log->write_to("=======================\nkilled $count genes\n");
$log->mail;


sub kill_gene {
    if($gene and $person and $remark) {
	$count++;
	#geneace kill
	my $geneObj = $ace->fetch('Gene', $gene);
	if($gene) {
	    my $ver = $geneObj->Version->name;
	    $ver++;
	    print ACE "\nGene : $gene\nVersion $ver\nHistory Version_change $ver now $person Event Killed\nDead\nRemark \"$remark\"\n";
	    print ACE "\nGene : $gene\n-D Sequence_name\n-D Method\n";
	    #nameserver kill
	    $db->kill_gene($gene) if $ns;
	}
	else {
	    $log->error("no such gene $gene\n");
	}
    }
    else {
	$log->error("missing info on $gene : $person : $remark\n");
    }
    undef $gene; undef $person ;undef $remark;
}

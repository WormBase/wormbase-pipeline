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

=head batch_merge.pl

=item Options:

  -old       Used the output format from the old nameserver

  -file	     file containing genes to merge <Mandatory>

  -override  This option forces the geneace ststus checks to be ignored for pushing 
             merges into the Nameserver when the changes have already been done in geneace
             via a batch adhoc submission. Cannot be used with the -load option

  -force     This option forces the .ace file to be loaded into geneace when there are Gene
             names involved as it's some what easier to resolve issues post merge instead of 
             overlooking when merges have failed because of this.

  -ns        Causes the merge to be done in the NS.

    FORMAT:

optional  EMAIL : mh6@sanger.ac.uk
optional   NAME : Michael Paulini
           DEAD : killed geneid for Bm9556 - WBGene00229817
           LIVE : retained geneid for Bm4144 - WBGene00224405
           USER : mh6 - WBPerson4055
optional   WARNING : Merged gene WBGene00229817 was dead

optional  EMAIL : pad@sanger.ac.uk
optional   NAME : Paul Davis
           DEAD : killed geneid for CBG11218 - WBGene00032378
           LIVE : retained geneid for CBG11131 - WBGene00032305
           USER : pad - WBPerson1983

optional  EMAIL : gw3@sanger.ac.uk
optional   NAME : Gareth Williams
           DEAD : killed geneid for F32E10.7 - WBGene00017994
           LIVE : retained geneid for F45E4.3 - WBGene00018468
           USER : gw3 - WBPerson4025

OLD:

GENE MERGE
USER : jspieth - WBPerson615
LIVE:retained geneID for CBN18483 - WBGene00157208
DEAD: killed geneID CBN20805 - WBGene00159530


    The blank line between entries is ESSENTIAL

  -debug     limits to specified user <Optional>
  -load      loads the resulting .ace file into geneace.

e.g. perl batch_merge.pl -file merger.txt


=cut

my ($help, $verbose, $store, $wormbase);
my ($file, $debug, $load, $old, $ns, $out, $force, $override);
GetOptions(
	   'file:s'     => \$file,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
	   'old'        => \$old,
	   'ns'         => \$ns,
	   'out:s'      => \$out,
	   'force'      => \$force,
	   'override'   => \$override,
	  ) or die;


my $species = 'elegans';
my $log;
if (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}

$wormbase = Wormbase->new(
                             -organism => $species, 
                             -debug => $debug, 
                            );


$log->write_to(sprintf("%s started at: %s\n",join( ' ',$0,@ARGV) ,`date +%H:%M:%S`));

my $db;
my $ecount;

if ($load && $override) {
   $log->log_and_die("OPTIONS ERROR: You cannot use the -load and -override option at the same time!\n");
}

my $database = "/nfs/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) killing genes in file [${file}]\n\n");


if ($ns) {
  $log->write_to("Contacting NameServer.....\n");
  $db = NameDB_handler->new($wormbase);
}

my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");

my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
$out ||= "batch_merge.ace";
my $output = $outdir.$out;

my %gene_versions; # remember the latest version used in all genes altered in case a gene is being merged in to more than once

##############################
# warn/notify on use of -load.
##############################
if (!defined $load) {$log->write_to("You have decided not to automatically load the output of this script\n\n");}
elsif (defined $load) { $log->write_to("Output has been scheduled for auto-loading.\n\n");}

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$output") or $log->log_and_die("cant write output: $!\n");
my($livegene, $deadgene, $user);
my $count=0;
while (<FILE>) {
    chomp;
    print "$_\n" if ($debug);
    unless (/\w/) {
	&merge_gene;
    }
    else { #gather info
	if (defined $old){
	    if   (/^LIVE:retained\s+geneID\s+for\s+\S+\s+-\s+(WBGene\d{8})/) { $livegene = $1; } 
	    elsif ((/^DEAD:\s+killed\s+geneID\s+\S+\s+-\s+(WBGene\d{8})/) || (/^DEAD:\s+killed\s+geneID\s+\s+-\s+(WBGene\d{8})/)) { $deadgene = $1; }
	    elsif (/^USER\s+:\s+\S+\s+-\s+(WBPerson\d+)/) { $user = $1;}
	    elsif (/^GENE MERGE/){} # ignore this line
	    else { $log->error("malformed line : $_\n") }
	}
	else {
	    if   (/\s+LIVE\s+:\s+retained\s+geneid\s+for\s+\S+\s+-\s+(WBGene\d{8})/) {
		$livegene = $1;
	    } 
	    elsif((/\s+DEAD\s+:\s+killed\s+geneid\s+for\s+\S+\s+-\s+(WBGene\d{8})/) || (/\s+DEAD\s+:\s+killed\s+geneid\s+for\s+\s+-\s+(WBGene\d{8})/)){ 
		$deadgene = $1; 
	    } 
	    elsif(/\s+USER\s+:\s+\S+\s+-\s+(WBPerson\d+)/) {
		$user = $1;
	    }
	    elsif(/\s+EMAIL/){} # ignore this line
	    elsif(/\s+NAME/){} # ignore this line
	    elsif(/\s+WARNING/){
		print "$_ before merging into $livegene \n"
	    }
	    else { 
		$log->error("malformed line : $_\n")
	    }
	}
    }
}

&merge_gene; # remember the last one!
close(ACE);
$log->write_to("$count genes in file to be merged\n\n");
$log->write_to("$count genes merged\n\n");
&load_data if ($load);
$log->write_to("Check $output file and load into geneace.\n") unless ($load);
$log->mail();
exit(0);

###############################################################################################

sub merge_gene {
  if($livegene and $deadgene and $user) {
    print "\nMerging $deadgene into $livegene for USER: $user\n\n";
    my $output = "";
    my $ok = 1; # error status

    # process LIVE gene
    my $livegeneObj = $ace->fetch('Gene', $livegene);
    if ($livegeneObj) {
      # is this a Live gene with no merge from the DEAD gene in the last Acquires_merge?
      my $status = $livegeneObj->Status->name;
      if ($debug) {print "LIVE GENE:$livegeneObj :  ".$livegeneObj->Status->name."\n";}
      if ($status ne 'Live') {
	  $log->error("ERROR: $livegene is not a Live gene\n");
	  $ok = 0;
      }
      # get the last acquires_merge
      my $acquires_merge;
      foreach my $acquires_merge_obj ($livegeneObj->at('Identity.History.Acquires_merge')) {
	if (defined $acquires_merge_obj) {
	  ($acquires_merge) = $acquires_merge_obj->row;
	}
      }
      unless ($override) {
      if (defined $acquires_merge && $acquires_merge eq $deadgene) {
	$log->error("Warning: $livegene has a tag saying it has already merged with $deadgene - this merge will not be done again\n");
	$ok = 0;
      }
    }

      # get the version
      my $ver;
      if (exists $gene_versions{$livegene}) {
	$ver = $gene_versions{$livegene};
      } else {
	$ver = $livegeneObj->Version->name;
      }
      $ver++;
      $gene_versions{$livegene} = $ver;

      $output .= "\nGene : $livegene\nVersion $ver\nHistory Version_change $ver now $user Event Acquires_merge $deadgene\nAcquires_merge $deadgene\n";

    } else {
      $log->error("ERROR: no such gene $livegene\n");
      $ok = 0;
    }

    # process DEAD gene
    my $deadgeneObj = $ace->fetch('Gene', $deadgene);
    if ($deadgeneObj) {
      # is this a Live gene with no merge into the LIVE gene in the last Merged_into?
      if ($debug) {print "DEAD GENE:$deadgeneObj :  ".$deadgeneObj->Status->name."\n";}
      unless ($override) {
	my $status = $deadgeneObj->Status->name;
	if ($status ne 'Live') {
	    $log->error("ERROR: $deadgene is not a Live gene\n");
	    $ok = 0;
	}
      }
      # get the last Merged_into tag
      my $merged_into;
      foreach my $merged_into_obj ($deadgeneObj->at('Identity.History.Merged_into')) {
	if (defined $merged_into_obj) {
	  ($merged_into) = $merged_into_obj->row;
	}
      }
      if (defined $merged_into && $merged_into eq $deadgene) {
	$log->error("Warning: $deadgene has a tag saying it has already merged into $livegene - this merge will not be done again\n");
	$ok = 0;
      }

      # get the CGC name
      my $dead_CGC_name = $deadgeneObj->CGC_name;
      if (defined $dead_CGC_name) {
	$log->error("Warning: $deadgene has a CGC name ($dead_CGC_name) - this merge needs checking as cgc names are involved. Check this with Tim. and only then load the .ace file\n");
	if (defined$load) {undef$load unless $force;}
	$ok = 1;
      }

      my $ver;
      if (exists $gene_versions{$deadgene}) {
	$ver = $gene_versions{$deadgene};
      } else {
	$ver = $deadgeneObj->Version->name;
      }
      $ver++;
      $gene_versions{$deadgene} = $ver;
      $output .= "\nGene : $deadgene\nVersion $ver\nHistory Version_change $ver now $user Event Merged_into $livegene\nMerged_into $livegene\nDead\n\n";
      
      #Stuff to be removed from the dead gene.
      $output .= "\nGene : $deadgene\n-D Map_info\n-D Sequence_name\n-D Allele\n-D Reference\n-D Ortholog\n-D Paralog\n-D Map_info\n-D method\n\n";

      # transfer operon connections.
      my $operon_connect = $deadgeneObj->at('Gene_info.Contained_in_operon');
      if (defined $operon_connect) {
	foreach my $operon_connect_list ($operon_connect->col) {
	  $output .= "\nGene : $livegene\nContained_in_operon $operon_connect_list\n";
	}
      }

      # transfer the Other_names
      my $dead_Other_names_col = $deadgeneObj->at('Identity.Name.Other_name');
      if (defined $dead_Other_names_col) {
	foreach my $dead_Other_names ($dead_Other_names_col->col) {
	  $output .= "\nGene : $livegene\nOther_name $dead_Other_names\n";
	}
      }	


      # transfer Alleles
      foreach my $Alleles ($deadgeneObj->at('Gene_info.Allele')) {
	if (defined $Alleles) {
	  $output .= "\nGene : $livegene\nAllele $Alleles\n";
	}
      }

      # transfer references
      foreach my $references ($deadgeneObj->at('Reference')) {
	if (defined $references) {
	  $output .= "\nGene : $livegene\nReference $references\n";
	}
      }
      # transfer Features
      foreach my $feature ($deadgeneObj->at('Associated_feature')) {
	  if (defined $feature) {
	      $output .= "\nGene : $livegene\nAssociated_feature $feature\n";
	  }
      }
      
      # transfer the Ortholog tags
      foreach my $dead_Orthologs ($deadgeneObj->at('Gene_info.Ortholog')) {
	if (defined $dead_Orthologs) {
	  my @row = $dead_Orthologs->row;
	  my $row0=$row[0]->name;
	  unless (defined $row[1]) {$log->log_and_die("\n\n!!!!!!!Missing ortholog species data in gene $deadgeneObj\n\n");}
	  my $row1=$row[1]->name;
	  unless (defined $row[2]) {$log->log_and_die("\n\n!!!!!!!Missing ortholog Evidence data in gene $deadgeneObj - $row0($row1)\n\n");}
	  my $row2=$row[2]->name;
	  my @col = $deadgeneObj->at("Gene_info.Ortholog.$row0.$row1.$row2")->col;
	  $row1 = '"' . $row1 . '"'; # add quotes to the species name
	  foreach my $col (@col) {
	    $output .= "\nGene : $livegene\nOrtholog $row0 $row1 $row2 $col\n";
	  }
	}
      }	
      
    } else {
      $log->error("ERROR: no such gene $deadgene\n");
      $ok = 0;
    }


    # we did this one successfully
    if ($ok) {
      print ACE $output;
      $count++;
      if ($ns) {
	#nameserver merge
	$log->write_to("NS->merge $livegene $deadgene\n");
	if ($debug) {print "merging $livegene $deadgene";}
	if (my $ids = $db->merge_genes($livegene, $deadgene)){
	  $log->write_to("Merge complete, $deadgene is DEAD and has been merged into gene $livegene\n");
	}
	else {
	  $log->write_to("ERROR: Failed to merge $deadgene into $livegene\n");
	}
      }
    }
  } elsif (!defined($livegene && $deadgene && $user)) {
    $log->write_to("Warning: additional blank line in input file has been ignored\n");
  } else {
    $log->error("missing info on $livegene : $deadgene : $user\n");
  }
  undef $livegene; undef $deadgene ;undef $user;
}

sub load_data {
# load information to $database if -load is specified
$wormbase->load_to_database("$database", "$output", 'batch_merge.pl', $log, undef, 1);
$log->write_to("Loaded $output into $database\n\n");
$wormbase->run_command("mv $output $backupsdir".$out. $wormbase->rundate. "\n"); #append date to filename when moving.
$log->write_to("Output file has been cleaned away like a good little fellow\n\n");
print "Finished!!!!\n";
}

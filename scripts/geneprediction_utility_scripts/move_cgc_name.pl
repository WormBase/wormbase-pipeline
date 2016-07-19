#!/usr/bin/env perl
# moves a cgc name from a gene to another and does the namedb and geneace interactions
# input file looks like:
# FROM_GENEID TO_GENEID
use Ace;
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use NameDB_handler;
use strict;

my ($acedb,$test,$file,$user,$species,$DB);
GetOptions(
           'database:s' => \$acedb,
           'user:s'     => \$user,
           'test'       => \$test,
           'file:s'     => \$file,
           'species:s'  => \$species,
          ) or die;


my $adb= Ace->connect(-path => $acedb)||die("cannot connect to acedb $acedb\n");
$species||='brugia';

if ($test) {
    $DB = 'test_wbgene_id;utlt-db;3307';
  } else {
    $DB = 'nameserver_live;web-wwwdb-core-02;3449';
}

my $db = NameDB_handler->new($DB,$user,$user);
$db->setDomain('Gene');

open IN,$file;
while(<IN>){
       chomp;
       my ($geneid1,$geneid2)=split; # from to

       my $gene1=$adb->fetch(Gene => $geneid1)||die("cannot fetch $geneid1\n");
       my $gene2=$adb->fetch(Gene => $geneid2)||die("cannot fetch $geneid2\n");

       my $gene1update = $db->delName("$gene1",CGC => $gene1->CGC_name);
       $db->_update_public_name("$gene1") if $gene1update;
       die("couldn't delete ${\$gene1->CGC_name} from $gene1\n") unless $gene1update;
       
       my $gene1Version=$gene1->Version+1;
       # remove crap
       print <<HERE;
Gene : $gene1
-D CGC_name
Public_name "${\$gene1->Sequence_name}"
Other_name "${\$gene1->CGC_name}" Curator_confirmed WBPerson4055
Version $gene1Version
Version_change $gene1Version now WBPerson4055 Event Name_change Other_name "${\$gene1->CGC_name}"
-D Gene_class

HERE

       # add crap
       my $gene2update = $db->add_name($gene2,$gene1->CGC_name,'CGC',$species);
       die("couldn't add CGC-name ${\$gene1->CGC_name} to $gene2\n") unless $gene2update;

       my $gene2Version=$gene2->Version+1;

       print <<HERE;
Gene : $gene2
CGC_name "${\$gene1->CGC_name}" Curator_confirmed WBPerson4055
Public_name "${\$gene1->CGC_name}"
Version $gene2Version
Gene_class ${\$gene1->Gene_class}
Version_change $gene2Version now WBPerson4055 Event Name_change CGC_name "${\$gene1->CGC_name}"

HERE
}

#!/usr/local/bin/perl -w


use strict;

use constant USERNAME => 'lstein';
use constant PASSWD   => undef;

use lib '../lib','./lib','./blib';
use NameDB;

my $db = NameDB->connect('name_test',USERNAME,PASSWD);
$db->initialize(1);  # start with a completely clean slate

$db->addDomain(Protein => 'P%010d');
$db->addDomain('Gene');
$db->addDomain('Sequence');

my @domains = $db->getDomains;
print join("\n",@domains),"\n";

$db->setDomain('Protein');
print $db->getDomain,"\n";
print $db->getTemplate,"\n";

$db->setDomain('Gene');
$db->setTemplate('G%010d','Gene');
print $db->getTemplate('Gene'),"\n";

my $result;

$db->setDomain('Protein');
$result  = $db->addNameType('CGC',1,1);
$result  = $db->addNameType('GenBank',1);
$result  = $db->addNameType('EMBL',0,0);
$result  = $db->addNameType('WormPep',0,1);
$result  = $db->addNameType('WTP',0,1);

$result  = $db->setDomain('Gene');
$result  = $db->addNameType('CGC',1,1);
$result  = $db->addNameType('GenBank',1,0);
$result  = $db->addNameType('Locus',0,1);
$result  = $db->addNameType('WTP',0,0);

$db->setDomain('Protein');
my $types = $db->getNameTypes;
print join "\n",@$types,"\n";
print join(",",$db->getNameTypeAttributes('CGC')),"\n";

$db->setDomain('Gene');
my $unc1 = $db->idCreate;
$result  = $db->addName($unc1,CGC => 'unc-1');
$result  = $db->addName($unc1,GenBank  => 'AB123');
$db->addName($unc1,WTP=>'I23Z31');
$db->addName($unc1,WTP => 'MJM32');

my $unc2 = $db->idCreate;
$result  = $db->addName($unc2,CGC => 'unc-2');
$db->addName($unc2,GenBank=>'xyz321');
$db->addName($unc2,WTP => 'I23Z29');
$db->addName($unc2,WTP => 'I23Z30');
$db->addName($unc2,WTP => 'I23Z31');

$db->setDomain('Protein');
my $ace = $db->idCreate;
$db->addName($ace,EMBL => 'ace');
$db->addName($ace,WormPep=>'wp123');
$db->addName($ace,WTP => 'I23Z31');

# we should get one and only one unc-1
$db->setDomain('Gene');
$result = $db->idGetByTypedName(CGC=>'unc-1');
print scalar @$result,"\n";

# we should get two WTP:I23Z31
$result = $db->idGetByTypedName(WTP => 'I23Z31');
print scalar @$result,"\n";

# we should get two wildcards
my @result = $db->idSearch('*');
print scalar @result,"\n";

# in the protein domain, we should get only one WTP:I23Z31
$db->setDomain('Protein');
$result = $db->idGetByTypedName(WTP => 'I23Z31');
print scalar @$result,"\n";

# search 'em all
@result = $db->idGetByAnyName('ace');
print join("\n",@result),"\n";

# test history
$db->setDomain('Gene');
$result = $db->idGetByTypedName(CGC=>'unc-1');
$unc1   = $result->[0];
print_history($db=>$unc1);

# exercise names
my @names = $db->idTypedNames($unc2,'WTP');
print join("\n",@names),"\n";

print_names($db,$unc2);

# the public name (as opposed to public ID) for unc-1 should be "unc-1"
print $db->idPublicName($unc1),"\n";

# the object should be alive
print $db->idLive($unc1),"\n";

# kill the object
$db->idKill($unc1);

# now the object should be dead
print $db->idLive($unc1),"\n";

# resurrect the object
$db->idResurrect($unc1);
print $db->idLive($unc1),"\n";

print_history($db=>$unc1);

# merge $unc1 into $unc2
print_names($db,$unc1);
print_names($db,$unc2);

$db->idMerge($unc1=>$unc2);
print $db->idLive($unc1),"\n";
print $db->idLive($unc2),"\n";
print_names($db,$unc1);
print_names($db,$unc2);

print_history($db=>$unc1);
print "\n";
print_history($db=>$unc2);

$db->setDomain('Protein');
($ace) = $db->idGetByTypedName(EMBL=>'ace');
my $ace2 = $db->idSplit($ace);
$db->addName($ace2,EMBL=>'ace-2');
print_history($db=>$ace);
print_history($db=>$ace2);

my @descendents = $db->idGetChildren($ace);
print join "\n",@descendents,"\n";

print $db->idGetUltimate($ace),"\n";

my @changed = $db->idChanged('2002-04-15');
print "@changed";

sub print_history {
  my ($db,$obj) = @_;
  local $^W = 0;
  for my $event (@{$db->idGetHistory($obj)}) {
    print join("\t",
              'version',
              @{$event}{qw(version event user date related_id name_type name_value)}
              ),"\n";
  }
print "\n";
}

sub print_names {
  my ($db,$id) = @_;
  my %typed_names = $db->idAllNames($id);
  print "$id:\n";
  for my $type (keys %typed_names) {
    my @values = ref $typed_names{$type} ? @{$typed_names{$type}} : $typed_names{$type};
    print "\t",$type,"=>",join(',',@values),"\n";
  }
}



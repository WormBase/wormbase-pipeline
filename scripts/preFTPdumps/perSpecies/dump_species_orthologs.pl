#!/usr/bin/env perl

use Time::localtime;
use IO::File;
use Storable;
use Getopt::Long;

use Wormbase;
use Log_files;

use strict;

my ($store,$debug,$test,$database,$species);
GetOptions(
       'store=s' => \$store,
       'debug=s' => \$debug,
       'test'    => \$test,
       'species=s'  => \$species,
       'database=s' => \$database,
)||die(@!);

my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $wormbase->autoace )||die Ace->error;

my $file = $wormbase->reports . '/'.
   join('.',$wormbase->gspecies_name,$wormbase->ncbi_bioproject,$wormbase->get_wormbase_version,'orthologs.txt');

my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");

print $of 
"# ${\$wormbase->long_name} orthologs\n".
"# WormBase version: " . $dbh->version . "\n".
'# Generated:'.get_date."\n".
'# File is in record format with records separated by "=\n"'."\n".
"#      Sample Record\n".
'#      WBGeneID \t PublicName \n'."\n".
'#      Species \t Ortholog \t MethodsUsedToAssignOrtholog \n'."\n".
'# BEGIN CONTENTS'."\n".
"=\n";

my $i = $dbh->fetch_many(-query=>"find Gene Species=\"${\$wormbase->long_name}\"");
while (my $gene = $i->next) {

    my %orthologs = ();   
    
    # Nematode orthologs
    foreach ($gene->Ortholog) {
	my $methods  = join('; ',map { "$_" } eval { $_->right(2)->col }) ;
	$orthologs{$_->Species} = [ $_,$methods ];
    }

    foreach ($gene->Ortholog_other) {
	my $methods  = join('; ',map { "$_" } eval { $_->right->col });
	$orthologs{$_->Species} = [ $_,$methods ];
    }

    print $of join("\t",$gene,$gene->Public_name),"\n";
    
    foreach (sort keys %orthologs) {
	my ($ortholog,$methods) = @{$orthologs{$_}};
	print $of join("\t",$_,$ortholog,$methods),"\n";
    }
    print $of "=\n";
}

$log->mail;
$of->close;


# yyyy-mm-dd
sub get_date{
 my $tm = localtime;
 my ($day,$month,$year)=($tm->mday,$tm->mon,$tm->year+1900);
 return "$year-$month-$day";
}


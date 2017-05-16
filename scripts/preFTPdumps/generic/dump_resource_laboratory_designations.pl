#!/usr/bin/env perl
use Getopt::Long;
use IO::File;
use Time::localtime;
use Storable;

use Wormbase;
use Log_files;

use strict;

my ($species,$store,$debug,$test,$database,$file);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'file=s'     => \$file,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)
}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $wormbase->autoace)||$log->log_and_die(Ace->error);

$file ||= $wormbase->ftp_site."/staging/releases/WS".$wormbase->get_wormbase_version.'/REPORTS/resource_laboratory_designations.txt';
my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");

my $labs = $dbh->fetch_many(-class=>'Laboratory',-name => '*',-filled=>1) or die $dbh->error;

print $of, "# WormBase Laboratory and Allele Designations\n";
print $of, "# WormBase, version " . $dbh->version . "\n";
print $of, "# Generated ",&get_date(),"\n";
print $of, join('/','#Lab Designation','Allele Designation','Lab Representative','Institute'),"\n";
while (my $obj = $labs->next) {
    my $representative = eval  { $obj->Representative->Full_name };
    my $laboratory     = $obj;
    my $allele         = $obj->Allele_designation;
    my ($institute)    = $obj->Mail;
    print $of,join("\t",$laboratory,$allele,$representative,$institute),"\n";    
}

$log->mail;
$of->close;

# yyyy-mm-dd
sub get_date{
 my $tm = localtime;
 my ($day,$month,$year)=($tm->mday,$tm->mon,$tm->year+1900);
 return "$year-$month-$day";
}



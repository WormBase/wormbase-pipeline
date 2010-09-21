#!/usr/local/bin/perl
# Script written to an ace file into a specified database.

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use File::Basename;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($debug,$file,$database,$store,$test,$wormbase);

GetOptions(
           'debug:s'       => \$debug,       # email output to this user pad
           'file:s'        => \$file,        # the input file to be loaded
	   'database:s'    => \$database ,   # database for data to be loaded, defaults to camace
	   'store:s'       => \$store, 
          );

##Usage:  perl ~pad/Scripts/Load2database.pl -database ~/DATABASES/camace/ -file ~/DATABASES/camace/DUMPS/dump_2010-06-23_A.1.ace -debug pad

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                             );
}

my $log = Log_files->make_build_log($wormbase);

# write out the Homol_data lines for the clones and superlinks used
if (!defined $database){
print 
my $database  = "/nfs/disk100/wormpub/DATABASES/camace";
}

# Main body of script
#---------------------

$log->write_to("Loading $file into $database\n");
$wormbase->load_to_database($database, $file, $debug, $log); #appropriate checks are made in the Wormbase.pm

print "Diaskeda same Poli\n"; #we had alot of fun#
print "είχαμε πολλή διασκέδαση\n"; #better spelling
my $mail = "${debug}\@sanger.ac.uk";
$log->mail($mail);
exit(0);

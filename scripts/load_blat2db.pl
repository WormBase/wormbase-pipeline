#!/software/bin/perl -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2014-08-04 14:38:21 $

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use Log_files;
use Storable;

my ($load,$delete,$dbdir,$help,$all,$wormbase,$test,$debug,$store,$species,$logfile);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
	    "all"     => \$all,
	    "dbdir=s" => \$dbdir,
	    "h"	      => \$help,
	    "debug=s" => \$debug,
            "logfile:s" => \$logfile,
	    "test"    => \$test,
	    "store:s" => \$store,
	    "species:s" => \$species,
	    );


if ( $store ) {
  $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
                             );
}

# establish log file.
my $log = $logfile ? Log_files->make_log($logfile, $debug) : Log_files->make_build_associated_log($wormbase);

my $tace =  $wormbase->tace;
my $blat_dir;
my $acefiles;
my $test_file = $wormbase->blat."/virtual_objects.".$wormbase->species.".blat.EST.".$wormbase->species.".ace";
$log->write_to("Loading BLAT data into $dbdir\n");
# has the build finished?? script needs to look at ~wormpub/BUILD/autoace/BLAT, if it doesnt exist use current_db/BLAT!!
if (!-e $test_file){
  $blat_dir = ($wormbase->basedir('elegans')."/elegans/BLAT");
  $acefiles = ($wormbase->basedir('elegans')."/elegans/acefiles");
  $log->write_to("The build must have finished you are now going to use $blat_dir\n\n");
  my $test_file2 = $blat_dir."/virtual_objects.".$wormbase->species.".blat.EST.".$wormbase->species.".ace";
  unless (-e $test_file2) {
    $log->log_and_die("Could find any BLAT data as $test_file and $test_file2 don't exist....abort\n");
  }
}  
elsif (-e $test_file){
  $log->write_to("The build is still there..... \n\n");
  $blat_dir = $wormbase->blat;
  $acefiles = $wormbase->acefiles;
}

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or all for both> -dbdir <path to your database>\n" if ($help);

#Can you get write access?
my $access = $wormbase->check_write_access($dbdir);
die "You do not have write access for $dbdir\n" if ($access eq "no");


my $dbname;
if ($dbdir =~ /\/(\w+)$/) {
    $dbname = $1;
}
else {
	die print STDERR "$dbdir not valid\n";
}
print "$dbname\n";
&delete($dbname) if ($delete || $all);
&load($dbname)   if ($load   || $all);


$log->mail();
exit(0);


#############################################

sub delete {

my $command=<<END;
query find Homol_data "BLAT_*"
kill
clear
query find Feature_data "Confirmed_intron*"
kill
clear
save
quit
END
    open (DB, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
    print DB $command;
    close DB;
    $log->write_to("BLAT Homold data and Confirmed_introns have been removed.\n");
}

sub load {
my $species = $wormbase->species;
my $command;
if ($dbname =~ /autoace/) {
$command=<<END;
pparse $blat_dir/virtual_objects.$species.blat.EST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.OST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.RST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.mRNA.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.nembase.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.washu.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.nematode.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.EST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.OST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.RST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.mRNA.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.embl.$species.ace
save
pparse $blat_dir/$species.blat.elegans_OST.ace
pparse $blat_dir/$species.blat.elegans_RST.ace
pparse $blat_dir/$species.blat.elegans_mRNA.ace 
pparse $blat_dir/$species.blat.embl.ace 
pparse $blat_dir/$species.blat.nembase.ace
pparse $blat_dir/$species.blat.washu.ace
pparse $blat_dir/$species.blat.nematode.ace
pparse $blat_dir/$species.good_introns.OST.ace
pparse $blat_dir/$species.good_introns.RST.ace
pparse $blat_dir/$species.good_introns.EST.ace 
pparse $blat_dir/$species.good_introns.mRNA.ace 
pparse $blat_dir/$species.good_introns.embl.ace 
save 
pparse $blat_dir/$species.blat.est.ace            
save
quit
END
}

else {
$command=<<END;
pparse $acefiles/chromlinks.ace
pparse $blat_dir/virtual_objects.$species.blat.EST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.OST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.RST.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.mRNA.$species.ace
pparse $blat_dir/virtual_objects.$species.blat.ncRNA.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.EST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.OST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.RST.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.mRNA.$species.ace
pparse $blat_dir/virtual_objects.$species.ci.ncRNA.$species.ace
save
pparse $blat_dir/$species.blat.${species}_EST.ace
pparse $blat_dir/$species.blat.${species}_OST.ace
pparse $blat_dir/$species.blat.${species}_RST.ace
pparse $blat_dir/$species.blat.${species}_mRNA.ace
pparse $blat_dir/$species.blat.${species}_ncRNA.ace
pparse $blat_dir/$species.good_introns.EST.ace 
pparse $blat_dir/$species.good_introns.OST.ace
pparse $blat_dir/$species.good_introns.RST.ace
pparse $blat_dir/$species.good_introns.mRNA.ace
pparse $blat_dir/$species.good_introns.ncRNA.ace
save 
quit
END
}
    open (LOAD, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
    print LOAD $command;
    close LOAD;
$log->write_to("BLAT Homold data and Confirmed_introns have been Updated.\n");
    
}


#!/software/bin/perl -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2008-03-10 13:28:00 $

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use Log_files;
use Storable;

my ($load,$delete,$dbdir,$help,$all,$wormbase,$test,$debug,$store);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
	    "all"     => \$all,
	    "dbdir=s" => \$dbdir,
	    "h"	      => \$help,
	    "debug=s" => \$debug,
	    "test"    => \$test,
	    "store:s" => \$store,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
                                                                                                                              
my $tace =  $wormbase->tace;
my $blat_dir;
my $acefiles;

# has the build finished?? script needs to look at ~wormpub/BUILD/autoace/BLAT, if it doesnt exist use current_db/BLAT!!
if (!-e $wormbase->blat."/virtual_objects.".$wormbase->species.".blat.EST.ace"){
  print "The build must have finished you are now going to use " .$wormbase->database('current')."/BLAT\n\n";
  #$log->write_to( "The build must have finished you are now going to load data from current_db\n");
  $blat_dir = ($wormbase->database('current')."/BLAT");
}  
elsif (-e $wormbase->blat."/virtual_objects.".$wormbase->species.".blat.EST.ace"){
  print "The build is still there..... \n\n";
  $blat_dir = $wormbase->blat;
  $acefiles = $wormbase->acefiles;
}

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or all for both> -dbdir <path to your database>\n" if ($help);

#Can you get write access?
my $access = $wormbase->check_write_access($dbdir);
die "You do not have write access for $dbdir\n" if ($access eq "no");


my $dbname;
if ($dbdir =~ /\/(\w+ace)/) {
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

my $command;
if ($dbname =~ /autoace/) {
$command=<<END;
pparse $blat_dir/virtual_objects.$wormbase->species.blat.EST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.OST.ace
 pparse $blat_dir/virtual_objects.$wormbase->species.blat.RST.ace
pparse $blat_dir/virtual_objects.$wormbase->species.blat.mRNA.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.nembase.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.washu.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.nematode.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.EST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.OST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.RST.ace
pparse $blat_dir/virtual_objects.$wormbase->species.ci.mRNA.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.embl.ace 
save
pparse $blat_dir/$wormbase->species.blat.elegans_OST.ace
pparse $blat_dir/$wormbase->species.blat.elegans_RST.ace
pparse $blat_dir/$wormbase->species.blat.elegans_mRNA.ace 
pparse $blat_dir/$wormbase->species.blat.embl.ace 
pparse $blat_dir/$wormbase->species.blat.nembase.ace
pparse $blat_dir/$wormbase->species.blat.washu.ace
pparse $blat_dir/$wormbase->species.blat.nematode.ace
pparse $blat_dir/$wormbase->species.good_introns.OST.ace
pparse $blat_dir/$wormbase->species.good_introns.RST.ace
pparse $blat_dir/$wormbase->species.good_introns.EST.ace 
pparse $blat_dir/$wormbase->species.good_introns.mRNA.ace 
pparse $blat_dir/$wormbase->species.good_introns.embl.ace 
save 
pparse $blat_dir/$wormbase->species.blat.est.ace            
save
quit
END
}
elsif ($dbname =~ /camace/) {
$command=<<END;
pparse $acefiles/chromlinks.ace
pparse $blat_dir/virtual_objects.$wormbase->species.blat.EST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.OST.ace
pparse $blat_dir/virtual_objects.$wormbase->species.blat.RST.ace
pparse $blat_dir/virtual_objects.$wormbase->species.blat.mRNA.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.blat.ncRNA.ace
pparse $blat_dir/virtual_objects.$wormbase->species.ci.EST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.OST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.RST.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.mRNA.ace 
pparse $blat_dir/virtual_objects.$wormbase->species.ci.ncRNA.ace
save
pparse $blat_dir/$wormbase->species.blat.elegans_EST.ace
pparse $blat_dir/$wormbase->species.blat.elegans_OST.ace
pparse $blat_dir/$wormbase->species.blat.elegans_RST.ace
pparse $blat_dir/$wormbase->species.blat.elegans_mRNA.ace
pparse $blat_dir/$wormbase->species.blat.elegans_ncRNA.ace
pparse $blat_dir/$wormbase->species.good_introns.EST.ace 
pparse $blat_dir/$wormbase->species.good_introns.OST.ace
pparse $blat_dir/$wormbase->species.good_introns.RST.ace
pparse $blat_dir/$wormbase->species.good_introns.mRNA.ace
pparse $blat_dir/$wormbase->species.good_introns.ncRNA.ace
save 
quit
END
}
elsif ($dbname =~ /stlace/) { #this is now defunct as the stlace files aren't created.
$command=<<END;
pparse $blat_dir/virtual_objects.stlace.blat.EST.ace 
pparse $blat_dir/virtual_objects.stlace.blat.OST.ace 
pparse $blat_dir/virtual_objects.stlace.blat.RST.ace 
pparse $blat_dir/virtual_objects.stlace.blat.mRNA.ace 
pparse $blat_dir/virtual_objects.stlace.blat.embl.ace 
pparse $blat_dir/virtual_objects.stlace.blat.nembase.ace 
pparse $blat_dir/virtual_objects.stlace.blat.washu.ace 
pparse $blat_dir/virtual_objects.stlace.blat.nematode.ace 
pparse $blat_dir/virtual_objects.stlace.ci.EST.ace 
pparse $blat_dir/virtual_objects.stlace.ci.OST.ace 
pparse $blat_dir/virtual_objects.stlace.ci.mRNA.ace 
save
pparse $blat_dir/stlace.blat.ost.ace 
pparse $blat_dir/stlace.blat.mrna.ace 
pparse $blat_dir/stlace.blat.embl.ace 
pparse $blat_dir/stlace.good_introns.est.ace 
pparse $blat_dir/stlace.good_introns.ost.ace 
pparse $blat_dir/stlace.good_introns.mrna.ace 
pparse $blat_dir/stlace.good_introns.embl.ace 
save 
pparse $blat_dir/stlace.blat.est.ace            
save
quit
END
}
    open (LOAD, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
    print LOAD $command;
    close LOAD;
$log->write_to("BLAT Homold data and Confirmed_introns have been Updated.\n");
    
}


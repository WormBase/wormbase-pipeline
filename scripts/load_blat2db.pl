#!/usr/local/bin/perl5.6.1 -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch
#
# 030507 : dl  : Update tace queries to restrict deleted/uploaded Homol_data and Feature data class objects
# 030507 : dl  : Upload OST data as part of script

# Last edited by: $Author: pad $
# Last edited on: $Date: 2004-09-07 12:25:21 $

use strict;
use Getopt::Long;
use lib "/wormsrv2/scripts/";
use Wormbase;

my $tace =  &tace;

my ($load,$delete,$dbdir,$help,$all);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
	    "all"     => \$all,
	    "dbdir=s" => \$dbdir,
	    "h"	      => \$help);

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or all for both> -dbdir <path to your database>\n" if ($help);

#Can you get write access?
my $access = &check_write_access($dbdir);
die "You do not have write access for $dbdir\n" if ($access eq "no");


my $dbname;
if ($dbdir =~ /\/(\w+ace)/) {
    $dbname = $1;
}
else {
	die print STDERR "$dbdir not valid\n";
}

&delete($dbname) if ($delete || $all);
&load($dbname)   if ($load   || $all);

exit(0);

#############################################

sub delete {

my $command=<<END;
query find Sequence "SUPERLINK*"
follow Homol_data
kill 
clear
query find Sequence "SUPERLINK*"
follow Feature_data
kill 
clear
query find Sequence "SUPERLINK*"
edit -D Homol_data
edit -D Feature_data
save
quit
END
    open (DB, "| $tace $dbdir |") || die "Couldn't open $dbdir\n";
    print DB $command;
    close DB;
    
}

sub load {

my $command;
if ($dbname =~ /autoace/) {
$command=<<END;
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.blat.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.blat.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.blat.nematode.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.ost.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.mrna.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.embl.ace 
save
pparse /wormsrv2/autoace/BLAT/autoace.blat.ost.ace            
pparse /wormsrv2/autoace/BLAT/autoace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/autoace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/autoace.nematode.ace
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.ost.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.est.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.mrna.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.embl.ace 
save 
pparse /wormsrv2/autoace/BLAT/autoace.blat.est.ace            
save
quit
END
}
elsif ($dbname =~ /camace/) {
$command=<<END;
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.blat.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.blat.ost.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.blat.nematode.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.ost.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.mrna.ace 
save
pparse /wormsrv2/autoace/BLAT/camace.blat.ost.ace            
pparse /wormsrv2/autoace/BLAT/camace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/camace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.est.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.ost.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.mrna.ace 
#pparse /wormsrv2/autoace/BLAT/camace.good_introns.embl.ace 
save 
pparse /wormsrv2/autoace/BLAT/camace.blat.est.ace            
save
quit
END
}
elsif ($dbname =~ /stlace/) {
$command=<<END;
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.blat.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.blat.ost.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.blat.nematode.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.est.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.ost.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.mrna.ace 
save
pparse /wormsrv2/autoace/BLAT/stlace.blat.ost.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.mrna.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.embl.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.est.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.ost.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.mrna.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.embl.ace 
save 
pparse /wormsrv2/autoace/BLAT/stlace.blat.est.ace            
save
quit
END
}
    open (LOAD, "| $tace $dbdir |") || die "Couldn't open $dbdir\n";
    print LOAD $command;
    close LOAD;
    
}


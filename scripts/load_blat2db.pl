#!/usr/local/bin/perl5.6.1 -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch
#
# 030507 : dl  : Update tace queries to restrict deleted/uploaded Homol_data and Feature data class objects
# 030507 : dl  : Upload OST data as part of script

# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2003-05-07 15:19:40 $

use strict;
use Getopt::Long;
use lib "/wormsrv2/scripts/";
use Wormbase;

my $tace =  &tace;

my ($load,$delete,$dbdir,$help);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
	    "dbdir=s" => \$dbdir,
	    "h"	      => \$help);

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database>\n" if ($help);

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

&delete($dbname) if ($delete);
&load($dbname)   if ($load);

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
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/autoace.blat.OST.ace            
pparse /wormsrv2/autoace/BLAT/autoace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/autoace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/autoace.NEMATODE.ace
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.OST.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.EST.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/autoace.good_introns.EMBL.ace 
save 
pparse /wormsrv2/autoace/BLAT/autoace.blat.EST.ace            
save
quit
END
}
elsif ($dbname =~ /camace/) {
$command=<<END;
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/camace.blat.OST.ace            
pparse /wormsrv2/autoace/BLAT/camace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/camace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.EST.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.OST.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.EMBL.ace 
save 
pparse /wormsrv2/autoace/BLAT/camace.blat.EST.ace            
save
quit
END
}
elsif ($dbname =~ /stlace/) {
$command=<<END;
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.OST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/stlace.blat.OST.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.EST.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.OST.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.EMBL.ace 
save 
pparse /wormsrv2/autoace/BLAT/stlace.blat.EST.ace            
save
quit
END
}
    open (LOAD, "| $tace $dbdir |") || die "Couldn't open $dbdir\n";
    print LOAD $command;
    close LOAD;
    
}


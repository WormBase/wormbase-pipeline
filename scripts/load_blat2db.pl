#!/usr/local/bin/perl5.6.0 -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch

use strict;
use Getopt::Long;

my $tace =  "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/tace";

my ($load,$delete,$dbdir,$help);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
		"dbdir=s" => \$dbdir,
		"h"	      => \$help,
);

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database>\n" if ($help);

my $dbname;
if ($dbdir =~ /\/(\w+ace)/) {
	$dbname = $1;
}
else {
	die print STDERR "$dbdir not valid\n";
}

&delete($dbname) if ($delete || !$load);
&load($dbname)   if ($load   || !$delete);



#############################################

sub delete {

my $command=<<END;
find Sequence
edit -D Homol_data
edit -D Feature_data
clear
find Homol_data
kill
clear
find Feature_data
kill
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
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.autoace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/autoace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/autoace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/autoace.blat.nematode.ace 
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
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.camace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/camace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/camace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/camace.blat.nematode.ace 
pparse /wormsrv2/autoace/BLAT/camace.good_introns.EST.ace 
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
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLAT_EMBL.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.BLATX_NEMATODE.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.EST.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/virtual_objects.stlace.ci.EMBL.ace 
save
pparse /wormsrv2/autoace/BLAT/stlace.blat.mRNA.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.EMBL.ace 
pparse /wormsrv2/autoace/BLAT/stlace.blat.nematode.ace 
pparse /wormsrv2/autoace/BLAT/stlace.good_introns.EST.ace 
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


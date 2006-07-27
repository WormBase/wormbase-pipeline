#!/usr/local/bin/perl5.6.1 -w
#
# usage load_blat2db.pl <-options: load/delete or none for both> -dbdir <path to your database> 
#
# deletes existent BLAT data and /or reads in BLAT data for a given database
# 19.02.02 Kerstin Jekosch
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2006-07-27 09:52:55 $
use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};
use Wormbase;

my ($load,$delete,$dbdir,$help,$all,$wormbase,$test,$debug);
GetOptions (
	    "load"    => \$load,
	    "delete"  => \$delete,
	    "all"     => \$all,
	    "dbdir=s" => \$dbdir,
	    "h"	      => \$help,
	    "debug=s" => \$debug);

$wormbase = Wormbase->new( -debug => $debug,
			     -test => $test,
			     );

my $tace =  $wormbase->tace;
my $blat_dir = "/nfs/disk100/wormpub/BUILD/autoace/BLAT";

print STDERR "Give the full path for the database you want to modify!\n" unless ($dbdir);
print STDERR "usage load_blat2db.pl <-options: load/delete or all for both> -dbdir <path to your database>\n" if ($help);

#Can you get write access?
my $access = $wormbase->check_write_access($dbdir);
die "You do not have write access for $dbdir\n" if ($access eq "no");


my $dbname;
#if ($dbdir =~ /\/(\w+ace)/) {
if ($dbdir =~ /\/(\w+cam)/) {
    $dbname = $1;
}
else {
	die print STDERR "$dbdir not valid\n";
}
print "$dbname\n";
&delete($dbname) if ($delete || $all);
&load($dbname)   if ($load   || $all);

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
    open (DB, "| $tace $dbdir |") || die "Couldn't open $dbdir\n";
    print DB $command;
    close DB;
    
}

sub load {

my $command;
if ($dbname =~ /autoace/) {
$command=<<END;
pparse $blat_dir/virtual_objects.autoace.blat.est.ace 
pparse $blat_dir/virtual_objects.autoace.blat.ost.ace 
pparse $blat_dir/virtual_objects.autoace.blat.mrna.ace 
pparse $blat_dir/virtual_objects.autoace.blat.embl.ace 
pparse $blat_dir/virtual_objects.autoace.blat.nembase.ace 
pparse $blat_dir/virtual_objects.autoace.blat.washu.ace 
pparse $blat_dir/virtual_objects.autoace.blat.nematode.ace 
pparse $blat_dir/virtual_objects.autoace.ci.est.ace 
pparse $blat_dir/virtual_objects.autoace.ci.ost.ace 
pparse $blat_dir/virtual_objects.autoace.ci.mrna.ace 
pparse $blat_dir/virtual_objects.autoace.ci.embl.ace 
save
pparse $blat_dir/autoace.blat.ost.ace            
pparse $blat_dir/autoace.blat.mrna.ace 
pparse $blat_dir/autoace.blat.embl.ace 
pparse $blat_dir/autoace.blat.nembase.ace
pparse $blat_dir/autoace.blat.washu.ace
pparse $blat_dir/autoace.blat.nematode.ace
pparse $blat_dir/autoace.good_introns.ost.ace 
pparse $blat_dir/autoace.good_introns.est.ace 
pparse $blat_dir/autoace.good_introns.mrna.ace 
pparse $blat_dir/autoace.good_introns.embl.ace 
save 
pparse $blat_dir/autoace.blat.est.ace            
save
quit
END
}
#elsif ($dbname =~ /camace/) {
elsif ($dbname =~ /PADcam/) {
$command=<<END;
pparse $blat_dir/virtual_objects.camace.blat.est.ace 
pparse $blat_dir/virtual_objects.camace.blat.ost.ace 
pparse $blat_dir/virtual_objects.camace.blat.mrna.ace 
pparse $blat_dir/virtual_objects.camace.blat.embl.ace 
pparse $blat_dir/virtual_objects.camace.blat.nembase.ace 
pparse $blat_dir/virtual_objects.camace.blat.washu.ace 
pparse $blat_dir/virtual_objects.camace.blat.nematode.ace 
pparse $blat_dir/virtual_objects.camace.ci.est.ace 
pparse $blat_dir/virtual_objects.camace.ci.ost.ace 
pparse $blat_dir/virtual_objects.camace.ci.mrna.ace 
save
pparse $blat_dir/camace.blat.ost.ace            
pparse $blat_dir/camace.blat.mrna.ace 
pparse $blat_dir/camace.blat.embl.ace 
pparse $blat_dir/camace.good_introns.est.ace 
pparse $blat_dir/camace.good_introns.ost.ace 
pparse $blat_dir/camace.good_introns.mrna.ace 
save 
pparse $blat_dir/camace.blat.est.ace            
save
quit
END
}
elsif ($dbname =~ /stlace/) {
$command=<<END;
pparse $blat_dir/virtual_objects.stlace.blat.est.ace 
pparse $blat_dir/virtual_objects.stlace.blat.ost.ace 
pparse $blat_dir/virtual_objects.stlace.blat.mrna.ace 
pparse $blat_dir/virtual_objects.stlace.blat.embl.ace 
pparse $blat_dir/virtual_objects.stlace.blat.nembase.ace 
pparse $blat_dir/virtual_objects.stlace.blat.washu.ace 
pparse $blat_dir/virtual_objects.stlace.blat.nematode.ace 
pparse $blat_dir/virtual_objects.stlace.ci.est.ace 
pparse $blat_dir/virtual_objects.stlace.ci.ost.ace 
pparse $blat_dir/virtual_objects.stlace.ci.mrna.ace 
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
    open (LOAD, "| $tace $dbdir |") || die "Couldn't open $dbdir\n";
    print LOAD $command;
    close LOAD;
    
}


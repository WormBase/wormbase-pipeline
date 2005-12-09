#! /usr/local/bin/perl5.8.0 -w

# load_interpolated_map.pl

# Use to load interpolated map positions generated for each build back to geneace 
# The applies to all locus connected to sequence and has no map position

# by Chao-Kung Chen [030625]

# Last updated on: $Date: 2005-12-09 13:55:16 $
# Last updated by: $Author: mt3 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;

my $user = `whoami`;
chomp $user;

if ($user ne "wormpub"){print "You need to be wormpub to run this script . . .\n\n"; exit(0)}

my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB/";

my $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; 
$current =~ s/NAME WS//; chomp $current;
$current++;

my $file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_WS$current*");

print "Uploading interpolated map position to geneace . . . \n";

my $command=<<END;
find locus * where interpolated_map_position
edit -D interpolated_map_position
pparse $file

save
quit
END

my $date = `date +%y%m%d`; chomp $date;
my $log = "/wormsrv2/logs/load_intp_map_to_geneace_from_WS$current.$date";
my $tace = &tace;

my $geneace_dir="/nfs/disk100/wormpub/DATABASES/geneace";
open (Load_GA,"| $tace -tsuser \"interpolated_map\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;




#!/usr/local/bin/perl5.8.0 -w

# load_interpolated_map_to_geneace.pl

# by Chao-Kung Chen

# Last updated on: $Date: 2003-08-01 10:41:39 $
# Last updated by: $Author: ck1 $


# Automatically update Geneace with interpolated map

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;

my $tace = &tace;  
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME //; chomp $current;
my $map = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_$current.*.ace");
my $ga = "/wormsrv1/geneace/";

print $map, "\n";
my $command=<<END;
pparse $map
save
quit
END

open (Load_GA,"| $tace -tsuser \"interpolated_map\" $ga") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;

exit(0);


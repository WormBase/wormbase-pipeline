#!/usr/local/bin/perl5.8.0 -w

# load_interpolated_map_to_geneace.pl

# by Chao-Kung Chen

# Last updated on: $Date: 2003-12-08 13:34:35 $
# Last updated by: $Author: ck1 $


# Automatically update Geneace with interpolated map

# Run this script AFTER the build is released, or 
# anytime during the build when get_interpolated_map step is finis#hed


use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;

my $tace = &tace;  
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME WS//; chomp $current;
my $autoace = $current+1; $autoace = "WS$autoace";

my $map = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_$autoace.*.ace");

my $command=<<END;

Find Locus * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map
save
quit

END


my $ga = "/wormsrv1/geneace/";

# upload when corresponding file exists
if ($map){
  open (Load_GA,"| $tace -tsuser \"interpolated_map\" $ga") || die "Failed to upload to Geneace";
  print Load_GA $command;
  close Load_GA;
}


exit(0);


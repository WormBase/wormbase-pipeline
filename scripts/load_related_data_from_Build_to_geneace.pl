#!/usr/local/bin/perl5.8.0 -w

# load_interpolated_map_to_geneace.pl

# by Chao-Kung Chen

# Last updated on: $Date: 2004-02-25 12:58:05 $
# Last updated by: $Author: ck1 $


# load Geneace related data from Build back to Geneace

# Run this script AFTER the build is released, or 
# anytime during the build when get_interpolated_map and update_inferred multi-pt data 


use strict;                    
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

#######################################################
# load Geneace related data from Build back to Geneace
#######################################################

my $tace = &tace;
my $geneace_dir = "/wormsrv1/geneace";
my $autoace = get_wormbase_version(); # only the digits   
my $log = "/wormsrv2/logs/load_WS".$autoace."_build_data_to_geneace";
open(LOG, ">$log") || die $!;
print LOG "\n\nLoading WS$autoace Geneace related data back to geneace\n";
print LOG "----------------------------------------------------------------\n";

# interpolated map data
my @map = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_WS$autoace.*.ace");

# corrected reverse physicals
my $rev_pyhs = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/rev_physical_update_$autoace");

# inferred multi-pt data
my $multi = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/inferred_multi_pt_obj_$autoace");

# updated person from Caltech
system("perl /wormsrv2/scripts/update_person.pl");

my $command;

if ($rev_pyhs && $multi){
  print LOG "\n\nLoading WS$autoace interpolated map data / corrected reverse physicals / inferred multi-pt data back to geneace . . .\n\n";
  $command=<<END;
Find Locus * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map[-1]
pparse $rev_pyhs
pparse $multi
save
quit
END
	      }
if ($rev_pyhs && !$multi){
  print LOG "\nLoading WS$autoace interpolated map data / corrected reverse physicals back to geneace . . .\n\n";
  $command=<<END;
Find Locus * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map[-1]
pparse $rev_pyhs
save
quit
END
}

if (!$rev_pyhs && $multi){
  print LOG "\nLoading WS$autoace interpolated map data / inferred multi-pt data back to geneace . . .\n\n";
$command=<<END;
Find Locus * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map[-1]
pparse $multi
save
quit
END
}

my $rundate = &rundate;
print LOG "\nLoading updated Person info from Caltech\n\n";
my @person_log = `cat /wormsrv2/logs/update_Person.$rundate`;
print LOG @person_log;

# upload when corresponding file exists
if ($map[-1]){
  open (Load_GA,"| $tace -tsuser \"interpolated_map\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
  print Load_GA $command;
  close Load_GA;
}
else {
  print LOG "\n\nInterpolated_map_file does not exist!\n";
}

mail_maintainer("Loading WS$autoace Geneace related data to Geneace", "ck1\@sanger.ac.uk", $log);

#!/usr/local/bin/perl5.8.0 -w

# load_interpolated_map_to_geneace.pl

# by Chao-Kung Chen

# Last updated on: $Date: 2004-05-28 12:56:05 $
# Last updated by: $Author: ck1 $


# load Geneace related data from Build back to Geneace

# RUN this script anytime during the build when get_interpolated_map and update_inferred multi-pt data are done


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

######################
# ----- globals -----
######################

my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");          # tace executable path
#my $geneace_dir = "/wormsrv1/geneace";
my $geneace_dir = "/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST";
my $autoace = get_wormbase_version(); # only the digits
my $log = "/wormsrv2/logs/load_WS".$autoace."_build_data_to_geneace";

open(LOG, ">$log") || die $!;
print LOG "\n\nLoading WS$autoace Geneace related data back to $geneace_dir\n";
print LOG "---------------------------------------------------------------------------------------------\n";


###################################################################
# ----- load Geneace related data from Build back to Geneace -----
###################################################################

# (1) interpolated map data
my $map = "/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_WS$autoace.ace";

# (2) corrected reverse physicals
my $rev_pyhs = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/rev_physical_update_WS$autoace");

# inferred multi-pt data

# (3) new multipt obj created for pseudo markers
my $multi = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/inferred_multi_pt_obj_WS$autoace");

# (4) existing multipt obj with updated flanking marker loci
my $multi_update = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/updated_multi_pt_flanking_loci_WS$autoace");


# (5) updated person from Caltech
system("perl /nfs/team71/worm/ck1/WORMBASE_CVS/scripts/GENEACE/update_person.pl");
#system("perl /wormsrv2/scripts/update_person.pl");

# (6) allele mapping data from build
my $allele_delete = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update_delete.WS$autoace.ace";
my $allele_update = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update.WS$autoace.ace";


################################################
# ----- short message about files to load -----
################################################

my $msg = <<MSG;
\nNOTE: The following 8 files should have been loaded.\n
(1) $map\n
(2) $rev_pyhs\n
(3) $multi\n
(4) $multi_update\n
(5) $allele_delete\n
(6) $allele_update\n
(7) /tmp/new_Person.ace\n
(8) /tmp/new_Person_name.ace\n
Depending on whether pseudo genetic markers are availabe to the build or if there are reverse physicals, files like (2) or (3) maybe missing. But the rest should always be available\n
MSG

print LOG $msg;

######################
# ----- Loading -----
######################

my  $command=<<END;
Find Gene * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map
pparse $rev_pyhs
pparse $multi
pparse $multi_update
pparse $allele_delete
pparse $allele_update

save
quit
END

my $rundate = &rundate;
print LOG "\nLoading updated Person info from Caltech\n\n";
my @person_log = `cat /wormsrv2/logs/update_Person.$rundate`;
print LOG @person_log;

# upload when corresponding file exists
open (Load_GA,"| $tace -tsuser \"interpolated_map\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;


###########################
# ----- email notice -----
###########################

mail_maintainer("Loading WS$autoace Geneace related data to Geneace - test only", "ALL", $log);

__END__

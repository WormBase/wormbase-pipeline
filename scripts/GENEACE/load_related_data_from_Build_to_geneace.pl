#!/usr/local/bin/perl5.8.0 -w

# load_interpolated_map_to_geneace.pl

# by Chao-Kung Chen

# Last updated on: $Date: 2004-06-09 12:50:41 $
# Last updated by: $Author: ck1 $


# load Geneace related data from Build back to Geneace

# RUN this script anytime during the build or after the build when get_interpolated_map and update_inferred multi-pt data are done


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

######################
# ----- globals -----
######################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){print "\nYour need to be wormpub to upload data to geneace\n"; exit 0 };

my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");          # tace executable path
my $build_ver = get_wormbase_version(); # only the digits

my $geneace_dir = "/wormsrv1/geneace";
my $autoace = "/wormsrv2/autoace";
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";

my $log = "/wormsrv2/logs/load_WS".$build_ver."_build_data_to_geneace";

open(LOG, ">$log") || die $!;
print LOG "\n\nLoading WS$build_ver Geneace related data back to $geneace_dir\n";
print LOG "---------------------------------------------------------------------------------------------\n";


##############################
# ----- preparing data -----
##############################

# (1) interpolated map data
my @map = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_map_to_geneace_WS$build_ver.*ace");
my $map = $map[-1];

# (2) corrected reverse physicals
my $rev_pyhs = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/rev_physical_update_WS$build_ver");

# inferred multi-pt data

# (3) new multipt obj created for pseudo markers
my $multi = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/inferred_multi_pt_obj_WS$build_ver");

# (4) existing multipt obj with updated flanking marker loci
my $multi_update = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/updated_multi_pt_flanking_loci_WS$build_ver");


# (5) updated geneace with person/person_name data from Caltech

######################################################################
# ----- retrieve  person/person_name datafrin autoace Or curr_db -----
######################################################################

my $person = "/tmp/new_Person.ace";
my $person_name = "/tmp/new_Person_name.ace";

my $dump_Person_Person_name_from_caltech=<<END;
find Person *
show -a -f $person
find Person_name *
show -a -f  $person_name
quit
END

# Two possibilities of doing ?Person/?Person_name update from Caltech:
# (1) duing the build:
# (2) after the build: retrieved from curr_db as Person/Person_name data will not be in autoace anymore

# doing update during the build
`echo "$dump_Person_Person_name_from_caltech" | $tace $autoace`;
`chmod 777 $person $person_name`;

# doing update after the build
if (! -e ($person && $person_name)){
  print LOG "Autoace Person and Perosn_name data are now in curr_db .. retrieving them from current_DB (WS$build_ver)...\n\n";
  `echo "$dump_Person_Person_name_from_caltech" | $tace $curr_db`;
  print LOG "Got Person and Perosn_name data from current_DB (WS$build_ver)...";
}

#######################################################################################
# ----- remove person/person_name data from Geneace : 
#       the value of "CGC_representative_for" is kept as geneace keeps the record -----
#######################################################################################

my $remove_Person_Person_name_from_ga=<<END;
find Person *
edit -D PostgreSQL_id
edit -D Name
edit -D Laboratory
edit -D Address
edit -D Comment
edit -D Tracking
edit -D Lineage
edit -D Publication

find Person_name *
edit -D Name
save
quit
END

`echo "$remove_Person_Person_name_from_ga" | $tace $geneace_dir`;


# (6) allele mapping data from build

my $allele_delete = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update_delete.WS$build_ver.ace";
my $allele_update = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update.WS$build_ver.ace";


################################################
# ----- short message about files to load -----
################################################

my $msg = <<MSG;
\n\nNOTE: The following 8 files should have been loaded.\n
(1) $map\n
(2) $rev_pyhs\n
(3) $multi\n
(4) $multi_update\n
(5) $allele_delete\n
(6) $allele_update\n
(7) $person\n
(8) $person_name\n
Depending on whether pseudo genetic markers are availabe to the build or if there are reverse physicals, files like (2) or (3) maybe missing. But the rest should always be available\n
MSG

print LOG $msg;

###################################################################
# ----- load Geneace related data from Build back to Geneace -----
###################################################################

my  $command=<<END;
Find Gene * where Interpolated_map_position
edit -D Interpolated_map_position
pparse $map
pparse $rev_pyhs
pparse $multi
pparse $multi_update
pparse $allele_delete
pparse $allele_update
pparse $person
pparse $person_name

save
quit
END

my $rundate = &rundate;

# upload when corresponding file exists
open (Load_GA,"| $tace -tsuser \"autoace_update\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;

#######################
# ----- clean up -----
#######################

# remove temp. person, person_name data
`rm -f $person $person_name`;


###########################
# ----- email notice -----
###########################

#mail_maintainer("Loading WS$build_ver Geneace related data to Geneace", "ck1\@sanger.ac.uk, krb\@sanger.ac.uk", $log);

mail_maintainer("Loading WS$build_ver Geneace related data to Geneace", "ck1\@sanger.ac.uk", $log);

__END__


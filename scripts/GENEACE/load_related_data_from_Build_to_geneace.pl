#!/usr/local/bin/perl5.8.0 -w
#
# load_interpolated_map_to_geneace.pl
#
# by Chao-Kung Chen
#
# loads Geneace related data from Build back to /wormsrv1/geneace
# RUN this script anytime during the build or after the build when get_interpolated_map 
# and update_inferred multi-pt data are done
#
# Last updated on: $Date: 2004-07-14 09:19:32 $
# Last updated by: $Author: krb $


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
# can use dumped Person class in /wormsrv2/wormbase/caltech/caltech_Person.ace

# First need to remove person/person_name data from /wormsrv1/geneace
# Not that the value of "CGC_representative_for" is kept as geneace keeps this record
# i.e. you can't delete *all* of the Person class from geneace

my $remove_Person_data_from_geneace=<<END;
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

`echo "$remove_Person_data_from_geneace" | $tace $geneace_dir`;

# new Person data will have been dumped from citace
my $person = "/wormsrv2/wormbase/caltech/caltech_Person.ace";

# (6) allele mapping data from build

my $allele_delete = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update_delete.WS$build_ver.ace";
my $allele_update = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update.WS$build_ver.ace";


################################################
# ----- short message about files to load -----
################################################

my $msg = <<MSG;
\n\nNOTE: The following 7 files should have been loaded.\n
(1) $map\n
(2) $rev_pyhs\n
(3) $multi\n
(4) $multi_update\n
(5) $allele_delete\n
(6) $allele_update\n
(7) $person\n

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
save
quit
END

# upload when corresponding file exists
open (Load_GA,"| $tace -tsuser \"autoace_update\" $geneace_dir >> $log") || die "Failed to upload to Geneace\n";
print Load_GA $command;
close Load_GA;


###########################
# ----- email notice -----
###########################

&mail_maintainer("Loading WS$build_ver Geneace related data to Geneace", "krb\@sanger.ac.uk", $log);

exit(0);

__END__


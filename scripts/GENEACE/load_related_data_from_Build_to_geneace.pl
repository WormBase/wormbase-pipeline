#!/usr/local/bin/perl5.8.0 -w
#
# load_related_data_from_Build_to_geneace.pl
#
# by Chao-Kung Chen
#
# loads Geneace related data from Build back to /nfs/disk100/wormpub/DATABASES/geneace
# RUN this script anytime during the build or after the build when get_interpolated_map 
# and update_inferred multi-pt data are done
#
# Last updated on: $Date: 2006-04-25 14:14:26 $
# Last updated by: $Author: ar2 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

######################
# ----- globals -----
######################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){print "\nYour need to be wormpub to upload data to geneace\n"; exit 0 };

my $wormbase = Wormbase->new();
my $tace = $wormbase->tace;          # tace executable path
my $release = $wormbase->get_wormbase_version_name(); # only the digits

my $geneace     = $wormbase->database('geneace');
my $autoace     = $wormbase->autoace;
my $curr_db     = $wormbase->database('current');


my $log = Log_files->make_build_log($wormbase);


##############################
# ----- preparing data -----
##############################

# (2) interpolated map data
$log->write_to("Loading interpolated map data\n");
#Generate the file from autoace first
my $int_map_pos = $wormbase->acefiles."/interpolated_map_position.ace";
my $command  = "query find Gene Interpolated_map_position\nshow -a -t Interpolated_map_position -f $int_map_pos\nquit";

open (INT_A,"| $tace ".$wormbase->autoace) or $log->log_and_die("Failed to get data from autoace\n");
print INT_A $command;
close INT_A;

#need to refresh all Interpolated map positions as they are all relative to each other.
$wormbase->load_to_database($geneace,$int_map_pos,'interpolated_map_positions_from_autoace');

# (1) Promoted map positions
$log->write_to("Loading pseudo map positions\n");
my $file = $wormbase->acefiles."/pseudo_map_positions.ace";

$wormbase->load_to_database($geneace, $file, 'pseudo_map_positions_from_autoace');


# (6) updated geneace with person/person_name data from Caltech
# can use dumped Person class in /wormsrv2/wormbase/caltech/caltech_Person.ace
$log->write_to("Updating person name information from caltech_Person.ace file\n");

# First need to remove person/person_name data from /nfs/disk100/wormpub/DATABASES/geneace
# Not that the value of "CGC_representative_for" is kept as geneace keeps this record
# i.e. you can't delete *all* of the Person class from geneace
$log->write_to("First removing old Person data\n");
$command=<<END;
find Person *
edit -D PostgreSQL_id
edit -D Name
edit -D Laboratory
edit -D Address
edit -D Comment
edit -D Tracking
edit -D Lineage
edit -D Publication
save
quit
END

open (Load_GA,"| $tace -tsuser \"person_update_from_autoace\" $geneace") || die "Failed to upload to Geneace\n";
print Load_GA $command;
close Load_GA;


# new Person data will have been dumped from citace
$log->write_to("Adding new person data\n");
my $person = $wormbase->acefiles."/primaries/caltech/caltech_Person.ace";
$wormbase->load_to_database($geneace, $person,"caltech_Person");


# new Person data will have been dumped from citace
#need to remove previous data
$command = "find Paper\nkill\nsave\nquit";
open (GA,"| $tace -tsuser \"remove_papers\" $geneace") or $log->log_and_die("Failed to edit to Geneace - paper removal\n");
print GA $command;
close GA;

#load the new stuff
$log->write_to("Adding new paper data\n");
my $paper = $wormbase->acefiles."/primaries/caltech/caltech_Paper.ace";
$wormbase->load_to_database($geneace, $paper,"caltech_Paper");

$log->mail();
exit(0);

__END__


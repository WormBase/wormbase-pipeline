#!/usr/local/bin/perl5.8.0 -w
#
# load_related_data_from_Build_to_geneace.pl
#
# by Chao-Kung Chen
#
# loads Geneace related data from Build back to /nfs/wormpub/DATABASES/geneace
# RUN this script anytime during the build or after the build when get_interpolated_map 
# and update_inferred multi-pt data are done
#
# Last updated on: $Date: 2011-04-18 09:03:17 $
# Last updated by: $Author: klh $


use strict;
use Getopt::Long;

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

my $db;
&GetOptions('db=s' => \$db);

######################
# ----- globals -----
######################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){print "\nYou need to be wormpub to upload data to geneace\n"; exit 0 };

my $wormbase = (defined $db)
    ? Wormbase->new(-autoace => $db)
    : Wormbase->new();
my $tace     = $wormbase->tace;          # tace executable path
my $release  = $wormbase->get_wormbase_version_name(); # only the digits

my $geneace     = $wormbase->database('geneace');

my $log = Log_files->make_build_log($wormbase);


##############################
# ----- preparing data -----
##############################
my $accept_large_differences = 1; 
my $command;

#
# Map data 1 = interpolated map data
#
$log->write_to("Loading interpolated map data\n");

my @int_map_files = glob($wormbase->acefiles . "/interpolated_gene_*.ace");
foreach my $mfile (@int_map_files) {
  $wormbase->load_to_database($geneace,$mfile,'interpolated_map_positions_from_autoace',$log, 0, $accept_large_differences);
}

#
# Map data 2 = Promoted map positions produced by interpolation manager
#
$log->write_to("Loading pseudo map positions\n");
my $file = $wormbase->acefiles."/pseudo_map_positions.ace";
$wormbase->load_to_database($geneace, $file, 'pseudo_map_positions_from_autoace', $log, 0, $accept_large_differences);

#
# Map data 3 = Genetic map fixes produced by interpolation manager
#
$log->write_to("Loading genetic map fixes\n");
$file = $wormbase->acefiles."/genetic_map_fixes.ace";
$wormbase->load_to_database($geneace, $file, 'genetic_map_fixes_from_autoace', $log, 0, $accept_large_differences);


#
# Update geneace with person/person_name data from Caltech
# 
my $person = $wormbase->acefiles."/primaries/citace/caltech_Person.ace";
if (-e $person) {

  $log->write_to("Updating person name information from caltech_Person.ace file\n");

  # 
  # First need to remove person/person_name data from geneace
  # Note that the value of "CGC_representative_for" is kept as geneace keeps this record
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

  #
  # new Person data will have been dumped from citace
  #
  $log->write_to("Adding new person data\n");

  $wormbase->load_to_database($geneace, $person,"caltech_Person",$log);
} else {
  $log->write_to("NOT updating person name information - could not find file $person\n");
}

#
# new Paper data will have been dumped from citace
#
my $paper = $wormbase->acefiles."/primaries/citace/caltech_Paper.ace";
if (-e $paper) {
  $log->write_to("Adding new paper data\n");

  $wormbase->load_to_database($geneace, $paper,"caltech_Paper",$log);
} else {
  $log->write_to("NOT updating Paper information - could not find file $paper\n");
}

$log->mail();
exit(0);

__END__


#!/usr/local/bin/perl5.8.0 -w

# update_Person_in_geneace.pl

# by Chao-Kung Chen [030728]

# Last updated on: $Date: 2005-12-09 13:55:16 $
# Last updated by: $Author: mt3 $


# Automatically update Geneace with Person/Person_name classes from autoace

use strict;                    
use lib "/wormsrv2/scripts/"; 
use Wormbase;

###########################
# variables
###########################

my $user = `whoami`;
chomp($user);

if ($user ne "wormpub"){print "You need to be wormpub to run this update!\n"; exit(0)}

my $tace = &tace;
my $autoace = "/wormsrv2/autoace";
my $ga = "/nfs/disk100/wormpub/DATABASES/geneace";

my $rundate = `date +%y%m%d`; chomp $rundate;
my $log = "/wormsrv2/logs/update_Person.$rundate";


my $dump_PI=<<END;
find Laboratory * 
show -a -T Representative -f /tmp/lab_PI.ace
quit
END

`echo "$dump_PI" | $tace $ga`;

my $remove_Person_Person_name_from_ga=<<END;
find Person * 
edit -D Name
edit -D CGC_representative_for
edit -D Laboratory
edit -D Address
edit -D Tracking
edit -D Publication
edit -D Comment
find Person_name *
edit -D Name
save
quit
END

`echo "$remove_Person_Person_name_from_ga" | $tace $ga`;

my $dump_Person_Person_name_from_autoace=<<END;
find Person * 
show -a -f /tmp/new_Person.ace
find Person_name *
show -a -f /tmp/new_Person_name.ace
quit
END

`echo "$dump_Person_Person_name_from_autoace" | $tace $autoace`;

my $parse_files=<<END;
pparse /tmp/lab_PI.ace 
pparse /tmp/new_Person.ace
pparse /tmp/new_Person_name.ace
save
quit
END

open (Load_GA,"| $tace -tsuser \"Person_update\" $ga > $log") || die "Failed to upload to Geneace";
print Load_GA $parse_files;
close Load_GA;

`rm -f /tmp/lab_PI.ace /tmp/new_Person.ace /tmp/new_Person_name.ace`;

print "Person and Person_name classes updated. Check /wormsrv2/logs/update_Person.$rundate to see if everything is ok!\n";

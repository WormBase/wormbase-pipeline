#!/usr/local/bin/perl5.8.0 -w

# update_Person_in_geneace.pl

# by Chao-Kung Chen [030728]

# Last updated on: $Date: 2004-05-19 08:58:37 $
# Last updated by: $Author: ck1 $


# Automatically update Geneace with Person/Person_name classes from autoace

use strict;

###########################
# variables
###########################

my $user = `whoami`;
chomp($user);

if ($user ne "wormpub"){print "You need to be wormpub to run this update!\n"; exit(0)}

my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");          # tace executable path
my $autoace = "/wormsrv2/autoace";
my $ga = "/wormsrv1/geneace";

my $rundate = `date +%y%m%d`; chomp $rundate;
my $log = "/wormsrv2/logs/update_Person.$rundate";


# updates geneace ?Person and ?Person_name classes from caltech via autoace

# the CGC_representative_for is kept as geneace keeps the record
my $remove_Person_Person_name_from_ga=<<END;
find Person *
edit -D Name
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

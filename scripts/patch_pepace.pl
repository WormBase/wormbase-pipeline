#!/usr/local/bin/perl -w


use strict;
use lib '/wormsrv2/scripts';
use Wormbase;

my $WP_number = &get_wormpep_version;
my $WP_previous = $WP_number -1;

my $tace = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";

my $command=<<END;
query Find Protein; Live
edit -D Live
pparse /wormsrv2/WORMPEP/wormpep${WP_number}/patch_wormpep.${WP_number}-${WP_previous}.ace
save
quit
END
 
open (TACE,"| $tace /wormsrv2/pepace") || die "Couldn't open pipe to tace\n";
print TACE $command;
close (TACE);

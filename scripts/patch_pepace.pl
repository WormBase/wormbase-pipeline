#!/usr/local/bin/perl -w


use strict;
use lib '/wormsrv2/scripts';
use Wormbase;

my $WS_number = &get_wormbase_version;
my $WS_previous = $WS_number -1;

my $tace = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";

my $command=<<END;
pparse /wormsrv2/WORMPEP/wormpep${WS_number}/patch_wormpep.${WS_number}-${WS_previous}.ace
save
quit
END
 
open (TACE,"| $tace /wormsrv2/pepace") || die "Couldn't open pipe to tace\n";
print TACE $command;
close (TACE);

#!/usr/local/bin/perl -w


use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use Ace;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;

my $log = "/wormsrv2/logs/$0.log.$rundate";
open(LOG,">$log")|| die "cant open $log";
print LOG "$0 started at ",`date`,"\n";
print LOG "=============================================\n";



my $WS_number = &get_wormbase_version;
my $WS_previous = $WS_number -1;

my $tace = "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/tace";

my $command=<<END;
pparse /wormsrv2/WORMPEP/wormpep${WS_number}/patch_wormpep.${WS_number}-${WS_previous}.ace
save
quit
END

print LOG "parsing patch_wormpep.${WS_number}-${WS_previous}.ace to /wormsrv2/pepace\n";
 
#open (TACE,"| $tace /wormsrv2/pepace") || die "Couldn't open pipe to tace\n";
#print TACE $command;
#close (TACE);

print LOG "parsing complete - about to check . . . . . ";

################################################# 
 ##############  test  pepace ##################
#################################################

my $query = <<END;
Find Protein; Live; !Corresponding_DNA
END


my $error;

my $db = Ace->connect('/wormsrv2/pepace') or die ("Could not connect with pepace\n");    
my @oddies = $db->fetch(-query=>"$query");
foreach my $pep(@oddies) {
    $error .= "$pep is Live but has no Corresponding_DNA\n";
  }

 $query = <<END;
Find Protein; Corresponding_DNA; !Live
END
@oddies = $db->fetch(-query=>"$query");
print LOG "\n\n";
foreach my $pep(@oddies) {
    $error .=  "$pep has Corresponding_DNA but is not Live\n";
  }

$db->close;

if (defined($error))
  {
    print LOG "\n\n$error\n\n";
  }
else{
  print LOG "ok!\n\n";}

print LOG "$0 finished at `date`\n";

close LOG;

$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer($0,$maintainer,$log);

#!/usr/local/bin/perl5.6.1 -w

# update_caltech.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-01-24 13:02:03 $
# Last updated by: $Author: ck1 $


# Automatically update Geneace with Erich's functional annotation update

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';

#######
# usage
#######

if (!$ARGV[0]){
  print "\nUsage: perl5.6.1 pheno_annotate.pl filename\n(First copied the file to /wormsrv1/geneace/ERICHS_DATA)\n\n";
  exit(0);
}

my $update_file = "$ARGV[0]";


#######################
# check user is wormpub
#######################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\nYou have to be wormpub to run this script!\n\n";
  exit(0);
}

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");


my $command=<<END;
find sequence * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/seq_TS_dump.ace

find locus * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/loci_TS_dump.ace

edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

save
quit
END

my $tace = &tace;

my $geneace_dir="/wormsrv1/geneace";

open (FH,"| $tace $geneace_dir ") || die "Failed to upload to Geneace";
print FH $command;
close FH;

print $update_file, "\n";


$command="pparse $update_file\nsave\nquit\n";

open (Load_GA,"| $tace $geneace_dir ") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;

if($!){
  print "##########################################\n";
  print "\nPhenotype annotation is now updated!\n\n";       
  print "\nIf everthing is OK, REMEMBER to remove\n"; 
  print "loci_TS_dump.ace and seq_TS_dump.ace\n";
  print "in /wormsrv1/geneace/ERICHS_DATA\n\n";
  print "##########################################\n\n"; 
}
else{
  print "##############################\n";
  print "Mission failed, Mr. Bond!\n";
  print "##############################\n\n";
}

__END__


#perl5.6.1 ~ck1/WORMBASE_CVS/scripts/update_caltech.pl /wormsrv1/geneace/ERICHS_DATA/annots-22jan2003.ace

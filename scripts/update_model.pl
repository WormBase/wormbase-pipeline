#!/usr/local/bin/perl5.8.0 -w
# 
# update_model.pl
#
#
# Script to update models.wrm in various directories for the build
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Getopt::Long;
use Log_files;

# ----- specify status of model
my $new;
GetOptions ("n|new" => \$new);

my $log = Log_files->make_build_log();

$log->write_to("Updating the models to start new build\n");
$log->write_to("\n* YOU SHOULD UPDATE PRIMARY DATABASE MODELS TOO *\n\n");
# ----- update models.wrm in various of build dir. under /wormsrv2 only when there are changes

if ($new){
  $log->write_to("Updating in . . \n");
  my @db = qw ("autoace" "geneace" "citace" "cshace" "camace" "stlace" "brigace");
  foreach(@db){
    $log->write_to("\t$_\n");
    system("cd /wormsrv2/$_/wspec/; cvs update");
  }
}

system("cd /wormsrv2/scripts/; cvs update");

print "\n\nmodel.wrm files upated for the build..\n\n";

$log->mail;

exit(0);

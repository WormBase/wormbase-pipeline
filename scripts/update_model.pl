#!/usr/local/bin/perl5.8.0 -w
# 
# update_model.pl
#
#
# Script to update models.wrm in various directories for the build
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-03-02 10:48:19 $


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Getopt::Long;

# ----- specify status of model
my $new;
GetOptions ("n|new" => \$new);

# ----- update models.wrm in various of build dir. under /wormsrv2 only when there are changes

if ($new){
  my @db = qw ("autoace" "geneace" "citace" "cshace" "camace" "stlace" "brigace");
  foreach(@db){
    system("cd /wormsrv2/$_/wspec/; cvs update");
  }
}

system("cd /wormsrv2/scripts/; cvs update");

print "\n\nmodel.wrm files upated for the build..\n\n";

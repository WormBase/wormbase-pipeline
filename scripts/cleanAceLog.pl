#!/usr/local/bin/perl
#
# RD 991015
#
#

use strict ;

(my $acedir = shift @ARGV) || die "Usage: cleanAceLog.pl <acedir>\n" ;

open (LOG, "$acedir/database/log.wrm") || die "failed to open $acedir/database/log.wrm\n" ;

while (<LOG>) {
    next if /^\s*$/ ;		# remove blank lines

    s/\/\/ // ;

    next if /\*\*\*\*\*\*\*\*\*/ ;
    next if /New start/ ;
    next if /Readlock manager : session/ ;

				# remove large object messages
    next if /Class \S+, object .*cells.$/ ;
    next if /^This is just a warning, acedb has no hard limits on the mumber of cells per object,/ ;
    next if /^but the performance degrades on very large objects, and it possible you should / ;
    next if /^reconsider your models design.  In particular, you may be overusing XREF or/ ;
    next if /^This message is issued once per offending class and per code run/ ;

				# remove redundant saving messages
    next if /Saving Session/ ;
    next if /Saving done/ ;

				# write access messages
    next if /Write access/ ;

				# disk extension
    next if /in partition block/ ;
    next if /Disk extend: Creating partition/ ;


    if (/error/ || /Error/ || /ERROR/) { print "* " ; } else { print "  " ; }
    print ;
}


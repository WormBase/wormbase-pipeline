#!/usr/local/bin/perl
#
# calculates cosmid coordinates from chrom coords
# usage: lookup_coords.pl -link <link name> -coord <number> -camace/stlace/autoace
#
# 28.2.2002 Kerstin Jekosch

use strict;
use Wormbase;
use Getopt::Long;
use Ace;

##################
# variables etc. #
##################

my ($link,$coord,$camace,$stlace,$autoace,$help);
GetOptions (
    "link=s"  => \$link,
    "coord=s" => \$coord,
    "camace"  => \$camace,
    "stlace"  => \$stlace,
    "autoace" => \$autoace,
    "-h"      => \$help,
);
if (($help) || (!$link) || (!$coord)) {
    print "usage: lookup_coords.pl -link <link name> -coord <number> -camace/stlace/autoace\n";
    exit(0);
}

my $dbdir;
$dbdir   = '/wormsrv1/camace'   if ($camace);
$dbdir   = '/wormsrv2/stlace'   if ($stlace);
$dbdir   = '/wormsrv2/autoace/' if ($autoace);
my $tace = &tace; 

###################
# be really fancy #
###################

if ($autoace) {
    print STDERR "This leads to /wormsrv2/autoace which might be in the build procedure right now.\n";
    print STDERR "Do you want to proceed? (y/n)?\n";
    my $choice = <STDIN>;
    if ($choice =~ /Y|y/) {
    }
    elsif ($choice =~ /N|n/) {
        print "Good bye...\n";
        exit(0);    
    }
    else {
        print "Wrong choice, bye :o)\n";
        exit(0);
    }
}

#######################
# open ace connection #
#######################

my $db      = Ace->connect(-path=>$dbdir,-program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
my $linkobj = $db->fetch(Sequence => "$link");

########
# Main #
########

my @subs = $linkobj->Subsequence(); 
my (%final);
foreach my $sub (@subs) {
    my ($name,$start,$end) = $sub->row();
    next if ($name =~ /\./); # ignore genes
    
    #################################
    # if subsequence is link itself #
    #################################
    
    if ($name =~ /LINK/) {
        print "$name\n";
        my $newlinkobj = $db->fetch(Sequence => "$name");
        my @newsubs = $newlinkobj->Subsequence(); 
        foreach my $newsub (@newsubs) {
            my ($newname,$newstart,$newend) = $newsub->row();
            next if ($newname =~ /\./); # ignore genes
            if (($coord >= $newstart) && ($coord <= $newend)) {
                $final{$newname} = $coord - $start - $newstart +2;
            }
            elsif (($coord <= $newstart) && ($coord >= $newend)) {
                $final{$newname} = $coord - $start - $newend +2;
            }
        }
    }
    
    ###############################
    # if subsequences are cosmids #
    ###############################
    
    else {
        if (($coord >= $start) && ($coord <= $end)) {
            $final{$name} = $coord - $start +1;
        }
        elsif (($coord <= $start) && ($coord >= $end)) {
            $final{$name} = $coord - $end +1;
        }
    }
}

############
# show off #
############

foreach my $name (sort keys %final) {
    print "$name\t$final{$name}\n";
}




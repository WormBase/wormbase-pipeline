#!/usr/local/bin/perl -w
#
# copies overlapcheck output files to /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks
# calculates the reported mistakes
# writes new html files and replaces the old ones   
#
# by Kerstin Jekosch
# 01/07/20
###################################################################################################

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
$|=1;
use Symbol 'gensym';

#####################
# variables and log #
#####################

my @chrom = qw ( I II III IV V X ); 
my $cvsversion = &get_cvs_version('/wormsrv2/copy2web.pl');
my $WSversion  = &get_wormbase_version;
my $maintainers = "dl1\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

my $logfile = "/wormsrv2/logs/copy2web.$WSversion.${rundate}.$$";
open (LOG,">$logfile") || die "Cannot open logfile $!\n";
LOG->autoflush();

print LOG "# copy2web\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "WormBase version : $WSversion\n";
print LOG "cvs version      : $cvsversion\n";
print LOG "\n";

my @output;
my @filenames = qw( overlapping_genes_cam overlapping_genes_stl EST_in_intron_cam EST_in_intron_stl repeat_in_exon_cam repeat_in_exon_stl );
 
############################################
# copy overlapcheck files to www directory #
############################################

print LOG "copying files from /wormsrv2/autoace/CHECKS/ to /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks\n"; 

system("mkdir /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion") && die "Cannot create web $WSversion dir $!\n";
system("mkdir /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks") && die "Cannot create web $WSversion/Checks dir $!\n";

foreach my $chrom (@chrom) {
    foreach my $file (@filenames) {
        print LOG "Copying file CHROMOSOME_$chrom.$file to /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks\n";
        system ("cp -f /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/")
            && die "Could not copy /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file $!\n";
    }
}

##################
# get linecounts #
##################

print LOG "get number of reports\n";

for (my $n = 0; $n < @filenames; $n++) {
    foreach my $chrom (@chrom) {
        print LOG "Calculating line nummbers for CHROMOSOME_$chrom.$filenames[$n]\n";
        my $line = `wc -l /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$filenames[$n]`;
        my ($new) = ($line =~ /(\d+)/);
        push @{$output[$n]}, $new;
    }
}    
 
####################### 
# make new html files #
####################### 

print LOG "Generating new html files in /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/\n";

my ($newnumber) = ($WSversion =~ /(\d+$)/);
my $oldnumber = $newnumber -1;
my $oldWSversion = "WS".$oldnumber;

for (my $m = 0; $m < @filenames; $m++) {    
    my $fh     = gensym();
    my $newfh  = gensym();
    my $file   = $filenames[$m];
    my $count = 0;
    
    open ($fh, "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$oldWSversion/Checks/$file.html") 
        || die "Cannot open old html file $!\n";
    open ($newfh, ">/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/$file.html")
        || die "Cannot open new html file $!\n";
    
    print LOG "Generating /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/$file.html\n";
    
    while (<$fh>) {
        if ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count < 6)) {
            my $old = $1;
            print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$output[$m][$count]." [$old]</B></TD>\n";
            $count++;
        }
        elsif ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count > 5)) {
            my $old = $1;
            my $sum;
            for (my $l = 0; $l < @chrom; $l++) {
                $sum += $output[$m][$l]; 
            }
            $count++;
            print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$sum." [$old]</B></TD>\n";
        }
        else {
            print $newfh $_;
        }
    }    
    close $fh    || die "Cannot close old html file $!\n";
    close $newfh || die "Cannot close new html file $!\n";
}

#################
# mail log file #
#################

close LOG;

open (mailLOG, "|/usr/bin/mailx -s \"copy2web Report: copy2web\" $maintainers ") || die "Cannot open mailLOG $!\n";
open (readLOG, "<$logfile") || die "Cannot open readLOG $!\n";
while (<readLOG>) {
    print mailLOG $_;
}
close readLOG || die "Cannot close readLOG $!\n";
close mailLOG || die "Cannot close mailLOG $!\n";

exit(0);

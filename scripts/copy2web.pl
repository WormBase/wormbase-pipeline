#!/usr/local/bin/perl -w

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
$|=1;
use Symbol 'gensym';

######################
# variables and logs #
######################

my @chrom = qw ( I II III IV V X ); 
my $cvsversion = &get_cvs_version('copy2web');
my $WSversion  = 'WS46'; #&get_wormbase_version;
my $maintainers = "kj2\@sanger.ac.uk";
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

my (@camgen,@stlgen,@camest,@stlest,@camrep,@stlrep,@output,@linecounts);
my @filenames   = qw( overlapping_genes_cam overlapping_genes_stl EST_in_intron_cam EST_in_intron_stl repeat_in_exon_cam repeat_in_exon_stl );
 
############################################
# copy overlapcheck files to www directory #
############################################

print LOG "copying files from /wormsrv2/autoace/CHECKS/ to /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks\n"; 

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

for (my $n = 0; $n < 6; $n++) {
    foreach my $chrom (@chrom) {
        print LOG "Calculating line nummbers for CHROMOSOME_$chrom.$filenames[$n]\n";
        my $line = `wc -l /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$filenames[$n]`;
        my ($new) = ($line =~ /(\d+)/);
        push @{$output[$n]}, $new;
    }
    print "see: $output[$n][0], $output[$n][1], $output[$n][2], $output[$n][3], $output[$n][4], $output[$n][5]\n"
}    
 
####################### 
# make new html files #
####################### 

print LOG "Generating new html files in /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/\n";

for (my $m = 0; $m < 6; $m++) {    
    my %handles;
    my $fh     = gensym();
    my $newfh  = gensym();
    my $file   = $filenames[$m];
    my $count = 0;
    
    open ($fh, "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/$file.html") 
        || die "Cannot open old html file $!\n";
    open ($newfh, ">/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/$file.html.new")
        || die "Cannot open new html file $!\n";
    
    print LOG "Generating /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/$WSversion/Checks/$file.html.new\n";
    
    while (<$fh>) {
        if ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count < 6)) {
            my $old = $1;
            print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$output[$m][$count]." [$old]</B></TD>\n";
            print "was $old is now $output[$m][$count]\n";
            $count++;
        }
        elsif ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count > 5)) {
            my $old = $1;
            my $sum;
            for (my $l = 0; $l < 6; $l++) {
                $sum += $output[$m][$l]; 
            }
            $count++;
            print "sum was $old is now $sum\n";
            print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$sum." [$old]</B></TD>\n";
        }
        else {
            print $newfh $_;
        }
    }    
    close $fh    || die "Cannot close old html file $!\n";
    close $newfh || die "Cannot close new html file $!\n";
}

#####################
# replace old files #
#####################

print LOG "Replacing old files\n"; 

foreach my $file (@filenames) {
    print LOG "Deleting /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/Checks/$file.html\n";
    unlink ("/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/Checks/$file.html") 
        || die "Cannot delete old file $file $!\n";
    print LOG "Renaming /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/Checks/$file.html.new\n\n";    
    rename ("/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/Checks/$file.html.new",
        "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/Checks/$file.html") 
        || die "Cannot rename old file $file $!\n";
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


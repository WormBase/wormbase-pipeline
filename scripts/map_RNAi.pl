#!/usr/local/bin/perl -w
#
# map_RNAi
# v 0.1
#
# Cronjob integrity check controls for generic ACEDB database.
#
# Usage: map_RNAi [-options]
#
# 010927 :  dl : modified output print line to include speech marks. this prevents any acefile
#                parsing problems
# 060801 :  kj : original version 
#
#############################################################################################


#############
# variables #
#############

$|=1;
BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}
use Bio::Seq;
use IO::Handle;
use Getopt::Std;
use Cwd;
use Ace;
use lib "/wormsrv2/scripts/";
use Wormbase;

##########################
# Script variables (run) #
##########################

my $maintainer = "dl1\@sanger.ac.uk kj2\@sanger.ac.uk krb\@sanger.ac.uk";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;
my %output      = ();
my %finaloutput = ();

########################
# command-line options #
########################

$opt_d="";   # Verbose debug mode
$opt_h="";   # Help/Usage page

getopts ('hd');
&usage(0) if ($opt_h);
my $debug = $opt_d;

#############
# Paths etc #
#############

my $tace = glob("~acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace");    # tace executable path
my $dbdir  = "/wormsrv2/autoace";                                  # Database path
my $gffdir = "/wormsrv2/autoace/CHROMOSOMES";
my @chromosomes = qw( I II III IV V X );
my $version = &get_cvs_version($0);
my $db_version = &get_wormbase_version_name;
 
################
# Open logfile #
################

my $log="/wormsrv2/logs/map_RNAi.$rundate";

open (LOG,">$log");
LOG->autoflush();

print LOG "# map_RNAi\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "# version        : $version\n";
print LOG "\n";

############################################
# get exons and RNAis out of the gff files #
############################################        
   
foreach $chromosome (@chromosomes) {

    my %exoncount = ();
    my %RNAicount = ();
    my %genes     = ();
    my %RNAi      = ();
    my %exon      = ();
        
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
        chomp;
        s/^\#.*//;
        next unless /\S/;
        @f = split /\t/;

	if ($f[1] eq "RNAi") {
	    my ($name) = ($f[8] =~ /\"(.*)\"$/);
            $RNAicount{$name}++;
            my $RNAiname = $name.".".$RNAicount{$name};
            $RNAi{$RNAiname} = [$f[3],$f[4]];
	}
	elsif (($f[2] eq "exon") && (($f[1] eq "curated") || ($f[1] eq "Pseudogene") || ($f[1] eq "provisional"))) {
            my ($name) = ($f[8] =~ /\"(\S+)\"/);
            $exoncount{$name}++;
            my $exonname = $name.".".$exoncount{$name};
            $exon{$exonname} = [$f[3],$f[4]];
        }
    }
    close GFF_in;
    
#########################   
# make exons into genes #
#########################
 
    foreach my $name (sort keys %exoncount) {
        my $v = $exoncount{$name};
        my $w = $name.".".$v;
        $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
    }
    
###################
# make indexlists #
###################
 
    my @exonlist= sort { $exon{$a}->[0]  <=> $exon{$b}->[0]  } keys %exon;
    my @genelist= sort { $genes{$a}->[0] <=> $genes{$b}->[0] } keys %genes;
    my @RNAilist = sort { $RNAi{$a}->[0]   <=> $RNAi{$b}->[0] || $a cmp $b } keys %RNAi;     

##########
# map it #
##########
 
    my $lastfail = 0;
    
    for (my $x = 0; $x < @RNAilist; $x++) {
        my $testRNAi   = $RNAilist[$x];
        my $RNAistart  = $RNAi{$testRNAi}->[0];
        my $RNAistop   = $RNAi{$testRNAi}->[1];
            
        for (my $y = 0; $y < @genelist; $y++) {
            my $testgene = $genelist[$y];
            my $genestart= $genes{$testgene}->[0];
            my $genestop = $genes{$testgene}->[1];
            
            if ($RNAistart > $genestop) {
                $lastfail = $y;
                next; 
            }
            
            elsif ($RNAistop < $genestart) {
                last; 
            }
            
            else {
                for (my $z = 1; $z <= $exoncount{$testgene}; $z++) {
                    my $exon_start = $exon{"$testgene.$z"}->[0];
                    my $exon_stop  = $exon{"$testgene.$z"}->[1];
                   
                    if ( not (($RNAistart > $exon_stop) || ($RNAistop < $exon_start))) {
                        my ($RNAi) = ($testRNAi =~ /(\S+)\.\d+$/);
                        push @{$output{$RNAi}}, $testgene;
                        print LOG "$RNAi mapped to $testgene\n";
                    }
                }
            }                
        }
    }
}

close LOG;

###################
# sort the output #
###################
 
foreach my $mess (sort keys %output) {
    @{$output{$mess}} = sort @{$output{$mess}};
    my $count = 0; 
    push @{$finaloutput{$mess}}, $output{$mess}->[0];
    for (my $m = 1; $m < (scalar @{$output{$mess}}); $m++) {
        if ($output{$mess}->[$count] ne $output{$mess}->[$m]) {
            push @{$finaloutput{$mess}}, $output{$mess}->[$m];
            $count = $m;    
        }    
    }
}

########################
# produce output files #
########################

open (OUT,    ">/wormsrv2/autoace/mappings/RNAi_mappings.$db_version");
open (OUTACE, ">/wormsrv2/autoace/mappings/RNAi_mappings.$db_version.ace");

foreach my $mapped (sort keys %finaloutput) {
    print OUT "$mapped\t@{$finaloutput{$mapped}}\n";
    for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
        print OUTACE "RNAi : $mapped\n";
        print OUTACE "-D Predicted_gene\n\n";
    }
} 

foreach my $mapped (sort keys %finaloutput) {
    print OUT "$mapped\t@{$finaloutput{$mapped}}\n";
    for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
        print OUTACE "RNAi : \"$mapped\"\n";
        print OUTACE "Predicted_gene \"$finaloutput{$mapped}->[$n]\"\n\n";
    }
} 


##############################
# read acefiles into autoace #
##############################

my $command =<<END;
pparse /wormsrv2/autoace/mappings/RNAi_mappings.$db_version.ace
save
quit
END

open (TACE,"| $tace -tsuser map_RNAi.pl $dbdir") || die "Couldn't open tace connection to $dbdir\n";
print TACE $command;
close (TACE);


exit(0);

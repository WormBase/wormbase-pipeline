#!/usr/local/bin/perl
#
# Check_DNA.pl v1.0
# dl
# 2000-04-26
#

########################################
# iniatialise                          #
########################################

$|=1;
#use IO::Handle;
#use Getopt::Std;
use lib '/wormsrv2/scripts';
use Wormbase;

my $gffdir    = "/wormsrv2/autoace/CHROMOSOMES";
my $agpdir    = "/wormsrv2/autoace/yellow_brick_road";
my $scriptdir = "/wormsrv2/scripts";

# prepare array of file names and sort names
@files = (
	  'CHROMOSOME_I.gff',
	  'CHROMOSOME_II.gff',
	  'CHROMOSOME_III.gff',
	  'CHROMOSOME_IV.gff',
	  'CHROMOSOME_V.gff',
	  'CHROMOSOME_X.gff',
	  );

@gff_files = sort @files; 
undef @files; 

############################################################
# loop through each GFF file                               #
############################################################

$debug = 1;

foreach $file (@gff_files) {
    
    next if ($file eq "");
    my ($chromosome) = $file =~ (/CHROMOSOME\_(\S+)\./);
    
    open (OUT, ">$agpdir/CHROMOSOME_$chromosome.clone_path.gff") || die "can't open output file '$agpdir/CHROMOSOME_$chromosome.clone_path.gff'\n";
    open (GFF, "<$gffdir/$file") || die "can't open gff file '$gffdir/$file'\n";
    while (<GFF>) {
	chomp;
	next if ($_ =~ /^\#/);
	($name,$method,$feature,$start,$stop,$score,$strand,$other,$name) = split (/\t/,$_);
	
	if (($method eq "Genomic_canonical") && ($feature eq "Sequence")) {
	    print OUT "$_\n";
	}
    }
    close GFF;
    close OUT;
    
    # modify to make the clone_acc lists
    system ("$scriptdir/GFF_with_accessions $agpdir/CHROMOSOME_$chromosome.clone_path.gff > $agpdir/CHROMOSOME_$chromosome.clone_acc.gff");

}




exit;


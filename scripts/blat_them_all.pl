#!/usr/bin/perl5.6.0
#
# wrapper to generate blat data by
# - getting sequence out of autoace
# - blating it against all ESTs/mRNAs
# - mapping it back to autoace
# - producing confirmed introns
# - producing virtual objects to put the data into
#
# 16.10.01 Kerstin Jekosch
# 17.10.01 kj: modified to get everything onto wormsrv2 and to include an mRNA and parasitic nemattode blatx option

use strict;
use Ace;
use Wormbase;
use Getopt::Std;
use IO::Handle;
use vars qw($opt_e $opt_m $opt_x $opt_b $opt_s);
$|=1;

$opt_e = ""; # run everything for ESTs
$opt_m = ""; # run everything for mRNAs
$opt_x = ""; # run everything for parasitic nematode ESTs
$opt_b = ""; # start with blating (autoace.fa, chromosome.ace already present)
$opt_s = ""; # start later with sorting/mapping (blat2ace.pl)
getopts('emxbs');
print "Choose between -e (ESTs), -m (mRNAs) or -x (blatx of nematode ESTs)\n" unless ($opt_m || $opt_e || $opt_x);

###############
# directories #
###############

my $dbdir  = '/wormsrv2/autoace';
my $bin    = '/wormsrv2/scripts';
my $query;
$query     = '/nfs/disk100/wormpub/analysis/ESTs/C.elegans_nematode_ESTs' if ($opt_e);
$query     = '/wormsrv2/mRNA_GENOME/NDB_mRNA.fasta'                       if ($opt_m);
$query     = '/wormsrv2/autoace/BLAT/nematode.fa'       if ($opt_x); 
my $seq    = '/wormsrv2/autoace/BLAT/autoace.fa';
my $chrom  = '/wormsrv2/autoace/BLAT/chromosome.ace';
my $blatex = '/nfs/disk100/wormpub/blat/blat';
my $blat   = '/wormsrv2/autoace/BLAT';
my $giface = '/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/giface';

#############
# LOG stuff #
#############

my $maintainer = "kj2\@sanger.ac.uk"; #dl1\@sanger.ac.uk krb\@sanger.ac.uk";
my $rundate    = `date +%y%m%d`;   chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;
my $version    = &get_cvs_version('/wormsrv2/scripts/blat_them_all.pl');
my $WS_version = &get_wormbase_version_name;

#my $logfile = "/wormsrv2/logs/blat_them_all.${WS_version}.$rundate.$$";
#system ("/bin/touch $logfile");
#open (LOG,">>$logfile") or die ("Could not create logfile\n");
#LOG->autoflush();
#open (STDOUT,">>$logfile");
#STDOUT->autoflush();
#open (STDERR,">>$logfile"); 
#STDERR->autoflush();

#print LOG "# blat_them_all\n\n";     
#print LOG "# version        : $version\n";
#print LOG "# run details    : $rundate $runtime\n";
#print LOG "\n";
#print LOG "WormBase version : ${WS_version}\n";
#print LOG "\n";
#print LOG "======================================================================\n";
#print LOG " -e : perform blat for ESTs\n"                     if ($opt_e);
#print LOG " -m : perform blat for mRNAs\n"                    if ($opt_m);
#print LOG " -x : perform blatx for parasitic nematode ESTs\n" if ($opt_x);
#print LOG "======================================================================\n";
#print LOG "\n";
#print LOG "Starting blat process .. \n\n";

###########################
# get data out of autoace #
###########################

unless (!$opt_b || $opt_s) {
	print "Getting sequence data out of autoace and putting it into $seq\n";
	print "Getting chromosome data out of autoace and putting it into $chrom\n";

	my $command1=<<EOF;
	find sequence "CH*"
	follow subsequence
	dna -f /wormsrv2/autoace/BLAT/autoace.fa
	clear
	find sequence "CH*"
	show -a -f /wormsrv2/autoace/BLAT/chromosome.ace
	quit
EOF

	system("echo '$command1' | $giface $dbdir");
}

############
# run blat #
############

unless ($opt_s) {
	if ($opt_e) {
		print "running blat and putting the result in $blat/est_out.psl\n";
		system("$blatex $seq $query $blat/est_out.psl") && die "Blat failed\n";
	}
	if ($opt_m) {
		print "running blat and putting the result in $blat/mRNA_out.psl\n";
		system("$blatex $seq $query $blat/mRNA_out.psl") && die "Blat failed\n";
	}
	if ($opt_x) {
		print "Getting nematode consensus sequences and putting them into /wormsrv2/autoace/BLAT/nematode.fa\n";
		my $command2=<<EOF;
		query find nematode_ESTs where Remark = "*Blaxter*"
		dna -f /wormsrv2/autoace/BLAT/nematode.fa
		quit
EOF
		system("echo '$command2' | $giface $dbdir ");
		print "running blat and putting the result in $blat/nematode_out.psl\n";
		system("$blatex $seq $query -t=dnax -q=dnax $blat/nematode_out.psl") && die "Blat failed\n";
	}
}

##################
# map to autoace #
##################

if ($opt_e) {
	print "Mapping blat data to autoace\n";
	print "Putting results to $blat/autoace.blat.EST.ace and $blat/autoace.ci.EST.ace\n";
	system("$bin/blat2ace.pl -a -i") && die "Mapping failed\n"; 
}
if ($opt_m) {
	print "Mapping blat data to autoace\n";
	print "Putting results to $blat/autoace.blat.mRNA.ace and $blat/autoace.ci.mRNA.ace\n";
	system("$bin/blat2ace.pl -a -i -m") && die "Mapping failed\n"; 
}
if ($opt_x) {
	print "Mapping blat data to autoace\n";
	print "Putting results to $blat/autoace.blat.nematode.ace\n";
	system("$bin/blat2ace.pl -a -x") && die "Mapping failed\n"; 
}

#############################
# produce confirmed introns #
#############################

if ($opt_e) {
	print "Producing confirmed introns in $blat/autoace.good_introns.EST.ace\n";
	system("$bin/confirm_introns.pl -a -e") && die "Intron confirmation failed\n";
}
if ($opt_m) {
	print "Producing confirmed introns in $blat/autoace.good_introns.mRNA.ace\n";
	system("$bin/confirm_introns.pl -a -m") && die "Intron confirmation failed\n";
}

#########################################
# produce files for the virtual objects #
#########################################

if ($opt_e) {
	print "Producing files for the virtual objects: $blat/virtual_objects.BLAT_EST.ace\n";
	system("$bin/superlinks.blat.pl $chrom > ! $blat/virtual_objects.BLAT_EST.ace")
		&& die "Producing virtual objects for blat failed\n";

	print "Producing files for the virtual objects: $blat/rawdata/virtual_objects.ci.EST.ace\n";
	system("$bin/superlinks.confirmed_introns.pl $chrom > ! $blat/virtual_objects.ci.EST.ace")
		&& die "Producing virtual objects for confirmed introns failed\n"; 
}
if ($opt_m) {
	print "Producing files for the virtual objects: $blat/virtual_objects.BLAT_mRNA.ace\n";
	system("$bin/superlinks.blat.pl -m $chrom > ! $blat/virtual_objects.BLAT_mRNA.ace")
		&& die "Producing virtual objects for blat failed\n";

	print "Producing files for the virtual objects: $blat/virtual_objects.ci.mRNA.ace\n";
	system("$bin/superlinks.confirmed_introns.pl -m $chrom > ! $blat/virtual_objects.ci.mRNA.ace")
		&& die "Producing virtual objects for confirmed introns failed\n"; 
}
if ($opt_x) {
	print "Producing files for the virtual objects: $blat/virtual_objects.BLATX_NEMATODE.ace\n";
	system("$bin/superlinks.blat.pl -x $chrom > ! $blat/virtual_objects.BLATX_NEMATODE.ace")
		&& die "Producing virtual objects for blat failed\n";
}

##########################
# give 'read in' message #
##########################

print "\n";
print "Read into autoace:\n";
if ($opt_e) {
	print "$blat/virtual_objects.BLAT_EST.ace\n";
	print "$blat/autoace.blat.EST.ace\n";
	print "$blat/virtual_objects.ci.EST.ace\n";
	print "$blat/autoace.good_introns.EST.ace\n";
}
if ($opt_m) {
	print "$blat/virtual_objects.BLAT_mRNA.ace\n";
	print "$blat/autoace.blat.mRNA.ace\n";
	print "$blat/virtual_objects.ci.mRNA.ace\n";
	print "$blat/autoace.good_introns.mRNA.ace\n";
}
if ($opt_x) {
	print "$blat/virtual_objects.BLATX_NEMATODE.ace\n";
	print "$blat/autoace.blat.nematode.ace\n";
}
print "Good bye :o)\n";


#&mail_maintainer("WormBase Report: blat_them_all ",$maintainer,$logfile);




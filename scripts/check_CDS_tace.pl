#!/usr/local/bin/perl
#
# check_CDS.pl v1.1
# dl
# 2000-04-26
#
# Aceperl script to check database gene models. 
#
###########################################################################################
#
# 000426 dl : PP version
# 000503 dl : Utilise tace rather than AcePerl for speed reasons
#           : screening camace can be done in 7 minutes rather than 1.5 hours
#
###########################################################################################

BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}
use Bio::Seq;
use Getopt::Std;
use Cwd;
$|=1;

##################################################
# command-line parsing                           #
##################################################

getopts ('csdwa:');

if ($opt_c) {$ace="1"; print "Using camace data ..\n";}
if ($opt_s) {$ace="2"; print "Using St Louis data ..\n";}
if ($opt_a) {$ace="3"; print "Using User-given  data ..\n";}
if ($opt_d) {$debug="1"; print "Using Debug/Verbose mode ..\n";}
if ($opt_w) {$html="1"; print "Using HTML mode ..\n";}

if ((!$opt_c)&&(!$opt_s)&&(!$opt_a)) {
    &usage;
    exit;
}

##################################################
# Connect with acedb database                    #
##################################################

my $stlacepath="/wormsrv2/stlace";
my $camacepath="/wormsrv2/camace";

if ($opt_a =~ /^(\~\w+)\//){
    $autoacepath=glob("$1");
    $autoacepath =~ s/\/tmp_mnt//;
    $FILENAME=$';
    $autoacepath="$autoacepath"."/"."$FILENAME";
 } elsif ($opt_a =~ /^(\w+)/) {
    $autoacepath="$CWD"."/"."$opt_p";
 } elsif ($opt_a =~ /\/\w+/) {
    $autoacepath=$opt_a;
 } else {
    &usage; 
}

if ($ace == 1) {
    $ENV{'ACEDB'}="$camacepath";
}
elsif ($ace == 2) {
    $ENV{'ACEDB'}="$stlacepath";
}
elsif ($ace == 3) {
    $ENV{'ACEDB'}="$autoacepath";
}

##################################################
# Arrays for storing CDS details                 #
##################################################

%CDS_len = "";
%CDS_chksum = "";
%source = "";
%protein = "";

##################################################
# Main Loop                                      #
##################################################

$exec=&tace;
$command=<<EOF;
query find Sequence where Method = \"curated\" 
show -a 
DNA
quit
EOF

open (textace, "echo '$command' | $exec  | ");
while (<textace>) {
    chomp;
#    print "$_\n";
   if (/Sequence : \"(\S+)\"/)  {
	$gene = $1;
	push (@genes, $gene);
    }
    if (/Source\s+\"(\S+)\"/)  {
	$source{$gene} = $1;
    }
    if (/Corresponding_protein\s+\"(\S+)\"/)  {
	$protein{$gene} = $1;
    }
    if (/acedb> >(\S+)/) {
	$CDS_pointer = $1;
	$dna_counter = 1;
	$CDS_seq = ">$1\n";
	next;
    }
    if (($dna_counter == 1) && ((/^>/) || (/object dumped/))) {
	$CDS_seq=~tr/a-z/A-Z/;
	$CDS_seq=~s/\>\w+//;
	$CDS_seq=~s/\W+//mg;
	if ($CDS_seq =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	    $CDS_seq=~s/[^ACGTUMRWSYKVHDBXN]//img;
	}
	$CDS_len{$CDS_pointer} = length($CDS_seq);
	$bioseq = Bio::Seq->new(-seq=>$CDS_seq,-id=>$CDS_pointer,-ffmt=>'Fasta',-type=>'Dna',);
	$CDS_chksum{$CDS_pointer}=$bioseq->GCG_checksum;
	undef $bioseq;

	(/>(\S+)/);
	$CDS_seq = ">$1\n";
	$CDS_pointer = $1;
	next;
    }
    if ($dna_counter == 1) {
	$CDS_seq .= $_;
	next;
    }
}
close (textace);

################################################## 
# output the tab-delimited flatfile              #
##################################################

foreach $gene (@genes) {
    if ($protein{$gene} eq "") {
	$protein{$gene} = "n/a     ";
    }
    print "$gene    \t$CDS_len{$gene}\t$CDS_chksum{$gene}\t$protein{$gene}\t$source{$gene}\n";
}

exit;

####################
# subs
####################

sub tace {
   local($prog);
   local($name); 
   $name=`uname -sr`;
    if    ($name=~/^SunOS/) {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
    elsif ($name=~/^IRIX/)  {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)   {($prog)=<~acedb/RELEASE.SUPPORTED/bin.ALPHA_4/giface>;}
    elsif ($name=~/^Linux/) {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}
    else {print STDERR "No known binary for $uname\n";exit;}
    return $prog;
}



##################################################
# usage subroutine                               #
##################################################

sub usage {
    print "\n\nUsage: $0 [-options] > outfile\n\n";
    print "Options:\n";
    print "-c   => Use Cambridge data  \n";
    print "-s   => Use St Louis data   \n";
    print "-a Pathname  => Use data in path Pathname  \n";
    print "-d   => Debug/Verbose mode\n\n";
}




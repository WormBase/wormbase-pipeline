#!/usr/local/bin/perl5.6.1 -w

# update_caltech.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-02-12 15:38:04 $
# Last updated by: $Author: ck1 $


# Automatically update Geneace with Erich's functional annotation update

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;

###################################################
# variables and command-line options with aliases # 
###################################################

my ($debug, $recipients);
$recipients ="ck1\@sanger.ac.uk, emsch\@its.caltech.edu, krb\@sanger.ac.uk";
my $tace = &tace;   # tace executable path

GetOptions ("d|debug=s"  => \$debug);

##################################################
# check user is wormpub otherwise script won't run
##################################################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\nYou have to be wormpub to run this script!\n\n";
  exit(0);
}

# if ($debug){$recipients ="$debug\@sanger.ac.uk"}

###############################
# touch logfile for run details
###############################

$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");
my $script = $0;

my $rundate = `date +%y%m%d`; chomp $rundate;

my $log = "/wormsrv2/logs/update_caltech.$rundate.$$";

#########################################
# check for new file to upate on FTP site
#########################################

my $caltech="/wormsrv1/geneace/ERICHS_DATA/FTP";
my $updatedir="/wormsrv1/geneace/ERICHS_DATA";
my @dates;

open(LOG, ">$log") || die "Can't write to log file $log\n";
print LOG "This file is generated automatically. If you have spotted any bug, please contact ck1\@sanger.ac.uk\n";
print LOG "--------------------------------------------------------------------------------------------------\n\n";

chdir $caltech;
system("echo '\$ caltech' | ftp -i caltech.wormbase.org") && print LOG "Failed to download file\n";

my @ftp_date=dataset($caltech);    # [0] is the date of the file , [1] is the filename
my @last_date=dataset($updatedir); # [0] is the date of the file , [1] is the filename

push(@dates, $ftp_date[0], $last_date[0]);

my @order = sort {$a <=> $b} @dates;

if ($last_date[0] == $order[0]){
  print LOG "UPDATE file $ftp_date[1] avilable on FTP...start to download\n\n";
  chdir $updatedir;
  system("mv $last_date[1] ARCHIVE/"); 
  system("echo '\$ caltech' | ftp -i caltech.wormbase.org") && print LOG "Failed to download file\n"; 
}
else {
  print LOG "No new update on FTP site\n";
  $recipients ="ck1\@sanger.ac.uk, krb\@sanger.ac.uk";
  mail_maintainer($script, $recipients, $log) unless $debug;
  exit(0);
}

sub dataset {
  my $dir = shift;
  my ($date, $name);
  opendir(DIR, $dir) || die "Can't read directory";
  my @dir=readdir DIR;
  splice(@dir, 0,2);
  closedir (DIR);
  foreach (@dir){
    if ($_ =~ /^annots-(\d+)(\w{3,3})(\d+)/){
      $name = $_;
      my $mon = $2;
      if ($mon eq "jan"){$mon = "01"}
      if ($mon eq "feb"){$mon = "02"}
      if ($mon eq "mar"){$mon = "03"}
      if ($mon eq "apr"){$mon = "04"}
      if ($mon eq "may"){$mon = "05"}
      if ($mon eq "jun"){$mon = "06"} 
      if ($mon eq "jul"){$mon = "07"}
      if ($mon eq "aug"){$mon = "08"}
      if ($mon eq "sep"){$mon = "09"}
      if ($mon eq "oct"){$mon = "10"}
      if ($mon eq "nov"){$mon = "11"}
      if ($mon eq "dec"){$mon = "12"}
      $date = $1.$mon.$3;
    }
  }
  return $date, $name; 
}

###############################################
# check if locus is now CGC_approved other_name
# if yes, replace by main_name in ace file
# if not, output to an file and check by hand
###############################################

my $ga_dir = "/wormsrv1/geneace";

my $locus_has_other_name=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name.def" quit
EOF

my @other_main=loci_as_other_name($locus_has_other_name, $ga_dir);

my %other_main=%{$other_main[0]};

my %main;
my @main=@{$other_main[1]}; foreach (@main){$main{$_}++};

my $change=0;
my ($non_cgc_other_name, $cgc_other_name);

$!=();

open (IN, "$updatedir/$ftp_date[1]");
print LOG "Can't read in update file $ftp_date[1] downloaded from FTP site\n" if ($!);

my $modify ="/wormsrv1/geneace/ERICHS_DATA/$ftp_date[1].modified";

open(OUT, ">$modify");
print LOG "Couldn't write modification to $modify" if ($!);

while(<IN>){
  if ($_ =~ /^Locus : \"(.+)\"/){

    ######################################################### 
    # check if locus is an other_name and also a locus object
    #########################################################

    if ($other_main{$1} && $main{$1}){    
      $non_cgc_other_name .= "$1 - this is a valid CGC gene name, but has also been used as an other name for @{$other_main{$1}}\n";
    }  
    if ($other_main{$1} && !$main{$1}){     # check if locus is an other_name and also not a locus object  
      $cgc_other_name .= "\n$1 is now an other_name of @{$other_main{$1}}";
      print OUT "Locus : \"@{$other_main{$1}}\"\n";
      $change++;
    }
    else {
      print OUT "$_";
    }
  }
  else {
    print OUT "$_";
  }
}
print LOG "Functional annotations are attached to the following genes which\n";
print LOG "*possibly* should be attached to different genes:\n\n";
print LOG "$non_cgc_other_name\n\n";

print LOG "Functional annotations are attached to the following gene names\n";
print LOG "which are no longer valid CGC names, and should be attached to new names as follows:\n";
print LOG "$cgc_other_name\n\n";
print LOG "There are $change change(s)\n";
print LOG "These have been changed in $ftp_date[1],\n";
print LOG "but Caltech needs to update on their side. THANKS.\n\n";

system("rm $caltech/*");
#=start
my $command=<<END;
find sequence * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/seq_TS_dump.ace
edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

find locus * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/loci_TS_dump.ace

edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

pparse $modify

save
quit
END

my $geneace_dir="/wormsrv1/geneace";

open (Load_GA,"| $tace -tsuser \"Functional_annotation\" $geneace_dir ") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;

if($!){
  print "##########################################\n";
  print "\nPhenotype annotation is now updated!\n\n";       
  print "\nIf everthing is OK, REMEMBER to remove\n"; 
  print "loci_TS_dump.ace and seq_TS_dump.ace\n";
  print "in /wormsrv1/geneace/ERICHS_DATA\n\n";
  print "##########################################\n\n"; 
}
else{
  print "#################################\n";
  print "Mission not 100% successful, Mr. Bond!\n";
  print "#################################\n\n";
}

system("mv $last_date[1].modified ARCHIVE/");

##################################################
# Mail log file
# Only mail to person running script in debug mode
##################################################

mail_maintainer($script, $recipients, $log) unless $debug;

#############
# subroutines
#############

sub loci_as_other_name {

  my ($def, $dir) = @_;
  my ($main, @main, $other_name, %other_main);
  
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//g; #print $other_name, "\n";
    #  print "$main -> $other_name\n";
      push(@{$other_main{$other_name}}, $main);
      push(@main, $main);
    }
  }
  return \%other_main, \@main
}



__END__




#!/usr/local/bin/perl5.8.0 -w

# update_caltech.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-06-17 16:17:22 $
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

my ($debug, $help, $recipients);
$recipients ="ck1\@sanger.ac.uk, emsch\@its.caltech.edu, krb\@sanger.ac.uk";
my $tace = &tace;   # tace executable path

GetOptions ("d|debug=s"  => \$debug,
            "h|help"  => \$help);

if ($help){
  system("perldoc update_caltech.pl");
  exit (0);
}

if ($debug){
  print "DEBUG : $debug\n";
  $recipients = "$debug\@sanger.ac.uk";
}
	

##################################################
# check user is wormpub otherwise script won't run
##################################################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\nYou have to be wormpub to run this script!\n\n";
  exit(0);
}

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

print $ftp_date[0], "#\n";
print $ftp_date[1], "#\n";
print $last_date[0], "##\n";
print $last_date[1], "##\n";

if ($last_date[0] != $ftp_date[0]){
  print LOG "UPDATE file $ftp_date[1] avilable on FTP\n\n";
  system("mv $ftp_date[1] ../");
}  
if ($last_date[0] == $ftp_date[0]){
  print LOG "No new update on FTP site\n";
  $recipients ="ck1\@sanger.ac.uk, krb\@sanger.ac.uk";  # notify ck1 & krb if no update
  mail_maintainer($script, $recipients, $log);
  exit(0);
}

########################################
# process locus / other_name information
########################################

my $ga_dir = "/wormsrv1/geneace";

my $locus_has_other_name_to_cds=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name_to_cds.def" quit
EOF
my $locus_to_CDS=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS.def" quit
EOF

my @other_main=CGC_loci_and_other_name($locus_has_other_name_to_cds, $locus_to_CDS, $ga_dir);

my %other_main=%{$other_main[0]};
my %locus_cds=%{$other_main[1]};


my $count=0;
my (@exceptions, %exceptions, @other_names, %other_names);

foreach (sort keys %other_main){
  my $locus = @{$other_main{$_}}->[0];
  my @cds = @{$other_main{$_}}->[1];

  if (exists $locus_cds{$_}){
     push(@exceptions, $_);
  }
  else {
    push(@other_names, $_);
  }
}

foreach (@exceptions){$exceptions{$_}++}
foreach (@other_names){$other_names{$_}++}

open (IN, "$updatedir/$ftp_date[1]");

my $modify ="$updatedir/$ftp_date[1].modified";

open(OUT, ">$modify");

print LOG "CURRENT LIST of main name/other name INFORMATION:\n";
print LOG "-------------------------------------------------\n";

while(<IN>){
  if ($_ =~ /^Locus : \"(.+)\"/){
   
    ######################################################### 
    # check if locus is an other_name and also a locus object
    #########################################################

    if ($exceptions{$1}){
      print $1, "\n";
      #$count++;
      my $locus = @{$other_main{$1}}->[0];
      my @cds = @{$other_main{$1}}->[1]; 
      if (@{$other_main{$1}}->[1]){
	print LOG "$1 (@{$locus_cds{$1}}) is different from $1 which is an other_name of $locus (@cds)\n";
	print LOG "Functional annotation is assigned to $1 (@{$locus_cds{$1}})\n\n";
	print OUT "Locus : \"$1\"\n";
      }
      else {
	print LOG "$1 (@{$locus_cds{$1}}) is different from $1 which is an other_name of $locus (no CDS connection yet)\n";
	print LOG "Functional annotation is assigned to $1 (@{$locus_cds{$1}})\n\n";
	print OUT "\nLocus : \"$1\"\n";
      }
    }

    ##################################################################### 
    # check if locus should become an other_name of a CGC locus object
    #####################################################################

    elsif ($other_names{$1}){
      $count++;
      my $locus = @{$other_main{$1}}->[0];
      my @cds = @{$other_main{$1}}->[1];
      if (@{$other_main{$1}}->[1]){
	print LOG "$1 is an other_name of $locus (@cds)\n";
	print LOG "Functional annotation is assigned to $locus\n\n";
	print OUT "\nLocus : \"$locus\"\n";
	print OUT "Other_name \"$1\"\n";
      }
      else {
	print LOG "$1 is an other_name of $locus (no CDS connection yet)\n";
	print LOG "Functional annotation is assigned to $locus\n\n";
	print OUT "\nLocus : \"$locus\"\n";
	print OUT "Other_name \"$1\"\n";
      }
    }
    else {
      print OUT "$_";
    }
  }
  else {
    print OUT "$_";
  }
}

my $diff_out = `diff $updatedir/$ftp_date[1] $modify`;
if (!$diff_out){ 
  print LOG "------------------------------------------\n\n";
  print LOG "EVERYTHING WENT OK\n\n";
  print LOG "THANKS.\n\n";
}
else {
  print LOG "------------------------------------------\n\n";
  print LOG "CHANGES MADE\n\n";  
  print LOG $diff_out, "\n";
  print LOG "------------------------------------------\n\n";
  print LOG "Please update accordingly\n";
  print LOG "THANKS.\n\n";
}

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

if (!$debug){
  open (Load_GA,"| $tace -tsuser \"Functional_annotation\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
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
    print "######################################\n";
    print "Mission not 100% successful, Mr. Bond!\n";
    print "######################################\n\n";
  }
  
  ######################################
  # move modified file/last update file to ARCHIVE folder
  ######################################
  
  chdir $updatedir;
  system("mv $last_date[1]* ARCHIVE/");
}  

##################################################
# Mail log file
# Only mail to person running script in debug mode
##################################################
  
mail_maintainer($script, $recipients, $log);
#mail_maintainer($script, "$debug\@sanger.ac.uk", $modify);
exit(0);


#############
# subroutines
#############

sub dataset {
  my $dir = shift;
  my ($date, $name);
  opendir(DIR, $dir) || die "Can't read directory";
  my @dir=readdir DIR;
  splice(@dir, 0,2);
  closedir (DIR);
  foreach (@dir){
    if ($_ =~ /^annots-(\d+)(\w{3,3})(\d+)\.ace$/){
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

sub CGC_loci_and_other_name {

  my ($def1, $def2, $dir) = @_;
  my ($main, $other_name, %other_main, $locus, $cds, %locus_cds);
  
  open (FH1, "echo '$def1' | tace $dir |") || die "Couldn't access geneace\n";
  while (<FH1>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//;
      if ($3){
	$cds = $3;
	push(@{$other_main{$other_name}}, $main, $cds);
      }
      else {
	push(@{$other_main{$other_name}}, $main);
      }
    }
  }
 
  open (FH2, "echo '$def2' | tace $dir |") || die "Couldn't access geneace\n";
  while (<FH2>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $locus = $1;
      $locus =~ s/\\//;
      $cds = $2;
      push(@{$locus_cds{$locus}}, $cds);
    }
  }
  return \%other_main, \%locus_cds;
}


__END__

=head2 NAME - update_caltech.pl  


            This script automatically checks for available functional annotation update at Caltech FTP site.
            It emails to ck1 and krb no matter update is available or not.
            When update is available log messages is also included in the mail notice. 
            Erich at Caltech will also receive mail notice automatically.

            The script also checks for loci that are other_names in Caltech update file. 
            If an other_name is still a separate object, it will be listed in the mail message 
            and need to be looked at by hand.
            If a locus has been merged to a main name, that locus will be changed to its 
            main name in the update ace file. 

=head3 <USAGE> 
 

=head2 Options: [d or debug] [h or help]

 
B<-debug:>   
            Debug mode, follow by user name, eg. -d ck1 or -debug ck1

            No email will be sent to anyone including the person running the script,
            because there are too MANY emails already. You just need to look for 
            "/wormsrv2/logs/update_caltech.rundate.pID" file. 


B<-help:>      Quess what!


=head3 <RUN update_caltech.pl>

            The script is now set to run on Friday mornings at 007.
            
           


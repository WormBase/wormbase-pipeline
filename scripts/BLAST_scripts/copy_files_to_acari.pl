#!/usr/local/bin/perl5.6.1 -w
#
# copy_files_to_acari.pl
#
# dl
#
# copy the latest dna and wormpep files from wormsrv2 to acari

use strict;
use lib "/wormsrv2/scripts/";  
use Wormbase;
use Getopt::Std;
use vars qw($opt_c $opt_w);
# $opt_w copy wormpep data
# $opt_c copy chromosomes

getopts ('cw');

my $maintainers = "All";
my $dir = "/nfs/acari/wormpipe/BlastDB";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $log = "/wormsrv2/logs/copy_file_to_acari.$rundate.$runtime";
open (LOG,">$log") or die "cant open $log\n";

print LOG "$0\n";
print LOG "started at ",`date`,"with options \n";
print LOG "c - copy chromosomes\n" if ($opt_c);
print LOG "w - create wormpep blast dbase\n" if ($opt_w);
print LOG "no options chosen - does nothing" unless ( $opt_c || $opt_w );
print LOG "\n";
print LOG "=============================================\n";
print LOG "\n";

my @CHROMOSOME = ('I','II','III','IV','V','X');
if ( $opt_c )
  {
    print "-------------------------------------------------------\nCopying new Chromosomal DNA\n\n";
    print LOG  "Removing old DNA files . . . ";
    &test_system_command(  system ("rm -f $dir/*dna") );

    
    print LOG "Copying the new DNA files from the FTP site . . . ";
    print LOG "\n";

    foreach my $chrom (@CHROMOSOME) {
      print LOG "Chromosome $chrom . . . ";
      &test_system_command( system ("scp -r wormpub\@wormsrv2:/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_$chrom.dna $dir/") );
      print LOG "\nProcess chromsome $chrom\n";
     # print LOG " \tunzipping the DNA  for chromosome $chrom . . .";
     # &test_system_command( system ("gzip -d $dir/CHROMOSOME_$chrom.dna.gz") );
      print LOG "\treplacing padding characters by N's in chromosome consensus DNA . .";
      &test_system_command( system ("cat $dir/CHROMOSOME_$chrom.dna | perl -n -e 's/\-/n/g;print' > $dir/tmp_file") );
      print LOG "\tmoving chromosome $chrom from tmp_file  . . ";
      &test_system_command(system ("mv -f $dir/tmp_file $dir/CHROMOSOME_$chrom.dna") );
      
      print LOG "chromosome $chrom complete\n\n";
    }

    print LOG "Finished copying DNA to acari\n\n-------------------------------------------------------\n";
  }

if( $opt_w )
  {
    print LOG "-------------------------------------------------------\nCopy over the latest wormpep relese\n\n";
    print "\nUpdating WORMPEP\n";
    our $WS_version = &get_wormbase_version;
    
    my $wp_file = "wormpep" . $WS_version . ".pep";
    my $wp_old  = "wormpep" . ($WS_version - 1) . ".pep";
    
    print LOG "copying wormpep version $WS_version . . ";
    &test_system_command( system ("scp wormpub\@wormsrv2:/wormsrv2/WORMPEP/wormpep$WS_version/wormpep$WS_version $dirBlastDB/") );
    
    print LOG "and make it blastable. . . ";
    print "making blastable database. . .\n ";
    #note which setdb this command uses
    &test_system_command( system ("/usr/local/pubseq/bin/setdb $dir/$wp_file") );
    
    print LOG "Removed old version of wormpep . . ";
    &test_system_command(  system ("rm -f $dir/${wp_old}*") );

    print LOG "Wormpep copied and ready for action\n\n-------------------------------------------------------\n";
    print "finished\n";
  }


print LOG "\n\nfinished at ",`date`,"\n";
close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = qw(Copy_files_to_acari.pl);
$maintainers = "ar2\@sanger.ac.uk";
&mail_maintainer ($name,$maintainers,$log);
#########################################

exit(0);

#this routine check wether the system command was successfull
sub test_system_command
  {
    my $system_test = shift;
    #NOTE: system calls return '0' if SUCCESSFUL
    if( $system_test == 0){
      print LOG "done\n";}
    else{
      print LOG "FAILED!\n\n";
    }
  }


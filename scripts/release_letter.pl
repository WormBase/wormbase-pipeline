#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# release_letter.pl                            
# 
# by Keith Bradnam                              
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-07-24 10:11:59 $         
#
#Generates a release letter at the end of build.
#Three subroutines are called during the build - 
#release_wormpep by make_wormpep
#release_composition by dump_chromosomes.pl
#release_databases by dbcomp
#These write to a file in autoace/RELEASE_LETTER and are incorperated in to the letter at the end. 
#This allows for overwriting during the build when errors are fixed and scripts rerun




#test the subs#############################
#&release_composition;
#&release_wormpep(19462,20639,1177);#$number_cds $number_total $number_alternate
#&release_databases;
###########################################

use strict;                                      # always use
use lib "/wormsrv2/scripts/";                    # Try to use Wormbase.pm where it useful to do so
use Wormbase;
use Getopt::Std;

# Try to keep different parts of code cleanly separated using comments...
#######################################
# command-line options                #
#######################################

use vars qw($opt_c $opt_d $opt_l);
# $opt_d do release_databases
# $opt_c do release_composition
# $opt_l write the letter out

getopts ('dcl');

##############
# variables  #                                                                   
##############

my $ver = &get_wormbase_version;
my $old_ver = $ver -1;

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainer  = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $date        = `date`;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/finish_build.$rundate";


&release_databases if defined($opt_d);
&release_composition if defined($opt_c);

if( defined($opt_l))
  {
    #make the release letter
    my $release_letter = "/wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver";
    open (RL,">$release_letter");
    
    print RL "New Release of acedb WS$ver, Wormpep$ver and Wormrna$ver $date\n======================================================================\n\n";
    print RL "This directory includes:\n
i)\tdatabase.WS$ver.*.tar.gz     -   compressed data for new release\n
ii)\tmodels.wrm.WS$ver           -   the latest database schema (also in above database files)\n
iii)\tCHROMOSOMES/subdir         -   contains 3 files (DNA, GFF & AGP per chromosome)\n
iv)\tWS$ver-WS$old_ver.dbcomp    -   log file reporting difference from last release\n
v)\twormpep$ver.tar.gz           -   full Wormpep distribution corresponding to WS$ver\n
vi)\twormrna$ver.tar.gz          -   latest WormRNA release containing non-coding RNA's in the genome\n
vii)\tconfirmed_genes.WS$ver.gz  -   DNA sequences of all genes confirmed by EST &/or cDNA\n\n\n";

    print RL "Release notes on the web:\n-------------------------\n";
    print RL "http://www.sanger.ac.uk/Projects/C_elegans/WORMBASE\n\n\n\n";
    
    my $release_dir = ("/wormsrv2/autoace/RELEASE_LETTERS");
    my @release_files = ("$release_dir/dbases","$release_dir/composition","$release_dir/wormpep");
    
    #include all the pre-generated reports
    my $file = shift(@release_files);
    while (defined($file))
      {
	open (READIN, "<$file") || die "cant open $file\n";
	while(<READIN>) {
	  print RL "$_";
	}
	close READIN;
	print RL "\n\n\n";
	$file = shift(@release_files);
      }
    
    
    # Get the GeneModel corrections
    my %cam_introns;
    $cam_introns{$ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/GFF/CHROMOSOME_*.check_intron_cam.gff | wc -l`;
    $cam_introns{$old_ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/GFF/CHROMOSOME_*.check_intron_cam.gff | wc -l`;
    $cam_introns{change} = $cam_introns {$ver} - $cam_introns {$old_ver};
    
    my %stl_introns;
    $stl_introns{$ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/GFF/CHROMOSOME_*.check_intron_stl.gff | wc -l`;
    $stl_introns{$old_ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/GFF/CHROMOSOME_*.check_intron_stl.gff | wc -l`;
    $stl_introns{change} = $stl_introns {$ver} - $stl_introns {$old_ver};
    
    print RL "GeneModel correction progress WS$old_ver -\> WS$ver\n-----------------------------------------\n
Confirmed introns not is a CDS gene model;\n\n\t\t+---------+--------+\n\t\t| Introns | Change |\n\t\t+---------+--------+\n";
    printf RL ("Cambridge\t|  %5d  |  %4d  |\n", $cam_introns{$ver},$cam_introns{change});
    printf RL ("St Louis \t|  %5d  |  %4d  |\n", $stl_introns{$ver},$stl_introns{change});
    print RL "\t\t+---------+--------+\n\n\n";
    
    
    #Members of known repeat families that overlap predicited exons
    my %cam_repeats;
    $cam_repeats{$ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/Checks/CHROMOSOME_*.repeat_in_exon_cam | wc -l`;
    $cam_repeats{$old_ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/Checks/CHROMOSOME_*.repeat_in_exon_cam | wc -l`;
    $cam_repeats{change} = $cam_repeats{$ver} - $cam_repeats{$old_ver};
    
    my %stl_repeats;
    $stl_repeats{$ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/Checks/CHROMOSOME_*.repeat_in_exon_stl | wc -l`;
    $stl_repeats{$old_ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/Checks/CHROMOSOME_*.repeat_in_exon_stl | wc -l`;
    $stl_repeats{change} = $stl_repeats{$ver} - $stl_repeats{$old_ver};
    
    print RL "Members of known repeat families that overlap predicted exons;\n\n\t\t+---------+--------+\n\t\t| Introns | Change |\n\t\t+---------+--------+\n";
    printf RL ("Cambridge\t|  %5d  |  %4d  |\n", $cam_repeats{$ver},$cam_repeats{change});
    printf RL ("St Louis \t|  %5d  |  %4d  |\n", $stl_repeats{$ver},$stl_repeats{change});
    print RL "\t\t+---------+--------+\n\n\n";
    
  

    #Synchronisation with GenBank / EMBL
    #grep ERROR /wormsrv2/autoace/yellow_brick_road/CHROMOSOME_X.agp_seq.log
    my @chromosomes = ("I","II","III","IV","V","X");
    my $csome = shift @chromosomes;
    print RL "Synchronisation with GenBank / EMBL:\n------------------------------------\n\n";
    my $check = 0;
    while ($csome)
      {
	my $errors = `grep ERROR /wormsrv2/autoace/yellow_brick_road/CHROMOSOME_$csome.agp_seq.log`;
	while( $errors =~ m/for\s(\p{IsUpper}\w+)/g ) {
	  print RL "CHROMOSOME_$csome\tsequence $1\n";
	  $check = 1;
	}
	$csome = shift @chromosomes;
      }
    if ($check == 0){
      print RL "No synchronisation issues\n\n";
    }
    
    print RL "----------------------------------------------------------------\n\n\n\n";
    print RL "New Data:\n---------\n\n\n";
    print RL "New Fixes:\n----------\n\n\n";
    print RL "Known Problems:\n--------------\n\n\n";
    print RL "Other Changes:\n--------------\n\n";
    print RL "Proposed Changes / Forthcoming Data:\n------------------------------------\n\n\n";
    
    
    print RL "Quick installation guide for UNIX/Linux systems\n-----------------------------------------------\n\n";
    print RL "1. Create a new directory to contain your copy of WormBase,\n\te.g. /users/yourname/wormbase\n\n";
    print RL "2. Unpack and untar all of the database.*.tar.gz files into\n\tthis directory. You will need approximately 2-3 Gb of disk space.\n\n";
    print RL "3. Obtain and install a suitable acedb binary for your system\n\t(available from www.acedb.org).\n\n";
    print RL "4. Use the acedb 'xace' program to open your database, e.g.\n\ttype 'xace /users/yourname/wormbase' at the command prompt.\n\n";
    print RL "5. See the acedb website for more information about acedb and\n\tusing xace.\n\n";
    
    
    print RL "____________  END _____________\n";
    
    print "DONT FORGET TO FILL IN THE LAST FEW FIELDS IN THE LETTER\n found at /wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver\n";
    
    my $name = "Release Letter for WS$ver";
    &mail_maintainer($name,$maintainer,$release_letter);
  }
# Always good to cleanly exit from your script
exit(0);


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item release_letter.pl  [-c -d -l]

=back

This script:
q
will kindly take away the pain of having to write a release letter at the end of the build
If the script is being run at the end of a build in which the 3 sub parts have been generated then use the B<-l> option.
Otherwise generate the sequence comparison and database comparison sections with the B<-c> and B<-d> options.

I<release_letter.pl MANDATORY arguments:> B<NONE>

I<script_template.pl  OPTIONAL arguments:>


B<-d> create the database comparison file.

B<-c> create the sequence composition comparison file.

B<-l> actually write the letter out.

=back
=over 4

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

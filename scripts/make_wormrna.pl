#!/usr/local/bin/perl5.6.0 -w
#
# make_wormrna.pl
# 
# Usage : make_wormrna.pl -r <release_number>
#
# Builds a wormrna data set from the current autoace database
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-04-25 10:57:25 $


#################################################################################
# variables                                                                     #
#################################################################################

$| = 1;
use strict;
use vars qw($opt_r $opt_d $opt_h);
use Getopt::Std;
use IO::Handle;
use Ace;
use Socket;
use lib '/wormsrv2/scripts/';
use Wormbase;

    
 #######################################
 # Script variables (run)              #
 #######################################

my $maintainer = "dl1\@sanger.ac.uk, krb\@sanger.ac.uk, kj2\@sanger.ac.uk";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;


 #######################################
 # command-line options                #
 #######################################

getopts ('dhr:');

 #######################################
 # sanity checks for the input         #
 #######################################

&error(0) if ($opt_h);
&error(1) if ((length $opt_r) == 0);
&error(2) if (($opt_r =~ /\D+/));

 #######################################
 # release data                        #
 #######################################

my $release = $opt_r; 
my $release_date     = &get_wormbase_release_date("long");
my $old_release = $release-1;
my $debug = $opt_d;

if ($debug) {$maintainer = "dl1\@sanger.ac.uk";}


my $dbdir     = "/wormsrv2/autoace";
my $wrdir     = "/wormsrv2/WORMRNA/wormrna$old_release";
my $new_wrdir = "/wormsrv2/WORMRNA/wormrna$release";
my $tace      = "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/tace";

$ENV{'ACEDB'} = $dbdir;

 ########################################
 # Open logfile                         #
 ########################################

my $log="/wormsrv2/logs/make_wormrna.$release.$rundate.$$";
open (LOG,">$log") || &error(3);
LOG->autoflush();

print LOG "# make_wormrna.pl\n";
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";
print LOG "Wormrna version  : wormrna$opt_r\n\n";
print LOG "=============================================\n";
print LOG "\n";

 ##################################################
 # Make new directory for current release         #
 ##################################################

mkdir ("$new_wrdir" , 0755) || &error(4);   

print LOG "# $runtime : making wormrna$release for $rundate\n\n";

 ###############################################
 # retrieve the desired RNA sequence objects   #
 ###############################################

$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "# $runtime : connect to primary database\n";

my $db = Ace->connect (-path => $dbdir, -program => $tace) || &error(5);
my @dotnames_1 = $db->fetch (-query => 'FIND Genome_sequence ; FOLLOW Subsequence ; where (Method = rna) OR (Method=tRNAscan-SE-1.11)');
my @dotnames_2 = $db->fetch (-query => 'FIND Sequence *LINK* ; FOLLOW Subsequence ; where (Method = rna) OR (Method=tRNAscan-SE-1.11)');

push (@dotnames_1 , @dotnames_2);
@dotnames_1 = sort @dotnames_1;
my $count = scalar(@dotnames_1)+scalar(@dotnames_2);
$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "=> " . scalar(@dotnames_1) . " RNA sequences\n";
print LOG "=> " . scalar(@dotnames_2) . " of which are attached to LINK objects\n";
print LOG "# $runtime : finished connection to database\n\n";

 ###########################################################################
 # get the rna sequence, write a rna.fasta file,
 ###########################################################################
 
$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "# $runtime : write wormrna.rna file\n\n";

my $obj = "";

open (DNA , ">$new_wrdir/wormrna$release.rna") || &error(6); 
my (%dot2num , @dotnames , @c_dotnames);

foreach my $dot (@dotnames_1) {    
  undef (my $dna);
  undef (my $locus);
  undef (my $brief_id);
  undef (my $method);
    
  print LOG "Extracting RNA sequence $dot\n";
  $obj = $db->fetch(Sequence=>"$dot");
  
  $locus = $obj->Locus_genomic_seq(1);
  if ((!defined ($locus)) || ($locus eq "")) {
    print LOG "No locus designation for $dot\n";
    undef ($locus);
  }
  
  $brief_id = $obj->Brief_identification(1);
  if ((!defined ($brief_id)) || ($brief_id eq "")) {
    print LOG "No Brief_id for $dot\n";
    undef ($brief_id);
  }
  
  $dna = $obj->asDNA();
  if ((!defined ($dna)) || ($dna eq "")) {
    print LOG "cannot extract dna sequence for $dot\n";
  }
  $dna =~ /^>(\S+)\s+(\w.*)/s ; my $dseq = $2 ; $dseq =~ tr/a-z/A-Z/ ; $dseq =~ tr /T/U/ ; $dseq =~ s/\s//g;

  my $rseq = &reformat($dseq);
  if (defined $locus) {
    print DNA ">$dot $brief_id locus:$locus\n$rseq";
  }
  else {
    print DNA ">$dot $brief_id\n$rseq";
  }
  
  $obj->DESTROY();
  
}   

close DNA;
chmod (0444 , "$new_wrdir/wormrna$release.rna") || print LOG "cannot chmod $new_wrdir/wormrna$release.rna\n";

$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "# $runtime : finished writing wormrna.rna file\n";

###########################################################################
# Create the associated README file          
###########################################################################

print LOG "# $runtime : creating README file\n\n";
open (README , ">$new_wrdir/README") || &error(7); 

my $readme = <<END;
WormRNA
-------
WormRNA is an additional database that accompanies the main WormBase 
database (see http://www.wormbase.org for more details) and simply comprises
the sequences of all known non-coding RNA molecules in the C. elegans genome.  

This release (WormRNA$release) corresponds to release WS$release of WormBase.
The accompanying file (wormrna$release.rna) contains $count RNA sequences in FASTA
format.

WormBase group, Sanger Institute
$release_date

END

print README "$readme";
close(README);

 ##############################
 # mail $maintainer report    #
 ##############################

close LOG;

open (mailLOG, "|/usr/bin/mailx -s \"WormBase Report: make_wormrna.pl\" $maintainer ");
open (readLOG, "<$log");
while (<readLOG>) {
    print mailLOG $_;
}
close readLOG;
close mailLOG;

 ##############################
 # hasta luego                #
 ##############################

exit(0);

#################################################################################
# Subroutines                                                                   #
#################################################################################


sub reformat {
    my $in_string = shift;
    my $out_string = "";

    my $string_len = length ($in_string);
    my $lines = int ($string_len / 60) ;

    for (my $i = 0; $i <= $lines; $i++) {
	$out_string = $out_string . substr($in_string,($i*60),60) . "\n";
    }
    return ($out_string);
}


 ##########################
 # run details            #
 ##########################

sub run_details {
    print "# make_wormrna.pl\n";
    print "# run details    : $rundate $runtime\n";
    print "\n";
    print "wormrna version  : wormrna$opt_r\n";
    print "Primary database : $dbdir\n\n";

    if ($opt_d) {
	print "Usage : make_wormrna.pl [-options]\n";
	print "=============================================\n";
	print " -r <int>     : release version number\n";
	print " -h           : help pages   \n";
	print " -d           : verbose (debug) mode\n";
	print "=============================================\n";
	print "\n";
    }
} # end of sub 'run details'

 ##########################
 # errors from the script #
 ##########################

sub error {
    my $error = shift;
    # Error 0 - help page
    if ($error == 0) {
        exec ('perldoc',$0);
        exit (0);
    }
    # Error  1 - no wormrna release number
    elsif ($error == 1) {
        # No wormrna release number file
	&run_details;
        print "=> No wormrna release number supplied\n\n";
        exit(0);
    }
    # Error  2 - invalid wormrna release number
    elsif ($error == 2) {
        # Invalid wormrna release number file
        &run_details;
        print "=> Invalid wormrna release number supplied.\n=> Release number must be an interger (e.g. 30)\n\n";
	exit(0);
    }
    # Error  3 - cannot open new wp.log file 
    elsif ($error == 3) {
        &run_details;
        print "=> Failed to create a new wp.log for wormrna release wormrna$release\n\n";
        $runtime = `date +%H:%M:%S`; chomp $runtime;
        print LOG "=> Exiting at $rundate $runtime\n";
        close LOG;
        &mail_maintainer("WormBase Report: make_wormrna.pl",$maintainer,$log);
    }
    # Error  4 - cannot create new wormrna directory 
    elsif ($error == 4) {
        &run_details;
	print "=> Failed to create a new directory for wormrna release wormrna$release\n\n";
        $runtime = `date +%H:%M:%S`; chomp $runtime;
	print LOG "=> Failed to create a new directory for wormrna release wormrna$release\n\n";
	print LOG "=> Exiting at $rundate $runtime\n";
        close LOG;
        &mail_maintainer("WormBase Report: make_wormrna.pl",$maintainer,$log);
    }
    # Error  5 - cannot connect to ACEDB database 
    elsif ($error == 5) {
        &run_details;
        print "=> Failed to connect to primary database 'dbdir'\n\n";
        $runtime = `date +%H:%M:%S`; chomp $runtime;
        print LOG "=> Exiting at $rundate $runtime\n";
        close LOG;
        &mail_maintainer("WormBase Report: make_wormrna.pl",$maintainer,$log);
    }
    # Error  6 - cannot open new wp.log file 
    elsif ($error == 6) {
        &run_details;
        print "=> Failed to create a new wormrna.rna for wormrna release wormrna$release\n\n";
        $runtime = `date +%H:%M:%S`; chomp $runtime;
        print LOG "=> Exiting at $rundate $runtime\n";
        close LOG;
        &mail_maintainer("WormBase Report: make_wormrna.pl",$maintainer,$log);
    }
    # Error  7 - cannot open new wp.log file 
    elsif ($error == 7) {
        &run_details;
        print "=> Failed to create README file for wormrna release wormrna$release\n\n";
        $runtime = `date +%H:%M:%S`; chomp $runtime;
        print LOG "=> Exiting at $rundate $runtime\n";
        close LOG;
        &mail_maintainer("WormBase Report: make_wormrna.pl",$maintainer,$log);
    }

    exit(1);
} # end of sub 'error'



__END__

=pod

=head2   NAME - make_wormrna.pl


=head1 USAGE

=over 4

=item make_wormrna.pl [-options]

=back

make_wormrna.pl will generate a rna data set from the autoace
database directory.

make_wormrna.pl mandatory arguments:

=over 4

=item -r, release number

=back

make_wormrna.pl OPTIONAL arguments:

=over 4

=item -h, Help page

=item -d, Verbose/Debug mode

=back

=head1 EXAMPLES:

=over 4

=item make_wormrna.pl -r 4

=back

Creates a new wormrna data set in the (new) /wormsrv2/WORMRNA/wormrna4 directory

=cut





















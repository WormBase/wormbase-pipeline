#!/usr/local/bin/perl5.8.0 -w
#
# make_wormrna.pl
# 
# Usage : make_wormrna.pl -r <release_number>
#
# Builds a wormrna data set from the current autoace database
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $


#################################################################################
# variables                                                                     #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use IO::Handle;
use Ace;
use Socket;

    
######################################
# variables and command-line options # 
######################################

my ($help, $debug, $release);
my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;
my $log;
my $errors      = 0;    # for tracking how many errors there are

GetOptions ("help"      => \$help,
	    "release=s" => \$release,
            "debug=s"   => \$debug);


# Display help if required
&usage("Help") if ($help);
&usage("releasename") if (!$release);
&usage("releasename") if ($release =~ /\D+/);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;


#######################################
# release data                        #
#######################################

my $release_date = &get_wormbase_release_date("long");
my $old_release  = $release-1;

my $dbdir     = "/wormsrv2/autoace";
my $wrdir     = "/wormsrv2/WORMRNA/wormrna$old_release";
my $new_wrdir = "/wormsrv2/WORMRNA/wormrna$release";
my $tace      = &tace;

$ENV{'ACEDB'} = $dbdir;


# Make new directory for current release      
mkdir ("$new_wrdir" , 0755) || die "Couldn't create $new_wrdir\n";   


###############################################
# retrieve the desired RNA sequence objects   #
###############################################

my $db = Ace->connect (-path => $dbdir, -program => $tace) || die "Couldn't connect to $dbdir\n";
# Get RNA genes, but not other Transcript objects
my @transcripts = $db->fetch (-query => 'FIND elegans_RNA_genes where NOT Ignore');

@transcripts = sort @transcripts;
my $count = scalar(@transcripts);
print LOG "Finding ". $count . " RNA sequences, writing output file\n\n";


###########################################################################
# get the rna sequence, write a rna.fasta file,
###########################################################################

open (DNA , ">$new_wrdir/wormrna$release.rna") || die "Couldn't write wormrna$release.rna file\n"; 
my (%dot2num , @dotnames , @c_dotnames);

foreach my $transcript (@transcripts) {    
  undef (my $dna);
  undef (my $cgc_name);
  undef (my $brief_id);
  undef (my $method);
  undef (my $gene);
    
  my $obj = $db->fetch(Transcript=>"$transcript");
  

  # Grab Brief_identification
  $brief_id = $obj->Brief_identification;
  if ((!defined ($brief_id)) || ($brief_id eq "")) {
    print LOG "ERROR: No Brief_id for $transcript\n";
    $errors++;
    undef ($brief_id);
  }
  
  $dna = $obj->asDNA();
  if ((!defined ($dna)) || ($dna eq "")) {
    print LOG "ERROR: cannot extract dna sequence for $transcript\n";
    $errors++;
    # can't include in WormRNA if there is no DNA!
    next; 
  }
  $dna =~ /^>(\S+)\s+(\w.*)/s; 
  my $dseq = $2; 
  $dseq =~ tr/a-z/A-Z/; 
  $dseq =~ tr /T/U/; 
  $dseq =~ s/\s//g;
  
  # Grab locus name if present
  $gene = $obj->Gene;
  $cgc_name = $gene->CGC_name if (defined $gene);
  

  my $rseq = &reformat($dseq);
  if ($cgc_name) {
    print DNA ">$transcript $brief_id locus:$cgc_name\n$rseq";
  }
  elsif($brief_id) {
    print DNA ">$transcript $brief_id\n$rseq";
  }
  else{
    print DNA ">$transcript\n$rseq";
  }
  
  $obj->DESTROY();
  
}   

close DNA;
chmod (0444 , "$new_wrdir/wormrna$release.rna") || print LOG "cannot chmod $new_wrdir/wormrna$release.rna\n";


###########################################################################
# Create the associated README file          
###########################################################################

open (README , ">$new_wrdir/README") || die "Couldn't creat README file\n"; 

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


##################
# Tidy up
##################

close(README);
print LOG "\n============================================\n";
print LOG "Finished at ",&runtime, "\n";
close(LOG);
$db->close;

if($errors == 0){
   &mail_maintainer ("make_wormrna.pl WS$release",$maintainers,$log);
}
elsif($errors == 1){
   &mail_maintainer ("make_wormrna.pl WS$release : $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer ("make_wormrna.pl WS$release : $errors ERRORS!!!",$maintainers,$log);
}

exit(0);




##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "Started at ",&runtime," \n";
  print LOG "=============================================\n";
  print LOG "\n";

}

#################################################################

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



#################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  # Error  2 - invalid wormrna release number
  elsif ($error eq "releasename") {
    # Invalid wormrna release number file
    print "=> Invalid/missing wormrna release number supplied.\n";
    print "=> Release number must be an interger (e.g. 30)\n\n";
    exit(0);
  }

}




__END__

=pod

=head2   NAME - make_wormrna.pl


=head1 USAGE

=over 4

=item make_wormrna.pl [-options]

=back

make_wormrna.pl will generate a rna data set from the autoace
database directory.  Finds RNA genes in Transcript class, but filters out 
any other Transcript objects which represent full length transcripts of
coding genes (i.e. where Method = Transcript)

make_wormrna.pl mandatory arguments:

=over 4

=item -release <release number>

=back

make_wormrna.pl OPTIONAL arguments:

=over 4

=item -help, Help page

=item -debug <username> = Verbose/Debug mode

=back

=head1 EXAMPLES:

=over 4

=item make_wormrna.pl -release 4

=back

Creates a new wormrna data set in the (new) /wormsrv2/WORMRNA/wormrna4 directory

=cut





















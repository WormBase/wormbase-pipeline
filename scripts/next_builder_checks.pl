#!/usr/local/bin/perl5.8.0 -w
# next_builder_checks.pl
#
# by Keith Bradnam
#
# A simple script to send a check list to the person who will be performing the next
# build to check the current build
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2009-12-23 11:09:29 $
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;
use File::Compare;

my $store;                                          # to specify saved commandline arguments
my $maintainers = "All";
my ($help,$debug, $species);
my ($clones, $pfam, $seq, $wormpep);


GetOptions(
	   "help"    => \$help,
	   "debug=s" => \$debug,
	   'store=s' => \$store,
	   'clones'  => \$clones,
	   'pfam'    => \$pfam,
	   'seq'     => \$seq,
	   'wormpep' => \$wormpep, 
	   'species:s'=>\$species
);

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $debug, -organism => $species ) }

##########################################
# Variables Part II (depending on $wb)    #
###########################################
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user
my $WS_current = $wb->get_wormbase_version;
my $wormbase = $wb;

# Use debug mode?
if ($debug) {
    $wb->debug($debug);
}

my $log = Log_files->make_build_log($wb);

$wormpep=$pfam=$seq=$clones=1 unless($wormpep or $pfam or $seq or $clones);
#only do some checks if elegans.
unless ($wb->species eq 'elegans') {
    $pfam=$seq=$clones=0;
}

$log->write_to("Checking ".$wb->full_name.": - ".$wb->orgdb."\n\n");
my $ace;
my $aceold;
#only connect once and if required.
if($wormpep or $clones or $pfam){ 
  $ace    = Ace->connect('-path' => $wormbase->orgdb);
  #$aceold = Ace->connect('-path' => $wormbase->database('currentdb'));
  $aceold = Ace->connect('-path' => $wormbase->database('WS209'));
}
if($clones) {
  $log->write_to("##################################\nChecking clones . .\n\n");
  my @clones = qw(C25A1 F56A3 C04H5 B0432 C07A9 F30H5 C10C6 B0545 C12D8 K04F1 C02C6 AH9);

  #expected no of objects per type
  my %data_types = ('Homol_data' => 11,
		    'Feature_data' => 3,
		   );

  foreach my $clone (@clones) {
    $log->write_to("checking clone $clone \n");	
    my $query = "find Homol_data $clone*";
    my @hd    = $ace->fetch(-query => $query);
    my @hdold = $aceold->fetch(-query => $query);
    
    &check_for_missing_data(\@hd, \@hdold, 'Homol_data');
    
    # check the blastx Homol_data
    my $hd;
    my $count = 0;
    foreach my $hd (@hd) {
      if ($hd->name =~ /wublastx/) {
	print $hd->name,"\n";
	$count++;
	if (defined $hd->Pep_homol(3) ) { # check for presence of a score.
	  #$log->write_to($hd->name." OK\n");
	} else {
	  $log->error("ERROR: ".$hd->name." missing data\n");
	}
      }
    }
    
    if($count < 11 ) {
      $log->error("ERROR: $clone has stuff missing\n");
    }
    
    #check Feature_data
    $query = "find Feature_data $clone*";
    @hd    = $ace->fetch(-query => $query);
    @hdold = $aceold->fetch(-query => $query);

    &check_for_missing_data(\@hd, \@hdold, 'Feature_data');

    $count = 0;
    foreach my $hd (@hd) {
      print $hd->name,"\n";
      $count++;
      if( scalar $hd->Feature(2)->row > 3 ) {
	#$log->write_to($hd->name." OK\n");
      } else {
	$log->error("ERROR: ".$hd->name." missing data\n");
      }
    }
    if($count < 3 ) {
      $log->error("ERROR: $clone has stuff missing\n");
    }
  }
}

if ($seq) {
    $log->write_to("\n##################################\nChecking sequence composition . .\n\n");
    #check composition.all for n's
    my $file = $wormbase->chromosomes ."/composition.all";
    undef $/;
    open (ALL,"<$file") or $log->log_and_die("cant open $file : $!\n");
    my $in = <ALL>;
    close ALL;
    $/ = "\n";
    $in =~ /n\s+(\d+)/;
    if($1 == 0){
	$log->write_to("no n's thanksfully\n");
    }else {
	$log->error("\n\nthere are n's in the genome!\n\n");
    }

    #check composition is same as start of build
    $wormbase->run_command("ls ".$wormbase->orgdb."/CHROMOSOMES/*.dna | grep -v masked |grep -v Mt| xargs composition > /tmp/comp", $log);
    if(compare($file,"/tmp/comp") == 0) { 
	$log->write_to("composition same as start of build\n\n");
    }else {
	$log->error("composition has changed during build!\n\n");
    }
    $wormbase->run_command("rm -f /tmp/comp", $log);
}

if($pfam){
    $log->write_to("\n##################################\nChecking PFAM motifs . .\n\n");
    #check PFAM motifs have title
    my $query = "query find motif PFAM* where !Title";
    my $no_tits = $ace->count("-query"=>$query);
    if($no_tits > 20) {
	$log->error("$no_tits PFAM domains are missing a Title\n");
    }else {
	$log->write_to("Only $no_tits PFAM domains are missing a Title\n");
    }
}

if($wormpep){
    $log->write_to("\n##################################\nChecking new proteins . .\n\n");
    #check that new wormpep entries have domains and blastp
    my $new_pep_file = $wormbase->wormpep."/new_entries.".$wormbase->get_wormbase_version_name;
    open (PEP,"<$new_pep_file") or $log->log_and_die("cant open $new_pep_file : $!\n");
    my @newpeps;
    while(<PEP>) {
	if(/>(\S+)/) {
	    push(@newpeps, $1);
	}
    }
    close PEP;
    my ($Pcount, $Mcount); #pephomol motifhomol
    foreach my $pep(@newpeps){
	my $pepObj = $ace->fetch('Protein' => $wormbase->wormpep_prefix.":$pep");
	$Pcount++ if (defined $pepObj->Pep_homol);
	#print STDERR $pepObj->name," P\n" unless(defined $pepObj->Pep_homol);
	$Mcount++ if (defined $pepObj->Motif_homol);
	#print STDERR $pepObj->name," M\n" unless(defined $pepObj->Motif_homol);
    }
    ($Pcount / scalar @newpeps < 0.5) ?
	$log->error("ERROR: more than third ($Pcount / ".scalar @newpeps.") of new proteins dont have Pep_homols\n") :
	$log->write_to("new proteins Pep_homols look ok\n");

    ($Mcount / scalar @newpeps < 0.3) ?
	$log->error("ERROR: only ($Mcount / ".scalar @newpeps.") of new proteins have Motif_homols\n") :
	$log->write_to("new proteins Motif_homols look ok\n");
}

$ace->close if(defined $ace);

$log->mail;
exit;



my $log_msg= <<'LOGMSG';
1) The following 12 clones are representative of the whole genome in that
they include one Sanger and one St. Louis clone for each chromosome.  Check
each clone to ensure that it contains BLAT data (EST and mRNA), BLAST data,
waba data, gene models, UTRs etc.  Also check for presence of tandem and inverted
repeats which have gone missing in the past

i)    C25A1
ii)   F56A3
iii)  C04H5
iv)   B0432
v)    C07A9
vi)   F30H5
vii)  C10C6
viii) B0545
ix)   C12D8
x)    K04F1
xi)   C02C6
xii)  AH9

2) If the genome sequence has changed, you should inspect the clones containing
those changes to see if there are any strange errors (e.g. duplicate sets of data
which are slightly out of sync.)

3) Check ~wormpub/BUILD/autoace/CHROMOSOMES/composition.all - are there any non-ATCGN
characters

4a) Check that the latest WormPep proteins have proper protein and motif homologies
This has been a problem in some builds where all new WormPep proteins have not got any
BLAST analyses.  Pick a few random Wormpep proteins and especially check that all of
the various blastp homologies are there (human, fly, worm, yeast etc.) and try to
check at least one protein from the ~wormpub/BUILD/WORMPEP/wormpepXXX/new_entries.WSXXX file

4b) Now that we have a curated set of brigpep, should do this periodically for
C. briggase protein objects too...these now have their own set of blastp hits

5) Check PFAM Motif objects have a title tag. It is a problem if there are more than about 20.

6) Run: 
  ls ~wormpub/BUILD/autoace/CHROMOSOMES/*.dna | grep -v masked |grep -v Mt| xargs composition
Make sure this is the same as it was at the start of the build:
  cat ~wormpub/BUILD/autoace/CHROMOSOMES/composition.all
Bad Homol objects can lead to errors esp when chromosome length has been reduced

Thats all...for now!  If you are satisfied the build is ok, please inform the person
building the database. Please continue to add to this list as appropriate.

========================================================================================
LOGMSG

$log->write_to($log_msg);

$log->mail("$maintainers", "Please check the ongoing build of WS${WS_current}");

exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

##################################################################

# check for missing data compared to currentdb, 
# (see Perl Cookbook p.126)

sub check_for_missing_data {

  my ($hd_aref, $hdold_aref, $data_name) = @_;    

  my %seen=();
  my @oldonly=();
  foreach my $hd (@{$hd_aref}) {$seen{$hd->name}=1} # get the name of objects in the clone in the Build
  foreach my $hdold (@{$hdold_aref}) {
    unless ($seen{$hdold->name}) {push(@oldonly, $hdold->name)} # compare them to the objects in the clone in currentDB
  }
  if (@oldonly) {
    foreach my $hd (@oldonly) {
      $log->error("ERROR: ".$hd->name." missing $data_name data compared to currentdb\n");
    }
  }
}

##################################################################

__END__

=pod

=head1 NAME - next_builder_checks.pl

=head1 USAGE

=over 4

=item next_builder_checks.pl --user <user>

=back

This script simply sends a list of check items to the next person doing the build.
They should 'sign off' on each build and hand back to the main person when they
are happy all is ok.

=item MANDATORY arguments: -user <valid unix username>

Needed to send email

=back

=over 4

=item OPTIONAL arguments: -help, -debug <user>


 
=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk


=cut

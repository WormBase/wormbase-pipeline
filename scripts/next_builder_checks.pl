#!/usr/local/bin/perl5.8.0 -w
# next_builder_checks.pl
#
# by Keith Bradnam
#
# A simple script to send a check list to the person who will be performing the next
# build to check the current build
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2012-09-24 09:06:49 $
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
my ($clones, $pfam, $seq, $wormpep, $test);


GetOptions(
	   "help"    => \$help,
	   "debug=s" => \$debug,
	   'store=s' => \$store,
	   'clones'  => \$clones,
	   'pfam'    => \$pfam,
	   'seq'     => \$seq,
	   'wormpep' => \$wormpep, 
	   'species:s'=>\$species,
	   'test'    => \$test,
);

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species ) }

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
if (
    $wb->species eq 'brenneri' ||
    $wb->species eq 'remanei' ||
    $wb->species eq 'japonica' ||
    $wb->species eq 'briggsae' ||
    $wb->species eq 'pristionchus') {
  $pfam=$seq=0;
} elsif ($wb->species ne 'elegans') {
  $log->log_and_die("\tSORRY, Don't know how to check this species: $species\n");
}

$log->write_to("Checking ".$wb->full_name.": - ".$wb->orgdb."\n\n");
my $ace;
my $aceold;
#only connect once and if required.
if($wormpep or $clones or $pfam){ 
  print "Connecting to Ace\n";
  $ace    = Ace->connect('-path' => $wormbase->orgdb);
  print "Connecting to currentdb Ace\n";
  $aceold = Ace->connect('-path' => $wormbase->database('current'));
}

# these clones are chosen because they all have protein matches to all of the protein databases blasted against
if($clones) {
  $log->write_to("##################################\nChecking clones . .\n\n");
  my @clones;
  if ($wb->species eq 'elegans') {
    @clones = qw(C25A1 F56A3 C04H5 B0432 C07A9 F30H5 C10C6 B0545 C12D8 K04F1 C02C6 AH9);
  } elsif ($wb->species eq 'brenneri') {
    @clones = qw(Cbre_Contig1 Cbre_Contig10 Cbre_Contig20 Cbre_Contig50 Cbre_Contig100 Cbre_Contig200  Cbre_Contig400 Cbre_Contig600 Cbre_Contig800);
  } elsif ($wb->species eq 'briggsae') {
    # briggae contains a mixture of data on chromosomes and supercontigs, so include both
    @clones = qw(cb25.fpc0002 cb25.fpc0011c cb25.fpc0081 cb25.fpc0143a chrI chrII);
  } elsif ($wb->species eq 'remanei') {
    @clones = qw(Crem_Contig0 Crem_Contig10 Crem_Contig15 Crem_Contig30 Crem_Contig100 Crem_Contig200 Crem_Contig300 Crem_Contig500 Crem_Contig800);
  } elsif ($wb->species eq 'japonica') {
    @clones = qw(Cjap.Contig0 Cjap.Contig10 Cjap.Contig15 Cjap.Contig30 Cjap.Contig100 Cjap.Contig200 Cjap.Contig300 Cjap.Contig500 Cjap.Contig800);
  } elsif ($wb->species eq 'pristionchus') {
    @clones = qw(Ppa_Contig0 Ppa_Contig10 Ppa_Contig15 Ppa_Contig30 Ppa_Contig100 Ppa_Contig200);
  }

  foreach my $clone (@clones) {
    $log->write_to("\n##################################\nchecking clone $clone\n");	
    my $query = "find Homol_data \"$clone:*\"";
    my @hd    = $ace->fetch(-query => $query);
    my @hdold = $aceold->fetch(-query => $query);

    # add in the BLAT Homol_data objects in brenneri and briggsae
    if ($wb->species eq 'brenneri' || 
	$wb->species eq 'briggsae' || 
	$wb->species eq 'remanei') {
      $query = "find Homol_data \"*:${clone}_*\"";
      push @hd, $ace->fetch(-query => $query);
      push @hdold, $aceold->fetch(-query => $query);
    }
    if ($wb->species eq 'briggsae') {
      $query = "find Homol_data \"*:${clone}\"";
      push @hd, $ace->fetch(-query => $query);
      push @hdold, $aceold->fetch(-query => $query);
    }
    
    &check_for_missing_data(\@hd, \@hdold, 'Homol_data', 'currentdb');
    
    # check the blastx Homol_data
    my $hd;
    my $count = 0;
    foreach my $hd (@hd) {
      if ($hd->name =~ /wublastx/) {
	print $hd->name,"\n";
	$count++;
	# check for presence of an alignment of one of this type of protein to the clone.
	if (defined $hd->Pep_homol(3) ) { 
	  #$log->write_to($hd->name." OK\n");
	} else {
	  $log->error("\tERROR: Homol_data ".$hd->name." does not contain any wublastx alignments\n");
	}
      }
    }
    
    # check for 11 wublastx Homol_data objects (fly, brenenri, briggsae,
    # human, japonica, pristionchus, remanei, slimSwissProt,
    # slimTrEmbl, worm, yeast)
    my @expected = qw(fly brenneri briggsae human japonica pristionchus remanei slimSwissProt slimTrEmbl worm yeast);
    &check_for_missing_data2(\@hd, \@expected, 'Feature_data', 'what is expected');

#    if($count < 11) {
#      $log->error("\tERROR: $clone has wublastx Homol_data objects missing\n");
#    }
    
    #check Feature_data
    $query = "find Feature_data \"$clone:*\"";
    @hd    = $ace->fetch(-query => $query);
    @hdold = $aceold->fetch(-query => $query);

    &check_for_missing_data(\@hd, \@hdold, 'Feature_data', 'currentdb');

    $count = 0;
    foreach my $hd (@hd) {
      if (! defined $hd->Feature(2)) {$log->write_to("Undefined object for ".$hd."\n");next}
      print $hd->name,"\n";
      $count++;
      # check the Feature_data line
      #$log->write_to("Testing line for ".$hd->name."\n");
      if (scalar $hd->Feature(2)->row >= 3 ) {
	#$log->write_to($hd->name." OK\n");
      } else {
	$log->error("\tERROR: ".$hd->name." missing clone-length data?\n");
      }
    }
    
    if ($wb->species eq 'briggsae') {
      @expected = qw(TRF Dust); # briggsae inverted feature_data is on the clones, not the chromosomes
    } else {
      @expected = qw(TRF Dust inverted);
    }

    &check_for_missing_data2(\@hd, \@expected, 'Feature_data', 'what is expected');

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
    if ($newpeps[0]){ # else you get a funky division by zero
     ($Pcount / scalar @newpeps < 0.5) ?
	$log->error("ERROR: more than third ($Pcount / ".scalar @newpeps.") of new proteins dont have Pep_homols\n") :
	$log->write_to("new proteins Pep_homols look ok\n");

     ($Mcount / scalar @newpeps < 0.3) ?
	$log->error("ERROR: only ($Mcount / ".scalar @newpeps.") of new proteins have Motif_homols\n") :
	$log->write_to("new proteins Motif_homols look ok\n");
    }
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

  my ($hd_aref, $hdold_aref, $data_name, $compared_to) = @_;    

  my %seen=();
  my @oldonly=();
  foreach my $hd (@{$hd_aref}) {$seen{$hd->name}=1} # get the name of objects in the clone in the Build
  foreach my $hdold (@{$hdold_aref}) {
    unless ($seen{$hdold->name}) {push(@oldonly, $hdold->name)} # compare them to the objects in the clone in currentDB
  }
  if (@oldonly) {
    foreach my $hd (@oldonly) {
      $log->error("\tERROR: $hd missing $data_name data compared to $compared_to\n");
    }
  }
}

##################################################################

# check for missing data compared to currentdb, 
# that compares the objects to a simple list of things that should match in a regexp

sub check_for_missing_data2 {

  my ($hd_aref, $list_aref, $data_name, $compared_to) = @_;    

  my %seen=();
  my @notseen=();
  foreach my $hd (@{$hd_aref}) {$seen{$hd->name}=1} # get the name of objects in the clone in the Build
  foreach my $expected (@{$list_aref}) {
    my $found=0;
    foreach my $seen (keys %seen) {
      if ($seen =~ /$expected/) {
	$found=1;
	last;
      }
    }
    if (!$found) {
      push @notseen, $expected;
    }
  }
  if (@notseen) {
    foreach my $hd (@notseen) {
      $log->error("\tERROR: missing $data_name '$hd' compared to $compared_to\n");
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

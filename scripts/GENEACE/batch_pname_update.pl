#!/software/bin/perl -w
#
# batch_pname_update.pl
# 
# This script has been written to automatically change the Public_name
# of all NameServerIDs specified in a text file.
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2015-03-31 10:16:20 $
#



use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use lib '/software/worm/lib/perl';
use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;

=pod

=head batch_pname_update.pl

=item Options:

  -file      file containing IDs and old/new names  <Mandatory>

    FORMAT: standard
    WBGene00012345 AH6.4 Sequence
    WBGene00012345 bla-4 CGC

    FORMAT: -newonly
    WBGene00012345 unc-20
    WBGene00054321 AH6.4

  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans
  -newonly   modifies what is necessary in the input file - it assumes it is only CGC names

e.g. perl pname_update.pl -species elegans -file variation_name_data

=cut




#connect to name server
my ($help, $debug, $verbose, $store, $wormbase, $species);
my ($file, $cgc, $sequence, $newonly);

GetOptions (
	    "file=s"     => \$file,
	    "debug=s"    => \$debug,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "cgc"        => \$cgc,
	    "sequence"   => \$sequence,
            "newonly"    => \$newonly,
            "seqforce"   => \$seqforce,
	   );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -organism => $species
			     );
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

# establish species
unless (defined $species){
  $species = "elegans";
  $log->write_to ("Species has been set to elegans as it was not specified on command line\n\n");
}

#ENFORCE command line options
die "-file option is mandatory\n" unless $file;
if (!$newonly) {
  die "-cgc is only available with -newonly" if ($cgc);
  die "-sequence is only available with -newonly" if ($sequence);
}

$db = NameDB_handler->new($wormbase);

open (DATA, "<$file");

while (<DATA>) {
  chomp;
  # if you only have the new names
  if ($newonly) {
      if ($_ =~/(WB\S+\d{8})\s+(\S+)/)	{
	  $log->write_to ("\nProcessing $_\n");
	  if (&check_name_exists($2)) {
	    $log->write_to ("\nERROR: $_ CGC name already exists\n");
            if ($seqforce) {
                    print "\nOnly adding Sequence name data\n";
                &addname($1, $2, $cgc, $public, $sequence)
            }
            
	  } else {
	    &addname($1, $2, 'CGC')
	  }
      }
      else {
	  $log->write_to ("\nERROR: $_ has a line format error\n");
      } 
      next;      
  }
  else {
      if ($_ =~/(WB\S+\d{8})\s+(\S+)\s+(\S+)/)	{
	$log->write_to ("\nProcessing $_\n");
	if (&check_name_exists($2)) {
	  $log->write_to ("\nERROR: $_ name already exists\n");
	} else {
	  &addname($1, $2, $3);
	}
      } else {
	$log->write_to ("\nERROR: $_ has a line format error\n");
      } 
      next;
  }
}

close DATA;
$log->mail;

sub check_name_exists {
  my $name = shift;
  my ($ID) = $db->idGetByAnyName($name);
  unless($ID) {  
    return undef;
  } else {
    $log->write_to ("Match: $name exists as the Public_name for $ID\n");
    return 1;
  }
}


sub addname {
  my $ID = shift;
  my $NEW = shift;
  my $TYPE = shift;

  if ($TYPE eq 'CGC') {
    $db->addName($ID, 'CGC', $NEW);
    $log->write_to ("$NEW added to $ID as CGC_name\n");

  } elsif ($TYPE eq 'Sequence') {
    $db->addName($ID,'Sequence', $NEW);
    $log->write_to ("$NEW added to $ID as Sequence_name\n");

  }
  
}

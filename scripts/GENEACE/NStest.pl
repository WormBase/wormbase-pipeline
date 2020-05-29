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

use lib '/software/worm/lib/perl';
use Wormbase;
use NameDB_handler;
use Storable;
use Getopt::Long;

=pod

=head batch_pname_update.pl

=item Options:

  -file      file containing IDs and old/new names  <Mandatory>

    FORMAT: standard
    WBVar00278423 ot582 Blah
    WBGene00012345 AH6.4 Blah

    FORMAT: -newonly
    WBGene00012345 unc-20
    WBGene00054321 AH6.4

  -domain    set domain to whatever (Gene Variation Feature)<Mandatory>
  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans
  -test      use the test nameserver  <Optional 4 testing>
  -user      username                 <Mandatory>
  -pass      password                 <Mandatory>
  -cgc       this also adds the name as a cgc name in the nameserver if used with the -cgc option
  -public    this only tries to add the name as a Public_name - used to update missing Public_names
  -sequence    this only tries to add the name as a Sequence_name - used to update missing Sequence_names
  -newonly   modifies what is necessary in the input file

e.g. perl pname_update.pl -user blah -pass blah -species elegans -test -file variation_name_data

=cut

my ($debug, $test, $help);

GetOptions (
            "debug=s"    => \$debug,
            "test"       => \$test,
	    "help"       => \$help,
#            "file=s"     => \$file,
#            "store:s"    => \$store,
#            "species:s"  => \$species,
#            "cgc"        => \$cgc,
#            "public"     => \$public,
#            "sequence"   => \$sequence,
#            "newonly"    => \$newonly,
           );


my $wormbase;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
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


$db = NameDB_handler->new($wormbase, $log);


if ($help) {
  $db->{'db'}->print_authentication_instructions();
} else {
  $db->test();
}


$db->close;

$log->mail();
print "Finished.\n";
exit(0);

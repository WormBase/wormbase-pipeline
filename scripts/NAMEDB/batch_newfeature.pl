#!/usr/bin/env perl
# create consecutive new features in the name server based on a number requested

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";

use lib $ENV{'CVS_DIR'};

use strict;
use NameDB_handler;
use Getopt::Long;
use Wormbase;

my($test,$request);
my ($help, $debug, $verbose, $store, $wormbase);

GetOptions (	"help"       => \$help,
		"debug=s"    => \$debug,
		"verbose"    => \$verbose,
		"store:s"    => \$store,
		"request:i"  => \$request, # how many features to create
           );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             );
}


die "give me a number of features that you want\n"  unless ($request =~ /^\d+$/);

my $first_id = undef;
my $last_id;

my $db = NameDB_handler->new($wormbase);

my ($new_ids, $batch) = $db->new_features($request);

print "Made $request new feature IDs\n";
print "First ID assigned:  ".$new_ids->[0]."\n";
print "Last ID assigned:   ".$new_ids->[-1]."\n";
print "Nameserver batch ID: $batch\n";

exit;

=head2  batch_newfeature.pl

Adds batches of feature ids to the nameserver

=over 4

batch_newfeature.pl -request 10

=end

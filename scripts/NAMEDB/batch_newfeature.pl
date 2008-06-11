use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans/';
use lib $ENV{'CVS_DIR'};

use strict;
use NameDB_handler;
use Getopt::Long;

my($test,$request,$USER);

GetOptions (	"test"       => \$test,
		"request:i"  => \$request,
		"user:s"     => \$USER
           );

die "give me a number of features that you want\n"  unless ($request =~ /^\d+$/);

#connect to name server and set domain to 'Feature'
my $DB    	= 'wbgene_id;mcs2a;3305';
$DB = 'test_'.$DB if $test;
my $DOMAIN  = 'Feature';
my $db = NameDB_handler->new($DB,$USER,$USER);
    $db->setDomain('Feature');

my $c = 0;
while( $c < $request ){
    my $id = $db->idCreate;
    $c++;
}

print "made $request new feature IDs\n";

exit;

=head2  batch_newfeature.pl

Adds batches of feature ids to the nameserver

=over 4

batch_newfeature.pl -user ar2 -request 10

=end

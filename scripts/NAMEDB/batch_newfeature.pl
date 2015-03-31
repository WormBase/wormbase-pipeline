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

my $first_id = undef;
my $last_id;

#connect to name server and set domain to 'Feature'
my $DB = 'wbgene_id;shap;3303';
$DB    = 'test_wbgene_id;utlt-db:3307' if $test;
my $DOMAIN  = 'Feature';
my $db = NameDB_handler->new($DB,$USER,$USER,"/nfs/WWWdev/SANGER_docs/data");
    $db->setDomain('Feature');

my $c = 0;
while( $c < $request ){
    my $id = $db->idCreate;
    print "$id\n";
    if (!defined $first_id) {$first_id = $id}
    $last_id = $id;
    $c++;
}

print "Made $request new feature IDs\n";
print "First ID assigned:  $first_id\n";
print "Last ID assigned:   $last_id\n";

exit;

=head2  batch_newfeature.pl

Adds batches of feature ids to the nameserver

=over 4

batch_newfeature.pl -user ar2 -request 10

=end

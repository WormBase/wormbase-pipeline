#!/software/bin/perl -w
                                  
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Log_files;
use Modules::CarryOver;

#my $store = '/nfs/disk100/wormpub/TEST_BUILD/autoace/Elegans.store';
#my $wormbase = retrieve( $store );
my $wormbase = Wormbase->new('-test' => 1,
							'-debug' => 'ar2');
$wormbase->debug('ar2');
my $log = Log_files->make_build_log($wormbase);

my $carry = CarryOver->new($wormbase,$log);
$carry->carry_wormpep();
$carry->carry_wormrna();

$log->write_to("done\n");
$log->mail;
exit;

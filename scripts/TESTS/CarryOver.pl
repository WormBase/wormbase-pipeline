#!/software/bin/perl -w
                                  
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Log_files;
use Modules::CarryOver;
use Getopt::Long;

my ($store, $old, $new);
my ($pep, $rna);
GetOptions (
		"old:i" => \$old,
		"new:i" => \$new,
		"store:s" => \$store,
		"pep"	  => \$pep,
		"rna"   => \$rna
	   );
	   
die("need store old(int) and new(int)\n") unless ($store and $old and $new);
my $wormbase = retrieve( $store );
$wormbase->debug('ar2');
my $log = Log_files->make_build_log($wormbase);

my $carry = CarryOver->new($wormbase,$log);
$carry->carry_wormpep($new, $old) if $pep;
$carry->carry_wormrna($new, $old) if $rna;

$log->write_to("done\n");
$log->mail;
exit;

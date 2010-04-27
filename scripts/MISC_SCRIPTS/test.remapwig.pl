use lib $ENV{'CVS_DIR'}."/Modules";
use Remap_wig;

my $mapper = Remap_wig->new('elegans');
my $file = '/nfs/team71/worm/ar2/Desktop/RNAPIIIP2_dyeswap_WS180.wig';
$mapper->remap_file ('wig'   => $file,
		     'genome' => 'I',
		     'out'   => $file."_remap",
		     'old'   => 140,
		     'new'   => 200,
		     );

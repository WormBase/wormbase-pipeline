#!/nfs/disk100/wormpub/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my $database;
my $dir;
my $tace;
my $recursive;
my ($test, $debug, $store);

GetOptions ('database:s' => \$database,
	    'dir:s'      => \$dir,
	    'recursive'  => \$recursive,
	    'debug:s'    => \$debug,
	    'test'       => \$test,
	    'store:s'    => \$store
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);
$tace = $wormbase->tace;
$database = $database ? $database : $wormbase->autoace;

$log->log_and_die("bad options") unless (-d $dir and -e $database);
&read_dir($dir);

my $command;
foreach my $file (@to_load) {
  next if( $file eq '.' or $file eq '..');
  if( substr($file,-3,3 ) eq "ace" ) {
    $command .= "pparse $file\n";
  }
}

$command .= "save\nquit\n";
print $command;
 
open (WRITEDB, "| $tace $database -tsuser reload |") or $log->log_and_die("Couldn't open pipe to $database\n");
print WRITEDB $command;
close WRITEDB;

$log->mail;
exit(0);


sub read_dir {
  my $dir = shift;
  opendir (DIR,$dir) or die "cant open that directory $dir\n";
  print "reading $dir\n";
  my @files = readdir DIR;
  foreach my $file ( @files ) {
    next if( $file eq '.' or $file eq '..');
    if (-d "$dir/$file"){
      &read_dir ("$dir/$file") if $recursive;
      print "passing $dir/$file for inclusion\n";
      next;
    }
    if( substr($file,-3,3 ) eq "ace" ) {
      push (@to_load, "$dir/$file");
      print "\t$file\n"
    }
  }
  close DIR;
}

__END__

=pod

=head2 NAME - load_files_in_dir.pl

=over 4

loads all *.ace files in the specified directory in to the specified database

can be recursive

=cut

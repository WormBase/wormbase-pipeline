#!/usr/bin/perl5.6.1 -w

use lib "/nfs/team71/worm/ar2/wormbase/scripts";
use Wormbase;

my $database = shift;
my $dir = shift;
my $tace = &tace;

opendir (DIR,$dir) or die "cant open that directory $dir\n";
my @files = readdir DIR;
close DIR;
my $command;

foreach my $file (@files) {
  if( substr($file,-3,3 ) eq "ace" ) {
    $command .= "pparse $dir/$file\n";
  }
}

$command .= "save\nquit\n";
 
open (WRITEDB, "| $tace $database |") || die "Couldn't open pipe to $database\n";
print WRITEDB $command;
close WRITEDB;

exit(0);

__END__

=pod

=head2 NAME - load_files_in_dir.pl

=over 4

loads all *.ace files in the specified directory in to the specified database

=cut

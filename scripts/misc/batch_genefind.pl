#!/usr/local/bin/perl -w

use strict;
use lib  $ENV{'CVS_DIR'};
use Getopt::Long;
use POSIX;
use File::Path;

my $config;       # text file to specify parameters
my $database;     # for coordinate info
my $chromosome;   # so that Coords converter can give clone details
my $window;
my $chunksize = 1000000;
my $out_dir;

GetOptions ( 
	    "database:s"    => \$database,
	    "chromosome:s"  => \$chromosome,
	    "out_dir:s"     => \$out_dir,
	    "config:s"      => \$config,
	    "window:s"      => \$window,
	    "chunk:s"       => \$chunksize
	   );

die unless ($chromosome and $config);
my @files = &split_files($config);

my $min = 0;
my $max = $min + $chunksize;
my $test_file = "CHUNK/$files[0]_${min}_${max}";
my $errdir = glob("~wormpub/BSUB_ERRORS");

while( -e $test_file ) {

  my $bsub = $ENV{'CVS_DIR'}."/misc/genefinder.pl -config $config -chromosome $chromosome -min $min -max $max";
  $bsub .= " -out_dir $out_dir" if $out_dir;
  my $error_file = "$errdir/$files[0]_${min}_${max}.err";
  print "$bsub\n";
  system("bsub -e $error_file $bsub");

  $min += $chunksize;
  $max += $chunksize;
  $test_file = "CHUNK/$files[0]_${min}_${max}";
}

exit(0);


sub split_files 
  {
    my $config = shift;
    open(CF,"<$config") or die "cant open $config :$!\n";
    my @all_files;

    my $parameter;
    my $ignore;
    while ( <CF> ) {
      next unless /\w/;
      next if (/\#/);
      chomp;
      s/\s//g;
      if ( />(\w+)/ ) {
	$ignore = 0;
      CASE:{
          ($1 eq "feature_gff") && do { last CASE; };
          ($1 eq "exclude_gff") && do { last CASE; };
	  $ignore = 1;
        }
      } else {
	push(@all_files,"$_") unless ($ignore == 1);
      }
    }

    my $split_file;
    foreach my $file ( @all_files ) {
      my @file_handles;
 
      open(FH, "<$file") or die "cant open $file\n";

      #create dir for split files
      mkpath("CHUNK",0777);

      while ( <FH> ) {
	my @data = split;
	my $file_index = floor( $data[3] / $chunksize);

	unless ( $file_handles[$file_index] ) {
	my $min = $file_index * $chunksize;
	my $max = $min + $chunksize;
	my $filename = "CHUNK/${file}_${min}_${max}";

	open ($file_handles[$file_index],">$filename") or die "cant open $filename: $!\n";
	
	}
	print {$file_handles[$file_index]} $_;
      }
      close FH;
    }

    return @all_files;
  }

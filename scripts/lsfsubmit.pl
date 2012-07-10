#!/usr/bin/env perl 

use Getopt::Long qw(:config pass_through);

my $bsubmem = 2;
my $queue = 'normal';

&GetOptions('m=s'   => \$bsubmem, 
            'q=s'   => \$queue,
    );

my $minus_M = $bsubmem * 1000 * 1000;
my $minus_R = sprintf("'select[mem>%d] rusage[mem=%d]'", $bsubmem * 1000, $bsubmem * 1000);

unshift @ARGV, ("bsub", "-q $queue", "-I", "-M $minus_M", , "-R $minus_R");

print "@ARGV\n";
my $exit_code = system(@ARGV);

if ($exit_code) {
  die "EXITED with code $exit_code\n";
}



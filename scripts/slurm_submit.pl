#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config pass_through);
use Modules::WormSlurm;
use Path::Class;

my $scratch = $ENV{BUILD_TMP};
my $mem = 2; # default memory allocation of 2Gb
my $time = '01:00:00'; # default time allocation of 1 hour
my $queue = 'production'; # use production queue as default
my $outfile;
my $errfile;
my $keep_out = 1;
my $keep_err = 1;

GetOptions(
    "mem|m=s"   => \$mem,
    "time|t=s"  => \$time,
    "queue|q=s" => \$queue,
    "outfile|o=s" => \$outfile,
    "errfile|e=s" => \$errfile,
    );

if (!defined $outfile) {
    $outfile = $scratch . '/' . $$ . '.out';
    $keep_out = 0;
}
if (!defined $errfile) {
    $errfile = $scratch . '/' .$$ . '.err';
    $keep_err = 0;
}


my $cmd = join(' ', @ARGV);
my $mem_mb = $mem * 1000;

my $job_id = WormSlurm::submit_job_and_wait($cmd, $queue, $mem_mb . 'm', $time, $outfile, $errfile);

if (-e $errfile) {
    my $err_fh = file($errfile)->openr;
    while (my $err_line = $err_fh->getline()) {
	print STDERR $err_line;
    }
    remove_hps_file($errfile) unless $keep_err;
}
if (-e $outfile) {
    my $out_fh = file($outfile)->openr;
    while (my $out_line = $out_fh->getline()) {
	print STDOUT $out_line;
    }
    remove_hps_file($outfile) unless $keep_out;
}

unless (defined $job_id && $job_id != 0) {
    print STDERR "ERROR: Job submission failed\n";
} else {
    my $exit_code = WormSlurm::get_exit_code($job_id);
    if (not defined $exit_code) {
	print STDERR "ERROR: Unable to determine exit code for submitted job (${job_id}) - investigate\n";
    } elsif ($exit_code) {
	print STDERR "ERROR: Job terminated abnormally with exit code ${exit_code}\n";
    } else {
	print STDERR "JOB SUCCESSFULLY COMPLETED\n";
    }
}

exit(0);


sub remove_hps_file {
    my $file = shift;

    # Need to remove file with sbatch as need to be on production queue to access /hps/nobackup
    WormSlurm::submit_job("rm $file", 'production', '200m', '0:00:10', '/dev/null', '/dev/null');

    return;
}

#!/usr/bin/env perl 

use Getopt::Long qw(:config pass_through);

my $bsubmem = 2;

&GetOptions('m=s'   => \$bsubmem, 
            'q=s'   => \$queue,
            'E=s'   => \$pre_exec,
	    'g'     => \$gpfs
    );

my $gpfsstr = '';
$gpfsstr = 'select[gpfs]' if ($gpfs); # http://www.ebi.ac.uk/systems-srv/public-wiki/index.php/EBI_Good_Computing_Guide#GPFS
my $minus_M = $bsubmem * 1000;
my $minus_R = sprintf("'select[mem>%d] rusage[mem=%d] %s'", $bsubmem * 1000, $bsubmem * 1000, $gpfsstr);

#unshift @ARGV, ("bsub", "-q $queue", "-I", "-M $minus_M", , "-R $minus_R");
my $command = "bsub ";
$command .=  "-q $queue " if defined $queue;
$command .= "-E '$pre_exec'" if defined $pre_exec;
$command .= "-I -M $minus_M -R $minus_R @ARGV";

my ($job_id) = open_command_line($command);

if (not $job_id) {
  print STDERR "ERROR: Could not get job_id from job, so unable to determine exit status - investigate\n";
} else {
  my $found = 0;
    
  while(not $found) {
    # pause; wait for LSF to register the job as complete
    sleep(10);
    open(my $blist, "bjobs -d $job_id |");
    while(<$blist>) {
      /^(\d+)/ and do {
        if ($1 eq $job_id) {
          $found = 1;
          last;
        }
      };
    }
  }
  
  my $error_code;

  #
  # bjobs -l sometimes takes a while to register the completion 
  #
  while(not defined $error_code) {
    open(my $bjobs, "bjobs -l $job_id |");
    while(<$bjobs>) {
      /Exited with exit code (\d+)/ and do {
        $error_code = $1;
        last;
      };
      /Done successfully/ and do {
        $error_code = 0;
        last;
      };
    }

    sleep(10);
  } 
  
  if (not defined $error_code) {
    print STDERR "ERROR: Unable to determine exit code for submitted job - investigate\n";
  } elsif ($error_code) {
    print STDERR "ERROR: job terminated abnormally (exit code $error_code)\n";
  } else {
    print STDERR "SUCCESSFULLY COMPLETE\n";
  }
}


sub open_command_line {
  my ( $bsub_command ) = @_;
  
  my $lsf = 0;
  
  if ( open( my $pipe, '-|' ) ) {
    while (<$pipe>) {
      if (/Job <(\d+)>/) {
        $lsf = $1;
      } 
      print;
    }
    
    return $lsf;
  } else {
    # We want STDERR and STDOUT merged for the bsub process
    # open STDERR, '>&STDOUT';
    # probably better to do with shell redirection as above can fail
    exec( $bsub_command . ' 2>&1' ) || die("Could not run bsub");
  }
} ## end sub open_command_line

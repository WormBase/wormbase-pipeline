=pod 

=head1 NAME

 Slurm

=head1 SYNOPSIS

 my $job_id = WormSlurm::submit_job($cmd, $queue, $memory, $time, $outfile, $errfile);

 my $job_id = WormSlurm::submit_job_and_wait($cmd, $queue, $memory, $time, $outfile, $errfile);

 my $exit_code = WormSlurm::get_exit_code($job_id);

 WormSlurm::wait_for_jobs(@job_ids);

 WormSlurm::submit_jobs_and_wait_for_all(\@cmds, \@queues, \@memory, \@times, \@outfiles, \@errfiles);

=head1 DESCRIPTION

 Provides functions for running an managing Slurm jobs

=head1 CONTACT

Mark Quinton-Tulloch  mqt@ebi.ac.uk


=head1 METHODS

=cut


package WormSlurm;


sub submit_job {
    my ($cmd, $queue, $memory, $time, $outfile, $errfile) = @_;

    return submit_job_with_name($cmd, $queue, $memory, $time, $outfile, $errfile, undef);
}

sub submit_job_with_name {
    my ($cmd, $queue, $memory, $time, $outfile, $errfile, $name) = @_;
    return submit_job_with_name_and_dependencies($cmd, $queue, $memory, $time, $outfile, $errfile, $name, undef);
}

sub submit_job_with_name_and_dependencies {
    my ($cmd, $queue, $memory, $time, $outfile, $errfile, $name, $dependencies) = @_;

    my $parameters = "--parsable --container=${queue} --mem=${memory} --time=${time} -e ${errfile} -o ${outfile}";
    if (defined $name) {
	$parameters .= " -J ${name}";
    }
    if (defined $dependencies) {
	my @dependency_strings;
	for my $dependency(@$dependencies) {
	    push @dependency_strings, 'afterany:' . $dependency;
	}
	$parameters .= ' -d ' . join(',', @dependency_strings);
    }

    my $sbatch = "sbatch $parameters --wrap=\"${cmd}\"";
    print STDERR $sbatch . "\n";
    my $job_id = `$sbatch`;
    chomp $job_id;
    print STDERR "Unexpected job ID $job_id\n" unless defined $job_id and $job_id != 0;
    return $job_id;
}

sub submit_job_and_wait {
    
    my ($cmd, $queue, $memory, $time, $outfile, $errfile) = @_;

    my $job_id = submit_job($cmd, $queue, $memory, $time, $outfile, $errfile);

    my $running = 1;
    system("sleep 1");
    while ($running) {
	$running = job_is_running($job_id);
	system("sleep 15");
    }

    return $job_id;
}

sub job_is_running {
    my $job_id = shift;

    my $status = `squeue --job $job_id --noheader | wc -l`;
    chomp $status;

    return $status;
}

sub wait_for_jobs {
    my @job_ids = @_;

    my $running = 1;
    system("sleep 1");
    while ($running) {
	my @still_running;
	for my $job_id (@job_ids) {
	    push @still_running, $job_id if job_is_running($job_id);
	}
	$running = scalar @still_running;
	@job_ids = @still_running;
	system("sleep 15");
    }

    return;
}

sub submit_jobs_and_wait_for_all {
    my ($cmds, $queue, $memory, $time, $outfile_prefix, $errfile_prefix) = @_;

    my %slurm_jobs;
    for (my $ix = 0; $ix++; $ix < @$cmds) {
	my $outfile = $outfile_prefix eq '/dev/null' ? $outfile_prefix : $outfile_prefix . '.' . $ix . '.out';
	my $errfile = $outfile_prefix eq '/dev/null' ? $errfile_prefix : $errfile_prefix . '.' . $ix . '.err';
	my $job_id = submit_job($cmds->[$ix], $queue, $memory, $outfile, $errfile);
	$slurm_jobs{$job_id} = $cmds->[$ix];
    }

    wait_for_jobs(keys $slurm_jobs);

    return \%slurm_jobs;
}

sub get_exit_code {
    my $job_id = shift;

    my $exit_code = 0;
    open(SACCT, "sacct -b -p -j $job_id |");
    while(<SACCT>) {
	chomp;
	next if $_ =~ /^JobID/;
	my @cols = split(/\|/, $_);
	if ($cols[1] eq 'FAILED') {
	    ($exit_code) = $cols[2] =~ /^(\d+):\d+$/;
	    last;
	}
	if ($cols[1] eq 'TIMEOUT' || $cols[1] eq 'OUT_OF_MEMORY') {
	    $exit_code = 1;
	    last;
	}
    }
    close(SACCT);

    return $exit_code;
}

1;

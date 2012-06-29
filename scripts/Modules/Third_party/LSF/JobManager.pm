package LSF::JobManager; $VERSION = 0.5;

use strict;
use warnings;
use base qw( LSF );
use IPC::Run qw( start pump finish );
use LSF::Job;

sub import{
    my ($self, %p) = @_;
    $p{RaiseError}  ||= 1;
    $p{PrintOutput} ||= 1;
    $p{PrintError}  ||= 1;
    $self->PrintOutput($p{PrintOutput}) if exists $p{PrintOutput};
    $self->PrintError ($p{PrintError} ) if exists $p{PrintError};
    $self->RaiseError ($p{RaiseError} ) if exists $p{RaiseError};
}

# create a new manager object. Store default parameters that will
# be used each time submit is called
sub new{
    my($type,@params) = @_;
    my $class = ref($type) || $type || "LSF::JobManager";
    return bless {-params => [@params] }, $class;
}

# submit the command line to LSF. Any new parameters are used in 
# addition to those parameters passed to new and override them.
sub submit{
    my $this = shift;
    my($cmd,@params);
    if( @_ == 1 ){
        $cmd = shift;
    }else{
        $cmd = pop;
        @params = @_;
    }
    my @new_params = new_flags( [ $this->params ], \@params );
    my $job = LSF::Job->submit(@new_params, $cmd);
    $this->{-jobs}->{"$job"} = $job if $job;
    return $job;
}

# wait for all the submitted LSF jobs in a blocking manner
# achieved by submitting a job with the -I flag (stay connected to 
# the terminal) and a dependancy expression that tests that all
# the previously submitted LSF jobs have ended (any exit status)
sub wait_all_children{
    my ($this,%options) = @_;
    my @jobs = $this->jobs;
    unless(@jobs){
        warn "No LSF::Job's in this LSF::Manager\n" if $this->PrintError;
        return;
    }
    my $dependancy;
    if(exists $options{depend}){
        $dependancy = $options{depend};
    }else{
        for(@jobs){
            $dependancy .= "ended($_)&&";
        }
        $dependancy =~ s/\&\&$//;
    }
    
    my $job = LSF::Job->submit_top( new_flags( [ $this->params ],['-I',-w => $dependancy])
                                              ,"echo LSF::JobManager finished waiting for " . scalar @jobs . " jobs" );
    
	$options{history} = 1 unless exists $options{history};
    $this->cache_history if $options{history};
    return;
}

sub cache_history{
    # batch run the history of each job b/c its much quicker this way.
    # because we may have more jobs than the arg list allows then we have a 
    # second syntax for calling LSF::JobHistory->new
    # if the first argument is and array reference then 
    # see JobHistory.pm for details but basically we use xargs
	my $this = shift;
    my @jobs = $this->jobs;
    my @history = LSF::JobHistory->new([],@jobs);
    my %jobs = %{ $this->{-jobs} };
    for my $hist (@history){
        $jobs{$hist->id}->{-cached_history} = $hist;
    }
    return;
}

# return an array of the parameters submitted to new
sub params{
    my $this = shift;
    $this->{-params} = [@_] if @_;
    return @{ $this->{-params} };
}

# return an array of the lsf jobs submitted to this point.
sub jobs{
    my $this = shift;
    my $jobs = $this->{-jobs};
    return values %$jobs;
}

# clear a job manager object of previously submitted jobs
sub clear{
    my $this = shift;
    my $jobs = delete $this->{-jobs};
    return values %$jobs;
}

# internal sub to parse a facile command line of flags 
# for name=value pairs
sub parse_flags{
    my @defaults = @_;
    my %hash;
    while(local $_ = shift @defaults){
        if(/^(-\w)(.*)/){
            my($flag,$value) = ($1,$2);
            if($value ne ''){
                $hash{$flag} = $value;
            }elsif($defaults[0] !~ /^-\w/){
                $hash{$flag} = shift @defaults;
            }else{
                $hash{$flag} = undef;
            }
        }
    }
    return ( %hash );
}

# internal routine to allow new flags to override old ones
sub new_flags{
    my @defaults = @{ $_[0] };
    my @new = @{ $_[1] };
    my %defaults = parse_flags(@defaults);
    my %new = parse_flags(@new);
    my %hash = ( %defaults, %new );
    my @returned;
    while( my($key,$val) = each %hash ){
        push @returned, $key;
        push @returned, $val if defined $val;
    }
    return @returned;
}

1;

__END__


=head1 NAME

LSF::JobManager - submit and wait for a set of LSF Jobs

=head1 SYNOPSIS

    use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
    use LSF::JobManager;

    my $m = LSF::JobManager->new(-q=>'small');

    my $job = $m->submit("echo hello");
    $m->submit("echo world");
    
    for my $job ($m->jobs){
        $job->top;
    }
    
    $m->wait_all_children( history => 1 );
    print "All children have completed!\n";

    for my $job ($m->jobs){ # much quicker if history is pre-cached
        print STDERR "$job exited non zero\n" if $job->history->exit_status != 0;
    }
    
    $m->clear; # clear out the job manager to reuse.

=head1 DESCRIPTION

C<LSF::JobManager> provides a simple mechanism to submit a set of command lines
to the LSF Batch system and then wait for them all to finish in a blocking
(efficient) manner. 

=head1 INHERITS FROM

B<LSF>

=head1 CONSTRUCTOR

=over 4

=item new ( [ ARGS ] )

$manager = LSF::JobManager->new(-q=>'small'
                               ,-m=>'mymachine');

Creates a new C<LSF::JobManager> object.

Any parameters are used as defaults passed to the submit method.

=back

=head1 METHODS

=over

=item $manager->submit( [ [ ARGS ] ], [CMD] )

Submits a command line to LSF. This is a wrapper around the LSF::Job->submit
call. The required argument is the command line to submit. Optional arguments
override the defaults given to C<new>. The submitted LSF::Job object is returned
on success, otherwise undef is returned, $? and $@ set. See C<LSF::Job>

=item $manager->wait_all_children( [ depend => DEPENDANCY, history => 1 ] )

Waits for all previously submitted LSF Jobs to complete in a blocking manner.

If a dependancy expression is provided then this is used in preference to 
an autogenerated dependancy based on the job id's of all the submitted jobs.
In conjunction with named jobs ( the -J flag ) this can be used to get around
the shell arg length limit which can be hit by extremely long dependancy 
expressions. See the L<bsub> documentation for details of job names and 
depenancy expressions.

If the history flag is used to pass a true value then after all jobs have 
completed the history of each job is fetched in batch. It is much quicker
this way than fetching each job history separately. The time taken is still 
significant however. If you do not intend to check the exit status or the
command line or any other property of the Job::History object then it is worth
passing a false value and avoiding the overhead.

=item $manager->params( [ [ PARAMS ] ] )

Sets the default submission parameters.

=item $manager->jobs()

Returns an array of the submitted LSF::Job objects.

=item $manager->clear()

Empties the job array so that the manager can be reused. 
Returns an array of the jobs removed.

=head1 HISTORY

The B<LSF::Batch> module on cpan didn't compile easily on all platforms i wanted.
The LSF API didn't seem very perlish either. As a quick fix I knocked these
modules together which wrap the LSF command line interface. It was enough for
my simple usage. Hopefully they work in a much more perly manner.

=head1 SEE ALSO

L<LSF>,
L<LSF::Job>

=head1 AUTHOR

Mark Southern (mark_southern@merck.com)

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)

=cut

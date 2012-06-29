package LSF::Job; $VERSION = 0.7;

use strict;
use warnings;
use base qw( LSF );
# sugar so that we can use job id's in strings
use overload '""' => sub{ $_[0]->{-id} };
use LSF::JobHistory;
use IPC::Run qw( start pump finish );

sub import{
    my ($self, %p) = @_;
    $p{RaiseError}  ||= 1;
    $p{PrintOutput} ||= 1;
    $p{PrintError}  ||= 1;
    $self->PrintOutput($p{PrintOutput}) if exists $p{PrintOutput};
    $self->PrintError ($p{PrintError} ) if exists $p{PrintError};
    $self->RaiseError ($p{RaiseError} ) if exists $p{RaiseError};
}

sub new{
    my($type, $id) = @_;
    my $class = ref($type) || $type || 'LSF::Job';
    unless( $id =~ /^\d+$/ ){
        my $msg = "Invalid Job <$id>";
        die $msg if $class->RaiseError;
        warn $msg, "\n" if $class->PrintError;
        return;
    }
    return bless {-id => $id}, $class;
}

sub jobs{
    my($self,@params) = @_;
    my($out,$err);
    IPC::Run::run ['bjobs',@params],\undef,\$out,\$err;
    print $out if $out && $self->PrintOutput;
    if($?){
        $@ = $err;
        if( $err =~ /No unfinished job found/i or $err =~ /is not found/i ){
            warn $err if $err && $self->PrintError;
            return wantarray ? () : 0;
        }else{
            die $err if $self->RaiseError;
        }
    }
    print $err if $err && $self->PrintError;
    my @rows = split(/\n/,$out);
    shift @rows; # first row contains column headings
    if( wantarray ){
        my @return;
        for (@rows){
            /^(\d+)/ && push @return, $self->new($1);
        }
        return @return;
    }
    return scalar @rows;
}

sub submit{
    my ($self,@params) = @_;
    my @output;
    @output = $self->do_it('bsub',@params);
    return unless @output;
    my $idx = 0;
    $idx = 1 if $self->LSF =~ /^4/;
    $output[$idx] =~ /Job <(\d+)>/;
    return $self->new($1);
}

sub submit_top{
	my $self = shift;
    $self->submit_pos(1,@_);
}

sub submit_bottom{
	my $self = shift;
    $self->submit_pos(0,@_);
}

sub submit_pos{
    my ($self,$pos,@params) = @_;
    my ($job,@output);
    my $h = start( ['bsub',@params], \undef, \$output[0], \$output[1] );
    my $idx = $self->LSF =~ /^4/;
    pump( $h );
    if( $output[$idx] =~ /Job <(\d+)>/ ){
        $job = $self->new($1);
        # do post processing, assuming that the command succeeded
        # because we got the jobid back
        $self->post_process(0,@output);
        # then reset the output and error
        $output[0] = $output[1] =  '';
        my $re = $self->RaiseError(0);
        if($pos){
	        $job->top; # would throw an error if its already running
        }else{
        	$job->bottom;
        }
        $self->RaiseError($re);
    }
    finish($h);
    #  print out the rest of the output and error
    $self->post_process($?,@output);
    return $job;
}

sub id         { "$_[0]" }

sub bottom     { my $self = shift; $self->do_it('bbot',     @_, "$self") }

sub checkpoint { my $self = shift; $self->do_it('bchkpnt',  @_, "$self") }

sub delete     { my $self = shift; $self->do_it('bdel',     @_, "$self") }

sub kill       { my $self = shift; $self->do_it('bkill',    @_, "$self") }

sub migrate    { my $self = shift; $self->do_it('bmig',     @_, "$self") }

sub modify     { my $self = shift; $self->do_it('bmod',     @_, "$self") }

sub peek       { my $self = shift; $self->do_it('bpeek',    @_, "$self") }

sub requeue    { my $self = shift; $self->do_it('brequeue', @_, "$self") }

sub restart    { my $self = shift; $self->do_it('brestart', @_, "$self") }

sub resume     { my $self = shift; $self->do_it('bresume',  @_, "$self") }

sub run        { my $self = shift; $self->do_it('brun',     @_, "$self") }

sub stop       { my $self = shift; $self->do_it('bstop',    @_, "$self") }

sub switch     { my $self = shift; $self->do_it('bswitch',  @_, "$self") }

sub top        { my $self = shift; $self->do_it('btop',     @_, "$self") }


sub history{ 
    my ($self) = @_;
    return $self->{-cached_history} if $self->{-cached_history};
    my ($hist) = LSF::JobHistory->new("$self"); 
    if( defined $hist->exit_status ){
        $self->{-cached_history} = $hist;
    }
    return $hist;
}

1;

__END__

=head1 NAME

LSF::Job - create and manipulate LSF jobs

=head1 SYNOPSIS

use LSF::Job;

use LSF::Job RaiseError => 0, PrintError => 1, PrintOutput => 0;

$job = LSF::Job->new(123456);

...

$job = LSF::Job->submit(-q => 'default'
                       ,-o => '/dev/null'
                       ,"echo hello");

$job2 = LSF::Job->submit(-q => 'default'
                        ,-o => '/home/logs/output.txt'
                        ,"echo world!");

@jobs = LSF::Job->jobs( -J => "/mygroup/*" );

$job2->modify(-w => "done($job)" );

$job2->del(-n => 1);

...

$job->top();

$job->bottom();

$exit = $job->history->exit_status;

... etc ...

=head1 DESCRIPTION

C<LSF::Job> is a wrapper arround the LSF b* commands used to submit and
manipulate jobs. for a description of how the LSF commands work see the 
man pages of:

    bsub bbot bchkpnt bkill bmig bmod brequeue brestart bresume bstop bswitch btop

=head1 INHERITS FROM

B<LSF>

=head1 CONSTRUCTORS

=over 4

=item new ( [NUM] )

$job = LSF::Job->new(123456);

Creates a new C<LSF::Job> object.

Required argument is a LSF jobid. This does not *have* to exist in the system
but would probably be a good idea!

=item submit ( [ [ARGS] ], [COMMAND_STRING] )

$job = LSF::Job->submit(-q => 'default'
                       ,-o => '/dev/null'
                       ,"echo hello");

Creates a new C<LSF::Job> object.

Arguments are the LSF parameters normally passed to 'bsub'.

Required parameter is the command line (as a string) that you want to execute.

=item jobs ( [ARGS] )

@jobs = LSF::Job->jobs( -J => "/mygroup/*" );

Creates an array of LSF::Job objects corresponding to jobs that match the query

Arguments are the LSF parameters normally passed to 'bjobs'.

=back

=head1 METHODS

=over

=item $job->id() (or object in string context)

C<id> returns the jobid of the LSF Job. This is particularly useful when
building up job interdependencies

=item $job->history

Returns a LSF::JobHistory object with information about the LSF job. 
See the LSF::JobHistory perldoc page.

=item $job->bottom

Moves the LSF job to the bottom of its queue. See the bbot man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->checkpoint

Checkpoints a checkpointable job. See the bchkpnt man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->delete( [ARGS] )

*** Deprecated in LSF 5.0. Use kill instead. ***

Deletes the LSF job from the system. See the bdel man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->kill

Kills the LSF job. See the bkill man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->migrate

Migrates the LSF job. See the bmigrate man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->modify( [ARGS] )

Modifies the LSF job. See the bmod man page.
Since the objects are overloaded to return the job id when used in string 
context this allows easy build up of job dependancies e.g.
Returns true on success, false on failure. Sets $? and $@;

$job3->modify(-w => "done($job1) && done($job2)" );

=item $job->restart

Restarts a checkpointed job. See the brestart man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->resume

Resumes a suspended job. See the bresume man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->run

Starts the LSF job now. See the brun man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->stop

Stops the LSF job. See the bstop man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->switch( [ARGS] )

Switches the LSF job between LSF queues. See the bswitch man page.
Returns true on success, false on failure. Sets $? and $@;

=item $job->top

Moves the LSF job to the top of its queue. See the btop man page.
Returns true on success, false on failure. Sets $? and $@;

=back

=head1 BUGS

The use of the '-l' flag of the LSF command lines can be considered a bug.
Using LSF job names with non alphabetic characters can also be considered a bug.
Otherwise, please report them.

=head1 HISTORY

The B<LSF::Batch> module on cpan didn't compile easily on all platforms i wanted.
The LSF API didn't seem very perlish either. As a quick fix I knocked these
modules together which wrap the LSF command line interface. It was enough for
my simple usage. Hopefully they work in a much more perly manner.

=head1 SEE ALSO

L<LSF>,
L<LSF::JobHistory>,
L<bsub>,
L<bhist>,
L<bswitch>,
L<bdel>,
L<bkill>,
L<bstop>,
L<bmod>,
L<btop>,
L<bbot>,
L<brun>

=head1 AUTHOR

Mark Southern (mark_southern@merck.com)

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)

=cut

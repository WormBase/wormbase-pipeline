package LSF::JobHistory; $VERSION = 0.3;

use strict;
use warnings;
use base qw( LSF );
use IPC::Run qw( run );
use Date::Manip;

sub import{
    my ($self, %p) = @_;
    $p{RaiseError}  ||= 1;
    $p{PrintOutput} ||= 1;
    $p{PrintError}  ||= 1;
    $self->PrintOutput($p{PrintOutput}) if exists $p{PrintOutput};
    $self->PrintError ($p{PrintError} ) if exists $p{PrintError};
    $self->RaiseError ($p{RaiseError} ) if exists $p{RaiseError};
}

sub AUTOLOAD{
    my $key = our $AUTOLOAD;
    $key =~ s/.*:://;
    my $self = shift;
    $self->{$key} = shift if @_;
    return $self->{$key};
}

sub new{
    my($self,@params) = @_;
    my $class = ref($self) || $self || __PACKAGE__;
    my @output;
    if( ref $params[0] eq 'ARRAY' ){
        my $extra = shift @params;
        my $in;
        $in .= "$_\n" for @params;
        my ($out,$err);
        IPC::Run::run ['xargs','bhist','-l', @$extra ], \$in,\$output[0],\$output[1];
        $self->post_process($?,@output);
    }else{
        my @output = $class->do_it('bhist','-n0','-l',@params);
    }
    return unless @output;
    my @jobhistory;
    my @histories = split(/^-+$/m,$output[0]);
    for my $job (@histories){ # each entry divided by ------ line
        $job =~ s/Summary of time .+//sg;
        $job =~ s/\n +//g;
        my $this = bless {}, $class;
        if( $job =~ /^Job <(\d+)>/m ){
            $this->id("$1");
        }
        if( $job =~ /^.+: Done successfully\./m ){
            $this->exit_status(0);
        }
        elsif( $job =~ /Exited with exit code (\d+)/ ){
            $this->exit_status("$1");
        }
        if( $job =~ /Execution Pid <(\d+)>/ ){
            $this->pid("$1");
        }
        if( $job =~ /Command <(.+)>\n/ ){
            $this->command("$1");
        }
        if( $job =~ /CWD <(.+?)>[,\s]/ ){
            ( my $cwd = $1 ) =~ s/\$(\w+)/$ENV{$1}/g;
            $this->cwd($cwd);
        }
        if( $job =~ /^(.+): Starting /m ){
            $this->started( UnixDate( "$1", "%d-%b-%y %T" ) );
        }
        if( $job =~ /^(.+): (Done|Exited)/m ){
            $this->ended( UnixDate( "$1", "%d-%b-%y %T" ) );
        }
        if( $job =~ /Dispatched to <(\w+)>;/ ){
            $this->host("$1");
        }
        push @jobhistory,$this;
    }
    return @jobhistory;
}

sub process{

}

1;

__END__

=head1 NAME

LSF::JobHistory - get historical information about LSF jobs.

=head1 SYNOPSIS

use LSF::JobHistory;

use LSF::JobHistory RaiseError => 0, PrintError => 1, PrintOutput => 0;

( $jobhistory ) = LSF::JobHistory->new( [ARGS] );

( $jobhistory ) = LSF::JobHistory->new( $job );

( $jobhistory ) = LSF::JobHistory->new( [JOBID] );

@jobhistory = LSF::JobHistory->new( -J => '/MyJobGroup/*');

( $jobhistory ) = LSF::JobHistory->new($job);

$jobhistory = $job->history;

... etc ...

$exit_status = $jobhistory->exit_status;

$pid = $jobhistory->pid;

$command = $jobhistory->command;

$cwd = $jobhistory->cwd;


=head1 DESCRIPTION

C<LSF::JobHistory> is a wrapper arround the LSF 'bhist' command used to obtain 
historical information about jobs. See the 'bhist' man page for more
information. This provides a more reliable way to obtain the exit status of
an LSF job than from the LSF::JobInfo object because the bhist command can
search all of the available LSF logs to find the information.

=head1 INHERITS FROM

B<LSF>

=head1 CONSTRUCTOR

=over 4

=item new( [ARGS] || [JOBID] || $job );

($jobhistory) = LSF::JobHistory->new(  [ARGS]
                                    || [JOBID]
                                    );

Creates a new C<LSF::JobHistory> object.

Arguments are the LSF parameters normally passed to 'bhist' or
a valid LSF jobid or LSF::Job object. The bhist command is automatically called 
with the -n 0 and -l flags. 

Returns an array of LSF::JobHistory objects. Of course if your argument to new is
a single jobid then you will get an array with one item. If you query for a 
number of jobs with the same name or path then you will get a list.
In scalar context returns the number of jobs that match that criteria.

=head1 BUGS

Please report them.
Otherwise... the parsing of the LSF output can fail if the job names have 
non-alphanumeric characters in them. You probably shouldn't do this anyway.

=head1 HISTORY

The B<LSF::Batch> module on cpan didn't compile easily on all platforms i wanted.
The LSF API didn't seem very perlish either. As a quick fix I knocked these
modules together which wrap the LSF command line interface. It was enough for
my simple usage. Hopefully they work in a much more perly manner.

=head1 SEE ALSO

L<LSF>,
L<LSF::Job>,
L<bjobs>

=head1 AUTHOR

Mark Southern (mark_southern@merck.com)

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)

=cut

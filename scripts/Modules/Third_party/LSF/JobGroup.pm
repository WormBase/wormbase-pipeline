package LSF::JobGroup; $VERSION = 0.3;

use strict;
use warnings;
use base qw( LSF );
use overload '""' => sub{ $_[0]->{-name} };
use LSF::Job;
use IPC::Run qw( run );

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
    my($type,$name) = @_;
    my $class = ref($type) || $type || 'LSF::JobGroup';
    unless( $name ){
        my $msg = "Invalid group name <$name>";
        die $msg if $class->RaiseError;
        warn $msg, "\n" if $class->PrintError;
        return;
    }
    return bless { -name => $name }, $class;
}

sub name { "$_[0]" }

sub exists{ 
    my($self) = @_;
    my $job = eval{ LSF::Job->submit(-J => "$self/", "echo dummy") };
    return 0 unless $job;
    $job->delete;
    return 1;
}

sub delete  { my $self = shift; $self->do_it('bgdel', @_,"$self") }
sub add     { my $self = shift; $self->do_it('bgadd', @_,"$self") }
sub hold    { my $self = shift; $self->do_it('bghold',@_,"$self") }
sub release { my $self = shift; $self->do_it('bgrel', @_,"$self") }

sub modify{
    my($self) = shift;
    my $newname;
    my $flag;
    for(my $i = 0; $i < @_; $i++){
        if( $_[$i] =~ /^-G(.*)/ ){
            if($1){
                $newname = $1;
            }else{
                $newname = $_[$i+1];
            }
        }
    }
    my $retval = $self->do_it('bgmod',@_,"$self");
    $self->{-name} = $newname if $newname && $retval;
    return $retval;
}

sub jobs{
    my $self = shift;
    return LSF::Job->jobs(-J => "$self/", @_);
}

1;

__END__

=head1 NAME

LSF::JobGroup - manipulate LSF job groups

=head1 SYNOPSIS

use LSF::JobGroup;

use LSF::JobGroup RaiseError => 0, PrintError => 1, PrintOutput => 0;

$jobgroup = LSF::JobGroup->new( [GROUP_NAME] );

...

$jobgroup->add( [ARGS] ) unless $jobgroup->exists;

$jobgroup->delete;

$jobgroup->hold;

$jobgroup->release;
...
$jobgroup->modify(-w => 'exited(/mygroup/,==0)' );
...
@jobs = $jobgroup->jobs('-r');
... etc ...

=head1 DESCRIPTION

C<LSF::JobGroup> is a wrapper arround the LSF b* commands used to manipulate job
groups. for a description of how the LSF commands work see the man pages of:

    bgadd bgdel bghold bgrel bgmod bjobs

=head1 INHERITS FROM

B<LSF>

=head1 CONSTRUCTOR

=over 4

=item new ( [NUM] )

$jobgroup = LSF::JobGroup->new('/MyGroup');

Creates a new C<LSF::JobGroup> object.

Required argument is a job group name. This can be a single group name or a
path, much like a filesystem path. This does not *have* to exist in the system
as new job groups can be created. Names should only contain alphanumeric 
characters plus '_' and '-'. Not only my code but also LSF job dependancy
expressions will fail if you attempt otherwise.

=back

=head1 METHODS

=over

=item $jobgroup->exists

C<id> returns 1 if the job group exists, 0 otherwise. The method attempts to
create the group and if it fails it examines the LSF output to see if the group
existed. I couldn't find a better test to use. Answers on a postcard...

=item $jobgroup->add

Adds a job group, or group path.
Returns true on success, false on failure. Sets $? and $@;

=item $jobgroup->delete

Deletes a job group.
Returns true on success, false on failure. Sets $? and $@;

=item $jobgroup->hold

Holds a LSF job group. All pending jobs will wait until the group is released.
Returns true on success, false on failure. Sets $? and $@;

=item $jobgroup->release

Releases a LSF job group. Pending jobs are free to run.
Returns true on success, false on failure. Sets $? and $@;

=item $job->modify([ARGS])

Modifies the LSF job group. For example, changing its name or its dependancy
expression. See the L<bgmod> man page.

$jobgroup->modify(-w => "started($jobgroup,==0)" );
$jobgroup->modify(-w => "ended($jobgroup,==0)" );

=item $jobgroup->jobs([[ARGS]])

Returns an list of LSF::Job objects of jobs contained within this job group.
Remember to use the '-r' flag if you want to include jobs in sub groups.

=back

=head1 BUGS

The use of the '-l' flag of the LSF command lines can be considered a bug.
Using group names with non alphabetic characters can also be considered a bug.
Otherwise please report them.

=head1 HISTORY

The B<LSF::Batch> module on cpan didn't compile easily on all platforms i wanted.
The LSF API didn't seem very perlish either. As a quick fix I knocked these
modules together which wrap the LSF command line interface. It was enough for
my simple usage. Hopefully they work in a much more perly manner.

=head1 SEE ALSO

L<LSF>,
L<LSF::Job>,
L<bgadd>,
L<bgdel>,
L<bghold>,
L<bgrel>,
L<bgmod>,
L<bjobs>

=head1 AUTHOR

Mark Southern (mark_southern@merck.com)

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)

=cut

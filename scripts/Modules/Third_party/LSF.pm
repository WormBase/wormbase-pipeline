package LSF; 

$VERSION = 0.9;

use strict;
use warnings;
use IPC::Run qw( run );

sub BEGIN {
    my ($out,$err);
    eval{ 
        run ['lsid','-V'],\undef,\$out,\$err; 
    };
    if($@){
        die "Cannot find LSF executables. \$PATH is $ENV{PATH}\n";
    }else{
        die "Cannot determine LSF version\n" unless $err =~ /^LSF ([^,]+)/m;
        $__PACKAGE__::LSF = $1;
    }
}

sub LSF { $__PACKAGE__::LSF }

# sugar to preload all the LSF modules
sub import {
    shift;
    my %p = @_;
    $p{RaiseError}  ||= 1;
    $p{PrintOutput} ||= 1;
    $p{PrintError}  ||= 1;
    my @modules = qw(Job JobHistory JobGroup Queue JobManager);
    my @code = map { "use LSF::$_ PrintOutput => $p{PrintOutput}, PrintError => $p{PrintError}, RaiseError => $p{RaiseError};\n"; 
                   } @modules;
    eval join('', @code);
    die $@ if $@;
}

# Class method to control whether we die on errors
sub RaiseError { shift->static('RaiseError',@_) }

# Class methods to control whether LSF output/error is printed.
sub PrintOutput { shift->static('PrintOutput',@_) }

sub PrintError { shift->static('PrintError',@_) }

sub static{
    my ($self,$var) = (shift,shift);
    my $class = ref($self) || $self;
    my $varname = "${class}::${var}";
    no strict 'refs';
    my $retval = $$varname;
    $$varname = shift if @_;
    return $retval;    
}

sub do_it{
    my($self,@params) = @_;
    my ($out,$err);
    run \@params,\undef,\$out,\$err;
    return $self->post_process($?,$out,$err);
}

sub post_process{
    my($self,$exit,$out,$err) = @_;
    print $out if $out && $self->PrintOutput;
    if($exit){
        die $err if $self->RaiseError;
        $@ = $err;
        warn $err if $err && $self->PrintError;
        return undef;
    }
    warn $err if $err && $self->PrintError;
    return ($out,$err);
}

1;

__END__

=head1 NAME

LSF - A perl API built on top of the LSF command line tools

=head1 SYNOPSIS

    use LSF;
    use LSF RaiseError => 1, PrintError => 1, PrintOutput => 1;

=head1 DESCRIPTION

NOTE: FOR THESE MODULES TO WORK IT IS ESSENTIAL THAT YOU INCLUDE THE LSF 
COMMAND LINES IN YOUR PATH.

This is the base class of the LSF suite of modules. 'use LSF' will also 
preload all of the LSF modules at one time. Currently this includes:

      LSF::Job
      LSF::JobHistory
      LSF::JobGroup
      LSF::Queue
      LSF::JobManager

Two error reporting strategies are available and can be set globally via the
'use LSF' statement or individually in each of the LSF modules. By setting the
'RaiseError' directive to true, or by using the RaiseError class method, the 
LSF modules will die on error, otherwise they will return false, setting $? to 
the exit value and $@ to the stderr of the LSF command line. Additionally
the printing of LSF command line stdout and stderr can be controlled via the 
'PrintOutput' and 'PrintError' directives or class methods of the same names.
Defaults are as above.

For more information on any of these modules, please see its respective
documentation.

=head1 CLASS METHODS

=over

=item LSF()

Returns the LSF version string

=item RaiseError( [ [ TRUE or FALSE ] ] )

Controls whether LSF command line errors will be thrown. The default is
FALSE. When called with no arguments returns the current value.

=item PrintError( [ [ TRUE or FALSE ] ] )

Controls printing to STDERR the stderr of the LSF command line. The default is
TRUE. When called with no arguments returns the current value.

=item PrintOutput( [ [ TRUE or FALSE ] ] )

Controls printing to STDOUT the stdout of the LSF command line. The default is
FALSE. When called with no arguments returns the current value.

=back

=head1 HISTORY

The B<LSF::Batch> module on cpan didn't compile easily on all platforms i wanted.
The LSF API didn't seem very perlish either. As a quick fix I knocked these
modules together which wrap the LSF command line interface. It was enough for
my simple usage. Hopefully they work in a much more perly manner.

=head1 SEE ALSO

L<LSF::Batch>,
L<LSF::Job>,
L<LSF::JobHistory>,
L<LSF::JobManager>,
L<LSF::JobGroup>,
L<LSF::Queue>,
L<bsub>,
L<bhist>,
L<bjobs>,
L<bswitch>,
L<bdel>,
L<bkill>,
L<bstop>,
L<bmod>,
L<btop>,
L<bbot>,
L<brun>,
L<bqueues>,
L<bgadd>,
L<bgdel>,
L<bgmod>,
L<bghold>,
L<bgrel>

=head1 AUTHOR

Mark Southern (mark_southern@merck.com)

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)

=cut

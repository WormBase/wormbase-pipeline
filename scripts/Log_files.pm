#=POD

#=head1 Log_files

#  Create a log file to write script loggin info to.

#=head2
#  SYNOPSIS

#  my $log = Log_files->make_build_log();
#  $log->write_to("This text will appear in log file\n");
#  $log->mail();

#=cut

package Log_files;

use lib $ENV{CVS_DIR};
use Wormbase;
use Carp;
use File::Spec;

=head2 make_log

    Title   :   make_log
    Usage   :   my $log = Log_files->make_log("logfile.log");
                print $log "This is a log file\n";
    Returns :   ref to logging object
    Args    :   filename to use (inlcuding path)

=cut

sub make_log {
    my $class     = shift;
    my $file_name = shift;
    my $debug     = shift;
    my $log;
    my ( $volume, $directories, $fname ) = File::Spec->splitpath($file_name);
    if ($directories) {
        unless ( -d "$directories" and -w "$directories" ) {
            unless ( mkdir( "$directories", 0771 ) ) {
                warn "cant create log dir $directories - using /tmp\n";
                $file_name = "/tmp/$fname";
            }
        }
        $file_name = "/tmp/$fname";
    }
    else {
        $file_name = "/tmp/$file_name";
    }

    open( $log, ">$file_name" ) or croak "cant open file $file_name : $!";

    my $self = {};
    $self->{"FILENAME"} = $file_name;
    $self->{"SCRIPT"}   = $file_name;
    $self->{"FH"}       = $log;
    $self->{"DEBUG"}    = $debug if $debug;
    bless( $self, $class );

    return $self;
}

=head2 make_build_log

Generate log file in the logs directory with WS version and processID appended to script name.

  Title   :   make_build_log
  Usage   :   my $log = Log_files->make_build_log("$debug");
              print $log "This is a log file\n";
  Returns :   filehandle ref
  Args    :   Optional - Debug_status

=cut

sub make_build_log {
    my $class    = shift;
    my $wormbase = shift;
    my $ver      = $wormbase->get_wormbase_version;
    my $filename;
    $0 =~ m/([^\/]*$)/ ? $filename = $0 : $filename =
      $1;    # get filename (sometimes $0 includes full path if not run from its dir )

    my $path     = $wormbase->logs;
    my $log_file = "$path/$filename" . ".WS${ver}." . $$;
    my $log;
    open( $log, ">$log_file" ) or croak "cant open file $log_file : $!";
    print $log "WS$ver Build script : $filename \n\n";
    print $log "Started at ", $wormbase->runtime, "\n";
    print $log "-----------------------------------\n\n";

    my $self = {};
    $self->{"FILENAME"} = $log_file;
    $self->{"SCRIPT"}   = $filename;
    $self->{"FH"}       = $log;
    $self->{"DEBUG"}    = $wormbase->debug;
    $self->{'wormbase'} = $wormbase;

    bless( $self, $class );

    return $self;
}

sub write_to {
    my $self = shift;
    return if $self->{'end'};
    my $string = shift;
    my $fh     = $self->{"FH"};
    print $fh "$string";
    print STDERR "$string" if $self->{'DEBUG'};    # write log to terminal if debug mode
}

sub mail {
    my ( $self, $recipient, $subject ) = @_;
    my $fh = $self->{"FH"};
    print $fh "\n\n-----------------------------------\n";
    print $fh "Finished at ", $self->{'wormbase'}->runtime, "\n";
    close $fh;

    $recipient = "All" unless $recipient;
    $recipient = $self->{"DEBUG"} if ( defined $self->{"DEBUG"} );
    my $file = $self->{"FILENAME"};

    # use subject line from calling script if specified or default to script name
    my $script;
    if ($subject) {
        $script = $subject;
    }
    else {
        $script = $self->{"SCRIPT"};
    }
    $self->{'wormbase'}->mail_maintainer( $script, $recipient, $file );
    $self->{'MAILED'} = 1;
}

sub end {
    my $self = shift;
    my $fh   = $self->{"FH"};
    close $fh;
    $self->{'end'} = 1;
}

# this will cause the log file to be mailed even if the calling script dies unexpectantly
sub DESTROY {
    my $self = shift;
    unless ( defined $self->{'MAILED'} ) {
        my $fh = $self->{"FH"};
        print $fh
"\nTHIS MAIL WAS NOT SENT BY THE SCRIPT THAT WAS RUN\nThe log file object was not \"mailed\", which may mean the script did not finish properly\n";
        $self->{'SCRIPT'} .= " FAILED";
        $self->mail;
    }
}

# this sub is to allow an error to be reported to screen and log file prior to script death
sub log_and_die {
    my $self   = shift;
    my $report = shift;

    if ($report) {
        $self->write_to("$report\n");
        print STDERR "$report\n";
    }
    die;
}

sub error {
    my $self = shift;
    $self->{'ERRORS'}++;
    return $self->{'ERRORS'};
}

sub report_errors {

    # if no errors have been reported return zero
    my $self = shift;
    if ( !defined( $self->{'ERRORS'} ) ) {
        return (0);
    }
    else {
        return $self->{'ERRORS'};
    }
}

# in some circumstances a script may want to read the log file eg make_autoace::check_make_autoace
sub get_file {
    my $self = shift;
    return $self->{'FILENAME'};
}

1;

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


use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");
use Wormbase;
use Carp;

=head2 make_log

    Title   :   make_log
    Usage   :   my $log = Log_files->make_log("logfile.log");
                print $log "This is a log file\n";
    Returns :   ref to logging object
    Args    :   filename to use (inlcuding path)

=cut

sub make_log
  {
    my $class = shift;
    my $file_name = shift;
    my $log;
    open($log,">$file_name") or croak "cant open file $file_name : $!";

    my $self = {};
    $self->{"FILENAME"} = $file_name;
    $self->{"SCRIPT"} = $file_name;
    $self->{"FH"} = $log;
    bless ($self,$class);

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

sub make_build_log
  {
    my $class = shift;
    my $debug = shift;
    my $ver = &Wormbase::get_wormbase_version;
    $ver = $debug unless $ver;  #if wormsrv2 not accessible then Wormbase module wont get version
    my $filename;
    $0 =~ /([^\/]*$)/ ? $filename = $0 :  $filename = $1 ; # get filename (sometimes $0 includes full path if not run from its dir )

    my $path = "/wormsrv2/autoace/logs";
    if( defined $debug ) {
      $path = "/tmp/logs" ;
      system("mkdir $path") unless (-e $path);
    }

    my $log_file = "$path/$filename".".$ver.".$$;
    my $log;
    open($log,">$log_file") or croak "cant open file $log_file : $!";
    print $log "WS$ver Build script : $filename \n";
    print $log "Started at ",&Wormbase::runtime,"\n";
    print $log "-----------------------------------\n\n";

    my $self = {};
    $self->{"FILENAME"} = $log_file;
    $self->{"SCRIPT"} = $filename;
    $self->{"FH"} = $log;

    bless($self, $class);

    return $self;
  }

sub write_to
  {
    my $self = shift;
    my $string = shift;
    my $fh = $self->{"FH"};
    print $fh "$string";
  }


sub mail
  {
   my ($self, $recipient) = @_;

   my $fh = $self->{"FH"};
   print $fh "\n\n-----------------------------------\n";
   print $fh "Finished at ",&Wormbase::runtime,"\n"; 
   close $fh;

   $recipient = "All" unless $recipient;
   my $file = $self->{"FILENAME"};
   my $script = $self->{"SCRIPT"};
   &Wormbase::mail_maintainer($script,$recipient,$file);
  }


1;

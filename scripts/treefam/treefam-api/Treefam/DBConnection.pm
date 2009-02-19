# Author: jkh1
# 2006-08-18
#

=head1 NAME

 Treefam::DBConnection

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection (
                                  -database => 'treefam_3',
                                  -host     => 'vegasrv.sanger.ac.uk',
		                  -user     => 'anonymous',
                                  -port     => 3308 );

 OR to read configuration file:
 my $dbc = Treefam::DBConnection->new(-file=>file);

 OR to use default configuration from Treefam::Config
 my $dbc = Treefam::DBConnection->new();


 my $tree_handle = $dbc->get_TreeHandle();
 my $tree = $tree_handle->get_by_gene($geneID);


=head1 DESCRIPTION

 Provides database connection and handles for
 Treefam objects.

=head1 CONTACT

 jkh1@sanger.ac.uk

=cut


package Treefam::DBConnection;

use strict;
use Carp;
use DBI;
use Treefam::Config;
use Treefam::TreeHandle;
use Treefam::FamilyHandle;
use Treefam::GeneHandle;

sub new {

  my $class = shift;

  my %dbc = @_;

  if (!defined($dbc{'-database'})) {
    if (defined($dbc{'-file'}) && -e $dbc{'-file'}) {
      # use configuration file if provided
      unless (my $return = do $dbc{'-file'}) {
	croak "couldn't parse $dbc{'-file'}: $@" if $@;
	croak "couldn't do $dbc{'-file'}: $!"    unless defined $return;
	croak "couldn't run $dbc{'-file'}"       unless $return;
      }
    }
    $dbc{'-database'} = $TFDBNAME;
    $dbc{'-host'}     = $TFDBHOST;
    $dbc{'-user'}     = $TFDBUSER;
    $dbc{'-password'} = $TFDBPASS;
    $dbc{'-port'}     = $TFDBPORT;
  }

  my $db = $dbc{'-database'};
  die "ERROR: No database selected\n" unless ($db);
  # check that database and API versions match
  unless ($dbc{'-database'}=~/\w_$APIVERSION/) {
    warn "\nWARNING: API version doesn't match database version for $dbc{'-database'}. Some things won't work.\n\n";
    sleep(1); # give time to read message
  }
  my $host = $dbc{'-host'};
  my $port = $dbc{'-port'} if defined($dbc{'-port'});
  my $user_name = $dbc{'-user'};
  my $password = $dbc{'-password'} || "";
  my $dsn="DBI:mysql:$db:$host:$port";

  my $self = {};
  bless ($self, $class);

  my $dbh;

  $SIG{ALRM} = sub { die "timeout"; };
  alarm(10); # time in seconds to wait
  eval {
    $dbh= DBI->connect ($dsn, $user_name, $password, {RaiseError=> 1, PrintError=> 0});
    warn "connecting to database $db on host $host\n";
  };
    if (!$dbh || $@) {
      croak "Could not connect to database $db on host $host: $DBI::errstr\n";
    }
  alarm(0);

  $self->{'database'} = $dbc{'-database'};
  $self->{'database_handle'} = $dbh;

  return $self;

}

sub get_Database {

  my $self = shift;
  return $self->{'database'};
}

sub get_DatabaseHandle {

  my $self = shift;
  $self->{'database_handle'} = shift if @_;
  return $self->{'database_handle'};

}

sub disconnect {

  my $self = shift;
  $self->{'database_handle'}->disconnect();
}

sub get_FamilyHandle {

  my $self = shift;
  return new Treefam::FamilyHandle($self);

}

sub get_TreeHandle {

  my $self = shift;
  return Treefam::TreeHandle->new($self);

}

sub get_GeneHandle {

  my $self = shift;
  return Treefam::GeneHandle->new($self);

}



1;

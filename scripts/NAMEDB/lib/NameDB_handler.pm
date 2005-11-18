#author ar2
package NameDB_handler;

use lib "/nfs/WWWdev/SANGER_docs/cgi-bin/Projects/C_elegans/lib";
use NameDB;
our @ISA = qw( NameDB );

sub new
  {
    my $class = shift;
    my $dsn = shift;
    my ($name,$password) = @_;
    my $self = NameDB->connect($dsn,$name,$password);

    bless ($self, $class);
    return $self;
  }

sub web  # set flag to specify web or script for error (dienice)
  {
    my $self = shift;
    my $set = shift;
    $self{'web'} = $set if( $set );
    return $self{'web'};
  }

sub printAllNames
  {
    my $self = shift;
    my $obj = shift;
    my %names = $self->idAllNames($obj);
    print "<br>Current gene names for <b>$obj</b> <br>";
    foreach (keys %names ) {
      print "$_ : ",join(" ",@{$names{$_}}) || $names{$_} ,"<br>";
    }
  }

sub add_name
  {
    my $self = shift;
    my $id = shift;
    my $name = shift;
    my $type = shift;
    eval {
      $self->addName($id,$type => $name);
    };
    if ($@) {
      $self->dienice("$@");
    }
  }

sub isoform_exists
  {
    my $self = shift;
    my $name = shift;
    my $type = shift;
    if ( $name =~ /^(\w+\.\d+)\w/) {
      $id = $self->idGetByTypedName($type,$1);
      if ( $id->[0] ) {
	print "adding $name as an isoform of ".$id->[0]."<BR>";
	$self->addTypeName($id->[0],"$type => $name");
      }
    }
  }

sub validate_name
  {
    my $self = shift;
    my $name = shift;
    my $type = shift;
    #is this a valid name tpye?
    my @types = $self->getNameTypes;
    if( grep {$_ eq $type} @types) {
      #check name structure matches format eg CDS = clone.no
      my %name_checks = ( "CDS" => '^\w+\.\d+\w?',
			  "CGC" => '[a-z]{3,4}-\d+'
			);
      unless( $name =~ /$name_checks{$type}/ ) {
	$self->dienice("$name is incorrect format for $type<br>");
      }
    }
    else {
      $self->dienice("$type is a invalid typename:<br>@types");
    }
  }

sub validate_id
  {
    my $self = shift;
    my $id = shift;
    unless ( $self->idExists($id) ) {
      $self->dienice("$id does not exist<br>");
    }
    unless( $self->idLive($id) ) {
      $self->dienice("$id is not live<br>");
    }
  }

sub check_pre_exists
  {
    my $self = shift;
    my $name = shift;
    my $type = shift;
    $id = $self->idGetByTypedName($type,$name);
    if(defined $id->[0] ){
      my $err = "$name already exists as ".$id->[0].":<br>";
      $self->dienice("$err");
    }
  }


sub make_new_obj
  {
    my $self = shift;
    my $name = shift;
    my $type = shift;
    $id = $self->create_named_object($type, $name);
    print "The new id for $type: $name is $id\n";
  }

sub dienice {
    my $self = shift;
    my($errmsg) = @_;
    if( $self->web ) {
      print "<h2>Error</h2>\n";
      print "<p>$errmsg</p>\n";
      print end_html;
    }
    else {
      print "$errmsg\n";
    }
    exit;
  }

sub print_history {
  my $self = shift;
  my $id = shift;

  $self->validate_id($id);

  my $history = $self->idGetHistory($id);
  print "$id<br>";
  foreach my $event (@{$history}) {
    print $event->{version}," ";
    print $event->{date}," ";
    print $event->{event}," ";
    print $event->{name_type}," " if (defined $event->{name_type});
    print $event->{name_value}," " if (defined $event->{name_value});
    print $event->{user},"<br>";
  }
}

1;

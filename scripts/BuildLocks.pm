package BuildLocks;

sub new
  {
    my $class = shift;
    my $self = {};
    $self->{'dir'}   = shift;

    my $existing_locks = "$self->{'dir'}"."/LOCKS";

    # load current status

    %{$self->{'locks'}} = ( "test.pl"      => 'A1',
			    "trythis.pl"   => 'A2',
			    "phase2.pl"    => 'B1',
			    "test_lock.pl" => 'B2',
			    "phase3.pl"    => 'C1'
			  );

    foreach my $lock ( keys %{$self->{'locks'}} ) {
      my $phase = substr($self->{'locks'}->{"$lock"},0,1);
      my $event = substr($self->{'locks'}->{"$lock"},1);
      $self->{'key_locks'}->{$self->{'locks'}->{"$lock"}} = $lock;
      $self->{'status'}->{"$phase"}->{"$event"} = 0;
    }

    bless( $self, $class);
    return $self;
  }

sub lock_this
  {
    my $self = shift;
    my $script = shift;

    if( $self->{'locks'}->{$script} ) {
      my $lock_file = "$self->{'locks'}->{$script}";
      system("touch $self->{'dir'}/$lock_file:$script");
      
      my $phase = substr($self->{'locks'}->{"$script"},0,1);
      my $event = substr($self->{'locks'}->{"$script"},1);

      $self->{'status'}->{"$phase"}->{"$event"} = 1;
    }
    else {
      warn "$script does not have a designated lock\n";
    }
  }

sub check_lock
  {
    my $self = shift;
    my $script = shift;
    $script = caller unless defined $script;

    my $phase = substr($self->{'locks'}->{"$script"},0,1);
    my $event = substr($self->{'locks'}->{"$script"},1);

    if( $self->{'status'}->{"$phase"}->{"$event"} == 1 ) {
      warn "already done this locking event\n";
    }
    else {
      return 0 if $phase eq "A"; # can always do this
      my $prev_phase = chr(ord($phase)-1);
      foreach $p_event ( keys %{$self->{'status'}->{"$prev_phase"}} ) {
	if( $self->{'status'}->{"$prev_phase"}->{"$p_event"} == 0 ) {
	  my $key = "$prev_phase"."$p_event";
	  print "$script attempting to run when lock $key - ",$self->{'key_locks'}->{"$key"}," not done\n";
	  return 1;
	}
      }
    }
    return 0;
  }

1;

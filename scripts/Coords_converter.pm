package Coords_converter;
use Carp;

sub invoke 
  {
    my $class = shift;
    my $self = {};
    open (SL,"</nfs/disk100/wormpub/DATABASES/TEST_DBs/ANT_matchingace/dump/superlink_coords") or croak "cant open superlinks";
      my $parent;
    while (<SL>) {
      if(/Sequence.*(CHROMOSOME_\w+)/) {
	$parent = $1;
      }
      elsif( /Subsequence\s+\"(SUPERLINK_\w+)\"\s+(\d+)\s+(\d+)/ ) {
	$self->{"$parent"}->{'SUPERLINK'}->{"$1"} = [$2,$3];
      }
    }

    undef $parent;
    open (CLONE,"</nfs/disk100/wormpub/DATABASES/TEST_DBs/ANT_matchingace/dump/clone_coords") or croak "cant open clones";
    while(<CLONE>) {
      if(/Sequence.*\"(SUPERLINK_\w+)/) {
	$parent = $1;
      }
      elsif( /Subsequence\s+\"(\w+)\"\s+(\d+)\s+(\d+)/ ){
	$self->{'SUPERLINK'}->{"$parent"}->{"$1"} = [$2,$3];
      }
    }
    %{$self->{"numerals"}} = ( "1" => "I",
			    "2" => "II",
			    "3" => "III",
			    "4" => "IV",
			    "5" => "V",
			    "X" => "X"
			  );
			
    bless( $self, $class);
    return $self;
  }

sub GetSuperlinkFromCoord
  {
    my $self = shift;
    my $chrom = shift;
    my $coord = shift;

    foreach my $slink ( keys %{$self->{"CHROMOSOME_$chrom"}->{'SUPERLINK'}} ) {
      if($self->{"CHROMOSOME_$chrom"}->{'SUPERLINK'}->{$slink}->[0] < $coord and 
	 $self->{"CHROMOSOME_$chrom"}->{'SUPERLINK'}->{$slink}->[1] > $coord
	) {
	return $slink;
      }
    }
    return 0;
  }

sub GetCloneFromCoord
  {
    my $self = shift;
    my $parent = shift;
    my $coord = shift;
    my $chrom;

    if( length("$parent") < 5 ) {
      # parent is a chromosome number so lets find the right slink
      $chrom = $parent;
      $parent = &GetSuperlinkFromCoord($self,"$parent", $coord);
      my $sl_start = $self->{"CHROMOSOME_$chrom"}->{'SUPERLINK'}->{"$parent"}->[0];
      $coord -= $sl_start;
    }

    # get data with SUPERLINK as parent
    $chrom = &_getChromFromSlink($self,"$parent") unless $chrom;

    foreach $clone (keys %{$self->{'SUPERLINK'}->{"$parent"}} ) {
      if( $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0] < $coord and
	  $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[1] > $coord
	){
	return $clone;
      }
    }
  }

sub _getChromFromSlink
  {
    my $self = shift;
    my $sl = shift;
    return 0 unless defined $sl;

    if( $sl =~ /SUPERLINK_CB_(\w+)/) {
      return $1;
    }
    elsif( $sl =~ /SUPERLINK_RW(\w+)/) {
      return $self->{'numerals'}->{"$1"};
    }
    else{
      return 0;
    }
  }

sub LocateSpan
  {
    my $self = shift;
    my $chrom = shift;
    my $x = shift;
    my $y = shift;

    my $strand = "+";
    my ($seq, $rel_x, $rel_y);  # coordinates relative to returned seq
    
    # set strand and make sure $x is always lower coord
    if( $y < $x ) {
      $strand = "-";
      &_swap(\$x, \$y);
    }

    my $x_slink = &GetSuperlinkFromCoord( $self, $chrom, $x); # need this whatever
    my $sl_start = $self->{"CHROMOSOME_$chrom"}->{'SUPERLINK'}->{"$x_slink"}->[0];

    # see if the span is within a clone
    my $x_clone = &GetCloneFromCoord( $self, $chrom, $x);
    my $y_clone = &GetCloneFromCoord( $self, $chrom, $y);
    if( $x_clone eq $y_clone ) {# span maps within clone
      my $clone_start = $self->{SUPERLINK}->{"$x_slink"}->{"$x_clone"}->[0];
      $rel_x = $x - $sl_start - $clone_start +1;
      $rel_y = $y - $sl_start - $clone_start +1;
      $seq = $x_clone;
    }
    else {
      # locate on Slink ;
      my $y_slink = &GetSuperlinkFromCoord( $self, $chrom, $y);
 
      if( $x_slink eq $y_slink ) { # span maps within superlink
	$rel_x = $x - $sl_start + 1;
	$rel_y = $y - $sl_start + 1;
	$seq = $x_slink;
      }
      else { # spans superlink so goes on to the chromosome
	$rel_x = $x; 
	$rel_y = $y;
	$seq = "CHROMOSOME_$chrom";
      }
    }
    if( $strand eq "-" ){&_swap(\$rel_x,\$rel_y)}
    return ($seq, $rel_x, $rel_y);
  }

sub _swap
  {
    my $x = shift;
    my $y = shift;
    my $tmp = $$x;
    $$x = $$y;
    $$y = $tmp;
  }


return 1;

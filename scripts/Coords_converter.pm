=pod 

=head1 NAME

Coords_converter

=head1 SYNOPSIS

my $coords       = Coords_converter->invoke

my $superlink    = $coords->GetSuperlinkFromCoord("I", 10000000)

my $clone        = $coords->GetCloneFromCoord( "I", 1000000 )

$clone           = $coords->GetCloneFromCoord( "SUPERLINK_RW1", 100000 )

my @clone_coords = $coords->LocateSpan("I",5239404,5271341 )

=head1 DESCRIPTION

This object can be used to get from chromosomal coordinates eg those in gff file to clone or superlink coordinates.
Can also give a clone or superlink from specified chromosome coordinates

=head1 CONTACT

Anthony  ar2@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut




package Coords_converter;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

use Carp;
use Wormbase;

=head2 invoke

    Title   :   invoke
    Usage   :   Coords_converter->invoke($database,1);
    Function:   Creates the object and loads in the data (generates fresh if requested)
    Returns :   ref to self
    Args    :   Database  (optional) - which database to use. Default is current_DB
                refresh   default is NULL - connect to the database and update coordinates

=cut


sub invoke 
  {
    my $class = shift;
    my $database = shift;
    my $refresh = shift;

    if( $database ) {
      croak "$database does not exist\n" unless( -d "$database" );
    }
    else {
      $database = glob("~wormpub/DATABASES/current_DB");
      undef $refresh;
    }

    unless( -e "$database/SL_coords.ace" and -e "$database/clone_coords.ace") {
      warn "\nno coordinate info was present - extracting new data\n"; 
      $refresh = 1;
    }

    if( $refresh ) {
      print "refreshing coordinates for $database\n";
      my $tace = &tace;
      my $SL_coords_file = "$database/SL_coords.ace";
      my $clone_coords_file = "$database/clone_coords.ace";

      my @command;
      $command[0] = "find sequence CHROM*\nshow -a subsequence -f ${SL_coords_file}\n";
      $command[1] = "clear\nfind sequence SUPER*\nshow -a subsequence -f ${clone_coords_file}\n";

      # temp fix as &tace isn't working in this script
      $tace = "/nfs/disk100/wormpub/ACEDB/bin_ALPHA/tace";


      open (ACE,"| $tace $database") or croak "cant open $database\n";
      foreach (@command) {
	print ACE $_ ;
      }
      close(ACE);

      # need to remove genes SMapped to slinks
      open( FH, "$clone_coords_file") or die "bugger";
      open( NEW,">$clone_coords_file.bk" ) or die "bugger2";
      while (<FH>) {
	print NEW unless /\./;
      }
      close NEW;
      close FH;

      system("mv -f $clone_coords_file.bk $clone_coords_file") and croak "cant mv $clone_coords_file.bk\n" ;
      system("chmod 777 $clone_coords_file") and carp "cant chmod on $clone_coords_file - This could cause problems in future" ;

    }

    my $self = {};
    open (SL,"<$database/SL_coords.ace") or croak "cant open superlink coordinate file ";
      my $parent;
    while (<SL>) {
      if(/Sequence.*(CHROMOSOME_\w+)/) {
	$parent = $1;
      }
      elsif( /Subsequence\s+\"(SUPERLINK_\w+)\"\s+(\d+)\s+(\d+)/ ) {
	$self->{"$parent"}->{'SUPERLINK'}->{"$1"} = [$2,$3];
	$self->{'SUPERLINK2CHROM'}->{"$1"} = "$parent";
      }
    }

    undef $parent;
    open (CLONE,"</$database/clone_coords.ace") or croak "cant open clone coordinate file";
    while(<CLONE>) {
      if(/Sequence.*\"(SUPERLINK_\w+)/) {
	$parent = $1;
      }
      elsif( /Subsequence\s+\"(\w+)\"\s+(\d+)\s+(\d+)/ ){
	$self->{'SUPERLINK'}->{"$parent"}->{"$1"} = [$2,$3];
	$self->{'CLONE2CHROM'}->{"$1"} = $self->{'SUPERLINK2CHROM'}->{"$parent"};
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

=head2 GetSuperlinkFromCoord

    Title   :   GetSuperlinkFromCoord
    Usage   :   my $superlink = $coords->GetSuperlinkFromCoord('III',4500000)
    Function:   Determine what superlink any base is in
    Returns :   superlink as string eg "SUPERLINK_CB_V"
    Args    :   chromosome as string eg "V",  coordinate as int

=cut

  sub GetSuperlinkFromCoord
  {
    my $self = shift;
    my $chrom = shift;
    my $coord = shift;

    $chrom = "CHROMOSOME_$chrom" unless $chrom =~ /CHROMOSOME/;

    foreach my $slink ( keys %{$self->{"$chrom"}->{'SUPERLINK'}} ) {
      if($self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] < $coord and
	 $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[1] > $coord
	) {
	return $slink;
      }
    }
    return 0;
  }



=head2 GetCloneFromCoord

    Title   :   GetCloneFromCoord
    Usage   :   my $clone = $coords->GetCloneFromCoord('III',1200000)
    Function:   Determine what clone any base is in
    Returns :   clone as string eg "AH6"
    Args    :   chromosome as string eg "III",  coordinate as int

=cut

sub GetCloneFromCoord
  {
    my $self = shift;
    my $parent = shift;
    my $coord = shift;
    my $chrom;

    unless( "$parent" =~ /SUPER/ ) {
      # parent is a chromosome so lets find the right slink
      $chrom = $parent;
      $chrom = "CHROMOSOME_$chrom" unless $chrom =~ /CHROMOSOME/;

      $parent = &GetSuperlinkFromCoord($self,"$parent", $coord);
      my $sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{"$parent"}->[0];
      $coord -= $sl_start;
    }

    # get data with SUPERLINK as parent
    $chrom = $self->{SUPERLINK2CHROM}->{"$parent"} unless $chrom;

    foreach $clone (keys %{$self->{'SUPERLINK'}->{"$parent"}} ) {
      if( $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0] < $coord and
	  $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[1] > $coord
	){
	return $clone;
      }
    }
  }


=head2 LocateSpan

    Title   :   LocateSpan
    Usage   :   my @data = $coords->LocateSpan("CHROMOSOME_I", 100000, 101000);
    Function:   Returns smallest sequence object containing the given coordinate ie clone, slink or csome and the coordinates relative to that sequence.
    Returns :   Clone as string eg "V"; coordinates relative to clone obj 
    Args    :   Sequence obj (CHROMOSOME or SUPERLINK)  and two coordinates within that eg ("CHROMOSOME_I", 3000, 3100)

=cut

sub LocateSpan
  {
    my $self = shift;
    my $chrom = shift;
    my $x = shift;
    my $y = shift;

    #this doesn't handle single letter style chromosome names


    # if anything other than chromsome or superlink is passed you get back what you put in
    unless( $chrom =~ /CHROMOSOME/ or $chrom =~ /SUPERLINK/ ) { 
      my ($superlink, $chromosome, $seq);
      $seq = $chrom;
      undef $chrom;
    SLINKS:
      foreach my $slink (keys %{$self->{'SUPERLINK'}} ) {
	foreach my $clone (keys %{$self->{SUPERLINK}->{$slink}} ) {

	  if( "$clone" eq "$seq" ) {
	    $chrom = $self->_getChromFromSlink("$slink");

	    # modify the starting coordinate
	    my $offset = $self->{"$chrom"}->{SUPERLINK}->{"$slink"}->[0] -1 ; # superlink base coords
	    $offset += $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[0] - 1; # clone base coords
	
	    $x += $offset;
	    $y += $offset;
	    last SLINKS;
	  } 
	}
      }
    }

    if( $chrom =~ /SUPERLINK/ ) {
      my $slink = $chrom;
      $chrom = $self->{SUPERLINK2CHROM}->{"$slink"};
      $x += $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] - 1;
      $y += $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] - 1;
    }

    my $strand = "+";
    my ($seq, $rel_x, $rel_y);  # coordinates relative to returned seq
    
    # set strand and make sure $x is always lower coord
    if( $y < $x ) {
      $strand = "-";
      $self->swap(\$x, \$y);
    }

    $chrom = "CHROMOSOME_$chrom" unless $chrom =~ /CHROMOSOME/;
    my $x_slink = &GetSuperlinkFromCoord( $self, $chrom, $x); # need this whatever
    my $sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{"$x_slink"}->[0];

    # see if the span is within a clone
    my $x_clone = &GetCloneFromCoord( $self, $chrom, $x);
    my $y_clone = &GetCloneFromCoord( $self, $chrom, $y);
    if( $x_clone eq $y_clone ) {# span maps within clone
      my $clone_start = $self->{SUPERLINK}->{"$x_slink"}->{"$x_clone"}->[0];
      $rel_x = $x - $sl_start - $clone_start + 2;  # add 2 -  1 for clone, 1 for slink numbering adjustment
      $rel_y = $y - $sl_start - $clone_start + 2;
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
	$seq = "$chrom";
      }
    }
    if( $strand eq "-" ){$self->swap(\$rel_x,\$rel_y)}
    return ($seq, $rel_x, $rel_y);
  }

=head2 swap

    Title   :   swap
    Usage   :   swap(\$x,\$y);
    Function:   swaps values of variables passed in
    Returns :   nothing
    Args    :   references to two coordinates

=cut

sub swap
  {
    my $self = shift;
    my $x = shift;
    my $y = shift;
    carp "undefined values passed to swap \n" unless ($$x and $$y);
    my $tmp = $$x;
    $$x = $$y;
    $$y = $tmp;
  }

sub get_Chrom_from_clone
  {
    my $self = shift;
    my $clone = shift; 

    my $chrom = $self->{"CLONE2CHROM"}->{$clone} ? $self->{"CLONE2CHROM"}->{$clone} : undef;
    carp "unknown clone \"$clone\" requesting chromosome info" unless $chrom;

    return $chrom;
  }


sub _getChromFromSlink
  {
    my $self = shift;
    my $sl = shift;

    return $self->{SUPERLINK2CHROM}->{"$sl"};
  }

return 1;

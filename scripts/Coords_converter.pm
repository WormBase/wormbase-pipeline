=pod 

=head1 NAME

Coords_converter

=head1 SYNOPSIS

my $coords       = Coords_converter->invoke( $database, 1, $wormbase )

my $superlink    = $coords->GetSuperlinkFromCoord( "I", 10000000 )

my $offset       = $coords->CloneOffset( "AH6" )

my $clone        = $coords->GetCloneFromCoord( "I", 1000000 )
my $clone        = $coords->GetCloneFromCoord( "SUPERLINK_RW1", 100000 )

my $clone = $coords->GetCloneFromCoord( 'III',1200000 )

my @clone_coords = $coords->LocateSpan( "I",5239404,5271341 )

my $superlink = $coords->get_SuperlinkFromClone( "AH6" );

my $chromosome = $coords->get_Chrom_from_clone( "AH6" );

my $chromosome = $coords->_getChromFromSlink( "SUPERLINK_CB_I" );

my $chromosome = $coords->seq_obj_to_chrom( "$seq" )

my ($chrom, $coord) = $coords->Coords_2chrom_coords( "AH6", 2343 );

my $sl_length = $coords->Superlink_length( "SUPERLINK_RWXL" );

my result = $coords->isa_chromosome($seq)

my result = $coords->isa_superlink($seq)

my result = $coords->isa_clone($seq)

=head1 DESCRIPTION

This object can be used to get from chromosomal coordinates eg those
in gff file to clone or superlink coordinates.  Can also give a clone
or superlink from specified chromosome coordinates

=head1 CONTACT

Anthony  ar2@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Coords_converter;

use strict;
use lib $ENV{CVS_DIR};

use Carp;
use Wormbase;


=head2 invoke

Title   :   invoke
Usage   :   Coords_converter->invoke($database, 1, $wormbase);
Function:   Creates the object and loads in the data (generates fresh if requested)
Returns :   ref to self
Args    :   Database  (optional) - which database to use. Default is current_DB
		refresh   default is NULL - connect to the database and update coordinates
		wormbase - wormbase object

=cut

sub invoke {
  my $class = shift;
  my $database = shift;
  my $refresh = shift;
  my $wormbase = shift;
  
  my $self = {};
  
  $self->{species} = $wormbase->species;
  $self->{chromosome_prefix} = $wormbase->chromosome_prefix ;
  my $species = $self->{species};
  
  if( $database ) {
    croak "$database does not exist\n" unless( -d "$database" );
  }
  else {
    # the sequence files and coords files in current_DB are elegans-specific
    # so use the appropriate species' database if this is not elegans
    if ($species eq 'elegans') {
      $database = $wormbase->database("current");
    } else {
      $database = $wormbase->database($species);
    }
    undef $refresh;
  }
  
  #print "using database: $database\n";
  #print "species: $species\n";
  
  my $SL_coords_file = "$database/SL_coords";
  my $clone_coords_file = "$database/clone_coords";
  
  unless( -e $SL_coords_file and -e $clone_coords_file) {
    warn "\nno coordinate info was present - extracting new data\n";
    $refresh = 1;
  }

  if ($refresh) {
    my $db = Ace->connect(-path => $database);
    
    my $ncbi_bioproj = $wormbase->ncbi_bioproject;

    my $spe_obj = $db->fetch(Species => $wormbase->full_name);
    my @assemblies = $spe_obj->Assembly;
    my @relevant_assembly;
    foreach my $assembly (@assemblies) {
      foreach my $dbxref ($assembly->Database) {
        next if $dbxref->name ne "NCBI_BioProject";
        next if $assembly->Dead;
        my $bp = $dbxref->right->right->name;
        
        if ($bp eq $ncbi_bioproj) {
          push @relevant_assembly, $assembly;
        }
      }
    }

    if (scalar(@relevant_assembly) != 1) {
      die "Bad data when trying to refresh Coords_converter data\n";
    }
    
    open(my $sl_coords, ">$SL_coords_file") or die "Could not open $SL_coords_file for writing\n";

    foreach my $seq ($relevant_assembly[0]->Sequences) {
      print $sl_coords "\nSequence : \"", $seq->name, "\"\n";
      my @segs;

      my @subseqs = $seq->Subsequence;
      if (scalar(@subseqs) <= 1) { 
        my $seqobj;
        if (not @subseqs) {
          $seqobj = $seq->fetch;
        } else {
          $seqobj = $subseqs[0]->fetch;
        }
        my $seq_len = $seqobj->at('DNA')->right->right;

        push @segs, [$seqobj->name, 1, $seq_len];
      } else {
        foreach my $subseq (@subseqs) {
	  #The Superlinks are still connected to their parent sequence but not mapped so need to exclude.
          push @segs, [$subseq->name, $subseq->right->name, $subseq->right->right->name] unless ($subseq->name =~ /^SUPER/);
        }
      }

      foreach my $seg (sort { $a->[1] <=> $b->[1] } @segs) {
        printf $sl_coords "Subsequence    \"%s\" %d %d\n", @$seg;
      }
    }
    close($sl_coords);
    $db->close();
    system("cp $SL_coords_file $clone_coords_file") and croak "cant cp $SL_coords_file\n";    
  }

  &_read_data_file($database, $species, $self);
  
  bless( $self, $class);
  return $self;
}


=head2 _read_data_file

Title   :   _read_data_file
Usage   :   _read_data_file
Function:   loads in the data from a data file
Returns :
Args    :   Database - which database to use. Species - species being used. $self

=cut

sub _read_data_file {

  my $database = shift;
  my $species = shift;
  my $self = shift;
  
  my $parent;
  
  # In all species (now including elegans) there is no superlink, there is just
  # contig and supercontig (or chromosome/supercontig).  Supercontigs will 
  # eventually becomesynonymous with chromosomes, so we treat them as such here.  
  # We also create dummy superlinks that are the same size as the
  # supercontigs. This makes it easy to re-use the subroutines that
  # were written for elegans and which expect there to be a
  # superlink layer in the data.
  
  # The SL_coords and clone_coords files are duplicates of each
  # other in non-elegans sepcies, so only read from one of them

  my $supercontig;
  my $contig;
  open (CLONE,"<$database/clone_coords") or croak "cant open clone coordinate file $database/clone_coords\t$!";
  while(<CLONE>) {
    if (/Sequence\s+:\s+\"(\S+)\"/) {
      $supercontig = $1;
      
    } elsif (/Subsequence\s+\"(\S+)\"\s+(\d+)\s+(\d+)/) {
      $contig = $1;
      
      # some species (e.g. briggsae) have the contigs in a mixture of orientations
      my $sense = '+';
      my ($x, $y) = ($2, $3);
      if( $y < $x ) {
        $sense = "-";
        swap($self, \$x, \$y);
        #print "reverse sense found for $species $supercontig contig $contig\n";
      }
      $self->{'SENSE'}->{$contig} = $sense; # store the sense of the contig
      $self->{'LENGTH'}->{$contig} = $y - $x + 1; # store the length of the contig
      
      # we wish to store the end of the dummy superlink as being the length of the supercontig
      # so get the largest contig position in the supercontig
      if (! exists $self->{$supercontig}->{'SUPERLINK'}->{$supercontig} ||
          $self->{$supercontig}->{'SUPERLINK'}->{$supercontig}->[1] < $y) {
        $self->{$supercontig}->{'SUPERLINK'}->{$supercontig} = [1,$y]; # a dummy superlink is the supercontig
        $self->{'LENGTH'}->{$supercontig} = $y; # store the length of the supercontig
      }
      $self->{'SUPERLINK2CHROM'}->{$supercontig} = $supercontig; # the dummy superlink is the supercontig
      $self->{'SUPERLINK'}->{$supercontig}->{$contig} = [$x,$y]; # the superlink of a contig is the supercontig
      $self->{'CLONE2CHROM'}->{$contig} = $supercontig;
      $self->{'CLONE2SUPERLINK'}->{$contig} = $supercontig;
    }
  }
  close (CLONE);
}


=head2 GetSuperlinkFromCoord

Title   :   GetSuperlinkFromCoord
Usage   :   my $superlink = $coords->GetSuperlinkFromCoord('III',4500000)
Function:   Determine what superlink any base is in
Returns :   superlink as string eg "SUPERLINK_CB_V"
Args    :   chromosome as string eg "V",  coordinate as int
Non-elegans : Simply returns the input supercontig name

=cut

sub GetSuperlinkFromCoord {
  my $self = shift;
  my $chrom = shift;
  my $coord = shift;
  
  
  if ($self->{chromosome_prefix} eq "CHROMOSOME_" &&  $chrom !~ /CHROMOSOME_/) {
    $chrom = "CHROMOSOME_$chrom";
  }
  
  foreach my $slink ( keys %{$self->{"$chrom"}->{'SUPERLINK'}} ) {
    if($self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] <= $coord and
       $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[1] >= $coord
       ) {
      return $slink;
    }
  }
  return 0;
}


=head2 CloneOffset

Title   :   CloneOffset
Usage   :   my ($chromosome, $offset) = $coords->CloneOffset("AH6")
Function:   Return the starting position in the chromosome of the specified clone or superlink
Returns :   1) chromosome that the clone is in
		2) chromosomal coordinate of the first base of the clone as int
Args    :   clone name as string
Non-elegans : Returns the starting position of the contig in the supercontig

=cut

sub CloneOffset {
  my $self = shift;
  my $clone = shift;
  
  my $superlink;
  
  if ( $self->isa_chromosome($clone) ) {return ($clone, 1);}
  
  if ( $self->isa_superlink($clone) ) {
    $superlink = $clone;
    
  } else {
    $superlink = $self->get_SuperlinkFromClone($clone);
  }
  
  my $chromosome = $self->_getChromFromSlink($superlink);
  my $sl_start = $self->{"$chromosome"}->{'SUPERLINK'}->{"$superlink"}->[0];
  
  if ( $self->isa_superlink($clone) ) {
    return ($chromosome, $sl_start);
  }
  
  my $clonecoord = $sl_start + $self->{'SUPERLINK'}->{"$superlink"}->{"$clone"}->[0] - 1;
  return ($chromosome, $clonecoord);
  
}


=head2 GetCloneFromCoord

Title   :   GetCloneFromCoord
Usage   :   my $clone = $coords->GetCloneFromCoord('III',1200000)
Function:   Determine what clone any base is in
Returns :   clone as string eg "AH6"
Args    :   chromosome as string eg "III",  coordinate as int
Non-elegans : returns the contig (or gap) at that position of the supercontig

=cut

sub GetCloneFromCoord {
  my $self = shift;
  my $parent = shift;
  my $coord = shift;
  
  if ($self->{'species'} eq 'elegans') {
    # only elegans has superlinks that are different to the chromosomes
    # only the mitochondrial superlink MTCE has a name that does not
    # start 'SUPER'
    
    unless( "$parent" =~ /SUPER/ || $parent eq 'MTCE') {
      # parent is a chromosome so lets find the right slink
      my $chrom = $parent;
      if ($self->{chromosome_prefix} eq "CHROMOSOME_" &&  $chrom !~ /CHROMOSOME_/) {
	$chrom = "CHROMOSOME_$chrom";
      }
      
      $parent = $self->GetSuperlinkFromCoord($chrom, $coord);
      my $sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{"$parent"}->[0];
      $coord = $coord - $sl_start + 1;	# so $coord is now the coordinate in the $parent superlink
    }
  }
  
  my @matching_clones;
  foreach my $clone (keys %{$self->{'SUPERLINK'}->{"$parent"}} ) {
    if( $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0] <= $coord and
	$self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[1] >= $coord
	){
      push @matching_clones, {
        clone => $clone,
        lower_bound => $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0],
      };
    }
  }
  @matching_clones = sort { $b->{lower_bound} <=> $a->{lower_bound} } @matching_clones;
  return $matching_clones[0]->{clone};
}

=head2 GetCloneFromCoords

Title   :   GetCloneFromCoords
Usage   :   my $clone = $coords->GetCloneFromCoords('III', 1200000, 1200100)
Function:   Determine what clone any pair of bases  is in
Returns :   clone as string eg "AH6", or "" if no clone covers both positions,
		superlink for coordinate 1,
		superlink for coordinate 2,
		Start pos of clone in superlink1
Args    :   chromosome as string eg "III",  coordinate as int, coordinate as int
Non-elegans : returns the contig (or gap) at those positions of the supercontig, or ""

=cut


sub GetCloneFromCoords {
  my $self = shift;
  my $parent = shift;
  my $coord1 = shift;
  my $coord2 = shift;
  
  my $parent2;
  my $sl_start;
  
  if ($coord1 < 1) {
    carp "Chromosomal coordinate of less than 1 (coord1 = $coord1) passed to GetCloneFromCoords\n";
    $coord1 = 1;
  }
  if ($coord2 < 1) {
    carp "Chromosomal coordinate of less than 1 (coord2 = $coord2) passed to GetCloneFromCoords\n";
    $coord2 = 1;
  }
  if ($self->isa_chromosome($parent) && $coord1 > $self->Superlink_length($parent)) {
    carp "Chromosomal coordinate greater than length of chromosome (coord1 = $coord1) passed to GetCloneFromCoords\n";
    $coord1 = 1;
  }
  if ($self->isa_chromosome($parent) && $coord2 > $self->Superlink_length($parent)) {
    carp "Chromosomal coordinate greater than length of chromosome (coord2 = $coord2) passed to GetCloneFromCoords\n";
    $coord2 = 1;
  }
  
  if ($self->{'species'} eq 'elegans') {
    # only elegans has superlinks that are different to the chromosomes
    # only the mitochondrial superlink MTCE has a name that does not
    # start 'SUPER'
    
    unless( "$parent" =~ /SUPER/ || $parent eq 'MTCE') {
      # parent is a chromosome so lets find the right slink
      my $chrom = $parent;
      if ($self->{chromosome_prefix} eq "CHROMOSOME_" &&  $chrom !~ /CHROMOSOME_/) {
	$chrom = "CHROMOSOME_$chrom";
      }
      
      # see if we can find a single superlink with both positions on
      
      my $got_parent = "";
      foreach my $slink ( keys %{$self->{"$chrom"}->{'SUPERLINK'}} ) {
	if(
           $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] <= $coord1 &&
           $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[1] >= $coord1 &&
           $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] <= $coord2 &&
           $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[1] >= $coord2
           ) {
          $got_parent = $slink;
	}
      }
      
      if ($got_parent) {
	$parent = $got_parent;
	$parent2 = $got_parent;
	$sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{"$parent"}->[0];
	#print "found coords in one slink $got_parent sl_start $sl_start\n";
	$coord1 = $coord1 - $sl_start + 1;	# so $coord1 is now the coordinate in the $parent superlink
	$coord2 = $coord2 - $sl_start + 1;	# so $coord2 is now the coordinate in the $parent superlink
        
      } else {		# didn't find a single superlink with both positions on
	$parent = $self->GetSuperlinkFromCoord($chrom, $coord1);
	$sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{"$parent"}->[0];
	$coord1 = $coord1 - $sl_start + 1;	# so $coord1 is now the coordinate in the $parent superlink
        
	$parent2 = $self->GetSuperlinkFromCoord($chrom, $coord2);
	if ($parent ne $parent2) {return ("", $parent, $parent2, $sl_start);}
	$coord2 = $coord2 - $sl_start + 1;	# so $coord2 is now the coordinate in the $parent superlink
      }
    }
  } else {			# non-elegans species have superlinks that are the same as supercontigs
    $parent2 = $parent;
    $sl_start = 1;
  }
  
  my @matched_clones;
  foreach my $clone (keys %{$self->{'SUPERLINK'}->{"$parent"}} ) {
    if(
       $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0] <= $coord1 &&
       $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[1] >= $coord1 &&
       $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0] <= $coord2 &&
       $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[1] >= $coord2
       ){
      push @matched_clones, {
        clone => $clone,
        lower_bound => $self->{'SUPERLINK'}->{"$parent"}->{"$clone"}->[0],
      };
    }
  }
  
  if (@matched_clones) {
    @matched_clones = sort { $b->{lower_bound} <=> $a->{lower_bound} } @matched_clones;
    return ($matched_clones[0]->{clone},$parent,$parent2,$sl_start);
  } else {
    return ("", $parent, $parent2, $sl_start);
  }
}

=head2 LocateSpan

Title   :   LocateSpan
Usage   :   my @data = $coords->LocateSpan("CHROMOSOME_I", 100000, 101000);
Function:   Returns smallest sequence object containing the given coordinates and the coordinates relative to that sequence.
Returns :   Clone as string eg "V"; coordinates relative to clone obj
Args    :   Sequence obj (CHROMOSOME or SUPERLINK) and two coordinates within that eg ("CHROMOSOME_I", 3000, 3100)
Non-elegans : Works as expected.

=cut

sub LocateSpan {
  my $self = shift;
  my $chrom = shift;
  my $x = shift;
  my $y = shift;
  
  # if a clone is passed (handles negative coords eg AH6, -500, 12000)
  unless( $self->isa_chromosome($chrom) || # test if $chrom is a chromosome
          $self->isa_superlink($chrom) ) { # test if $chrom is a superlink
    my ($superlink, $chromosome, $seq);
    $seq = $chrom;
    undef $chrom;
    
  SLINKS:
    foreach my $slink (keys %{$self->{SUPERLINK}} ) {
      foreach my $clone (keys %{$self->{SUPERLINK}->{$slink}} ) {
        
        if( "$clone" eq "$seq" ) {
          $chrom = $self->_getChromFromSlink($slink);
          
          # modify the starting coordinate
          my $offset = $self->{$chrom}->{SUPERLINK}->{$slink}->[0] - 1 ; # superlink base coords
          $offset += $self->{SUPERLINK}->{$slink}->{$clone}->[0] - 1; # clone base coords
          
          if (!exists $self->{SENSE}->{$clone} || $self->{SENSE}->{$clone} eq '+') {
            $x += $offset;
            $y += $offset;
          } else {
            my $seqlen = $self->{LENGTH}->{$clone};
            $x = $offset + $seqlen - $x + 1;
            $y = $offset + $seqlen - $y + 1;
            ($x, $y) = ($y, $x);
          }
          last SLINKS;
        }
      }
    }
  }

  if( $self->isa_superlink($chrom) ) { # test if $chrom is a superlink
    my $slink = $chrom;
    $chrom = $self->{SUPERLINK2CHROM}->{$slink};
    $x += $self->{$chrom}->{'SUPERLINK'}->{$slink}->[0] - 1;
    $y += $self->{$chrom}->{'SUPERLINK'}->{$slink}->[0] - 1;
  }
  
  my ($seq, $rel_x, $rel_y);  # coordinates relative to returned seq
   
  # see if the span is within a clone
  my ($clone, $x_slink, $y_slink, $sl_start) = $self->GetCloneFromCoords($chrom, $x, $y);

  if( $clone ) { # span maps within clone
    $seq = $clone;
    my $clone_start = $self->{SUPERLINK}->{$x_slink}->{$clone}->[0];

    $rel_x = $x - $sl_start - $clone_start + 2;  # add 2 -  1 for clone, 1 for slink numbering adjustment
    $rel_y = $y - $sl_start - $clone_start + 2;

    if (exists $self->{SENSE}->{$clone} and $self->{SENSE}->{$clone} eq '-') {
      my $seqlen = $self->{LENGTH}->{$clone};
      $rel_x = $seqlen - $rel_x + 1;
      $rel_y = $seqlen - $rel_y + 1;
    } 

  } else {
    # locate on Slink
    if( $x_slink eq $y_slink ) { # span maps within superlink
      $rel_x = $x - $sl_start + 1;
      $rel_y = $y - $sl_start + 1;
      $seq = $x_slink;
    } else { # spans superlink so goes on to the chromosome
      $rel_x = $x; 
      $rel_y = $y;
      $seq = "$chrom";
    }
  }
  
  return ($seq, $rel_x, $rel_y);
}


=head2 LocateSpanUp

Title   :   LocateSpanUp
Usage   :   my @data = $coords->LocateSpan("AH6", 25, 100);
Function:   Returns the toplevel sequence coords of the given low(er) level piece
Returns :   ($chromosome, $start, $end)
Args    :   Sequence obj (CHROMOSOME or SUPERLINK) and two coordinates within that eg ("CHROMOSOME_I", 3000, 3100)
Non-elegans : Works as expected.

=cut

sub LocateSpanUp {
  my $self = shift;
  my $seqid = shift;
  my $x = shift;
  my $y = shift;
  
  $x = 1 if not defined $x;
  $y = $self->{LENGTH}->{$seqid} if not defined $y;

  if ( $self->isa_superlink( $seqid ) ) {
    my $chrom = $self->{SUPERLINK2CHROM}->{$seqid};

    $x += $self->{$chrom}->{'SUPERLINK'}->{$seqid}->[0] - 1;
    $y += $self->{$chrom}->{'SUPERLINK'}->{$seqid}->[0] - 1;

    $seqid = $chrom;
  } elsif ( not $self->isa_chromosome($seqid) ) {
    my ($superlink, $chromosome, $seq);
    
    SLINKS:
    foreach my $slink (keys %{$self->{SUPERLINK}} ) {
      foreach my $clone (keys %{$self->{SUPERLINK}->{$slink}} ) {
        
        if( "$clone" eq "$seqid" ) {
          my $chrom = $self->_getChromFromSlink($slink);
          
          # modify the starting coordinate
          my $offset = $self->{$chrom}->{SUPERLINK}->{$slink}->[0] - 1 ; # superlink base coords
          $offset += $self->{SUPERLINK}->{$slink}->{$clone}->[0] - 1; # clone base coords
          
          if (!exists $self->{SENSE}->{$clone} || $self->{SENSE}->{$clone} eq '+') {
            $x += $offset;
            $y += $offset;
          } else {
            my $seqlen = $self->{LENGTH}->{$clone};
            $x = $offset + $seqlen - $x + 1;
            $y = $offset + $seqlen - $y + 1;
          }
          $seqid = $chrom;
          last SLINKS;
        }
      }
    }
  }

  return ($seqid, $x, $y);
}


=head2 isa_chromosome

Title   :   isa_chromosome
Usage   :   $coords->isa_chromosome($seq)
Function:   tests if $seq is a chromosome - assumes that the chromosome name has any prefix (e.g. CHROMOSOME_) added
Returns :   true if $seq is a chromosome
Args    :   $seq - name of a sequence object
Non-elegans : Returns true if seq is a supercontig name
=cut

sub isa_chromosome {
  my $self = shift;
  my $seq = shift;
  if (!exists $self->{$seq} || !defined $self->{$seq}) {return 0;}
  return exists $self->{$seq}->{SUPERLINK};
}

=head2 isa_superlink

Title   :   isa_superlink
Usage   :   $coords->isa_superlink($seq)
Function:   tests if $seq is a superlink
Returns :   true if $seq is a superlink
Args    :   $seq - name of a sequence object
Non-elegans : Returns true if seq is a supercontig name (supercontigs are the same as superlinks)
=cut

sub isa_superlink {
  my $self = shift;
  my $seq = shift;
  return exists $self->{SUPERLINK}->{$seq};
}

=head2 isa_clone

Title   :   isa_clone
Usage   :   $coords->isa_clone($seq)
Function:   tests if $seq is a clone
Returns :   true if $seq is a clone
Args    :   $seq - name of a sequence object
Non-elegans : Returns true if seq is a contig
=cut

sub isa_clone {
  my $self = shift;
  my $seq = shift;
  return exists $self->{CLONE2CHROM}->{$seq};
}

=head2 swap

Title   :   swap
Usage   :   $coords->swap(\$x,\$y);
Function:   swaps values of variables passed in
Returns :   nothing
Args    :   references to two coordinates
Non-elegans : Works as expected
=cut

sub swap {
  my $self = shift;
  my $x = shift;
  my $y = shift;
  carp "undefined values passed to swap \n" unless ($$x and $$y);
  my $tmp = $$x;
  $$x = $$y;
  $$y = $tmp;
}

=head2 get_SuperlinkFromClone

Title   :   get_SuperlinkFromClone
Usage   :   my $superlink = $coords->get_SuperlinkFromClone("AH6");
Function:   get the superlink for a given clone
Returns :   superlink for the given clone eq "SUPERLINK_CB_II"
Args    :   clone name as string
Non-elegans : Returns the supercontig name of the input contig

=cut

sub get_SuperlinkFromClone {
  my $self = shift;
  my $clone = shift;
  
  my $superlink = $self->{"CLONE2SUPERLINK"}->{$clone} ? $self->{"CLONE2SUPERLINK"}->{$clone} : undef;
  carp "unknown clone \"$clone\" requesting superlink info" unless $superlink;
  
  return $superlink;
}

=head2 get_Chrom_from_clone

Title   :   get_Chrom_from_clone
Usage   :   my $chromosome = $coords->get_Chrom_from_clone("AH6");
Function:   get chromosome for given clone
Returns :   chromosome for given clone eq "CHROMOSOME_I"
Args    :   clone name as string
Non-elegans : Returns the supercontig name for a given contig

=cut

sub get_Chrom_from_clone {
  my $self = shift;
  my $clone = shift;
  
  my $chrom = $self->{"CLONE2CHROM"}->{$clone} ? $self->{"CLONE2CHROM"}->{$clone} : undef;
  carp "unknown clone \"$clone\" requesting chromosome info" unless $chrom;
  
  return $chrom;
}

=head2 _getChromFromSlink

Title   :   _getChromFromSlink
Usage   :   my $chromosome = $self->_getChromFromSlink("SUPERLINK_CB_I");
Returns :   chromosome for given clone eg "CHROMOSOME_I"
Args    :   superlink as string
Non-elegans : Returns the input supercontig because a dummy superlink is the same size and name as its supercontig

=cut

sub _getChromFromSlink {
  my $self = shift;
  my $sl = shift;
  
  return $self->{SUPERLINK2CHROM}->{"$sl"};
}

=head2 

Title   :   seq_obj_to_chrom
Usage   :   my $chromosome = $coords->seq_obj_to_chrom("$seq")
Function:   find the chromosome when you dont know what type of sequence object youre dealing with
Returns :   Chromosome for given $seq eg "CHROMOSOME_I"
Args    :   any sequence obj and string
Non-elegans : Returns the supercontig name for any (supercontig or contig) input object name

=cut

sub seq_obj_to_chrom {
  my ($self, $seq) = @_;
  my $chromosome = 0;
  
  # find which chromosome passed seq obj is on
  if ( $self->isa_chromosome($seq) ) {
    $chromosome = $seq;
  } elsif ( $self->isa_superlink($seq) ) {
    $chromosome = $self->_getChromFromSlink($seq);
  } elsif ( $self->isa_clone($seq) ) {
    $chromosome = $self->get_Chrom_from_clone($seq);
  } else {
    carp "unknown object '$seq' passed to Coords_converter seq_obj_to_chrom()\n";
  }

  return $chromosome;

}

=head2 Coords_2chrom_coords

Title   :   Coords_2chrom_coords
Usage   :   my ($chrom, $coord) = $coords->Coords_2chrom_coords("AH6", 2343);
Function:   convert non-chromosomal coordinates to chromosome relative coords
Returns :   chromosome as string, coordinate relative to that chromosome
Args    :   any sequence obj as string
		coordinate relative to seq obj
Non-elegans : Returns the supercontig and coord in the supercontig of a position in a contig

=cut

sub Coords_2chrom_coords {
  my ($self, $seq, $coord) = @_;
  
  my ($chrom, $offset) = $self->CloneOffset($seq);
  
  if (!exists $self->{SENSE}->{$seq} || $self->{SENSE}->{$seq} eq '+') {
    $coord += $offset - 1;
  } else {			# some contigs in briggsae are in the reverse orientation
    my $seqlen = $self->{LENGTH}->{$seq};
    $coord = $offset + $seqlen - $coord;
  }
  
  return ($chrom, $coord);
}

=head2 Chrom2slink_coords

Title   :   Chrom2slink_coords
Usage   :   my ($chrom, $coord) = $coords->Chrom2slink_coords("III", 2343, 2409);
Function:   convert chromosomal coordinates to slink relative coords
Returns :   slink as string, coordinate1 relative to that slink, coordinate2 relative to that slink, 
Args    :   chromosome as string
	    chromosome coordinate1
	    chromosome coordinate2
Non-elegans : Not sure - not tested this

=cut

sub Chrom2slink_coords {
  my ($self, $chrom, $coord1, $coord2) = @_;
  
  if ($self->{chromosome_prefix} eq "CHROMOSOME_" &&  $chrom !~ /CHROMOSOME_/) {
    $chrom = "CHROMOSOME_$chrom";
  }
  
  foreach my $slink ( keys %{$self->{"$chrom"}->{'SUPERLINK'}} ) {
    if($self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0] <= $coord1 and
       $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[1] >= $coord1
      ) {
	my $sl_start = $self->{"$chrom"}->{'SUPERLINK'}->{$slink}->[0];
	$coord1 = $coord1 - $sl_start + 1;	# so $coord1 is now the coordinate in the $parent superlink
	$coord2 = $coord2 - $sl_start + 1;	# so $coord2 is now the coordinate in the $parent superlink
	return ($slink, $coord1, $coord2);
    }
  }
}

=head2 Superlink_length

Title   :   Superlink_length
Usage   :   my $sl_length = $coords->Superlink_length("SUPERLINK_RWXL");
Function:   returns the length of the input superlink or clone or chromosome
Returns :   length of superlink
Args    :   name of superlink
Non-elegans : returns the length of the input contig or supercontig

=cut

sub Superlink_length {
  my ($self, $seq) = @_;
  my $len = undef;
  
  $len = $self->{LENGTH}->{$seq};
  
  return $len;
}

=head2 get_clones_in_superlink

Title   :   get_clones_in_superlink
Usage   :   $coords->get_clones_in_superlink($superlink)
Function:   returns a list of clones and their positions in the superlink
Returns :   returns a list of [clone, start_pos, end_pos] ordered by start_pos, sense ('+' or '-')
Args    :   name of superlink
Non-elegans : returns a list of contigs and their positions in the supercontig
=cut

sub get_clones_in_superlink {
  my ($self, $superlink) = @_;
  my @result;
  my @clones = keys %{$self->{SUPERLINK}->{$superlink}};
  if (!@clones) {croak "$superlink is not a valid superlink\n";}
  my @sorted_clones = sort {$self->{SUPERLINK}->{$superlink}->{$a}->[1] <=> $self->{SUPERLINK}->{$superlink}->{$b}->[1]} @clones;
  foreach my $clone (@sorted_clones) {
    push @result, [($clone, $self->{SUPERLINK}->{$superlink}->{$clone}->[0], $self->{SUPERLINK}->{$superlink}->{$clone}->[1]), $self->{'SENSE'}->{$clone}];
  }
  return @result;
}

=head2 get_superlinks_in_chromosome

Title   :   get_superlinks_in_chromosome
Usage   :   $coords->get_superlinks_in_chromosome($chromosome)
Function:   returns a list of superlinks and their positions in the chromosome
Returns :   returns a list of [superlink, start_pos, end_pos] ordered by start_pos
Args    :   name of chromosome
Non-elegans : returns a list of supercontigs and their positions in the supercontig
		i.e. simply returns the name of the input supercontig and 1 and the length of the supercontig
=cut

sub get_superlinks_in_chromosome {
  my ($self, $chrom) = @_;
  my @result;
  my @sl = keys %{$self->{$chrom}->{SUPERLINK}};
  if (!@sl) {croak "$chrom is not a valid chromosome\n";}
  my @sorted_sl = sort {$self->{$chrom}->{SUPERLINK}->{$a}->[1] <=> $self->{$chrom}->{SUPERLINK}->{$b}->[1]} @sl;
  foreach my $sl (@sorted_sl) {
    push @result, [($sl, $self->{$chrom}->{SUPERLINK}->{$sl}->[0], $self->{$chrom}->{SUPERLINK}->{$sl}->[1])];
  }
  return @result;
}


return 1;

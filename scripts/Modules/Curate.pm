# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-06-09 14:32:44 $      

=pod

=head1 NAME

 Curate.pm

=head1 SYNOPSIS

  use Curate;

  # initialise
  my $Curate = Curate->new($ace, $seq_obj, $wormbase, $log);

  # force the command to succeed, even if it fails normal sanity checking
  $self->set_force(1);
  # turn force off again
  $self->set_force(0);

  # run one command line
  $self->cmd($command_line);
  
  # get the next available sequence name 
  $self->next_name();
  $self->next_name($clone); # clone is required in elegans

  # get the CGC-name, if any
  $cgcname = $self->cgcname($seqname);

  # Making a simple prediction change
  $self->replace($class, $old_seqname, $new_seqname); # replace old sequence name with new sequence name

  # Adding a new Isoform
  $self->new_isoform($class, $existing_seqname, $new_seqname); # add new sequence name to the isoforms of the existing sequence name

  # Splitting a gene - it is your responsibility to ensure that any new structure containing the domains for a CGC-named gene are first in the list of new seqnames
  @new_seqnames = $self->split($class, $existing_seqname, $new1_seqname, $new2_seqname, ...); # split existing sequence name to two or more new structures

  # Merging a gene - CGC-names are checked to ensure that there are no CGC-names in the genes to die
  $self->merge($class, $live_seqname, $die1_seqname, $die2_seqname, ..., $new_seqname); merge genes to die into the gene to live using the structure of $new

  # Changing a CDS to a Pseudogene
  $self->CDS2Pseudogene($cds_seqname, $pseudogene_seqname);

  # Changing a CDS to a Transcript
  $self->CDS2Transcript($cds_seqname, $transcript_seqame, $type)

  # Changing a CDS to a Transposon_CDS
  $self->CDS2Transposon_CDS($cds_seqname, $tran_cds_seqname)

  # Deleting a gene, cds, isoform, whatever
  $self->delete($class, $seqname);

  # Making a new locus: CDS, Pseudogene, Transposon, non-coding Transcript
  $created_seqname = $self->create($class, $seqname)

  # Just make a new history object
  $self->just_make_history($class, $seqname);

  # Report any Gene_name for an object
  $self->check_gene_name($class, $seqname);




=head1 DESCRIPTION

Methods for curating gene models.

=head1 CONTACT

Gary gw3@ebi.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _


=cut


package Curate;

use strict;                                      
use Getopt::Long;
use Carp;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;
use DBI;
use DBD::mysql;
use FileHandle;
use Coords_converter;


=head2 

    Title   :   new
    Usage   :   my $Curate = Curate->new($out, $ace, $seq_obj, $wormbase, $log);
    Function:   initialises the object
    Returns :   RNASeq object;
    Args    :   

=cut


sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  
  $self->{'out'} = shift;       # output file handle
  $self->{'ace'} = shift;       # acedb db handle
  $self->{'seq_obj'} = shift;   # Sequence_extract object
  $self->{'wormbase'} = shift;  # wormbase object
  $self->{'log'} = shift;       # log object
  
  my $person;
  my $user = $ENV{'USER'};
  $self->{user} = $user;
  if ($user eq 'pad') {
    $person = 'WBPerson1983';
  } elsif ($user eq 'gw3') {
    $person = 'WBPerson4025';
  } elsif ($user eq 'mh6') {
    $person = 'WBPerson4055';
  } else {
    die "ERROR Unknown user: $user\nPlease set up WBPerson for $user in Curate.pm\n";
  }
  $self->{'person'} = $person;
  
  $self->{Isoform_details} = ();
  
  $self->{species} = $self->{wormbase}->species;

  if ($self->{species} eq 'elegans' || $self->{species} eq 'sratti') {
    $self->{next_cds_id} = ();
  } else {
    $self->{next_cds_id} = '';
  }

  $self->{force} = 0; # sanity checks fail normally, unless this is set to 1
  $self->{force_flag} = 0; # this is set to 1 if an error that requires forcing is found

  # +++ want to see if the OUT file has been moved to a directory for old OUT files
  # +++ i.e. there is no OUT file in existence
  # +++ otherwise we could get in trouble when we check to see what
  # +++ the next available sequence name to use is if we have not yet read
  # +++ in the last OUT file which may have assigned a sequence name
  # +++ already
  
  return $self;
}

######################################
# flag to force the command to succeed even if it fails the normal checks
sub set_force {
  my ($self, $force) = @_;
  $self->{force} = $force;
}

######################################
# check for -force and set the 'force_flag' if -force is not set
# extra checks: 
# 'elegans' - only raise the flag if species is elegans
sub force {
  my ($self, $extra_check) = @_;

  if (defined $extra_check && $extra_check eq 'elegans' && $self->{species} ne 'elegans') {return}
  if (!$self->{force}) {
    $self->{force_flag} = 1; # raise the flag
  }  


}

######################################
# check to see if the force_flag has been set
sub check_force_flag {
  my ($self) = @_;

  if ($self->{force_flag}) {
    die "\nUse -force if this is correct\n\n";
  }
}

######################################
# Making a simple prediction change
# $self->replace_cmd($class, $old_seqname, $new_seqname); # replace one old sequence name with one new sequence name

sub replace_cmd {
  my ($self, $class, $old_seqname, @new_seqnames) = @_;

  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The object to be replaced was not specified.\n";}
  if (@new_seqnames == ()) {push @new_seqnames, 'temp_gene_1'}
  $self->cmd("REPLACE CLASS $class EXISTING $old_seqname MODEL $new_seqnames[0]");
}

######################################
# Adding a new Isoform to  an existing seqname
# $self->isoform_cmd($class, $existing_seqname, $new_seqname); # add new sequence name(s) to the isoforms of the existing sequence name

sub isoform_cmd {
  my ($self, $class, $old_seqname, @new_seqnames) = @_;

  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  if (@new_seqnames == ()) {push @new_seqnames, 'temp_gene_1'}
  $self->cmd("NEW_ISOFORM CLASS $class EXISTING $old_seqname MODELS @new_seqnames");
}
######################################
# create a new gene
# $self->newgene_cmd($class, $new_seqname);

sub newgene_cmd {
  my ($self, $class, @new_seqnames) = @_;

  if (@new_seqnames == ()) {push @new_seqnames, 'temp_gene_1'}
  $self->cmd("NEW_GENE CLASS $class MODELS @new_seqnames");
}
######################################
# Splitting a gene - it is your responsibility to ensure that any new
# structure containing the domains for a CGC-named gene are first in
# the list of new seqnames
# if you wish to add several isoforms in a new locus, then split it first, then add isoforms to the locus

# @new_seqnames = $self->split_cmd($class, $existing_seqname, ($new1_seqname, $new2_sqname, ...); # split existing sequence name to two or more new structures

sub split_cmd {
  my ($self, $class, $old_seqname, @new_seqnames) = @_;
  
  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object to be split was not specified.\n";}
  my $splitstr = '';
  my $found=0;
  foreach my $model (@new_seqnames) {
    if (!defined $model || $model eq '') {
      if ($found) {die "ERROR In split($class, $old_seqname, new_seqnames) have got two or more undefined new_seqnames - they can't all be 'temp_gene_1'!\n"}
      $model = 'temp_gene_1';
      $found = 1;
    }
    $splitstr .= "SPLITGROUP $model ";
  }
  $self->cmd("SPLIT CLASS $class EXISTING $old_seqname $splitstr");
}
######################################
# Merging a gene - CGC-names are checked to ensure that there are no CGC-names in the genes to die
# existing seqnames - array ref of seqnames to merge (the first one will Live, the others will Die - a check if made for CGC-names and Live is reassigned if one is found)
# new_seqnames - array of isoforms of the new locus
# $self->merge_cmd($class, [$die1_seqname, $die2_seqname, ...,], [$new_isoform1, $new_isoform2]); merge genes to die into the gene to live using the structure of $new

sub merge_cmd {
  my ($self, $class, $existing_seqnames, @new_seqnames) = @_;
  
  if (!defined $existing_seqnames || $existing_seqnames eq '') {die "ERROR The existing sequence names to be merged were not specified.\n";}
  if (@new_seqnames == ()) {push @new_seqnames, 'temp_gene_1'}
  $self->cmd("MERGE CLASS $class DIE @{$existing_seqnames} MODELS @new_seqnames");
}

######################################
# Changing something to a CDS
# $self->to_cds_cmd($class, $cds_seqname);
sub to_cds_cmd {
  my ($self, $class, $old_seqname) = @_;
  
  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("CHANGE_CLASS CLASS $class NEWCLASS CDS EXISTING $old_seqname");
}


######################################
# Changing something to a Pseudogene
# $self->to_pseudogene_cmd($class, $cds_seqname)
sub to_pseudogene_cmd {
  my ($self, $class, $old_seqname) = @_;
  
  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("CHANGE_CLASS CLASS $class NEWCLASS Pseudogene EXISTING $old_seqname");
}


######################################
# Changing something to a Transcript
# $self->to_transcript_cmd($class, $cds_seqname, $type)
sub to_transcript_cmd {
  my ($self, $class, $old_seqname) = @_;
  
  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("CHANGE_CLASS CLASS $class NEWCLASS Transcript EXISTING $old_seqname");
}


######################################
# Changing something to a Transposon
# $self->make_transposon($class, $cds_familyname)
sub make_transposon_cmd {
  my ($self, $class, $old_seqname, $new_familyname) = @_;
  
  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  if (!defined $new_familyname || $new_familyname eq '') {die "ERROR The new Transposon Family was not specified.\n";}
  $self->cmd("MAKE_TRANSPOSON CLASS $class EXISTING $old_seqname FAMILY $new_familyname");
}


######################################
# Deleting a gene, cds, isoform, whatever
# $self->delete_cmd($class, $old_seqname);
sub delete_cmd {
die "ERROR NOT IMPLEMENTED YET\n";
}

######################################
# set the 'Last_reviewed' tag in the object
# $self->last_reviewed_cmd($class, $old_seqname);
sub last_reviewed_cmd {
  my ($self, $class, $old_seqname) = @_;

  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("LAST_REVIEWED CLASS $class EXISTING $old_seqname");
}
######################################
# Just make a history object
# $self->just_make_history_cmd($class, $old_seqname);
sub just_make_history_cmd {
  my ($self, $class, $old_seqname) = @_;

  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("JUST_MAKE_HISTORY CLASS $class EXISTING $old_seqname");
}

######################################
# Check Gene_name
# $self->check_gene_name_cmd($class, $old_seqname);
sub check_gene_name_cmd {
  my ($self, $class, $old_seqname) = @_;

  if (!defined $old_seqname || $old_seqname eq '') {die "ERROR The existing object was not specified.\n";}
  $self->cmd("CHECK_GENE_NAME CLASS $class EXISTING $old_seqname");
}

######################################
# execute one command line
# $self->cmd($command_line);
sub cmd {
  my ($self, $line) = @_;
  chomp($line);
  my @line = split /\s+/, $line;

  my $fh = $self->{out};
  print $fh "\n\n//\n// $line\n//\n\n";
  
  my %line = $self->parse_line(@line);
  
  # DELETE
  if ($line[0] =~ /DELETE/) {
    $self->delete_obj(%line);
    
    # KEEP
  } elsif ($line[0] =~ /KEEP/) {
    $self->keep_cds(%line);
    
    # MAKE_TRANSPOSON
  } elsif ($line[0] =~ /MAKE_TRANSPOSON/) {
    $self->make_transposon(%line);
    
    # MERGE
  } elsif ($line[0] =~ /MERGE/) {
    $self->merge_gene(%line);
    
    # NEW_GENE
  } elsif ($line[0] =~ /NEW_GENE/) {
    $self->new_gene(%line);
    
    # NEW_ISOFORM
  } elsif ($line[0] =~ /NEW_ISOFORM/) {
    $self->new_isoform(%line);
    
    # REPLACE
  } elsif ($line[0] =~ /REPLACE/) {
    $self->replace_cds(%line);
    
    # SPLIT
  } elsif ($line[0] =~ /SPLIT/) {
    $self->split_gene(%line);

    # CHANGE_CLASS
  } elsif ($line[0] =~ /CHANGE_CLASS/) {
    $self->change_class(%line);

    # LAST_REVIEWED
  } elsif ($line[0] =~ /LAST_REVIEWED/) {
    $self->last_reviewed(%line);

    # MAKE_HISTORY
  } elsif ($line[0] =~ /JUST_MAKE_HISTORY/) {
    $self->just_make_history(%line);

    # CHECK_GENE_NAME
  } elsif ($line[0] =~ /CHECK_GENE_NAME/) {
    $self->check_gene_name(%line);
  }
}

######################################
# make the FMAP display a location
# returns 0 if OK
sub goto_location {
  my ($self, $seqobj, $x, $y, $sense, $zoomout) = @_;

  if (!defined $zoomout) {
    $zoomout = 500;
  }
  $x -= $zoomout;
  $y += $zoomout;
  
  my $reversed_command = "";
  if (defined $sense && $sense eq '-') {
    $reversed_command = "; seqactions -rev_comp";
  }
  my $command = "xremote -remote 'gif seqget $seqobj -coords $x $y; seqdisplay $reversed_command'";
  
  my $return_status = system($command);
  if ( ( $return_status >> 8 ) != 0 ) {
    carp "X11 connection appears to be lost";
  }
  return $return_status;
}


######################################
# parse an ace file into the FMAP
# returns 0 if OK

sub load_ace {
  my ($self, $output) = @_;

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    carp "X11 connection appears to be lost";
  } else {
    system("mv -f $output $output.old");
  }
  return $return_status
}


######################################
# convert a Gene's sequence name to the Gene ID
# SeqName2Gene

sub SeqName2Gene {
  my ($self, $seqname) = @_;
  if ($seqname =~ /^WBGene\d+/) {return $seqname}
  
  my $genename_obj = $self->{ace}->fetch(Gene_name => "$seqname");
  if (defined $genename_obj) {
    return $genename_obj->Sequence_name_for->name;
  } else {
    # currently non-elegans species Gene's do not have Sequence_name data
    my @cdses = $self->{ace}->fetch(-query => "find CDS ${seqname}* WHERE Method = \"curated\"");
    my @pseuds = $self->{ace}->fetch(-query=>"find Pseudogene ${seqname}*");
    my @trans = $self->{ace}->fetch(-query=>"find Transcript ${seqname}* WHERE Method != \"Coding_transcript\"");
    foreach my $thing (@cdses, @pseuds, @trans) {
      if ($thing =~ /^${seqname}$/ || $thing =~ /^${seqname}[a-z]$/) {
	my $gene = $thing->Gene;
	if (!defined $gene) {die "The Gene ID is not defined for $seqname\n"}
	return $gene;
      }
    }
    
  }
  return undef;
}

######################################
# convert a Gene ID to the Gene's sequence name
# Gene2SeqName

sub Gene2SeqName {
  my ($self, $gene) = @_;
  if ($gene !~ /^WBGene\d+/) {return $gene}

  my $gene_obj = $self->{ace}->fetch(Gene => "$gene");
  if (!defined $gene_obj) {die "ERROR In Gene2SeqName: $gene was not found in the database\n"}
  my $seq_obj = $gene_obj->Sequence_name;
  if (!defined $gene_obj) {die "ERROR In Gene2SeqName: the tag Sequence_name was not populated for Gene $gene in the database\n"}
  my $seqname = $seq_obj->name;
  if (defined $seqname) {
    return $seqname;
  } else {

    # currently non-elegans species Gene's do not have Sequence_name data

    my @structures;
    my @corr_cds = $gene_obj->Corresponding_CDS;
    foreach my $struc (@corr_cds) {
      my $method = $struc->Method;
      if ($method eq 'curated') {push @structures, $struc}
    }

    # include Non_coding_transcripts
    my @corr_trans = $gene_obj->Corresponding_transcript;
    foreach my $tran (@corr_trans) {
      my $method = $tran->Method;
      if ($method ne 'Coding_transcript') {push @structures, $tran}
    }

    # this might be a Pseudogene
    my @corr_ps = $gene_obj->Corresponding_pseudogene;
    foreach my $pseud (@corr_ps) {
      my $method = $pseud->Method;
      if ($method eq 'Pseudogene') {push @structures, $pseud}
    }


    my $isoform_letter;
    foreach my $corr_cds (@structures) {
      my $letter;
      ($seqname, $letter) = $self->GeneSeqname($corr_cds);
      if (defined $seqname) {return $seqname}
    }

  }

}

######################################
# get the next available sequence name 
# $self->next_name($clone); # clone is required in elegans and sratti - so we get this from the model data
sub Next_CDS_ID {
  my ($self, $clone)  = @_; # $clone is required to get the clone name (in elegans) or the toplevel name (in sratti) or in other species it should be undef

  my %sratti_prefix = (
		       SRAE_chr1 => 'SRAE_1',
		       SRAE_chr2 => 'SRAE_2',
		       SRAE_chrX => 'SRAE_X',
		       SRAE_chrM => 'SRAE_M',
		      );

  if ($self->{species} eq 'elegans') {
    if (!defined $clone) {die "ERROR \$clone is not defined and this species is elegans\n";}
    if (!exists $self->{next_cds_id}{$clone}) {$self->Find_Next_CDS_ID($clone)}
    my $next_id = $self->{next_cds_id}{$clone};
    my ($number) = ($self->{next_cds_id}{$clone} =~ /^${clone}\.(\d+)$/);
    $number++;
    $self->{next_cds_id}{$clone} = $clone . '.' . $number;
    return $next_id;

  } elsif ($self->{species} eq 'sratti') {
    if (!defined $clone) {die "ERROR \$clone is not defined and this species is sratti\n";}
    if (!exists $self->{next_cds_id}{$clone}) {$self->Find_Next_CDS_ID($clone)}
    my $sratti_clone_start = substr($clone, 9);
    my $prefix;
    if (exists $sratti_prefix{$sratti_clone_start}) {
      $prefix = $sratti_prefix{$sratti_clone_start};
    } else {
      $prefix = 'SRAE_0';
    }
    return $self->{next_cds_id}{$prefix} += 100;

  } else {
    if ($self->{next_cds_id} eq '') {$self->Find_Next_CDS_ID()}
    return $self->{next_cds_id}++;
  }

}

######################################
# Populate the store for the next available sequence name ID
sub Find_Next_CDS_ID {
  my ($self, $clone) = @_; # $clone is either the clone name (in elegans) or the toplevel name (in sratti) or undef

  my %prefix = (
		briggsae    => 'CBG',
		remanei     => 'CRE',
		brenneri    => 'CBN',
		japonica    => 'CJA',
		pristiochus => 'PPA',
		brugia      => 'Bm',
		ovolvulus   => 'OVOC',
		sratti      => 'SRAE_', # followed by a digit or a 'X' then the number ending '00'
	       );

  my %sratti_prefix = (
		       SRAE_chr1 => 'SRAE_1',
		       SRAE_chr2 => 'SRAE_2',
		       SRAE_chrX => 'SRAE_X',
		      );

  my $maxnumber = 0;
  if ($self->{species} eq 'elegans') {
    my @cdses = $self->{ace}->fetch(-query => "find CDS ${clone}.* WHERE Method = \"curated\" OR Method = \"Transposon_CDS\"");
    my @pseuds = $self->{ace}->fetch(-query=>"find Pseudogene ${clone}.*");
    my @trans = $self->{ace}->fetch(-query=>"find Transcript ${clone}.*");
    foreach my $thing (@cdses, @pseuds, @trans) {
      if ($thing =~ /^${clone}\.(\d+)$/ || $thing =~ /^${clone}\.(\d+)[a-z]$/) {
	my $number = $1;
	if ($number > $maxnumber) {$maxnumber = $number}
      }
    }
    $maxnumber++;
    # check it has never been used as a CDS
    while ($self->check_history_seqname("${clone}.${maxnumber}")) {$maxnumber++}
    $self->{next_cds_id}{$clone} = $clone . '.' . $maxnumber;


  } elsif ($self->{species} eq 'sratti') {
    my $sratti_clone_start = substr($clone, 9);
    my $prefix;
    if (exists $sratti_prefix{$sratti_clone_start}) {
      $prefix = $sratti_prefix{$sratti_clone_start};
    } else {
      $prefix = 'SRAE_0';
    }
    my @cdses = $self->{ace}->fetch(-query => "find CDS ${prefix}* WHERE Method = \"curated\"");
    my @pseuds = $self->{ace}->fetch(-query=>"find Pseudogene ${prefix}*");
    my @trans = $self->{ace}->fetch(-query=>"find Transcript ${prefix}*");
    foreach my $thing (@cdses, @pseuds, @trans) {
      if ($thing =~ /^${prefix}(\d+)$/ || $thing =~ /^${prefix}(\d+)[a-z]$/) {
	my $number = $1;
	if ($number > $maxnumber) {$maxnumber = $number}
      }
    }
    $maxnumber += 100;
    # check it has never been used as a CDS
    while ($self->check_history_seqname("${prefix}${maxnumber}")) {$maxnumber += 100}
    $self->{next_cds_id}{$prefix} = $prefix . $maxnumber;


  } else {
    my $prefix = $prefix{$self->{species}};
    my @cdses = $self->{ace}->fetch(-query => "find CDS ${prefix}* WHERE Method = \"curated\"");
    my @pseuds = $self->{ace}->fetch(-query=>"find Pseudogene ${prefix}*");
    my @trans = $self->{ace}->fetch(-query=>"find Transcript ${prefix}*");
    foreach my $thing (@cdses, @pseuds, @trans) {
      if ($thing =~ /^${prefix}(\d+)$/ || $thing =~ /^${prefix}(\d+)[a-z]$/) {
	my $number = $1;
	if ($number > $maxnumber) {$maxnumber = $number}
      }
    }
    $maxnumber++;
    # check it has never been used as a CDS
    while ($self->check_history_seqname("${prefix}${maxnumber}")) {$maxnumber++}
    $self->{next_cds_id} = $prefix . $maxnumber;
  }

}


######################################
# check the wormpep history file to see if the sequence name has been used in the past for a CDS
# return 1 if it is found
sub check_history_seqname {
  my ($self, $seqname) = @_;
  my $result=0;

  my $version = $self->get_history_version($self->{wormbase}->database('current'));
  my $pepname = $self->{wormbase}->pepdir_prefix . 'pep';
  if (!-e "/nfs/wormpub/BUILD/WORMPEP/${pepname}$version/${pepname}.history${version}") {$version--;}
  open (PEP, "< /nfs/wormpub/BUILD/WORMPEP/${pepname}$version/${pepname}.history${version}") || die "ERROR Can't find /nfs/wormpub/BUILD/WORMPEP/${pepname}$version/${pepname}.history${version}\n";
  while (my $line = <PEP>) {
    my @f = split /\s+/, $line;
    if ($f[0] =~ /^${seqname}[a-z]*$/) {return 1;}
  }
  close(PEP);
  return $result;
}

######################################
######################################
######################################
######################################
######################################
######################################
######################################

######################################
# change() from original program
# reads a file of commands and executes them

sub  change {
  my ($self) = @_;
  
  my $batch_split_file = "batch_split.dat";
  my $batch_merge_file = "batch_merge.dat";
  my $batch_new_file = "batch_new.dat";
  my $batch_kill_file = "batch_kill.dat";
  my $cds_to_isoforms = "cds_to_isoforms.dat";
  open (BSPLIT, ">$batch_split_file") || die "ERROR Cant open file\n";
  open (BMERGE, ">$batch_merge_file") || die "ERROR Cant open file\n";
  open (BNEW, ">$batch_new_file") || die "ERROR Cant open file\n";
  open (CDS2ISO, ">$cds_to_isoforms") || die "ERROR Cant open file $cds_to_isoforms\n";
  
  my @commands = $self->read_commands();
  $self->execute_commands(@commands);
  
  close(BSPLIT);
  close(BMERGE);
  close(BNEW);
  close(CDS2ISO);
  
  print "\nNow sort out the few isoforms that need to be deleted and un-isoform their siblings if needed\n";
  print "\n";
  
  print "\nNow run:\n";
  print "batch_merge.pl -user gw3 -pass gw3 -ns -debug gw3 -file batch_merge.dat (-test)\n";
  print "output file: /nfs/wormpub/DATABASES/geneace/NAMEDB_Files/batch_merge.ace\n";
  print "No association between the merged gene and the CDS objects needs to be made - it is already in the curation db\n";
  print "\n";
  
  print "batch_split.pl -user gw3 -pass gw3 -ns -debug gw3 -file batch_split.dat (-test)\n";
  print "add the Gene tag by doing some editing of the ace output database/NAMEDB_Files/batch_split.ace\n";
  print "output file: /nfs/wormpub/DATABASES/geneace/NAMEDB_Files/batch_split.ace\n";
  print "\n";
  
  print "rm new_sneak_file\n";
  print "newgene.pl -input batch_new.dat -pseudoace new_sneak_file -species $self->{wormbase}->species -bio CDS -who WBPerson4025 -user gw3 -password gw3 -namedb -debug gw3 (-test)\n";
  print "output is /nfs/wormpub/DATABASES/geneace/NAMEDB_Files/newgene_input.ace\n";
  print "the new_sneak_file gets appended to!\n";
  print "add the Gene tag by reading in the new_sneak_file output of newgene.pl\n";
  print "The sneak file looks like:\n";
  print "\n";
  print "Object : CBG29763\n";
  print "Gene WBGene00270325\n";
  print "\n";
  
  print "\n";
  print "Then read the output files into geneace and briggsae_curation\n";
  print "\n";
  
}

######################################
# read command list
sub read_commands {
  my ($self) = @_;
  
  my @commands;
  
  my $commands = "commands.dat";
  open (COMM, "<$commands") || die "ERROR Can't read $commands";
  while (my $line = <COMM>) {
    push @commands, $line
  }
  close (COMM);
  return @commands;
}

######################################
# execute multiple commands
sub execute_commands {
  
  my ($self, @commands) = @_;
  
  foreach my $line (@commands) {
    
    $self->cmd($line);
  }
}

######################################
# get the data from the line
sub parse_line {
  my ($self, @line) = @_;
  my %line;
  
  my %keyword = (
		 'MODEL' => 1,
		 'MODELS' => 1,
		 'SPLITGROUP' => 1,   # this is for specifying groups of models
		 'EXISTING' => 1,
		 'DIE' => 1,          # this is the only circumstance when there should be more than one existing gene in a command (so this is 'EXISTINGS')
		 'CLASS' => 1,
		 'NEWCLASS' => 1,
		 'FAMILY' => 1,
		);

  my $key;
  my $splitgroup = -1;
  foreach my $word (@line) {
    if (keys %line == 0) { # the first word of the command is the VERB
      $line{VERB} = $word
    } else {
      if (exists $keyword{$word}) { # if word is a keyword
	$key = $word;
	if ($key eq 'SPLITGROUP') {$splitgroup++}
      } else {
	if ($key eq 'MODELS' || $key eq 'DIE') {
	  push @{$line{$key}}, $word;

	} elsif ($key eq 'SPLITGROUP') {
	  push @{$line{$key}[$splitgroup]}, $word;

	} else {
	  $line{$key} = $word;
	}
      }
    }
  }

  return %line;
}
######################################
######################################
######################################
######################################
######################################
# KEEP CLASS CDS CBG10381 MODEL briggsae_GG1698|c0_g1_i1-cb25.fpc2310b-1b.extended
# KEEP/REPLACE and NEW_ISOFORM on a gene should occur in that order

sub keep_cds {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my $model = $line{MODEL};

  if (defined $model && $model ne '') {
    $self->Add_remark($class, $existing, "The structure of this $class has been confirmed by the structure '$model'.");
  } else {
    $self->Add_remark($class, $existing, "The structure of this $class has been reviewed by a curator and deemed satisfactory.");
  }
  $self->Last_reviewed($class, $existing);

}
######################################
# REPLACE CLASS <CLASS> EXISTING <EXISTING-SEQ-NAME> MODEL <NEW-SEQ-NAME>
# REPLACE CLASS CDS EXISTING CBG02835 MODEL briggsae_GG6414|c50_g3_i1-chrII-1
# KEEP/REPLACE and NEW_ISOFORM on a gene should occur in that order

sub replace_cds {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my $model = $line{MODEL};

  if ($self->structure_comparison_sanity_check($class, $existing, $model)) {print "WARNING $existing and $model have the same structure.\n"; ; $self->force()}

  $self->structure_update_sanity_check($class, $existing, $model);

  my $history = $self->Make_history($class, $existing);

  $self->Delete_existing_structure($class, $existing);

  $self->Update_CDS_structure($class, $existing, $model);
  
  $self->Add_remark($class, $existing, 'This $class has been updated based on the published alternate intron splicing.' );
#  $self->Add_remark($class, $existing, "This $class has been updated based on the RNAseq data.");
  $self->Last_reviewed($class, $existing);

  $self->Print_CDS_to_Isoform($existing, $existing);

}
######################################
# add a new isoform to an existing structure
# CLASS <CLASS> EXISTING <EXISTING-SEQ-NAME> MODELS <NEW-SEQ-NAME>
# NEW_ISOFORM CLASS CDS EXISTING CBG02835 MODELS briggsae_GG5841|c0_g2_i1-chrII-1
# KEEP/REPLACE and NEW_ISOFORM on a gene should occur in that order

sub new_isoform {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my @models = @{$line{MODELS}};

  my $gene = $self->SeqName2Gene($existing);

  foreach my $model (@models) {

    $self->structure_update_sanity_check($class, $existing, $model);

    # get CDS ID for Gene
    my $update = 1;
    my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($gene, $update);

    if ($single) {
      if ($self->structure_comparison_sanity_check($class, $existing, $model)) {print "WARNING $existing and $model have the same structure.\n"; $self->force()}
      $self->Make_CDS_into_isoform($class, $cds, 'a');
      $self->Print_CDS_to_Isoform($cds, $cds.'a');
      $self->Last_reviewed($class, $cds.'a');
    }
  
    # check existing isoforms to see if we are trying to make a duplicate structure
    pop @used_letters; # get rid of last letter as we have not made it yet and anyway don't wish to compare the new structure to it.
    foreach my $iso (@used_letters) {
      my $iso_cds = $cds . $iso;
      if ($self->structure_comparison_sanity_check($class, $iso_cds, $model)) {print "WARNING $iso_cds and $model have the same structure.\n"; $self->force()}
    }

    my $cdsl = $cds; # CDS name with letter
    $cdsl .= $isoform_letter;
    
    $self->Update_CDS_structure($class, $cdsl, $model);
    
    $self->Add_stub($class, $cdsl);
    
    $self->Add_Gene($class, $cdsl, $gene);
    
    #$self->Add_remark($class, $cdsl, "This structure has been created based on the RNASeq data.");
    $self->Add_remark($class, $cdsl, 'This structure has been created based on the published alternate intron splicing.');

    $self->Last_reviewed($class, $cdsl);
    
    $self->Add_isoform_tag($class, $cdsl);
    
    $self->Print_CDS_to_Isoform($cds, $cdsl);
  }
}
######################################
# DELETE GENE WBGene00032898 CDS CBG11852a

sub delete_obj {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my $history = $self->Make_history($class, $existing);

  $self->Add_remark($class, $history, "This structure does not appear to be supported by the available evidence. It has been retired.");
  $self->Delete($class, $existing);

  # Check to see if there are other CDS isoforms on this Gene - if not we must manually kill the Gene - for only 6 instances it is not worth putting in the effort to automate it
  # may have to Change Isoform tag in other isoforms
  print "Have Deleted $existing - Please Kill $existing in the Nameserver and/or check on the Isoform status of other CDS on this gene - maybe there is just one isoform that needs an Isoform letter and tag removed\n";

}
######################################
# MAKE_TRANSPOSON CLASS CDS EXISTING WBGene00032898 FAMILY LINE2A_CB
# convert the CDS to Transposon_CDS and make a new Transposon object with the same span

sub make_transposon {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my $family = $line{FAMILY};

  my $gene = $self->SeqName2Gene($existing);

  # check to see if the existing gene has isoforms - die if so
  my $update = 0; # don't update
  my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($gene, $update);
  if (!$single) {die "ERROR Can't make a Transposon object from $existing because it has isoforms - are you sure this is correct?\n"}

  $self->Change_to_Transposon_CDS($class, $cds);
  $self->Last_reviewed($class, $cds);
  $self->Add_remark($class, $cds, "Converted this Transposon as it has a strong similarity to a retrotransposase ($family).");

  # find ID of next available Transposon
  my $transposon_id = $self->Find_Next_Transposon_ID();

  # set Corresponding_transposon in the CDS
  $self->set_Tag($class, $cds, 'Corresponding_transposon', $transposon_id);

  # make stub of Transposon object
  $self->Add_stub('Transposon', $transposon_id);

  # set Gene
  $self->set_Tag('Transposon', $transposon_id, 'Gene', $gene);

  # set Member_of to FAMILY
  $self->set_Tag('Transposon', $transposon_id, 'Member_of', $family);

  # set Remark
  $self->Add_remark('Transposon', $transposon_id, "Transposon span entered to represent a member of the $family transposon family.");

  # get start and end of CDS or Pseudogene object in Sequence
  # set start and end of Transposon in Sequence
  my ($sequence, $start, $end) = $self->get_Sequence($class, $cds);
  $self->set_Tag('Sequence', $sequence, 'Transposon', $transposon_id, $start, $end);
  
  print "In the Nameserver:\nChange class of $existing from $class to Transposon\n";

}

######################################
# NEW_GENE MODELS briggsae_GG6688|c25_g1_i1-chrII-1b

sub new_gene {
  my ($self, %line) = @_;

  # Create new gene ID

  my $class = $line{CLASS};
  my @models = @{$line{MODELS}};
  my $sequence = $self->get_Clone($class, $models[0]);
  my $cds = $self->Next_CDS_ID($sequence);

  my @new_cds;

  my $cdsl = $cds;
  my $isoform_letter = 'a';
  if (scalar @models > 1) {
    $cdsl = $cds . $isoform_letter;
  }


  foreach my $model (@models) {

    $self->Update_CDS_structure($class, $cdsl, $model);

    $self->Add_stub($class, $cdsl);

    if (scalar @models > 1) {
      $self->Add_isoform_tag($class, $cdsl)
    }
    
    $self->Add_remark($class, $cdsl, "This Gene has been created based on the RNASeq data.");


    $self->Last_reviewed($class, $cdsl);

    $self->Print_CDS_to_Isoform($cds, $cdsl);

    $isoform_letter++;
    $cdsl = $cds . $isoform_letter;

  }

  # write data for newgene.pl
  # rm new_sneak_file
  # newgene.pl -input batch_new_file -pseudoace new_sneak_file -species $species -bio CDS -who WBPerson4025 -user gw3 -password gw3 -namedb -debug gw3 (-test)
  # output is $database/NAMEDB_Files/newgene_input.ace
  # the new_sneak_file gets appended to!
  # run newgene.pl after the commands have been parsed
  # add the Gene tag by reading in the new_sneak_file output of newgene.pl
  print BNEW "$cds\n";
  my $species = $self->{wormbase}->species;
  print "New Gene seq-name $cds Species $species\n";

}
######################################
# SPLIT CLASS CDS EXISTING CBG11852 SPLITGROUP briggsae_GG6582|c0_g1_i3-chrII-1 SPLITGROUP briggsae_GG6553|c6_g1_i3-chrII-1 briggsae_GG6553|c6_g1_i1-chrII-1 

sub split_gene {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};
  my @splitgroups = @{$line{SPLITGROUP}};

  my $gene = $self->SeqName2Gene($existing);

  # Make History for each old CDS in this gene before we split it
  my $update = 0; # don't update
  my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($gene, $update);
  my $Parent_CDS = $cds;

  pop @used_letters;
  my %unused;
  if ($single) {
    $self->Make_history($class, $cds);
    $self->Delete_existing_structure($class, $cds);

  } else {
    foreach my $letter (@used_letters) {
      my $cdsl = $cds . $letter;
      $self->Make_history($class, $cdsl);
      $unused{$cdsl} = 1;
      $self->Delete_existing_structure($class, $cdsl);
    }
  }

  my $CGC_name = $self->get_CGC_name($existing, $gene);
  if (defined $CGC_name) {print "WARNING: The existing gene has a CGC_name: $CGC_name.\nThe half inheriting the CGC_name and Gene ID must be specified as the first sequence name.\n\n"; $self->force()}

  my $Parent_splitgroup_count=0; # splitgroup to inherit the original gene and CGC_name, if any (the longest resulting gene)
  my $splitgroup_count = 0;

  $splitgroup_count = 0;
  foreach my $splitgroup (@splitgroups) {

    if ($splitgroup_count == 0) { # assume the curator has done due diligence to put the half to inherit as the first one
      # use existing $cds - set above
      $cds = $Parent_CDS;
      # use existing $gene ID from the old gene that is being split
    } else {
      # get new $cds
      my $sequence = $self->get_Clone($class, $splitgroup->[0]);
      $cds = $self->Next_CDS_ID($sequence);
      # new gene ID required from batch_splitgene.pl - see below
    }

    my $isoform_letter = 'a';
    foreach my $model (@{$splitgroup}) {

      $self->structure_update_sanity_check($class, $existing, $model);

      my $cdsl = $cds; # CDS name with letter
      if (scalar @{$splitgroup} > 1) {
	$cdsl .= $isoform_letter;
	if ($splitgroup_count == 0) { # re-using isoform IDs from old gene being split
	  delete $unused{$cdsl};
	}
      }

      $self->Update_CDS_structure($class, $cdsl, $model);
      
      $self->Add_stub($class, $cdsl);

      if (scalar @{$splitgroup} > 1) {
	$self->Add_isoform_tag($class, $cdsl);
      }

      if ($splitgroup_count == 0) {
	$self->Add_Gene($class, $cdsl, $gene);
      }

      $self->Add_remark($class, $cdsl, "This gene has been split based on the RNASeq data.");
      $self->Last_reviewed($class, $cdsl);

      $self->Print_CDS_to_Isoform($cds, $cdsl);

      $isoform_letter++;
    }

    # write the data for batch_split.pl
    # batch_splitgene.pl -user gw3 -pass gw3 -ns -debug gw3 -file batch_split_file (-test)
    # run batch_splitgene.pl after the commands have been parsed
    # add the Gene tag by doing some editing of the ace output $database/NAMEDB_Files/batch_split.ace
    
    # format
    #splitgene.pl -old WBGene00238231 -new OVOC13433 -who 2062 -id WBGene00255445 -load -species ovolvulus
    #splitgene.pl -old WBGene00238817 -new OVOC13435 -who 2062 -id WBGene00255447 -load -species ovolvulus
    #splitgene.pl -old WBGene00239023 -new OVOC13436 -who 2062 -id WBGene00255448 -load -species ovolvulus
    if ($splitgroup_count != 0) {
      print BSPLIT "splitgene.pl -old $gene -new $cds -who WBPerson4025 -id WBGene00000000 -load -species $self->{wormbase}->species\n";
      my $species = $self->{wormbase}->species;
      print "In the Nameserver:\nSplit Gene Existing $gene New seq-name $cds Species $species\n";
    }

    $splitgroup_count++;
  }

  # if any CDS objects in the original gene have not been reused in
  # the split structure updates, then delete them - they already have
  # History objects
  foreach my $cdsl (keys %unused) {
    $self->Delete($class, $cdsl);
  }


}
######################################
# MERGE CLASS CDS DIE CBG00001 CBG010101 MODELS briggsae_GG3798|c9_g2_i1-chrI-1

sub merge_gene {
  my ($self, %line) = @_;

  my $update = 0; # don't update

  my $class = $line{CLASS};
  my @existing = @{$line{DIE}};
  my @models = @{$line{MODELS}};

  my %unused;
  my %single;

  my $Live;
  my @Die;
  my @genes;
  foreach my $old_seqname (@existing) {
    my $gene = $self->SeqName2Gene($old_seqname);
    push @genes, $gene;

    my $CGC_name = $self->get_CGC_name($old_seqname, $gene);
    if (defined $CGC_name) {
      #print "Found CGC_name $CGC_name for $gene in merge\n";
      if (defined $Live) {print "WARNING: Both these genes: '$gene' and '$Live' have CGC_names - take action to move the names.\n\n"; $self->force()}
      $Live = $gene;
    } else {
      push @Die, $gene
    }
  }
  if (!defined $Live) {$Live = shift @Die} # use the first gene if we don't have a CGC_name


  # write the data file for batch_merge.pl
  # batch_mergegene.pl -user gw3 -pass gw3 -ns -debug gw3 -file batch_merge_file (-test)
  # run batch_mergegene.pl after the commands have been parsed
  foreach my $Die (@Die) {
    print BMERGE " LIVE : retained geneid for thing - $Live\n";
    print BMERGE " DEAD : killed geneid for thing - $Die\n";
    print BMERGE " USER : gw3 - WBPerson4025\n";
    print BMERGE "\n";
    print "In the Nameserver:\nMerge gene to Live: $Live and gene to Die: $Die\n";
  }


  # Make History for each old CDS in each gene before we merge it
  foreach my $gene (@genes) {
    my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($gene, $update);
    $single{$gene} = $single;
    
    if ($single{$gene}) {
      $self->Make_history($class, $cds);
      $unused{$gene}{$cds} = 1;

    } else {
      foreach my $letter (@used_letters) {
	my $cdsl = $cds . $letter;
	$self->Make_history($class, $cdsl);
	$unused{$gene}{$cds . $letter} = 1;
      }
    }
  }

  # now make new structures in the Live gene - this is the gene we decide will live after the merge
  if ($single{$Live}) {
    if (scalar @models == 1) { # 1 model, 1 cds isoform
      my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($Live, $update);

      $self->Delete_existing_structure($class, $cds);    
      $self->Update_CDS_structure($class, $cds, $models[0]);
      $self->Add_Gene($class, $cds, $Live);

      $self->Add_remark($class, $cds, "The set of genes: @genes has been merged and this structure is based on the RNASeq data.");
      $self->Last_reviewed($class, $cds);
      
      $self->Print_CDS_to_Isoform($cds, $cds);

      delete $unused{$Live}{$cds};

    } else { # many models, 1 cds isoform
      my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($Live, $update);

      $self->Delete_existing_structure($class, $cds);    

      my $next_letter = 'a';
      foreach my $model (@models) {
	my $cdsl = $cds . $next_letter;

	if ($next_letter eq 'a') { # use existing cds
	  $self->Update_CDS_structure($class, $cds, $model);
	  $self->Add_Gene($class, $cds, $Live);
	  $self->Make_CDS_into_isoform($class, $cds, 'a');
	  $self->Print_CDS_to_Isoform($cds, $cds.'a');

	  delete $unused{$Live}{$cds};

	} else { # create new cds
	  $self->Update_CDS_structure($class, $cdsl, $model);
	  $self->Add_Gene($class, $cdsl, $Live);
	  $self->Add_stub($class, $cdsl);
	  $self->Add_isoform_tag($class, $cdsl);
	  $self->Print_CDS_to_Isoform($cds, $cdsl);

	  delete $unused{$Live}{$cdsl};
	}

	$self->Add_remark($class, $cdsl, "The set of genes: @genes has been merged and this structure is based on the RNASeq data.");
	$self->Last_reviewed($class, $cdsl);
	  
	$next_letter++;
      }
    }

  } else { 
    if (scalar @models == 1) { # 1 model, many cds isoforms
      my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = &Get_next_isoform_ID($Live, $update);
      my $cdsl = $cds . $used_letters[0];

      $self->Delete_existing_structure($class, $cdsl);    
      $self->Update_CDS_structure($class, $cdsl, $models[0]);
      $self->Add_Gene($class, $cdsl, $Live);

      $self->Make_isoform_into_CDS($class, $cdsl, $cds);

      $self->Add_remark($class, $cds, "The set of genes: @genes has been merged and this structure is based on the RNASeq data.");
      $self->Last_reviewed($class, $cds);
      $self->Print_CDS_to_Isoform($cds, $cds);

      delete $unused{$Live}{$cds . $used_letters[0]};

    } else { # many models, many cds isoforms
      my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($Live, $update);
      my $next_letter = 'a';
      foreach my $model (@models) {
	my $cdsl = $cds . $next_letter;
	
	# in case the cds isoform doesn't exist (more models than isoforms)
	$self->Add_stub($class, $cdsl); 
	$self->Add_isoform_tag($class, $cdsl);

	$self->Delete_existing_structure($class, $cdsl);    
	$self->Update_CDS_structure($class, $cdsl, $model);
	$self->Add_Gene($class, $cdsl, $Live);

	$self->Add_remark($class, $cdsl, "The set of genes: @genes has been merged and this structure is based on the RNASeq data.");
	$self->Last_reviewed($class, $cdsl);
	$self->Print_CDS_to_Isoform($cds, $cdsl);

	delete $unused{$Live}{$cds . $next_letter};

	$next_letter++;
      }
    }
  }

  # if any CDS objects in the original gene have not been reused in
  # the merge structure updates, then delete them - they already have
  # History objects
  foreach my $gene (keys %unused) {
    foreach my $cds (keys %{$unused{$gene}}) {
      $self->Delete($class, $cds);	
    }
  }

  #print "gene to keep alive: $genes[0]\n";
  #print "genes to die: ";
  #for (my $i=1; $i < scalar @genes; $i++) {
  #  print "$genes[$i] ";
  #}
  #print "\n";

}
######################################
# CHANGE_CLASS CLASS $class NEWCLASS PSEUDOGENE EXISTING $old_seqname
# change a structure to a different class
# The general process for changing a CDS to a Pseudogene is two steps:
#   CHANGE_CLASS of the CDS structure to Pseudogene
#   REPLACE the structure with the updated pseudogene structure

sub change_class {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $new_class = $line{NEWCLASS};
  my $existing = $line{EXISTING};


  if ($class eq $new_class) {die "ERROR CLASS and NEW_CLASS are the same: '$class'\n";}

  my $gene = $self->SeqName2Gene($existing);
  my $history_name = $self->Make_history($class, $existing);
  $self->Add_remark($class, $history_name, "This has been converted to a $new_class based on the RNASeq data.");

  $self->Add_stub($new_class, $existing);
  $self->Add_Gene($new_class, $existing, $gene);
  $self->Copy_db_remark($class, $existing, $new_class, $existing);
  $self->Copy_remarks($class, $existing, $new_class, $existing); 
  $self->Add_remark($new_class, $existing, "This has been converted to a $new_class based on the RNASeq data.");
  $self->Last_reviewed($new_class, $existing);
  $self->Copy_structure_to_new_class($class, $existing, $new_class);
  $self->Delete($class, $existing);

  print "In the Nameserver:\nChange class of $existing from $class to $new_class\n";

}

######################################
# LAST_REVIEWED CLASS $class NEWCLASS EXISTING $old_seqname
# Set the 'Last_reviewed' tag for the object and all of its sister isoforms

sub last_reviewed {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};

  # does the specified gene structure exist?
  my $gene = $self->SeqName2Gene($existing);
  if (!defined $gene) {die "ERROR Can't set the Last_reviewed tag for $existing - it is not attached to a Gene\n"}

  # get CDS ID for Gene
  my $update = 0; # don't make a new isoform
  my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = $self->Get_next_isoform_ID($gene, $update);

  if ($single) {
    $self->Add_remark($class, $cds, "This $class has been inspected and looks satisfactory.");
    $self->Last_reviewed($class, $cds);
    print "$cds has been reviewed\n";

  } else {
    foreach my $letter (@used_letters) {
      my $cdsl = $cds . $letter;
      $self->Add_remark($class, $cdsl, "This $class has been inspected and looks satisfactory.");
      $self->Last_reviewed($class, $cdsl);
      print "$cdsl has been reviewed\n";
    }
  }

}

######################################
# JUST_MAKE_HISTORY CLASS $class NEWCLASS EXISTING $old_seqname
# Simply make a History object for the specified object

sub just_make_history {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};

  my $gene = $self->SeqName2Gene($existing);
  if (!defined $gene) {die "ERROR Can't make a history object from $existing - it is not attached to a Gene\n"}
  $self->Make_history($class, $existing);

}

######################################
# CHECK_GENE_NAME CLASS $class NEWCLASS EXISTING $old_seqname
# Report any Gene_name belonging to the Gene of the specified object

sub check_gene_name {
  my ($self, %line) = @_;

  my $class = $line{CLASS};
  my $existing = $line{EXISTING};

  my $gene = $self->SeqName2Gene($existing);
  if (!defined $gene) {die "ERROR Can't check Gene_name for $existing - it is not attached to a Gene\n"}
  my $CGC_name = $self->get_CGC_name($existing, $gene);
  if (defined $CGC_name && $CGC_name ne '') {
    print "The $class object $existing has a Gene_name: '$CGC_name'\n";
  } else {
    print "The $class object $existing has not got a Gene_name\n";
  }
}

######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
# add a Remark tag
sub Add_remark {
  my ($self, $class, $cds, $text) = @_;

  my $date = $self->date();
  my $fh = $self->{out};
  print $fh "\n// Add_remark\n";
  print $fh "$class : $cds\n";
  print $fh "Remark \"[$date $ENV{USER}] $text\" Curator_confirmed $self->{person}\n";
  print $fh "Remark \"[$date $ENV{USER}] $text\" From_analysis RNASeq\n";
  print $fh "Remark \"[$date $ENV{USER}] $text\" Paper_evidence WBPaper00049972\n";
  print $fh "\n";
}

######################################
sub date {
  my ($self) = @_;
  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
}
######################################
# add the Last_reviewed tag
sub Last_reviewed {
  my ($self, $class, $cds) = @_;

  my $fh = $self->{out};
  print $fh "\n// Last_reviewed\n";
  print $fh "$class : $cds\n";
  print $fh "Last_reviewed now $self->{person}\n";
  print $fh "\n";

}

######################################
sub Change_to_Transposon_CDS {
  my ($self, $class, $cds) = @_;
  
  my $fh = $self->{out};
  print $fh "\n// Change_to_Transposon_CDS\n";
  print $fh "$class : $cds\n";
  
  if ($class eq 'CDS') {
    print $fh "Method Transposon_CDS\n";

  } elsif ($class eq 'Pseudogene') {
    print $fh "Method Transposon_Pseudogene\n";
    print $fh "Type Transposon_pseudogene\n";
  } else {
    die "ERROR $cds is not a CDS or a PSeudogene. I don't know how to set this up as a Transposon!\n";
  }    
  print $fh "\n";
  }
######################################
# make a history object
sub Make_history {
  my ($self, $class, $existing) = @_;

  my $clone_tag; # where it lives under the Sequence object tags
  my $new_method; # the new method to give it
  my $transcript_type;
  my $transcript_type_value1;
  my $transcript_type_value2;
  my $pseudogene_type;

  my $date = $self->date();
  my $version = $self->get_history_version($self->{wormbase}->database('current'));

  my $wormpep_prefix  = $self->{wormbase}->wormpep_prefix;
  $wormpep_prefix= lc($wormpep_prefix);
  # japonica uses 'jp' not 'ja' to name the history objects - blame Phil for this!
  if ($wormpep_prefix eq 'ja') {$wormpep_prefix = 'jp'}
  if ($wormpep_prefix eq 'cn') {$wormpep_prefix = 'np'}
  my $history = "$existing:${wormpep_prefix}$version";

  #print "enter CDS to make history object \n";
  return unless $existing;

  my $obj = $self->{ace}->fetch($class => "$existing");
  if ($class eq 'CDS') {
    $new_method = 'history';
    if ($obj->Method->name ne "curated" ) {
      die "ERROR $existing is not curated| I only do curated CDS histories!\n";
    }
  } elsif ($class eq 'Pseudogene') {
    $new_method = 'history_pseudogene';
    $pseudogene_type = $obj->Type->name;

  } elsif ($class eq 'Transcript') {
    $new_method = 'history_transcript';
    if ($obj->Method->name eq "Coding_transcript" ) {
      die "ERROR I don't do Coding_transcript Transcript histories!\n";
    }
    ($transcript_type, $transcript_type_value1, $transcript_type_value2) = $obj->Transcript->row(0);
  }
  $clone_tag = $self->CloneTag($class);
  
  if (!defined $obj) {
    die "ERROR Invalid $class. $existing is not a valid $class name";
  }

  if ($self->{ace}->fetch('CDS' => $history) || $self->{ace}->fetch('Pseudogene' => $history) || $self->{ace}->fetch('Transcript' => $history)) {
    warn "The history object '$history' already exists - not making this again\n";
    return $history;
  }
    
  my $species = $obj->Species->name;
  my $gene = $obj->Gene->name;
  my $lab = $obj->From_laboratory->name;
  my $seq = $obj->Sequence->name;
  my $brief = $obj->Brief_identification;
  my $brief_identification = $brief->name if ($brief);

  # parent clone coords
  my $clone = $obj->Sequence;
  my @clone_entry;
  my $clone_tag = $self->CloneTag($class);
  @clone_entry = $clone->${clone_tag};

  my $start;
  my $end;
  foreach my $CDS ( @clone_entry ) {
    next unless ($CDS->name eq "$existing");
    $start = $CDS->right->name;
    $end = $CDS->right->right->name;
    last;
  }

  my $fh = $self->{out};
  print $fh "\n// Make_history\n";
  print $fh "Sequence : $seq\n";
  print $fh "$clone_tag \"$history\" $start $end\n";

  print $fh "\n$class : $history\n";

  foreach ($obj->Source_exons) {
    my ($start, $end) = $_->row(0);
    print $fh "Source_exons ",$start->name," ",$end->name,"\n";
  }

  print $fh "Sequence $seq\n";
  print $fh "From_laboratory $lab\n";
  print $fh "Gene_history $gene\n" if $gene;
  print $fh "Species \"$species\"\n" if $species;
  print $fh "Evidence Curator_confirmed $self->{person}\n";
  print $fh "Remark \"[$date $ENV{USER}] Automatically generated history object.\"\n";
  print $fh "Method $new_method\n";
  print $fh "Brief_identification \"$brief_identification\"\n" if ($brief_identification);

  print $fh "CDS\n" if ($class eq 'CDS');
  print $fh "CDS_predicted_by ${\$obj->CDS_predicted_by}\n" if ($class eq 'CDS' && $obj->CDS_predicted_by);
 
  print $fh "$pseudogene_type" if ($pseudogene_type);

  print $fh "$transcript_type" if ($transcript_type);
  print $fh " \"",$transcript_type_value1->name,"\"" if ($transcript_type_value1);
  print $fh " \"",$transcript_type_value2->name,"\"" if ($transcript_type_value2);
  print $fh "\n" if ($transcript_type);

  print $fh "\n";

  $self->Copy_remarks($class, $existing, $class, $history);

  return "$history";
}

######################################
# Return which tag the class lives under in the Sequence object

sub CloneTag {
  my ($self, $class) = @_;
  my %ct = ('CDS'        => 'CDS_child',
	    'Pseudogene' => 'Pseudogene',
	    'Transcript' => 'Transcript',
	   );
  return $ct{$class};
}

######################################
# if the Build has started, get the Build version
# else use the version in the reference database (usually CurrentDB)
sub get_history_version {
  my ($self, $rdatabase) = @_;

  my $WS_version = `grep "NAME WS" $rdatabase/wspec/database.wrm`;
  chomp($WS_version);
  $WS_version =~ s/.*WS//;
  die "ERROR no version\n" unless $WS_version =~ /\d+/;
  
  # check to see if the Build has started - if so then use the Build's version
  my $build_file = $self->{wormbase}->autoace . "/wspec/database.wrm";
  if (-e $build_file) {
    my $Build_WS_version = `grep "NAME WS" $build_file`;
    chomp($Build_WS_version);
    $Build_WS_version =~ s/.*WS//;
    if (defined $Build_WS_version ) {$WS_version = $Build_WS_version }
  }

  return($WS_version); 
}
######################################
# remove the old Source exons and the Sequence location
sub Delete_existing_structure {
  my ($self, $class, $cds) = @_;

  my ($sequence, $start, $end) = $self->get_Sequence($class, $cds);

  my $fh = $self->{out};
  print $fh "\n// Delete_existing_structure\n";
  print $fh "$class : $cds\n";
  print $fh "-D Source_exons\n";
  print $fh "-D Sequence\n";
  print $fh "-D Start_not_found\n" if ($class eq 'CDS');
  print $fh "-D End_not_found\n" if ($class eq 'CDS');
  print $fh "\n";

  my $clone_tag = $self->CloneTag($class);

  print $fh "\n";
  print $fh "Sequence : $sequence\n";
  print $fh "-D $clone_tag $cds\n";
  print $fh "\n";

}

######################################
# create the exon structure of the CDS and its location
# Source_exons start end
# Sequence ?Sequence
# Sequence CDS_child ?CDS start end
# Start_not_found
# End_not_found
  
sub Update_CDS_structure {
  my ($self, $class, $cds, $source_model) = @_;

  my ($sequence, $start, $end) = $self->get_Sequence($class, $source_model);
  my ($starts, $ends) = $self->get_Source_exons($class, $source_model); 

  my ($start_not_found, $end_not_found, $premature_stops) = $self->start_and_stop_sanity_check($class, $source_model, $sequence, $start, $end, $starts, $ends);
  if ($premature_stops && $class eq 'CDS') {print "WARNING: $source_model has premature STOP codons\n\n"; $self->force('elegans');}
  if ($start_not_found && $class eq 'CDS') {print "WARNING: $source_model has START_NOT_FOUND\n\n";  $self->force('elegans');}
  if ($end_not_found && $class eq 'CDS') {print "WARNING: $source_model has END_NOT_FOUND\n\n";  $self->force('elegans');}

  my $clone_tag = $self->CloneTag($class);

  # +++ add check and tag for small artificial introns Properties Artificial_intron Low_quality_sequence for ENA dumping

  my $fh = $self->{out};
  print $fh "\n// Update_CDS_structure\n";
  print $fh "Sequence : $sequence\n";
  print $fh "$clone_tag $cds $start $end\n";
  print $fh "\n";
  print $fh "$class : $cds\n";
  print $fh "Sequence $sequence\n";
  for (my $i = 0; $i < scalar @{$starts}; $i++) {
    print $fh "Source_exons $starts->[$i] $ends->[$i]\n";
  }
  if ($start_not_found && $class eq 'CDS') {
    print $fh "START_NOT_FOUND 1\n";
  }
  if ($end_not_found && $class eq 'CDS') {
    print $fh "END_NOT_FOUND\n";
  }
  print $fh "\n";

}
######################################
# copy the exon structure and location from one class to an object of a new_class but with the same sequence_name
  
sub Copy_structure_to_new_class {
  my ($self, $class, $existing, $new_class) = @_;

  my ($sequence, $start, $end) = $self->get_Sequence($class, $existing);
  my ($starts, $ends) = $self->get_Source_exons($class, $existing); 

  my ($start_not_found, $end_not_found, $premature_stops) = $self->start_and_stop_sanity_check($class, $existing, $sequence, $start, $end, $starts, $ends);

  my $clone_tag = $self->CloneTag($new_class);

  my $fh = $self->{out};
  print $fh "\n// Copy_structure_to_new_class\n";
  print $fh "Sequence : $sequence\n";
  print $fh "$clone_tag $existing $start $end\n";
  print $fh "\n";
  print $fh "$new_class : $existing\n";
  print $fh "Sequence $sequence\n";
  for (my $i = 0; $i < scalar @{$starts}; $i++) {
    print $fh "Source_exons $starts->[$i] $ends->[$i]\n";
  }
  if ($start_not_found && $new_class eq 'CDS') {
    print $fh "START_NOT_FOUND 1\n";
  }
  if ($end_not_found && $new_class eq 'CDS') {
    print $fh "END_NOT_FOUND\n";
  }
  print $fh "\n";

}
######################################
# copy DB_Remark text
#
sub Copy_db_remark {
  my ($self, $class, $existing, $new_class, $new_seqname) = @_;

  my $obj = $self->{ace}->fetch($class => "$existing");
  my $db_remark = $obj->DB_remark;

  my $fh = $self->{out};
  if (defined $db_remark && $db_remark ne '') {
    print $fh "\n// Copy_db_remark\n";
    print $fh "$new_class : $new_seqname\n";
    print $fh "DB_remark \"$db_remark\"\n";
    print $fh "\n";
  }

}

######################################
# copy Remark text
#
sub Copy_remarks {
  my ($self, $class, $existing, $new_class, $new_seqname) = @_;

  my $obj = $self->{ace}->fetch($class => "$existing");

  my $fh = $self->{out};
  print $fh "\n// Copy_remarks\n";
  print $fh "$new_class : $new_seqname\n";

  foreach ($obj->Remark) {
    my ($remark, $evidence, $evidence_value1, $evidence_value2) = $_->row(0);
    print $fh "Remark \"", $remark->name ,"\"";
    print $fh " \"", $evidence->name,"\"" if ($evidence);
    print $fh " \"", $evidence_value1->name,"\"" if ($evidence_value1);
    print $fh " \"", $evidence_value2->name,"\"" if ($evidence_value2);
    print $fh "\n";
  }
  print $fh "\n";

}

######################################
#  my ($sequence, $start, $end) = &get_Sequence($cds);

sub get_Sequence {
  my ($self, $class, $cds) = @_;

  my $cds_obj = $self->{ace}->fetch($class => "$cds");
  if (! defined  $cds_obj) {die "ERROR Can't find $class $cds\n"}
  my $cds_sequence_obj = $cds_obj->Sequence;
  if (!defined $cds_sequence_obj) {die "ERROR Can't find a Sequence for $class $cds - Does it exist? Did you save your session?\n"}
  my $cds_sequence_name = $cds_sequence_obj->name;

  my $clone_tag = $self->CloneTag($class);

  my @clone_locations = $cds_sequence_obj->${clone_tag};

  my $cds_start;
  my $cds_end;

  foreach my $cds_location ( @clone_locations ) {
    next unless ($cds_location->name eq $cds);
    $cds_start = $cds_location->right->name;
    $cds_end = $cds_location->right->right->name;
    last;
  }
  return ($cds_sequence_name, $cds_start, $cds_end);
}
######################################
# get the clone that the 5' end of the structure is on

# this is used when creating a new ID for a locus becuase in elegans
# the Sequence can span two clones and so is promoted to the top level
# chromosome sequence, but we don't name loci after chromosomes in
# elegans, so we want just the Clone name.

sub get_Clone {
  my ($self, $class, $cds) = @_;

  my ($sequence, $start, $end) = $self->get_Sequence($class, $cds);
  if ($sequence =~ /^CHROMOSOME/) {
    my $coords = Coords_converter->invoke($self->{wormbase}->autoace, 0, $self->{wormbase});
    $sequence = $coords->GetCloneFromCoords($sequence, $start, $start+1);
    if ($sequence =~ /^CHROMOSOME/) {
      $sequence = $coords->GetCloneFromCoords($sequence, $end-1, $end);
      if ($sequence =~ /^CHROMOSOME/) {
	die "ERROR The ends of this structure both appear to be on the top-level CHROMOSOME sequence\n";
      }
    }
  }

  return $sequence;
}

######################################
#  my ($starts, $ends) = &get_Source_exons($class, $source_model); 

sub get_Source_exons {
  my ($self, $class, $cds) = @_;

  my $cds_obj = $self->{ace}->fetch($class => "$cds");

  my @starts;
  my @ends;

  foreach ($cds_obj->Source_exons) {
    my ($start,$end) = $_->row(0);
    push @starts, $start->name;
    push @ends, $end->name;
  }
  return(\@starts, \@ends);
}

######################################
# get sequence of new model and do sanity check on start and end
sub start_and_stop_sanity_check {
  my ($self, $class, $source_model, $sequence, $start, $end, $starts, $ends) = @_;
  my $start_not_found=0;
  my $end_not_found=0;
  my $premature_stops=0;

  my $cds_seq = $self->get_DNA($class, $source_model);

  if ($cds_seq !~ /^atg/) {
    $start_not_found = 1;
  }

  my @stops;
  foreach my $stop ('taa', 'tga', 'tag') {
    my $st_offset = 0;
    my $st_result = index($cds_seq, $stop, $st_offset);
    while ($st_result != -1) {
      my $frame = $st_result % 3;
      if ($frame == 0) {
	push @stops, $st_result; # 0 is the first base
      }
      $st_offset = $st_result + 1;
      $st_result = index($cds_seq, $stop, $st_offset);
    }
  }
  
  # check the last codon is at the end of the Sequence span of the structure
  my $last_stop = $stops[$#stops];
  if (!defined $last_stop || (length $cds_seq) - 3 != $last_stop) {
    $end_not_found = 1;
    if (scalar @stops > 0) {$premature_stops = 1;}
  } else {
    if (scalar @stops > 1) {$premature_stops = 1;}
  }

  return ($start_not_found, $end_not_found, $premature_stops);
}
######################################
# get compare properties of an existing and a new model and do sanity checks
sub structure_update_sanity_check {
  my ($self, $class, $existing, $model) = @_;

  # check they are the same sense indicating that the right model is being used
  my ($existing_sequence, $cds_start, $cds_end) = $self->get_Sequence($class, $existing);
  my ($model_sequence, $start, $end) = $self->get_Sequence($class, $model);
  my $cds_sense = '+';
  if ($cds_start > $cds_end) {$cds_sense = '-'}
  my $sense = '+';
  if ($start > $end) {$sense = '-'}
  if ($cds_sense ne $sense) {
    print "WARNING: The sense of the existing structure and the new model are different: '$existing' ($cds_sense) and '$model' ($sense)\n\n";
    $self->force();
  }
  # check to warn if the Sequence tag is changed indicating that the right model is being used
  if ($existing_sequence ne $model_sequence) {
    print "WARNING: The locations of the existing structure and the new model are on different sequences: '$existing' ($existing_sequence) and '$model' ($model_sequence)\n\n";
    $self->force();
  }
  # check they overlap indicating that the right model is being used
  my $overlap = 0;
  if ($sense eq '+') {
    if ($cds_start <= $end && $cds_end >= $start) {$overlap = 1}
  } else {
    if ($cds_end <= $start && $cds_start >= $end) {$overlap = 1}
  }
  if (!$overlap) {print "WARNING: The existing structure '$existing' ($cds_start - $cds_end) and new model '$model' ($start - $end) do not overlap\n\n"; $self->force(); }

}
######################################
# compare the structures of two objects
# checks the number and length of the exons
# allows the Sequence to be different in case one object is on a clone and the other is on a Chromosome
# returns true if they look like thay are the same things
sub structure_comparison_sanity_check {
  my ($self, $class, $existing, $model) = @_;

  # are they the same named object?
  if ($existing eq $model) {
    return 1;
  }

  # check they are the same sense
  my ($existing_sequence, $cds_start, $cds_end) = $self->get_Sequence($class, $existing);
  my ($model_sequence, $start, $end) = $self->get_Sequence($class, $model);
  my $cds_sense = '+';
  if ($cds_start > $cds_end) {$cds_sense = '-'}
  my $sense = '+';
  if ($start > $end) {$sense = '-'}
  if ($cds_sense ne $sense) {
    return 0;
  }
  
  # check the span of the objects
  my $cds_span = $cds_start - $cds_end;
  my $model_span = $start - $end;
  if ($cds_span != $model_span) {
    return 0;
  }

  # check sequence
  my $cds_seq = $self->get_DNA($class, $existing);
  my $model_seq = $self->get_DNA($class, $model);
  if ($cds_seq eq $model_seq) {
    return 1;
  } else {
    return 0;
  }

}

######################################

sub get_DNA {
  my ($self, $class, $cds) = @_;

  my $cds_obj = $self->{ace}->fetch($class => $cds);
  my $cds_seq = $cds_obj->asDNA();

  $cds_seq =~ s/^>\S+//;      # remove the title line
  $cds_seq =~ s/\n//g;          # remove newlines
  $cds_seq = lc $cds_seq;

  return $cds_seq;
}
######################################
# my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters) = Get_next_isoform_ID($gene, $update);
# if $update is false then the existing details are returned ($single is true if the existing CDS is not an Isoform)
# if $update is true then details are updated for a new isoform ($single is not changed - it is only true if the old CDS was not an Isoform)
sub Get_next_isoform_ID {
  my ($self, $gene, $update) = @_;

  my ($cds, $isoform_letter, $number_of_existing_isoforms, $single, @used_letters);

  if ( exists $self->{Isoform_details}{$gene}) {
    $cds = $self->{Isoform_details}{$gene}{cds};
    if ($update) {$self->{Isoform_details}{$gene}{isoform_letter}++}
    $isoform_letter = $self->{Isoform_details}{$gene}{isoform_letter};
    if ($update) {$self->{Isoform_details}{$gene}{number_of_existing_isoforms}++}
    $number_of_existing_isoforms = $self->{Isoform_details}{$gene}{number_of_existing_isoforms};
    if ($update) {$self->{Isoform_details}{$gene}{single} = 0}
    $single = $self->{Isoform_details}{$gene}{single};
    if ($update) {push @{$self->{Isoform_details}{$gene}{used_letters}}, $self->{Isoform_details}{$gene}{isoform_letter}}

  } else {

    my $gene_obj = $self->{ace}->fetch('GENE' => "$gene");
    my @structures;
    my @corr_cds = $gene_obj->Corresponding_CDS;
    foreach my $struc (@corr_cds) {
      my $method = $struc->Method;
      if ($method eq 'curated') {push @structures, $struc}
    }

    # include Non_coding_transcripts
    my @corr_trans = $gene_obj->Corresponding_transcript;
    foreach my $tran (@corr_trans) {
      my $method = $tran->Method;
      if ($method ne 'Coding_transcript') {push @structures, $tran}
    }

    # this might be a Pseudogene
    my @corr_ps = $gene_obj->Corresponding_pseudogene;
    foreach my $pseud (@corr_ps) {
      my $method = $pseud->Method;
      if ($method eq 'Pseudogene') {push @structures, $pseud}
    }


    $isoform_letter = '';
    $number_of_existing_isoforms = 0;
    foreach my $corr_cds (@structures) {
      my ($this_cds, $letter) = $self->GeneSeqname($corr_cds);
      if (defined $this_cds) { # something like 'briggsae_non_coding_transcript_isoformer_WBGene00025804_6602' is Method=curated, but doesn't match the CDS regex
	$cds = $this_cds;
	$number_of_existing_isoforms++;
	if ($letter gt $isoform_letter) {$isoform_letter = $letter}
	push @{$self->{Isoform_details}{$gene}{used_letters}}, $letter;
      }
    }

    if ($isoform_letter eq '') {
      $self->{Isoform_details}{$gene}{single} = 1; # this CDS is currently on its own and so is not an isoform and currently has no letter at the end of its ID
    } else {
      $self->{Isoform_details}{$gene}{single} = 0;
    }
    if ($update) {
      if ($isoform_letter eq '') {$isoform_letter = 'a'}
      $isoform_letter++;
      $number_of_existing_isoforms++;
      push @{$self->{Isoform_details}{$gene}{used_letters}}, $isoform_letter;
    }
  
    # and store for next time
    $self->{Isoform_details}{$gene}{cds} = $cds;
    $self->{Isoform_details}{$gene}{isoform_letter} = $isoform_letter;
    $self->{Isoform_details}{$gene}{number_of_existing_isoforms} = $number_of_existing_isoforms;
  }

  return ($cds, $isoform_letter, $number_of_existing_isoforms, $self->{Isoform_details}{$gene}{single}, @{$self->{Isoform_details}{$gene}{used_letters}});
}
######################################
# strip off isoform letter from a sequence name to give a gene sequence name and an isoform letter
#  my ($gene_seqname, letter) = $self->GeneSeqname($seqname);

sub GeneSeqname {
  my ($self, $seqname) = @_;

  my $seqname_regex = $self->{wormbase}->seq_name_regex();
  my $isoform_letter = '';
  my $gene_seqname;

  ($gene_seqname, $isoform_letter) = ($seqname =~ /($seqname_regex)(.*)$/);

  
  return ($gene_seqname, $isoform_letter);
}

######################################
# &Add_stub($class, $cds);
# species, lab, Method, etc details for a typical object

sub Add_stub {
  my ($self, $class, $seqname) = @_;
  
  my $long_name = $self->{wormbase}->long_name();
  my $fh = $self->{out};
  
  print $fh "\n// Add_stub\n";
  print $fh "$class : $seqname\n";
  print $fh "From_laboratory HX\n";
  print $fh "Species \"$long_name\"\n";
  
  if ($class eq 'CDS') {
    print $fh "CDS\n";
    print $fh "Method curated\n";
    print $fh "\n";
  } elsif ($class eq 'Pseudogene') {
    print $fh "Type Coding_pseudogene\n";
    print $fh "Method Pseudogene\n";
    print $fh "Brief_identification \"Pseudogene.\"\n";
    print $fh "\n";
  } elsif ($class eq 'Transcript') {
    print $fh "Transcript ncRNA\n"; # if you need something different, then you will have to edit the object
    print $fh "Method ncRNA\n";
    print $fh "\n";
  } elsif ($class eq 'Transposon') {
    print $fh "Method Transposon\n";
    print $fh "\n";
  } else {
    die "ERROR Unknown class: $class\n";
  }

}

######################################
# &Add_Gene($class, $cds, $gene);

sub Add_Gene {
  my ($self, $class, $cds, $gene) = @_;
  
  my $fh = $self->{out};
  print $fh "\n// Add_Gene\n";
  print $fh "$class : $cds\n";
  print $fh "Gene $gene\n";
  print $fh "\n";
  
}

######################################
# &Add_isoform_tag($class, $cdsl);

sub Add_isoform_tag {
  my ($self, $class, $cds) = @_;

  my $fh = $self->{out};
  print $fh "\n// Add_isoform_tag\n";
  print $fh "$class : $cds\n";
  print $fh "Isoform Curator_confirmed $self->{person}\n";
  print $fh "\n";

}

######################################
# &Make_CDS_into_isoform($class, $cds, 'a');

sub Make_CDS_into_isoform {
  my ($self, $class, $cds, $isoform_letter) = @_;
  
  $self->Add_isoform_tag($class, $cds);

  my $cdsl = $cds . $isoform_letter;
  
  # rename with isoform_letter
  my $fh = $self->{out};
  print $fh "\n// Make_CDS_into_isoform\n";
  print $fh "-R $class $cds $cdsl\n";
  print $fh "\n";

}

######################################
# set arbitraty tag to arbitray value(s)

sub set_Tag {
  my ($self, $class, $cds, $tag, @value) = @_;

  my $fh = $self->{out};
  print $fh "\n// set_Tag\n";
  print $fh "$class : $cds\n";
  print $fh "$tag @value\n";
  print $fh "\n";

  
}
######################################
# &Make_isoform_into_CDS($class, $cdsl, $cds);

sub Make_isoform_into_CDS {
  my ($self, $class, $isoform, $cds) = @_;
  
  # rename without isoform_letter
  my $fh = $self->{out};
  print $fh "\n// Make_isoform_into_CDS\n";
  print $fh "-R $class $isoform $cds\n";
  print $fh "\n";
  print $fh "$class : $cds \n";
  print $fh "-D Isoform\n";
  print $fh "\n";

}

######################################
# &Delete($class, $cds);

sub Delete {
  my ($self, $class, $cds) = @_;
  
  my $fh = $self->{out};
  print $fh "\n// Delete\n";
  print $fh "-D $class $cds\n";
  print $fh "\n";

}
######################################
# $transposon_id = $self->Find_Next_Transposon_ID();
# Warning! This finds the next available Transposon ID in the curation database - there may be higher IDs used in other databases! Manual checking is required!
sub Find_Next_Transposon_ID {
  my ($self) = @_;

  my $maxnumber = "WBTransposon00000000";

  my @trans = $self->{ace}->fetch(-query => "find Transposon");
  foreach my $thing (@trans) {
    if ($thing->name gt $maxnumber) {$maxnumber = $thing->name}
  }
  $maxnumber++;
  return $maxnumber;
}
######################################
# test for a CGC name
# look in the Gene class object first, then try the Nameserver
sub get_CGC_name {

  my ($self, $cds, $gene) = @_;

  if (!defined $gene) {die "ERROR In get_CGC_name: Gene ID not specified\n";}
  my $gene_obj = $self->{ace}->fetch('GENE' => "$gene");
  if (!defined $gene_obj) {die "ERROR In get_CGC_name: Can't find a Gene object called '$gene'\n";}
  my $CGC_name = $gene_obj->CGC_name;

  my $species = $self->{wormbase}->species;
  if (!defined $CGC_name && $species ne 'elegans') {
    $CGC_name = $self->cgcname($cds);
  }

  return $CGC_name;
}


######################################
# get the CGC-name, if any, from the Nameserver
# $cgcname = $self->cgcname($seqname);
  
sub cgcname {
  my ($self, $seqname) = @_;
  my $cgcname = '';
  
  my ($gene_seqname, $letter) = $self->GeneSeqname($seqname);

  # MYSQL CONFIG VARIABLES
  my ($port,$host,$database,);
  $host = 'web-wwwdb-core-02';
  $database = 'nameserver_live';
  $port = 3449;
  my $dsn = "dbi:mysql:$database:$host:$port" or die "ERROR Unable to connect: $DBI::errstr\n";

  #                           dsn, user,          pw
  my $connect = DBI->connect($dsn, $self->{user}, $self->{user});

  my (%storedIDs,$myquery, $execute,);

  # try to pull out CGC names (i.e. name_type_id = 1) linked with the specified sequence name (name_type_id = 3)
  my $myquery = "select object_name from secondary_identifier where name_type_id = 1 AND object_id in (select object_id from secondary_identifier where name_type_id = 3 AND object_name = '${gene_seqname}')";
  my $query_handle = $connect->prepare($myquery);
  $query_handle->execute(); 
  $query_handle->bind_columns(\$cgcname);
  while($query_handle->fetch()) {
    return $cgcname;
  }

}

######################################
# &Print_CDS_to_Isoform($cds, $cdsl);
# used for linking genes made by the nameserver to isoforms and for general stats

sub Print_CDS_to_Isoform {
  my ($self, $cds, $cdsl) = @_;

  print CDS2ISO "$cds\t$cdsl\n";

}


######################################

1;

######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################



# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - Curate.pm

=head1 USAGE

=over 4

=item  $Curate = Curate->new($db, $seqextract, $wb, $log);

=back

Module to do curation


=back

=over 4

=item -h, Help

=back

=over 4
 
No help.

=back


=head1 REQUIREMENTS

=over 4

=item An acedb database.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut

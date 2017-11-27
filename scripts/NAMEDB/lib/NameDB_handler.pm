package NameDB_handler;
use Carp;
use NameDB;
use warnings;
use strict;

#author ar2

=head1 

NameDB_handler

This is an API layer above the NameDB.pm module specifically for handling GeneID maninpulation

=item Synopsis

my $DOMAIN = 'Gene';
my $db = NameDB_handler->new($DB,$USER,$PASS);
	
$db || return(undef);
$db->setDomain($DOMAIN);

my $gene_id = $db->merge_genes($id_to_keep, $id_to_kill);
my $split   = $db->split_genes($cds_name, $name_type,$gene_id);

=cut

#use base NameDB;
our @ISA = qw( NameDB );

sub new
  {
    my $class = shift;
    my $dsn = shift;
    my ($name,$password, $path) = @_;
    my $db = NameDB->connect($dsn,$name,$password);

    bless ($db, $class);

    $path||="$ENV{CVS_DIR}/NAMEDB";

    #read in clone list to validate CDS names with
    my $clone_file = "$path/clonelist";
	
    #untaint file
    unless( $clone_file =~ m/^(.+)$/ ) {
      $db->dienice("tainted file\n");
    }
    $clone_file = $1;
	
    #read clones in
    open(CLONE,"<$clone_file") or $db->dienice("cant open $clone_file\n");
    my %clones;
    while (<CLONE>) {
      chomp;
      $db->{clones}->{uc($_)} = 1;
    }
    close CLONE;
    return $db;
  }

sub web				# set flag to specify web or script for error (dienice)
  {
    my $db = shift;
    my $set = shift;
    $db->{'web'} = $set if( $set );
    return $db->{'web'};
  }

sub printAllNames
  {
    my $db = shift;
    my $obj = shift;
    my $names = $db->idAllNames($obj);
    print "<br>Current gene names for <b>$obj</b> <br>";
    foreach (keys %{$names} ) {
      my $name_str;
      if (ref $names->{"$_"} eq 'ARRAY') {
	$name_str =  join(" ",@{$names->{"$_"}})
      } else {
	$name_str = $names->{$_};
      }
      print "$_ : $name_str<br>";
    }
  }


sub add_name {
  my ($db, $gene_id, $name, $type, $species) = @_;
  unless ($db and $gene_id and $name and $type) {
    $db->dienice("bad parameters");
    return undef;
  }
  
  #confirm GeneID and name are valid
  $db->validate_name($name, $type, $species)    or return undef;
  $db->validate_id($gene_id)          or return undef;
  $db->check_pre_exists($name, $type) or return undef;
  if ( $db->addName($gene_id,$type => $name) ) {
    $db->_update_public_name($gene_id, $species);
    return $name;  
  } else {
    $db->dienice("$name not added to $gene_id");
    return undef;
  }
}


# adding this for use by gene name curator to add lists of non-standard gene names.
# NOT for use with web interface.
sub force_name {
  my ($db, $gene_id, $name, $type, $species) = @_;
  unless ($db and $gene_id and $name and $type) {
    $db->dienice("bad parameters");
    return undef;
  }
  
  #confirm GeneID is valid
  $db->validate_id($gene_id)          or return undef;
  $db->check_pre_exists($name, $type) or return undef;
  if ( $db->addName($gene_id,$type => $name) ) {
    $db->_update_public_name($gene_id, $species);
    return $name;  
  } else {
    $db->dienice("$name not added to $gene_id");
    return undef;
  }
}




sub validate_name {
  my $db = shift;
  my $name = shift;
  my $type = shift;
  my $species = shift;
  #is this a valid name tpye?
  my @types = $db->getNameTypes;
  if ( grep {$_ eq $type} @types) {    		
    #check name structure matches format eg CDS = clone.no
    my $name_checks = {
		       'elegans' => { 
				     "CGC" => '^[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
				     "Sequence" => '(^([A-Z0-9_cel]+)\.\d+$)|(^([A-Z0-9_cel]+)\.t\d+$)',
				     "Public_name" => '^[a-z]{3,4}-[1-9]\d*(\.\d+)?$|^([A-Z0-9_]+)\.\d+$',
				    },
		       'briggsae' => {
				      "CGC" => '(^Cbr-[a-z21]{3,4}-[1-9]\d*(\.\d+)?)|(^Cbr-[a-z21]{3,4}\([a-z]+\d+\)?$)',
				      "Sequence" => '^CBG\d{5}$',
				      "Public_name" => '(^Cbr-[a-z21]{3,4}-[1-9]\d*(\.\d+)?)|(^Cbr-[a-z21]{3,4}\([a-z]+\d+\)?$)|^CBG\d{5}$',
				     },

		       'remanei' => {
				     "CGC" => '^Cre-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
				     "Sequence" => '^CRE\d{5}$',
				     "Public_name" => '^Cre-[a-z]{3,4}-[1-9]\d*(\.\d+)?$|^CRE\d{5}$',
				    },
		       'brenneri' => {
				      "CGC" => '^Cbn-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
				      "Sequence" => '^CBN\d{5}$',
				      "Public_name" => '^Cbn-[a-z]{3,4}-[1-9]\d*(\.\d+)?$|^CBN\d{5}$',
				     },
                        'pristionchus' => {
				     "CGC" => '^Ppa-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
				     "Sequence" => '^PPA\d{5}$',
				     "Public_name" => '^Ppa-[a-z]{3,4}-[1-9]\d*(\.\d+)?$|^PPA\d{5}$',
				    },
                        'japonica' => {
				     "CGC" => '^Cjp-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
				     "Sequence" => '^CJA\d{5}$',
				     "Public_name" => '^Cjp-[a-z]{3,4}-[1-9]\d*(\.\d+)?$|^CJA\d{5}$',
				    },
		        'brugia'   =>{
			             'CGC' => '^Bma-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
			             'Sequence' => '^Bm\d+$',
			             'Public_name' => '^Bma-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$|^Bm\d+$'	     
			            },
                        'ovolvulus' => {
                                     'CGC' => '^Ovo-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
                                     'Sequence' => 'OVOC\d+$',
                                     'Public_name' => '^Ovo-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$|^OVOC\d+$'
                                     },
                        'sratti' => {
                                     'CGC' => '^Sra-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
                                     'Sequence' => 'SRAE_[\dXM]\d+$',
                                     'Public_name' => '^Sra-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$|^SRAE_[\dXM]\d+$'
                                     },
		      };

    unless ( $name =~ /$name_checks->{$species}->{$type}/ ) {
      $db->dienice("$name is incorrect format for $species $type match ".$name_checks->{$species}->{$type});
      return undef;
    }
#    if (($species eq 'elegans') and ($type eq 'Sequence') and !(defined $db->{'clones'}->{uc($2)}) ) {
#      $db->dienice("$name isnt on a valid clone $2");
#      return undef;
#    }
  } else {
    $db->dienice("$type is a invalid typename: @types");
    return undef;
  }
  return $name;
}

sub validate_id
  {
    my $db = shift;
    my $id = shift;
    unless ( $db->idExists($id) ) {
      $db->dienice("$id does not exist");
      return undef ;
    }
    unless( $db->idLive($id) == 1) {
      $db->dienice("$id is not live");
      return undef;
    }
    return $id;
  }

sub check_pre_exists
  {
    my( $db, $name, $type) = @_;
    my $id = $db->idGetByTypedName($type => $name); #must be typed name as CDS and Sequence often have same name
    if (defined $id->[0] ) {
      $db->dienice("$name already exists as ".$id->[0]);
      return undef;
    }
    return $name;
  }


sub make_new_obj
  {
    my $db = shift;
    my $name = shift;
    my $type = shift;
    my $id = $db->create_named_object($type, $name);
    return $id;
  }

sub dienice {
  my $db = shift;
  my($errmsg) = @_;
  if ( $db->web ) {
    print "<h2>Error</h2>\n";
    print "<p>$errmsg</p>\n";
  } else {
    croak "$errmsg\n";
    return $errmsg;
  }
}

sub print_history {
  my $db = shift;
  my $id = shift;

  $db->validate_id($id);

  my $history = $db->idGetHistory($id);
  print "$id<br>";
  foreach my $event (@{$history}) {
    print $event->{version}," ";
    print $event->{date}," ";
    print $event->{event}," ";
    print $event->{name_type}," " if (defined $event->{name_type});
    print $event->{name_value}," " if (defined $event->{name_value});
    print $event->{user},"<br>";
  }
  return;
}


sub remove_all_names   {
  my $db = shift;
  my $id = shift;
  my %names = $db->idAllNames($id);
  foreach my $name_type (keys %names) {
    if ( ref $names{$name_type} eq 'ARRAY' ) {
      foreach my $name ( @{$names{$name_type}} ) {
	print "$name_type : $name\n";
	$db->delName($id,$name_type, $name);
      }
    } else {
      #print "$name_type : $names{$name_type}\n";
      $db->delName($id,$name_type, $names{$name_type});
    }
  }
  return;
}

 
sub merge_genes {
  my ($db, $gene_gene, $merge_gene) = @_;
  
  my ($gene_id, $merge_id);
  $gene_id  = ($db->idGetByAnyName($gene_gene)->[0]  or $gene_gene); #allow use of any name or gene_id
  $merge_id = ($db->idGetByAnyName($merge_gene)->[0] or $merge_gene); #allow use of any name or gene_id
  $db->validate_id($gene_id) or return undef;
  $db->validate_id($merge_id) or return undef;
  
  if ( $gene_id eq $merge_id) {
    $db->dienice("FAILED : $gene_gene and $merge_gene are the same!");
    return undef;
  }
  
  #enforce retention of CGC named gene
  my $names1 = $db->idAllNames($gene_id);
  my $names2 = $db->idAllNames($merge_id);
  
  unless( ($db->user eq 'pad') or ($db->user eq 'gw3' )){
    if ( $$names2{'CGC'} ) {
      if ($$names1{'CGC'} ) {
	#both genes have CGC name - confirm with geneace 
	$db->dienice("FAILED: Both genes have CGC names".$$names1{'CGC'}." and ".$$names2{'CGC'}.".  The correct course of action should be determined by the CGC admin (pad)<br>
					Please contact Geneace curator to resolve this");
      } else {
	#gene being eaten has a CGC name and eater doesnt
	$db->dienice("FAILED: $merge_gene has a CGC name ".$$names2{'CGC'}." and should probably be retained");
      }
      return undef;
    }
  }
  # warn that a gene with a CGC name has been killed
  if ( $$names2{'CGC'} ) { 
    print "$merge_id had a CGC name : $$names2{'CGC'}<br>";
  }
  #always remove names from merged id
  $db->remove_all_names($merge_id);
  #if this is a merger between a CGC gene and a sequence we need to transfer the Seq & CDS names to
  # the CGC gene id.
  unless ( $$names1{'Sequence'} ) {
    $db->addName($gene_id,'Sequence' => $$names2{'Sequence'}) if ($$names2{'Sequence'});
  }
  #do the merge
  if ($db->idMerge($merge_id,$gene_id)) {
    return ([$gene_id,$merge_id]);
  } else {
    $db->dienice("FAILED : merge failed");
    return undef;
  }
}


sub split_gene {
  my ($db,$name,$type,$gene,$species) = @_;
  unless ($type && $name && $gene){
    $db->dienice("bad parameters");
    return undef;
  }
  my $gene_id  = ($db->idGetByAnyName($gene)->[0]  or $gene); #allow use of any name or gene_id

  $db->validate_id($gene_id) or return undef;
  $db->validate_name($name, $type, $species) or return undef;
  $db->check_pre_exists($name, $type) or return undef;
  my $id = $db->idSplit($gene_id); ## hmm, does this work?
  if ( $id =~ /WBGene/ ) {
    $db->add_name($id, $name, $type, $species);
    return "$id";
  } else {
    $db->dienice("FAILED: spit gene $gene_id failed"); #error msg
    return undef;
  }
}


sub kill_gene {
  my ($db, $gene_id) = @_;
  $db->validate_id($gene_id) or return undef;
  if ($db->idKill($gene_id)) {
    return $gene_id;
  } else {	
    $db->dienice("cant kill GeneID $gene_id");
    return undef;
  }
}

sub remove_name {
  my ($db, $gene_id, $type, $name, $species) = @_;
  unless ($db and $gene_id and $type and $name) {
    $db->dienice("bad parameters");
    return undef;
  }
		
  $db->validate_id($gene_id) or return undef;
  $db->validate_name($name, $type, $species ) or return undef;
  my ($exist_id) = $db->idGetByTypedName($type,$name);
  if ( !$exist_id or ("$exist_id" ne "$gene_id") ) {
    $db->dienice("$name is not a name of $gene_id");
    return undef;
  }
  print "<br>GeneID = $gene_id<br>Type = $type<br>Name = $name<br>Species = $species<br>";
  if ( $db->delName($gene_id,$type,$name) ) {
    #need to update the public name - change to sequence name if CGC name removed
    if ($type eq 'CGC') {
      my $seq_name = $db->idTypedNames($gene_id,'Sequence');
      $db->add_name($gene_id,$seq_name->[0],'Public_name', $species);
    }
    return $name;
  }
}

sub new_gene {
  my ($db, $name, $type, $species) = @_;
	
  $db->validate_name($name, $type, $species)    or return undef;
  $db->check_pre_exists($name, $type) or return undef;

  my $id = $db->make_new_obj($name, $type);
  $db->_update_public_name($id);
  return $id ? $id : undef;
}
 
sub _update_public_name {
  my ($db, $id) = @_;
  my $names = $db->idAllNames($id);
  my $public = $names->{'Public_name'};
  my $new_public;
  if (!$public ) {
    if ( $names->{'CGC'} ) {
      $new_public = $names->{'CGC'};
    } elsif ( $names->{'Sequence'} ) {
      $new_public = $names->{'Sequence'};
    } else {
      $db->dienice("$id has no CGC or Sequence name");
      return undef;
      }
  } elsif ($names->{'CGC'} ) {
    if ($public ne $names->{'CGC'}) {
      $new_public = $names->{'CGC'};
    }
    } elsif ($names->{'Sequence'} and ($public ne $names->{'Sequence'})) {
      $new_public = $names->{'Sequence'};
    }
  $db->addName($id,'Public_name' => $new_public) if $new_public;
  return $new_public;
}


1;

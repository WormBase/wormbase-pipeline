#!/usr/local/bin/perl5.8.0 -w

# parse_sofa.pl

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-19 11:59:03 $

use strict;

#-----------------------
# parse SO definition
# ----------------------

my (%id_def, $id, $def_version);

open(IN, "/nfs/team71/worm/ck1/TMP/sofa\.definition") || die $!;

while(<IN>){
  my $def =();

  if ($_ =~ /version: \$Revision: 1.1 $/){$def_version = $1}
  if ($_ =~ /id: (.+)/){$id = $1}
  if ($_ =~ /definition: (.+)/){
    my @def =();
    @def = split(//, $1);
    foreach my $e (@def){
      $e = "\\". $e  if $e eq "\"";  # process quote " found in definition for Acedb
      $def .= $e;
    }
  }
  $id_def{$id} = $def if $id && $def;
}

#-----------------------
# parse SO ontology
# ----------------------

my ($so_term, $so_term_id, %hierarchy, $parellel_so_term, $parellel_so_term_id,
    $part_term, $part_term_id, $part_of_term_id, $part_of_term, $is_a_term, $is_a_term_id,
    $derived_term, $derived_term_id, $derived_from_term, $derived_from_term_id, @synonyms);

my ($level, @items, %level, $onto_version, $type);

open(IN, "/nfs/team71/worm/ck1/TMP/sofa.ontology") || die $!;

while(<IN>){

  chomp;
  my $part_of;

  if ($_ =~ /version: \$Revision: 1.1 $/){
    $onto_version = $1;
  }

  if ($_ =~ /^(\s{2,})(.+)/){

    $so_term = (); $so_term_id = ();                     # Is_a at the beginning
    $part_term = (); $part_term_id = ();                 # Part_of at beginning
    $derived_term = (); $derived_term_id = ();           # Derived_from at beginning
    $part_of_term = (); $part_term_id = ();              # Part_of in the middle 
    $is_a_term = (); $is_a_term_id = ();                 # Is_a in the middle 
    $derived_from_term = (); $derived_from_term_id = (); # Derived_from in the middle 

    $level = length($1);
    @items = split(/ ; /, $2);

    for (my $i = 0; $i < scalar @items; $i++){

      # parse SO term with @is_a@ at beginning of a line
      if ($items[0] =~ /^\@is_a\@(.+)/){
        $so_term = $1;
	
	if ($so_term eq "junction" or $so_term eq "region" or $so_term eq "sequence_variant"){
	  $type = $so_term;
	}
      }

      # parse SO term with @part_of@ at beginning of a line
      if ($items[0] =~ /^\@part_of\@(.+)/){
        $part_term = $1;
      }

       # parse SO term with @derived_from@ at beginning of a line
      if ($items[0] =~ /^\@derived_from\@(.+)/){
        $derived_term = $1;
      }

      # parse SO id for @is_a@ term at beginning of line
      if ($items[$i] =~ /(SO:\d+)/ && $so_term && !($part_of_term || $is_a_term || $derived_from_term) ){
        $so_term_id = $1;

	# assign version numbers of ontology and definition to only the top level terms
	push(@{$hierarchy{$so_term_id}{'Version'}}, $onto_version, $def_version) if $level == 2;

	# assign the one of the top 3 sequence features to each of their child terms
	$hierarchy{$so_term_id}{'Feature'} = $type if $so_term ne "located_sequence_feature";

	sort_out_hierarchy($so_term, $so_term_id, "Is");
      }

      # parse SO id for @part_of@ term at beginning of line
      if ($items[$i] =~ /(SO:\d+)/ && $part_term && !($part_of_term || $is_a_term || $derived_from_term) ){
        $part_term_id = $1;

	sort_out_hierarchy($part_term, $part_term_id, "Part");
      }

      # parse SO id for @derived_from@ term at beginning of line
      if ($items[$i] =~ /(SO:\d+)/ && $derived_term && !($part_of_term || $is_a_term || $derived_from_term) ){
        $derived_term_id = $1;

	sort_out_hierarchy($derived_term, $derived_term_id, "Derived");
      }

      # parse SO term of @part_of@ or @is_a@ or @derived_from@ in the middle of a line
      if ($items[$i] =~ /\@part_of\@ (.+)/ && ($so_term || $part_term || $derived_term) ){
        $part_of_term = $1;
      }
      elsif ($items[$i] =~ /\@is_a\@ (.+)/ && ($so_term || $part_term || $derived_term) ){
        $is_a_term = $1;
      }
      elsif ($items[$i] =~ /\@derived_from\@ (.+)/ && ($so_term || $part_term || $derived_term) ){
        $derived_from_term = $1;
      }

      # parse SO id after @part_of@ or @is_a@ or @derived_from@ in the middle of a line
      if ( $items[$i] =~ /(SO:\d+)/ && ($so_term || $part_term || $derived_term) &&
         ($part_of_term || $is_a_term || $derived_from_term) ){

	$part_of_term_id = $1 if $part_of_term;
	$is_a_term_id = $1 if $is_a_term;
	$derived_from_term_id = $1 if $derived_from_term;
	
	if ($so_term){
	  parse_inner_hierarchy($so_term_id);
        }

	if ($part_term){
	  parse_inner_hierarchy($part_term_id);
        }
	if ($derived_term){
	  parse_inner_hierarchy($derived_term_id);
        }

	$is_a_term = (); $is_a_term_id = ();
	$part_of_term = (); $part_of_term_id = ();
	$derived_from_term = (); $derived_from_term_id = ();
      }

      # parse synonym for SO term
      if ($items[$i] =~ /synonym:(.+)/){
        push(@{$hierarchy{$so_term_id}{'Synonym'}}, $1) if $so_term_id;
	push(@{$hierarchy{$part_term_id}{'Synonym'}}, $1) if $part_term_id;
	push(@{$hierarchy{$derived_term_id}{'Synonym'}}, $1) if $derived_term_id;
      }
    }
  }
}

# write SO_term ace file
open(ACE, ">sofa.ace") || die $!;

my $count_D =0;

foreach (keys %hierarchy){

  print ACE "\n\nSO_term : \"$_\"\n";
  print ACE "Term \"$hierarchy{$_}{'Term'}\"\n" if exists $hierarchy{$_}{'Term'};

  if (exists $id_def{$_}){
    print ACE "Definition \"$id_def{$_}\"\n";
  }
  else {
    $count_D++;  # output missing definitions in sofa.definition flatfile (if any)
  }
  if (exists $hierarchy{$_}{'Feature'}){
    print ACE "$hierarchy{$_}{'Feature'}\n";
  }
  if (exists $hierarchy{$_}{'Synonym'}){
    foreach my $e (@{$hierarchy{$_}{'Synonym'}}){
      print ACE "Synonym \"$e\"\n";
    }
  }

  if (exists $hierarchy{$_}{'Child_Is'}){
    foreach my $e (@{$hierarchy{$_}{'Child_Is'}}){
      print ACE "Is \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Child_Part'}){
    foreach my $e (@{$hierarchy{$_}{'Child_Part'}}){
      print ACE "Part \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Child_Derived'}){
    foreach my $e (@{$hierarchy{$_}{'Child_Derived'}}){
      print ACE "Derived \"$e\"\n";
    }
  }

  if (exists $hierarchy{$_}{'Parent_Is_a'}){
    foreach my $e (@{$hierarchy{$_}{'Parent_Is_a'}}){
      print ACE "Is_a \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Parent_Part_of'}){
    foreach my $e (@{$hierarchy{$_}{'Parent_Part_of'}}){
      print ACE "Part_of \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Parent_Derived_from'}){
    foreach my $e (@{$hierarchy{$_}{'Parent_Derived_from'}}){
      print ACE "Derived_from \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Ancestor'}){
    foreach my $e (@{$hierarchy{$_}{'Ancestor'}}){
      print ACE "Ancestor \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Descendent'}){
    foreach my $e (@{$hierarchy{$_}{'Descendent'}}){
      print ACE "Descendent \"$e\"\n";
    }
  }
  if (exists $hierarchy{$_}{'Version'}){
    print ACE "sofa_ontology \"$hierarchy{$_}{'Version'}->[0]\"\n";
    print ACE "sofa_definition  \"$hierarchy{$_}{'Version'}->[1]\"\n";
  }
}

print "Of ", scalar keys %hierarchy, " SO terms, $count_D do not have definition\n";


##########################
# s u b r o u t i n e s
##########################

sub sort_out_hierarchy {
  my ($term, $term_id, $type) = @_;

  $hierarchy{$term_id}{'Term'} = $term;

  #assign SO term id per level of depth
  $level{$level} = $term_id;

  if ( exists $level{$level-1} ){ # up one level
    foreach (my $i = 2; $i < $level; $i++){
      push(@{$hierarchy{$term_id}{'Ancestor'}}, $level{$i} );
      push(@{$hierarchy{$level{$i}}{'Descendent'}}, $term_id );
    }
    push(@{$hierarchy{$term_id}{'Parent_Is_a'}}, $level{$level-1} )                 if $type eq "Is";
    push(@{$hierarchy{$level{$level-1}}{'Child_Is'}}, $term_id )                    if $type eq "Is";
    push(@{$hierarchy{$part_term_id}{'Parent_Part_of'}}, $level{$level-1} )         if $type eq "Part";
    push(@{$hierarchy{$level{$level-1}}{'Child_Part'}}, $part_term_id)              if $type eq "Part";
    push(@{$hierarchy{$derived_term_id}{'Parent_Derived_from'}}, $level{$level-1} ) if $type eq "Derived";
    push(@{$hierarchy{$level{$level-1}}{'Child_Derived'}}, $derived_term_id)        if $type eq "Derived";
  }
}

sub parse_inner_hierarchy {
  my $term_id = shift;
  push(@{$hierarchy{$is_a_term_id}{'Child_Is'}}, $term_id)                    if $is_a_term_id;
  push(@{$hierarchy{$part_of_term_id}{'Child_Part'}}, $term_id)               if $part_of_term_id;
  push(@{$hierarchy{$derived_from_term_id}{'Child_Derived'}},$term_id)        if $derived_from_term_id; 
  push(@{$hierarchy{$term_id}{'Parent_Is_a'}}, $is_a_term_id)                 if $is_a_term_id;
  push(@{$hierarchy{$term_id}{'Parent_Part_of'}}, $part_of_term_id)           if $part_of_term_id;
  push(@{$hierarchy{$term_id}{'Parent_Derived_from'}}, $derived_from_term_id) if $derived_from_term_id;
}




__END__


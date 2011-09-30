#!/usr/bin/env perl 

# parse_SO_obo.pl

########################################################
# This script parses a SO obo file into Ace objects. It
# can be run either to parse the full SO (not normal), 
# or the sequence_feature subset (more normal), formerly
# known as the "SOFA"
#######################################################

use strict;
use Getopt::Long;

my (
    $in_term,
    $current_term,
    $so_version,
    $full_so,
    %terms);

&GetOptions('soversion=s' => \$so_version,
            'fullso'      => \$full_so,
            );


while(<>) {
  /^\[(\S+)\]/ and do {
    $in_term = ($1 eq 'Term') ? 1 : 0;        
    next;
  };

  next if not $in_term;

  /^id:\s*(\S+)/ and do {
    $current_term = $1;
    $terms{$current_term}->{acc} = $current_term;
    next;
  };

  /^name:\s*(\S+)/ and do {
    $terms{$current_term}->{name} = $1;
    next;
  };

  /^def:\s*\"(.+)\"\s+\[\S+\]\s*$/ and do {
    $terms{$current_term}->{def} = $1;
    next;
  };

  /^synonym:\s*\"(.+)\"\s+\S+\s+\[.*\]/ and do {
    push @{$terms{$current_term}->{synonyms}}, $1;
    next;
  };

  /^is_a:\s*(\S+)/ and do {
    $terms{$current_term}->{isa}->{$1} = 1;
    next;
  };

  /^relationship:\s*(\S+)\s+(\S+)/ and do {
    my ($rel_type, $rel_to) = ($1, $2);

    if ($rel_type eq 'part_of' or 
        $rel_type eq 'member_of' or
        $rel_type eq 'derives_from') {
      $terms{$current_term}->{$rel_type}->{$rel_to} = 1;
    }
    next;
  };

  /^is_obsolete:\s*true/ and do {
    $terms{$current_term}->{is_obsolete} = 1;
    next;
  };
}

foreach my $k (keys %terms) {
  next if $terms{$k}->{is_obsolete};

  if (not exists $terms{$k}->{isa}) {
    $terms{$k}->{is_root} = 1;
  }
}

foreach my $acc (sort keys %terms) {

  # in order to derive the Ancestor/Descendent relationship, 
  # we must traverse the tree downwards
  next if $terms{$acc}->{is_obsolete};

  my (%ancestors, $root_ancestor, $root_ancestor_child);
  my @list = ($acc);

  do {
    my %next_level_up;

    foreach my $k (@list) {
      
      if (exists $terms{$k}->{isa}) {
        foreach my $ot (keys %{$terms{$k}->{isa}}) {
          $ancestors{$ot} = 1;
          $next_level_up{$ot} = 1;

          if (exists $terms{$ot}->{is_root}) {
            # immediate parent of current term is root, so current term is one
            # level below
            $root_ancestor_child = $k;
            $root_ancestor = $ot;
          }
        }
      }
    }

    @list = sort keys %next_level_up;
  } while (@list);

  if (keys %ancestors) {
    $terms{$acc}->{ancestors} = \%ancestors;
  }
  if (defined $root_ancestor_child) {
    $terms{$acc}->{root_ancestor_child} = $root_ancestor_child;
  }
  if (defined $root_ancestor) {
    $terms{$acc}->{root_ancestor} = $root_ancestor;
  }

  &print_term($terms{$acc}, \*STDOUT);
}



sub print_term {
  my ($term, $fh) = @_;

  return if not exists $term->{root_ancestor};

  return if $terms{$term->{root_ancestor}}->{name} ne 'sequence_feature' and not $full_so;

  printf $fh "\nSO_term : \"%s\"\n", $term->{acc};
  printf $fh "SO_name \"%s\"\n", $term->{name};
  if ($term->{def}) {
    printf $fh "SO_definition \"%s\"\n", $term->{def};
  }
  printf $fh "Located_sequence_feature\t%s\n", ucfirst($terms{$term->{root_ancestor_child}}->{name});
  if ($term->{synonyms}) {
    foreach my $syn (@{$term->{synonyms}}) {
      printf $fh "SO_synonym\t\"%s\"\n", $syn;
    }
  }
  if ($term->{isa}) {
    foreach my $isa (sort keys %{$term->{isa}}) {
      next if $terms{$isa}->{is_root};
      printf $fh "Is_a\t\"%s\"\n", $isa;
    }
  }
  if ($term->{part_of}) {
    foreach my $po (sort keys %{$term->{part_of}}) {
      next if $terms{$po}->{is_root};
      printf $fh "Part_of\t\"%s\"\n", $po;
    }
  }
  if ($term->{member_of}) {
    foreach my $mo (sort keys %{$term->{member_of}}) {
      next if $terms{$mo}->{is_root};
      printf $fh "Member_of\t\"%s\"\n", $mo;
    }
  }
  if ($term->{derived_from}) {
    foreach my $df (sort keys %{$term->{derived_from}}) {
      next if $terms{$df}->{is_root};
      printf $fh "Derived_from\t\"%s\"\n", $df;
    }
  }
  if (exists $term->{ancestors}) {
    foreach my $an (sort keys %{$term->{ancestors}}) {
      next if $terms{$an}->{is_root};
      printf $fh "Ancestor\t\"%s\"\n", $an;
    }
  }

  if ($so_version) {
    my @non_root_parents = grep { not $terms{$_}->{is_root} } keys %{$term->{isa}};
    if (not @non_root_parents) {
      printf $fh "SO_version \"$so_version\"\n";
    } 
  }
}

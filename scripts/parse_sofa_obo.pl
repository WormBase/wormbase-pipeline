#!/software/bin/perl

use strict;
use Getopt::Long;

my ($sofa_version, 
    $in_term,
    $current_term,
    %terms);

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

my $root_term;
foreach my $k (keys %terms) {
  next if $terms{$k}->{is_obsolete};

  if (not exists $terms{$k}->{isa}) {
    if (not defined $root_term) {
      $terms{$k}->{is_root} = 1;
      $root_term = $k;
    } else {
      die "Should not be two root terms: $k, $root_term\n";
    }
  }
}
if (not defined $root_term) {
  die "Could not find any terms with no ancestors\n";
}


foreach my $acc (sort keys %terms) {

  # in order to derive the Ancestor/Descendent relationship, 
  # we must traverse the tree downwards
  next if $terms{$acc}->{is_obsolete};

  my (%ancestors, $root_child_ancestor);
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
            $root_child_ancestor = $k;
          }
        }
      }
    }

    @list = sort keys %next_level_up;
  } while (@list);

  if (keys %ancestors) {
    $terms{$acc}->{ancestors} = \%ancestors;
  }
  if (defined $root_child_ancestor) {
    $terms{$acc}->{root_child_ancestor} = $root_child_ancestor;;
  } 

  &print_term($terms{$acc}, \*STDOUT);
}



sub print_term {
  my ($term, $fh) = @_;

  return if $term->{acc} eq $root_term;

  printf $fh "\nSO_term : \"%s\"\n", $term->{acc};
  printf $fh "SO_name \"%s\"\n", $term->{name};
  if ($term->{def}) {
    printf $fh "SO_definition \"%s\"\n", $term->{def};
  }
  printf $fh "Located_sequence_feature\t%s\n", ucfirst($terms{$term->{root_child_ancestor}}->{name});
  if ($term->{synonyms}) {
    foreach my $syn (@{$term->{synonyms}}) {
      printf $fh "Synonym\t\"%s\"\n", $syn;
    }
  }
  if ($term->{isa}) {
    foreach my $isa (sort keys %{$term->{isa}}) {
      next if $isa eq $root_term;
      printf $fh "Is_a\t\"%s\"\n", $isa;
    }
  }
  if ($term->{part_of}) {
    foreach my $po (sort keys %{$term->{part_of}}) {
      next if $po eq $root_term;
      printf $fh "Part_of\t\"%s\"\n", $po;
    }
  }
  if ($term->{member_of}) {
    foreach my $mo (sort keys %{$term->{member_of}}) {
      next if $mo eq $root_term;
      printf $fh "Member_of\t\"%s\"\n", $mo;
    }
  }
  if ($term->{derived_from}) {
    foreach my $df (sort keys %{$term->{derived_from}}) {
      next if $df eq $root_term;
      printf $fh "Derived_from\t\"%s\"\n", $df;
    }
  }
  if (exists $term->{ancestors}) {
    foreach my $an (sort keys %{$term->{ancestors}}) {
      next if $an eq $root_term;
      printf $fh "Ancestor\t\"%s\"\n", $an;
    }
  }

}

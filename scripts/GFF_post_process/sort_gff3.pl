#!/usr/bin/env perl

#
# 2-pass algorithm:
#  - pass 1 : write a temporary file where all parent features come before their children; whilst
#    doing this pass, keep track of file locations of every feature with an ID
#  - pass 2 : parse the just-written temp file, and every time we encounter a feature with
#    and ID, write all features with that ID. This ensures that split features are grouped 
#


use File::Temp qw/ tempfile /;

# Keeps track which features have been serialized already (based on ID attribute).
my $serialized_ids = {};

# Mapping of feature IDs to their parent feature IDs (value: list; based on ID and Parent attributes).
my $parent_rels = {};

# Growing/shrinking set of features whose serialization is "on hold" because their parent feature has not been serialized yet.
my $feature_stash = {};

# Counter used to make up on-the-fly feature IDs for features that do not have one (i.e., lacking an ID attribute).
my $faux_feature_id = 0;

# Index of feature locations in the semi-sorted temporary file (see "$semisorted_gff3").
my $feature_index = {};

# Temporary file that is used to write out features that have been sorted by parent/dependent relationship, but which
# have not been grouped together by feature ID yet.
my $semisorted_gff3 = tempfile();

# First pass:
#   - write out into a temporary file
#   - serialize features so that "Parent" features come before its dependents
#   - create an index of serialized feature IDs
my $offset = 0;
my $last_seq;

while(<>) {
  chomp;
  my $line = $_;
  
  if ($line =~ /^#/) {
    print $semisorted_gff3 "$line\n";
  } elsif ($line =~ /^([a-zA-Z0-9.:^*$@!+_?-|]|%[0-9A-F]{2})+\s/) {
    my ($seq) = $line =~ /^(\S+)/;
    my ($id) = $line =~ /ID=([^;]+)/;
    my ($parent) = $line =~ /Parent=([^;]+)/;
    
    if (defined $last_seq and $seq ne $last_seq) {
      # Since we have now seen all the features for $last_seq, we can
      # compress the index 
      foreach my $idxid (keys %{$feature_index->{$last_seq}}) {
        if (scalar(@{$feature_index->{$last_seq}->{$idxid}}) == 1) {
          delete $feature_index->{$last_seq}->{$idxid};
        }
      }
      $serialized_ids = {};
    }
    $last_seq = $seq;

    if (not defined $parent) {
      # No parents. Just print the feature and record its ID.
      if (defined $id) {
        $serialized_ids->{$id} = 1;
        index_feature($feature_index, $seq, $id, $offset);
      }
      print $semisorted_gff3 "$line\n";
    } else {
      my @parents = split(/,/, $parent);
      my $parents_serialized = 0;
      foreach my $parent (@parents) {
        $parents_serialized++ if $serialized_ids->{$parent};
      }
      if ($parents_serialized == scalar(@parents)) {
        # Parents already have been serialized, so just write out the feature.
        if (defined $id) {
          $serialized_ids->{$id} = 1;
          index_feature($feature_index, $seq, $id, $offset);
        }
        print $semisorted_gff3 "$line\n";
      } else {
        if (defined $id) {
          # Stash the current input line for future serialization.
          $parent_rels->{$id} = \@parents;
          $feature_stash->{$id} = $line;
        } else {
          # No ID that we could keep track of. Create a temporary ID for the feature.
          $parent_rels->{"% $faux_feature_id"} = \@parents;
          $feature_stash->{"% $faux_feature_id"} = $line;
          $faux_feature_id++;
        }
      }
    }
    
    my $feature_serialized;
    do {
      $feature_serialized = undef;
      my @serialized_features = ();
      while(my ($feature_id, $value) = each %$feature_stash) {
        next unless $value;
        next if $serialized_ids->{$feature_id};
        my $parents_serialized = 0;
        my $parents = $parent_rels->{$feature_id};
        foreach my $parent (@$parents) {
          $parents_serialized++ if $serialized_ids->{$parent};
        }
        if ($parents_serialized == scalar(@$parents)) {
          # Parents have been serialized. Get stashed feature out too.
          $offset = tell($semisorted_gff3);
          index_feature($feature_index, $seq, $feature_id, $offset);
          print $semisorted_gff3 $feature_stash->{$feature_id} . "\n";
          $serialized_ids->{$feature_id} = 1;
          push(@serialized_features, $feature_id);
          $feature_serialized = 1;
        }
      }
      foreach my $feature_id (@serialized_features) {
        delete $feature_stash->{$feature_id};
      }
    } while($feature_serialized);
  } else {
    # Some other content. Just print it.
    print $semisorted_gff3 "$line\n";
  }
  
  $offset = tell($semisorted_gff3);
}

# Second pass:
#   - go through semi-sorted features in the temporary file
#   - if an indexed feature is encountered, output all its other occurrences (by feature ID) immediately as well
seek($semisorted_gff3, 0, 0);
while(<$semisorted_gff3>) {
  chomp;
  my $line = $_;
  
  if ($line =~ /^#/) {
    print "$line\n";
  } elsif ($line =~ /^([a-zA-Z0-9.:^*$@!+_?-|]|%[0-9A-F]{2})+\s/) {
    my ($seq) = $line =~ /^(\S+)/;
    my ($id) = $line =~ /ID=([^;]+)/;
    
    if (defined $id and exists $feature_index->{$seq}->{$id}) {
      # Print the current line, save the read offset, print the indexed lines, 
      # restore the read offset, mark that we serialized the feature, and proceed...
      my $current_position = tell($semisorted_gff3);
      my $indexed_positions = $feature_index->{$seq}->{$id};
      foreach my $indexed_position (@$indexed_positions) {
        seek($semisorted_gff3, $indexed_position, 0);
        my $indexed_line = <$semisorted_gff3>;
        print $indexed_line;
      }
      seek($semisorted_gff3, $current_position, 0);
      $feature_index->{$seq}->{$id} = [];
    } else {
      # Just print the line at hand.
      print "$line\n";
    }
  } else {
    # Again, some other content. Just print it.
    print "$line\n";
  }
}



#########################################

sub index_feature {
  my ($feature_index, $feature_seq, $feature_id, $offset) = @_;
  
  # Do not index features with faux IDs.
  return if ($feature_id =~ /^% /);
  
  push @{$feature_index->{$feature_seq}->{$feature_id}}, $offset;
}

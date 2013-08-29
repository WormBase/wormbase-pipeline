#!/usr/bin/perl

use Data::Dumper;

use File::Temp qw/ tempfile /;

sub index_feature {
    my ($file_handle, $feature_index, $feature_id, $offset) = @_;

    # Do not index features with faux IDs.
    return if ($feature_id =~ /^% /);

    my $index = $feature_index->{$feature_id};
    $index = () unless $index;
    push(@$index, $offset);
    $feature_index->{$feature_id} = $index;
}

# Keeps track which features have been serialized already (based on ID attribute).
my $serialized_ids = {};

# Mapping of feature IDs to their parent feature IDs (value: list; based on ID and Parent attributes).
my $parent_rels = {};

# Growing/shrinking set of features whose serialization is "on hold" because their parent feature has not been serialized yet.
my $feature_stash = {};

# Counter used to make up on-the-fly feature IDs for features that do not have one (i.e., lacking an ID attribute).
my $faux_feature_id = 0;

# Indicates whether we moved into the "FASTA" section of a GFF3 file, which means that all features have been seen already.
my $fasta_annex = undef;

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
while(<>) {
    chomp;
    my $line = $_;

    if ($fasta_annex) {
        # We are in the FASTA section. Just output the sequences as they float by...
        print $semisorted_gff3 "$line\n";
    } elsif ($line =~ /^#/) {
        $fasta_annex = 1 if ($line =~ /^##FASTA/); # When the FASTA section starts, mark this, because every following line can be just printed.
        print $semisorted_gff3 "$line\n";
    } elsif ($line =~ /^([a-zA-Z0-9.:^*$@!+_?-|]|%[0-9A-F]{2})+\s/) {
        my @ids = $line =~ /ID=([^;]+)/g;
        my @parents = $line =~ /Parent=([^;]+)/g;

        unless (@parents) {
            # No parents. Just print the feature and record its ID.
            if (@ids && scalar(@ids) == 1) {
                $serialized_ids->{$ids[0]} = 1;
                index_feature($semisorted_gff3, $feature_index, $ids[0], $offset);
            }
            print $semisorted_gff3 "$line\n";
        } else {
            @parents = split(/,/, $parents[0]);
            my $parents_serialized = 0;
            foreach my $parent (@parents) {
                $parents_serialized++ if $serialized_ids->{$parent};
            }
            if ($parents_serialized == scalar(@parents)) {
                # Parents already have been serialized, so just write out the feature.
                if (@ids && scalar(@ids) == 1) {
                    $serialized_ids->{$ids[0]} = 1;
                    index_feature($semisorted_gff3, $feature_index, $ids[0], $offset);
                }
                print $semisorted_gff3 "$line\n";
            } else {
                if (@ids && scalar(@ids) == 1) {
                    # Stash the current input line for future serialization.
                    $parent_rels->{$ids[0]} = \@parents;
                    $feature_stash->{$ids[0]} = $line;
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
                    index_feature($semisorted_gff3, $feature_index, $feature_id, $offset);
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
$fasta_annex = undef;
seek($semisorted_gff3, 0, 0);
while(<$semisorted_gff3>) {
    chomp;
    my $line = $_;

    if ($fasta_annex) {
        # Again, we are in the FASTA section. Just output the sequences as they float by...
        print "$line\n";
    } elsif ($line =~ /^#/) {
        $fasta_annex = 1 if ($line =~ /^##FASTA/); # Once more, when the FASTA section starts, mark this, because every following line can be just printed.
        print "$line\n";
    } elsif ($line =~ /^([a-zA-Z0-9.:^*$@!+_?-|]|%[0-9A-F]{2})+\s/) {
        my @ids = $line =~ /ID=([^;]+)/g;

        if (@ids && scalar(@ids) == 1) {
            # Only go ahead when the feature has not been serialized yet (based on ID).
            if ($feature_index->{@ids[0]}) {
                # Print the current line, save the read offset, print the indexed lines, restore the read offset, mark that we serialized the feature, and proceed...
                my $current_position = tell($semisorted_gff3);
                my $indexed_positions = $feature_index->{@ids[0]};
                foreach my $indexed_position (@$indexed_positions) {
                    seek($semisorted_gff3, $indexed_position, 0);
                    my $indexed_line = <$semisorted_gff3>;
                    print $indexed_line;
                }
                seek($semisorted_gff3, $current_position, 0);
                delete $feature_index->{@ids[0]};
            }
        } else {
            # Just print the line at hand.
            print "$line\n";
        }
    } else {
        # Again, some other content. Just print it.
        print "$line\n";
    }
}


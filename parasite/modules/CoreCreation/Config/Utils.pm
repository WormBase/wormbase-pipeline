package CoreCreation::Config::Utils;

=pod

=head1 SYNOPSIS

Miscellaneous functions for handling configuration during the ParaSite core database creation process.

Mainly created as a placeholder for validation/verifications functions that we could really do with adding.
No attempt yet to rationalize with CoreCreation::Config.

=head1 Functions
   
=cut

use strict;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(parasite_data_id parasite_data_id_is_valid flatten_hash);

use constant PARASITE_DATA_FILTER      => qr/^[a-z\d]+_[a-z\d]+_[a-z\d]+\d+$/;
use constant BIOPROJECT_FILTER         => qr/^[a-z]+\d+$/;
use constant MISSING_METADATA_PATTERN  => '?';

use constant ARRAY_REF_TYPE      => ref([]);
use constant HASH_REF_TYPE       => ref({});

=head2 parasite_data_id

Creates conventional ParaSite identifier for data of the form C<genus_species_bioproject>

=head3 Arguments

Two strings:  organism name, and bioproject identifier.

=head3 Return values.

Three strings: identifier that follows the ParaSite convention, the genus and the species.
(You probably only want the ID, but the genus and species are returned in case you want to
check how the organism name was parsed.)

=cut
sub parasite_data_id {
   my $this = (caller(0))[3];
   my $organism   = shift() || "$this must be passed an organism name";
   my $bioproject = shift() || "$this must be passed a bioproject identifier";
   
   $organism =~ s/sp. //;
   $organism =~ tr/ /_/;
   $organism = lc $organism;
   my ($genus, $species) = split (/_/, $organism);
   
   $bioproject = lc $bioproject;

   my $new_id  = parasite_data_id_is_valid( join("_", $genus, $species, $bioproject) )
               || die "couldn't derive a ParaSite data id from organism \"$organism\" and bioproject \"$bioproject\"";

   return($new_id, $genus, $species);
}

=head2 parasite_data_id_is_valid

Tests if s string conforms to the convention for a ParaSite data identifier, i.e. C<genus_species_bioproject>

=head3 Arguments

A string that is alleged to be an identifier that follows the ParaSite convention.

=head3 Return values.

Returns the string passed, if it follows the ParaSite convention; otherwise undef.

=cut
sub parasite_data_id_is_valid {
   my $this = (caller(0))[3];
   my $id = shift() || "$this must be passed an identifier";
   return($id =~ PARASITE_DATA_FILTER ? $id : undef);
}


=head2 flatten_hash

Function to flatten a hash.  Can be useful for quick & dirty checking of
values without needing code to recurse through an arbitrary structure.

Flattens any nested hashes or arrays by creating a top level key based on 
joining the sequence of hash keys/array indexes, with '/' as separator.

Any references found that aren't hashes or arrays are just converted into string 
representation

=head3 Arguments

A reference to a hash.

=head3 Return values.

A reference to the (now flattend) hash

=cut
sub flatten_hash {
   my $this = (caller(0))[3];
   my $hash_ref = shift() || "$this must be passed a hash name";
   
   # do..while iterates through all top level hash values until no references are found
   my $nested_ref;
   my $max_levels = 100;
   do {
      $nested_ref = 0;
      --$max_levels < 0 && die "insane number of nested levels in hash";
      foreach my $this_key (keys %{$hash_ref}) {
         my $this_value = $hash_ref->{$this_key};
         if( HASH_REF_TYPE eq ref($this_value) ) {
            map {$hash_ref->{join('/',$this_key,$_)} = $this_value->{$_}} (keys %{$this_value});
            delete $hash_ref->{$this_key};
            ++$nested_ref;
         } elsif (ARRAY_REF_TYPE eq ref($this_value) ) {
            my $i=0;
            map {$hash_ref->{join('/',$this_key,$i++)} = $_} (@{$this_value});
            delete $hash_ref->{$this_key};
            ++$nested_ref;
         } elsif ( ref($this_value) ) {
            # some other reference
            $hash_ref->{$this_key} = "$this_value";
            ++$nested_ref;
         }
      }
   } while ($nested_ref);
   
   return($hash_ref);
}

1;

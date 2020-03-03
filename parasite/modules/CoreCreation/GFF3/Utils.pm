package CoreCreation::GFF3::Utils;

=pod

=head1 SYNOPSIS

Miscellaneous functions for handling GFF3 during the ParaSite core database creation process.

   use CoreCreation::GFF3::Utils qw(validate substitute);
   use Try::Tiny; # merely a suggestion

   try {
      validate( 'foo.gff3' );
   } catch {
      die "foo.gff3 is not valid GFF3: $_";
   };

   my $new_gff3 = CoreCreation::GFF3::Utils::munge_line_by_line(
                     -file                => 'foo.gff3',
                     -strip_comments      => 1,
                     -strip_fasta         => 1,
                     -set_source          => 'WormBase_imported',
                     -handler             => sub{ my_transform_function(@_) },
                     );
                     
=head1 Dependencies

=over

=item Genometools (gt) excuteable in your $PATH

=item CPAN packages:

=over

=item File::Slurp

=item Storable

=item Sub::Identify

=item Try::Tiny

=item URI::Escape

=back

=back

=head1 Functions
   
=cut

use strict;

use File::Slurp;
use Storable;
use Sub::Identify;
use Try::Tiny;
use URI::Escape;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(validate substitute);

use constant GFF3_SEQID_FIELD    => 0;
use constant GFF3_SOURCE_FIELD   => 1;
use constant GFF3_TYPE_FIELD     => 2;
use constant GFF3_ATTR_FIELD     => 8;
use constant GFF3_LAST_FIELD     => 8;

use constant GFF3_VALIDATION_CMD => 'gt gff3validator';

use constant ARRAY_REF_TYPE      => ref([]);
use constant HASH_REF_TYPE       => ref({});

=head2 validate

Validates GFF3 file.  Dies if the validation cannot be executed.

=head3 Arguments

Name of GFF3 file.

=head3 Return values.

A flag indicating success (true) or failure (false) of the validation,
followed by a (hopefully empty) list of validation error messages.

=cut
sub validate {
   my $this = (caller(0))[3];
   my $gff3_file = shift() || die "$this must be passed a GFF3 file name";

   -e $gff3_file || die "$gff3_file does not exist";
   -s $gff3_file || die "$gff3_file is empty";
   open(GT, GFF3_VALIDATION_CMD." $gff3_file 2>&1 |") || die "can't execute ".GFF3_VALIDATION_CMD.": $!";
   my @validation_errors = map {chomp; $_} <GT>;
   # exit status is false on validation error *or* command error
   my $success = close(GT);
   # $! defined on command error
   my $msg = $!;
   if($success) {
      @validation_errors = ();
   } elsif($msg) {
      die "error running ".GFF3_VALIDATION_CMD.": $msg";
   }

   return($success, @validation_errors);
}

=head2 munge_line_by_line

Iterates through GFF3 file, line by line, and makes changes defined by
the parameters passed.  The (possibly modified) GFF3 lines are returned
as a list.  Note that lines returned have been chomp'ed.

=head3 Arguments

This function requies named parameters.

=over

=item -file [mandatory]

Name of the GFF3 file

=item -strip_comments

Flag. If true, all comment lines are removed from the GFF3. These lines are
deleted prior to any further processing.

=item -strip_fasta

Flag. If true, if a FASTA directive (##FASTA) is found then that line and
all remaining lines are deleted.

=item -strip_seqid

List.  Delete all features with a seqid in the list.  Alternatively pass
a file name; the file will be read to provide the list

=item -seqid_ncbi_to_submitter

String.  This should be a tab-separated file with submitter's
seqid values in field 2 and NCBI seqid in field 3. (This is envisaged to
be the C<.seq_region_synonyms.tsv> file created by C<prepare_data_folder.pl>).
If the seqid matches the NCBI value, it will be changed to the
submitter value.

=item -set_source

String.  The source (field 2) of each feature is set to this value.

=item -change_type

Hash of params.  If defined, changes a feature type.
Params:

=over

=item -from

Current feature type

=item -to

New feature type

=back

Can also past a list of hashes to change more than one feature type.

=item -name_gene_from_id

Hash of params.  If defined, creates a Name attribute from the ID attribute.
Params:

=over

=item -types

List of types for which new Names are to be created.  Defaults to all types.

=item -prefix, -suffix

Prefix and/or suffix that will be added to the ID.strip_seqid
If neither of theseare passed, the Name will be a copy of the ID.

=item -separator

Separator used after prefix/before suffix; defaults to '.'

=back

=item -copy_CDS_to_exon
1
Flag.  If true, CDS features are duplicated as additional exon features with
the ID tweaked accordingly.

=item -handler

Reference to a handler function. This is called for every feature line. 
The handler is passed a list of 9 fields; the first 8 fields are strings, the
9th deserializes attributes as a hash.  Attribute values are unescaped (%09
converted to a space, etc.), and multiple values are represented as a list.

The handler should return a list of 9 fields in the same form, which will be
turned into a string and substituted in place of the original line.  Note
attribute values should I<not> be escaped, and multiple values should be
in the form on a list (i.e. I<not> as a single string of comma-separated
values).

=back

The arguments cause GFF3 features to be modified are applied in the following order:

   -seqid_ncbi_to_submitter
   -set_source
   -change_type
   -name_gene_from_id
   -handler
   -copy_CDS_to_exon

=head3 Return values.

A reference to a list of new GFF3 lines.

=cut
sub munge_line_by_line {
   my $this = (caller(0))[3];
   my %params = HASH_REF_TYPE eq ref(@_[0]) ? %{@_[0]} : @_;
   
   my $gff3_file           = $params{-file} || die "$this must be passed a GFF3 file name (-file)";
   my $strip_comments      = $params{-strip_comments};
   my $strip_fasta         = $params{-strip_fasta};
   my $strip_seqid         = $params{-strip_seqid};
   my $seqid_ncbi_to_submitter = $params{-seqid_ncbi_to_submitter};
   my $set_source          = $params{-set_source};
   my $change_type         = $params{-change_type};
   my $name_gene_from_id   = $params{-name_gene_from_id};
   my $copy_CDS_to_exon    = $params{-copy_CDS_to_exon};
   my $handler             = $params{-handler};

   # if removing seqids, build hash for lookup of seqids to be removed
   if($strip_seqid) {
      if(ARRAY_REF_TYPE ne ref($strip_seqid)) {
         die "$this parameter -strip_seqid must be a list reference, or a file name"
            unless -e $strip_seqid;
         $strip_seqid = [ File::Slurp::read_file($strip_seqid, chomp=>1) ];
      }
      # convert to hash for fast lookup
      my %hash;
      @hash{@{$strip_seqid}} = (1) x @{$strip_seqid};
      $strip_seqid = \%hash;
   }
   
   # if changing NCBI seqid to submitter's values, build lookup
   my %ncbi_to_submitter_seqid = ();
   if($seqid_ncbi_to_submitter) {
      foreach my $line ( File::Slurp::read_file($seqid_ncbi_to_submitter, chomp=>1) ) {
         my @f = split(/\t/, $line, 4);
         my $sub  = $f[1] || die "seq region synonyms file has blank/missing submitter seqid value";
         my $ncbi = $f[2] || die "seq region synonyms file has blank/missing submitter NCBI value";
         $ncbi_to_submitter_seqid{$ncbi} = $sub;
      }
   }
   
   # validate -change-type parameters
   if($change_type) {
      die "$this parameter -change_type must be a hash reference, or list of hashes"
         if (HASH_REF_TYPE ne ref($change_type) && ARRAY_REF_TYPE ne ref($change_type))
         || (ARRAY_REF_TYPE eq ref($change_type) && grep(HASH_REF_TYPE ne ref($_), @{$change_type}));
      die "$this parameter -change_type hash(es) must include '-from' and '-to'"
         if grep(!exists $_->{-from} || !exists $_->{-to}, (ARRAY_REF_TYPE eq ref($change_type) ? @{$change_type} : $change_type));
   }
   
   die "$this parameter -name_gene_from_id must be a hash reference"
      if $name_gene_from_id && HASH_REF_TYPE ne ref($name_gene_from_id);

   die "$this parameter -handler must be a code reference"
      if $handler && ref(sub{}) ne ref($handler);
   
   my $status = '';
   my $num_lines = 0;
   my $num_headers = 0;
   my $num_comments = 0;
   my $num_features = 0;
   my $num_lines_fasta = 0;
   my $reading_fasta = 0;
   my @new_gff3 = ();
   open(GFF3, $gff3_file) || die die "can't read $gff3_file: $!";
   while(my $line = <GFF3>) {
      chomp($line);
      
      if($line =~ m/^##\s*FASTA\s*$/i) {
         # $line is a FASTA directive
         my @fasta_lines = map {chomp; $_} <GFF3>;
         push(@new_gff3,$line,@fasta_lines) unless $strip_fasta;
         $num_lines_fasta += 1 + scalar(@fasta_lines);
      } elsif($line =~ m/^##[^#]/) {
         # $line is a header
         push(@new_gff3,$line);
         ++$num_headers;
      } elsif($line =~ m/^#/) {
         # $line is a comment
         push(@new_gff3,$line) unless $strip_comments;
         ++$num_comments;
      } elsif($strip_seqid && exists $strip_seqid->{ (split(/\t/,$line,2))[0] }) {
         # line has a seqid indicating it should be skipped
         ++$num_features;
      } else {
         # $line is a feature
         my @fields = split(/\t/,$line);
         GFF3_LAST_FIELD == $#fields || die "$gff3_file line $. should be a GFF3 feature but it doesn't have 9 fields:\n$line";
         
         
         if($seqid_ncbi_to_submitter) {
            if( exists $ncbi_to_submitter_seqid{$fields[GFF3_SEQID_FIELD]} ) {
               $fields[GFF3_SEQID_FIELD] = $ncbi_to_submitter_seqid{$fields[GFF3_SEQID_FIELD]};
            }
         }
         if($set_source) {
            $fields[GFF3_SOURCE_FIELD] = $set_source;
         }
         if($change_type) {
            foreach my $change ( ARRAY_REF_TYPE eq ref($change_type) ? @{$change_type} : $change_type ) {
               $fields[GFF3_TYPE_FIELD] = $change->{-to} if $fields[GFF3_TYPE_FIELD] eq $change->{-from};
            }
         }
         if($name_gene_from_id) {
            if( !exists $name_gene_from_id->{-types} || grep($_ eq $fields[GFF3_TYPE_FIELD], @{$name_gene_from_id->{-types}}) ) {
               my $attr = _deserialize_attributes($fields[GFF3_ATTR_FIELD]);
               my $sep = $name_gene_from_id->{-separator} // '.';
               $attr->{Name} = exists $name_gene_from_id->{-prefix} ? $name_gene_from_id->{-prefix}.$sep : '';
               $attr->{Name} .= $attr->{ID};
               $attr->{Name} .= $sep.$name_gene_from_id->{-suffix} if exists $name_gene_from_id->{-suffix};
               $fields[GFF3_ATTR_FIELD] = _serialize_attributes($attr);
            }
         }
         if($handler) {
            # pass first 8 fields to handler as strings, and the attributes (field 9) as a hash
            my @new_fields = $handler->(@fields[0..(GFF3_ATTR_FIELD-1)],_deserialize_attributes($fields[GFF3_ATTR_FIELD]));
            GFF3_LAST_FIELD == $#new_fields || die "handler function ".(Sub::Identify::sub_name($handler))." did not return 9 fields";
            # serialize attributes
            $new_fields[GFF3_ATTR_FIELD] = _serialize_attributes($new_fields[GFF3_ATTR_FIELD]);
            @fields = @new_fields;
         }
         
         # finished modifying
         push(@new_gff3,join("\t",@fields));
         
         # extra feature(s) to add?
         if('CDS' eq $fields[GFF3_TYPE_FIELD] && $copy_CDS_to_exon) {
            try {
               push(@new_gff3,join("\t",_change_feature_type(@fields,'CDS','exon')));
            } catch {
               die "error at $gff3_file line $. whilst creating exon feature from CDS feature: $_";
            };
         }
         
         ++$num_features;
      }
      
      ++$num_lines;
      print STDERR "\b" x length($status);
      $status = "$num_headers headers $num_comments comments $num_features features $num_lines_fasta lines of FASTA";
      print STDERR $status;
   }
   print STDERR "\n";
   close(GFF3) || die "error whilst reading $gff3_file: $!";

   return(\@new_gff3);
}

# pass list of fields from a feature line, the old type, the new type, and optionally a new parent ID
# changes the type to the new type; changes the ID accordingly; and changes the parent ID (of provided)
sub _change_feature_type {
   my $this = (caller(0))[3];
   my @fields;
   for(my $i=0; $i<=GFF3_LAST_FIELD; ++$i) {
      push(@fields, shift()) || die "$this was passed $i GFF3 fields (expected ".(GFF3_LAST_FIELD+1).")";
   }
   my $old_type      = shift() || die "$this must be passed a value for the \"type\" field";
   my $new_type      = shift() || die "$this must be passed a value for the \"type\" field";
   my $new_parent_id = shift();

   # verify and then change feature type
   $old_type eq $fields[GFF3_TYPE_FIELD] || die "$this was asked to modify a $old_type feature but was passed a ".$fields[GFF3_TYPE_FIELD]." feature";
   $fields[GFF3_TYPE_FIELD] = $new_type;
   
   my $attr = _deserialize_attributes($fields[GFF3_ATTR_FIELD]);

   # change parent ID, if required
   if($new_parent_id) {
      exists $attr->{Parent} || die "$this was asked to change the Parent attribute of a feature which has no parent";
      $attr->{Parent} = $new_parent_id;
   }
   
   # change ID attribute; requires some  guesswork
   # first look for type within the ID, and change that
   unless( $attr->{ID} =~ s/$old_type/$new_type/i ) {
      # if that didn't work, take parent ID (if there is one) and add suffix based on new type
      if(exists $attr->{Parent}) {
         $attr->{ID} = join('.',$attr->{Parent},lc($new_type));
      } 
      # failing all else, just suffix existing ID with type
      else {
          $attr->{ID} = join('.',$attr->{ID},lc($new_type));
      }
   }
   
   $fields[GFF3_ATTR_FIELD] = _serialize_attributes($attr);
   
   return(@fields);
}

# pass attributes string (i.e. 9th field of a feature line)
# returns hash of attributes
sub _deserialize_attributes {
   my $this = (caller(0))[3];
   my $attr_string = shift() || "$this must be passed attributes string";
   
   my $attr = {};
   foreach my $attr_kev (split(/;/, $attr_string)) {
      # attributes are KEV (key/encoded value) pairs
      my($k,$ev) = split(/=/,$attr_kev,2);
      exists $attr->{$k} && die "multiple instances of attribute $k";
      # encoded value is actually comma-separated list (commas inside of values must be escaped)
      # => split on commas, and then decode
      my @v = map(URI::Escape::uri_unescape($_), split(/,/,$ev));
      if( defined $v[1] ) {
         $attr->{$k} = \@v;
      } else {
         $attr->{$k} = $v[0];
      }
   }
   
   return($attr);
}

# pass attributes as a hash
# returns attributes as a string (usable for 9th field of a feature line)
sub _serialize_attributes {
   my $this = (caller(0))[3];
   my $hashref = shift();
   die "$this must be passed attributes hash" unless $hashref && HASH_REF_TYPE eq ref($hashref);
   # want to be able to modify a local copy
   my $attr = Storable::dclone($hashref);
   
   my @attributes = ();
   # Name, ID and Parent are added to the string in order
   foreach my $k (grep(exists $attr->{$_}, qw(Name ID Parent))) {
      (ARRAY_REF_TYPE eq ref($attr->{$k})) && die "multiple values have been created for non-repeating attribute $k";
      push(@attributes, join('=',$k,URI::Escape::uri_escape($attr->{$k})));
      delete $attr->{$k};
   }
   foreach my $k (sort keys %{$attr}) {
      # if the value is a list, it should be rendered as a sting of encoded values
      # separated by *unencoded* commas
      my $ev   = (ARRAY_REF_TYPE eq ref($attr->{$k}))
               ? join(',', map(URI::Escape::uri_escape($_), @{$attr->{$k}}))
               : URI::Escape::uri_escape($attr->{$k});
      push(@attributes, join('=',$k,$ev));
   }
   
   return( join(';',@attributes) );
}


1;

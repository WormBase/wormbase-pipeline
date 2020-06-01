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
use constant CODE_REF_TYPE       => ref(sub{});

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

=item -gff3_input [mandatory]

Name of the GFF3 file I<or> GFF3 as a list.   List should be one GFF3 line
in each element, no newline char(s) -- i.e. the same format as the modified GFF3
returned by this function.

=item -file I<Deprecated>

Synonym for C<gff3_input>

=item -strip_comments

Flag. If true, all comment lines are removed from the GFF3. These lines are
deleted prior to any further processing.

=item -strip_headers

Flag. If true, all header lines (I<except> the C<##gff-version>) are removed
from the GFF3. These lines are deleted prior to any further processing.

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

I<Deprecated>:  use C<-copy_feature_to_new_type(-from=>'CDS' -to=>'exon')
instead

Flag.  If true, CDS features are duplicated as additional exon features with
the ID tweaked accordingly.

=item -copy_feature_to_new_type

Hash of params.  If defined, creates a copies feature to new feature
with a different type

Params:

=over

=item -from

Current feature type (mandatory)

=item -to

New feature type (mandatory)

=item -sources

List of sources.  Only features from one of these sources will be copied
(optional)

=item -matches

List of attribute values.  Only features with the relevant attribute (C<ID>
by default; can be set with C<-match_attr>) value in this list are copied
(optional).

=item -match_attr

The attribute that is matched against the values in C<-matches>; defaults to
C<ID>.

=back

=item -handler

Reference to a function. This is called for every feature line. 
The handler is passed a list of 9 fields; the first 8 fields are strings, the
9th deserializes attributes as a hash.  Attribute values are unescaped (%09
converted to a space, etc.), and multiple values are represented as a list.

The handler must return one of three things:

1- undef to indicate that the feature should
be deleted

2 - list of 9 fields in the same form, which will be
turned into a string and substituted in place of the original line.  Note
attribute values should I<not> be escaped, and multiple values should be
in the form on a list (i.e. I<not> as a single string of comma-separated
values).

3 - a list reference; each item in this list should be a list of 9 fields
as describe in option 2.   This allows the handler to create new features to
be added to the GFF.

=item -reader

Reference to a function.  This is like a read-only version of C<-handler>;
features are passed just the same waym, but nothing is returned, so no changes
can be made to the GFF.  Use this if you want to read the data for a "look ahead".

The reader function is called before any changes are made to features, but if
you want to look at all the features prior to making chaanges, you'll be
making a separate call to munge_line_by_line() anyway.

=back

=head3 Order of the operations

The arguments cause GFF3 features to be modified are applied in the following order:

   -seqid_ncbi_to_submitter
   -set_source
   -change_type
   -name_gene_from_id
   -handler
   -copy_CDS_to_exon (deprecated)
   -copy_feature_to_new_type

Currently there is no way to change this order.  The workaround is to call
munge_line_by_line() repeatedly, though you need to write the GFF3 file after
each call.   At the moment this seems acceptable as the requirement for
changing the order is likely to be an edge case.
   
=head3 Return values.

A reference to a list of new GFF3 lines.

=cut

sub munge_line_by_line {
   my $this = (caller(0))[3];
   my %params = HASH_REF_TYPE eq ref(@_[0]) ? %{@_[0]} : @_;
   
   my $gff3_input          = $params{-gff3_input}
                             || $params{-file} # DEPRECATED
                             || die "$this must be passed GFF3 (-gff3_input)";
   my $strip_headers       = $params{-strip_headers};
   my $strip_comments      = $params{-strip_comments};
   my $strip_fasta         = $params{-strip_fasta};
   my $strip_seqid         = $params{-strip_seqid};
   my $seqid_ncbi_to_submitter = $params{-seqid_ncbi_to_submitter};
   my $set_source          = $params{-set_source};
   my $change_type         = $params{-change_type};
   my $name_gene_from_id   = $params{-name_gene_from_id};
   my $copy_CDS_to_exon    = $params{-copy_CDS_to_exon}; # DEPRECATED
   my $copy_feature_to_new_type = $params{-copy_feature_to_new_type};
   my $reader              = $params{-reader};
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
      
   if( $copy_feature_to_new_type ) {
      die "$this parameter -copy_feature_to_new_type must be a hash reference"
         if HASH_REF_TYPE ne ref($copy_feature_to_new_type);
      die "$this parameter -copy_feature_to_new_type must pass a hash with keys '-from' and '-to'"
         if !exists $copy_feature_to_new_type->{-from} || !exists $copy_feature_to_new_type->{-to};
      die "$this parameter -matches must ne an ARRAY reference"
         if exists $copy_feature_to_new_type->{-matches} && ARRAY_REF_TYPE ne ref($copy_feature_to_new_type->{-matches});
      die "$this parameter -sources must ne an ARRAY reference"
         if exists $copy_feature_to_new_type->{-sources} && ARRAY_REF_TYPE ne ref($copy_feature_to_new_type->{-sources});
   }
   
   die "$this parameter -reader must be a code reference"
      if $reader && CODE_REF_TYPE ne ref($reader);

   die "$this parameter -handler must be a code reference"
      if $handler && CODE_REF_TYPE ne ref($handler);
   
   my $status = '';
   my $num_lines = 0;
   my $num_headers = 0;
   my $num_comments = 0;
   my $num_features = 0;
   my $num_lines_fasta = 0;
   my $reading_fasta = 0;
   my @new_gff3 = ();
   
   my @input; # array of GFF lines (no newline char(s))
   if( ref($gff3_input) && ARRAY_REF_TYPE eq ref($gff3_input) ) {
      @input = @{ $gff3_input };
   } else {
      @input = File::Slurp::read_file($gff3_input, chomp=>1);
   }
   
   foreach my $line (@input) {
      
      if($line =~ m/^##\s*FASTA\s*$/i) {
         # $line is a FASTA directive
         my @fasta_lines = map {chomp; $_} <GFF3>;
         push(@new_gff3,$line,@fasta_lines) unless $strip_fasta;
         $num_lines_fasta += 1 + scalar(@fasta_lines);
      } elsif($line =~ m/^##[^#]/) {
         # $line is a header
         push(@new_gff3,$line) unless $strip_headers && $line !~ m/^##gff\-version/;
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
         GFF3_LAST_FIELD == $#fields || die "input GFF3 line $. should be a GFF3 feature but it doesn't have 9 fields:\n$line";
         
         if($reader) {
            # pass first 8 fields to handler as strings, and the attributes (field 9) as a hash
            $reader->(@fields[0..(GFF3_ATTR_FIELD-1)],_deserialize_attributes($fields[GFF3_ATTR_FIELD]));
         }

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
            my @r = $handler->(@fields[0..(GFF3_ATTR_FIELD-1)],_deserialize_attributes($fields[GFF3_ATTR_FIELD]));
            if( defined $r[0] ){
               my $new_features; # list of lists
               if( ARRAY_REF_TYPE eq ref($r[0]) ) {
                  # if $r[0] references a list, everything in that list should be a list ref
                  foreach my $thing ( @{$r[0]} ) {
                     ARRAY_REF_TYPE eq ref($thing) || die "handler function ".(Sub::Identify::sub_name($handler))." should have returned a list of lists, but it returned a list containing $thing";
                  }
                  $new_features = $r[0];
               } 
               # handler retuned a simple list of 9 fields; convert to list-of-lists structure
               else {
                  $new_features = [\@r];
               }
               foreach my $new_feature ( @{$new_features} ) {
                  my @new_fields = @{$new_feature};
                  # every feature should be a list of 9 fields
                  GFF3_LAST_FIELD == $#new_fields || die "handler function ".(Sub::Identify::sub_name($handler))." did not return 9 fields @new_fields";
                  # serialize attributes
                  $new_fields[GFF3_ATTR_FIELD] = _serialize_attributes($new_fields[GFF3_ATTR_FIELD]);
                  push(@new_gff3,join("\t",@new_fields));
               }
               
            }
         }
         # no handler is defined => push fields onto new GFF
         else {
            # undef indicates that the feature should be deleted
            if( defined $fields[0] ) {
               push(@new_gff3,join("\t",@fields));
            }
         }
         
         # extra feature(s) to add?

         # DEPRECATED
         if('CDS' eq $fields[GFF3_TYPE_FIELD] && $copy_CDS_to_exon) {
            try {
               push(@new_gff3,join("\t",_change_feature_type(@fields,'CDS','exon')));
            } catch {
               die "error at input GFF3 line $. whilst creating exon feature from CDS feature: $_";
            };
         }
         if( $copy_feature_to_new_type && $copy_feature_to_new_type->{-from} eq $fields[GFF3_TYPE_FIELD] ) {
            my $copy_this_feature = 1;
            
            # if a source match is required, is this one of the matching sources?
            if( exists $copy_feature_to_new_type->{-sources} ) {
               # now conditional on source match
               $copy_this_feature = 0;
               if( grep( $_ eq @fields[GFF3_SOURCE_FIELD], @{$copy_feature_to_new_type->{-sources}} ) ) {
                  $copy_this_feature = 1;
               }
            }
            
            # if OK so far, and if an attribute match is required, is does this have a matching attribute??
            if( $copy_this_feature && exists $copy_feature_to_new_type->{-matches} ) {
               # now also conditional on attribute match
               $copy_this_feature = 0;
               my $attr_to_match = $copy_feature_to_new_type->{-match_attr} || 'ID';
               my $id = _deserialize_attributes($fields[GFF3_ATTR_FIELD])->{$attr_to_match} || die "feature has no ID attribute";
               if( grep( $_ eq $id, @{$copy_feature_to_new_type->{-matches}} ) ) {
                  $copy_this_feature = 1;
               }
            }
            
            # all OK: this feature is to be copied
            if( $copy_this_feature ) {
               try {
                  push( @new_gff3,
                        join( "\t",
                              _change_feature_type(@fields,
                                                   $copy_feature_to_new_type->{-from},
                                                   $copy_feature_to_new_type->{-to}
                                                   )
                              )
                        );
               } catch {
                  die "error at input GFF3 line $. whilst creating exon feature from CDS feature: $_";
               };
            }
         }
         
         ++$num_features;
      }
      
      ++$num_lines;
      print STDERR "\b" x length($status);
      $status = "$num_headers headers;  $num_comments comments;  $num_features features;  $num_lines_fasta lines of FASTA";
      print STDERR $status;
   }
   print STDERR "\n";

   return(\@new_gff3);
}


=head2 change_transcripts_to_noncoding

Finds designated transcripts in a GFF3 file, changes them to C<pseudogenic_transcript>,
and removes their child features.

What we're claiming biologically when we do this:

=over

=item * 
these transcripts belong to genes that have no function

=item *
they look like protein-coding genes but are not (because they have stops)

=item *
they are probably pseudogenes arising from retrotansposition, i.e. copying an re-integration
back into somewhere else in the genome, followed by degradation

=back

The (presumably modified) GFF3 lines are returned as a list.  Note that lines returned
have been chomp'ed.

=head3 Arguments

This function requies named parameters.

=over

=item -gff3_input [mandatory]

Name of the GFF3 file I<or> GFF3 as a list.   List should be one GFF3 line
in each element, no newline char(s) -- i.e. the same format as the modified GFF3
returned by this function.

=item -file I<Deprecated>

Synonym for C<gff3_input>

=item -transcripts

Transcript IDs that are to be changed.  I<Either> the name of a file containing
a list (one per line) I<or> an array of transcript IDs

=item -types

Optional list of types of transcripts.  Defaults to 'mRNA'.

=item -new_type

New type for the transcripts.  Defaults to 'pseudogenic_transcript'.

=head3 Note...

This was shamelessly ripped off from C<change_transcripts_to_noncoding.pl>.  Like that
script, this function will be used after a healthcheck has reported stop codons within
transcripts.  The (arguably flawed?) I<rationale> for adding it to this package
is to try and get the entire GFF3 munge process into the C<to_our_gff.pl> scripts
so it is repeatable and documented.

Within this package, is a separate function rather than a c<munge_line_by_line> option
because it removes features; it's cleaner to do this as a separate pass.  It's slower but
that's not significant.

=cut

sub change_transcripts_to_noncoding {
   my $this = (caller(0))[3];
   my %params = HASH_REF_TYPE eq ref(@_[0]) ? %{@_[0]} : @_;
   
   my $gff3_input          = $params{-gff3_input}
                             || $params{-file}  # DEPRECATED
                             || die "$this must be passed GFF3 (-gff3_input)";
   my $transcript_IDs      = $params{-transcripts} || die "$this must be passed a list of transcript IDs (file name or list) (-transcripts)";
   my $transcript_types    = $params{-types}       || ['mRNA'];
   my $new_transcript_type = $params{-new_type}    || 'pseudogenic_transcript';
   
   ARRAY_REF_TYPE eq ref($transcript_types) || die "When -types is passed to $this it must be an ARRAY ref ";

   # grab the list of "target" transcript IDs from the list/file provided
   my %transcripts_to_noncoding = ();
   {
      my $num_id = 0;
      my @tmp = ();
      if( ARRAY_REF_TYPE eq ref($transcript_IDs) ) {
         @tmp = @{$transcript_IDs};
      } else {
         @tmp = File::Slurp::read_file($transcript_IDs, chomp=>1);
      }
      foreach my $id (@tmp) {
         ++$num_id unless exists $transcripts_to_noncoding{$id};
         $transcripts_to_noncoding{$id}++;
      }
      die "$transcript_IDs contained no transcript IDs" unless $num_id;
   }
   my @new_gff3      = ();
   my $num_lines     = 0;
   my $num_noncoding = 0;
   my $num_removed   = 0;
   my $status        = '';

   my @input; # array of GFF lines (no newline char(s))
   if( ref($gff3_input) && ARRAY_REF_TYPE eq ref($gff3_input) ) {
      @input = @{ $gff3_input };
   } else {
      @input = File::Slurp::read_file($gff3_input, chomp=>1);
   }
   
   foreach my $line (@input) {
      ++$num_lines;
      if($line =~ m/^##\s*FASTA\s*$/i) {
         # $line is a FASTA directive
         my @fasta_lines = map {chomp; $_} <GFF3>;
         push(@new_gff3,$line,@fasta_lines);
      } elsif($line =~ m/^#/) {
         # $line is a header or comment
         push(@new_gff3,$line);
      } else {
         my @fields = split(/\t/,$line);
         GFF3_LAST_FIELD == $#fields || die "input GFF3 line $. should be a GFF3 feature but it doesn't have 9 fields:\n$line";
         my $attr = _deserialize_attributes($fields[GFF3_ATTR_FIELD]);
         # targetted transcripts changed into pseudogenic_transcript
         if( grep($_ eq $fields[GFF3_TYPE_FIELD], @{$transcript_types}) && exists $transcripts_to_noncoding{$attr->{ID}} ) {
            ++$num_noncoding;
            $fields[GFF3_TYPE_FIELD] = $new_transcript_type;
            push(@new_gff3,join("\t",@fields));
         }
         # children of targetted transcripts removed
         elsif( exists $transcripts_to_noncoding{$attr->{Parent}} ) {
            ++$num_removed;
         }
         # everything else unchanged
         else {
            push(@new_gff3,join("\t",@fields));
         }
         
         print STDERR "\b" x length($status);
         $status = "$num_lines lines;  $num_noncoding ".join('/',@{$transcript_types})." features changed to $new_transcript_type;  $num_removed child features removed";
         print STDERR $status;
      }
   }
   print STDERR "\n";

   return(\@new_gff3);
}

# pass list of fields from a feature line, the old type, the new type, and optionally a new parent ID
# changes the type to the new type; changes the ID accordingly; and changes the parent ID (if provided)
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
   unless( $attr->{ID} =~ s/\b$old_type\b/$new_type/i ) {
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

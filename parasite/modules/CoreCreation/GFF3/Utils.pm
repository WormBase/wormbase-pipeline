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

=item Sub::Identify

=item Try::Tiny

=item URI::Escape

=back

=head1 Functions
   
=cut

use strict;

use Sub::Identify;
use Try::Tiny;
use URI::Escape;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(validate substitute);

use constant GFF3_SOURCE_FIELD   => 1;
use constant GFF3_TYPE_FIELD     => 2;
use constant GFF3_ATTR_FIELD     => 8;
use constant GFF3_LAST_FIELD     => 8;

use constant GFF3_VALIDATION_CMD => 'gt gff3validator';

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

Flag. If true, all comments lines are removed from the GFF3. These lines are
deleted prior to any further processing.

=item -strip_fasta

Flag. If true, if a FASTA directive (##FASTA) is found then that line and
all remaining lines are deleted.

=item -set_source

String.  The source (field 2) of each feature is set to this value.

=item -name_gene_from_id

Hash of params.  If defined, creates a Name attribute from the ID attribute.
Params:

=over

=item -types

List of types for which new Names are to be created.  Defaults to all types.

=item -prefix, -suffix

Prefix and/or suffix that will be added to the ID.
If neither of theseare passed, the Name will be a copy of the ID.

=item -separator

Separator used after prefix/before suffix; defaults to '.'

=back

=item -copy_CDS_to_exon

Flag.  If true, CDS features are duplicated as additional exon features with
the ID tweaked accordingly.

=item -handler

Reference to a handler function. This is called for every feature line. 
The handler is passed a list of 9 fields; it should return a list of 9
fields, which will be turned into a string andsubstituted in place of
the original line.

=back

Name of GFF3 file, code reference.

=head3 Return values.

A reference to a list of new GFF3 lines.

=cut
sub munge_line_by_line {
   my $this = (caller(0))[3];
   my %params = ref({}) eq ref(@_[0]) ? %{@_[0]} : @_;
   
   my $gff3_file           = $params{-file} || die "$this must be passed a GFF3 file name (-file)";
   my $strip_comments      = $params{-strip_comments};
   my $strip_fasta         = $params{-strip_fasta};
   my $name_gene_from_id   = $params{-name_gene_from_id};
   my $set_source          = $params{-set_source};
   my $copy_CDS_to_exon    = $params{-copy_CDS_to_exon};
   my $handler             = $params{-handler};
   
   die "$this parameter -name_gene_from_id must be a hash reference"
      if $name_gene_from_id && ref({}) ne ref($name_gene_from_id);

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
      } else {
         # $line is a feature
         my @fields = split(/\t/,$line);
         GFF3_LAST_FIELD == $#fields || die "$gff3_file line $. should be a GFF3 feature but it doesn't have 9 fields:\n$line";
         
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
         if($set_source) {
            $fields[GFF3_SOURCE_FIELD] = $set_source;
         }
         if($handler) {
            # convert attributes to hash when passing to handler
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
      my $v = URI::Escape::uri_unescape($ev);
      if( exists $attr->{$k} ) {
         # repeating attribute: add as a list
         unless( ref([]) eq ref($attr->{$k}) ) {
            # currently a scalar, so make into a list containing existing value
            $attr->{$k} = [$attr->{$k}];
         }
         # push new value onto the list
         push( @{$attr->{$k}}, $v );
      } else {
         $attr->{$k} = $v;
      }
   }
   
   return($attr);
}

# pass attributes as a hash
# returns attributes as a string (usable for 9th field of a feature line)
sub _serialize_attributes {
   my $this = (caller(0))[3];
   my $attr = shift();
   die "$this must be passed attributes hash" unless $attr && ref({}) eq ref($attr);
   
   my $attr_string = '';
   foreach my $k (sort keys %{$attr}) {
      if(ref([]) eq ref($attr->{$k})) {
         # attributes that are lists are rensdered as repeating KEV pairs
         foreach my $v (@{$attr->{$k}}) {
            $attr_string .= ';' if $attr_string;
            $attr_string .= join('=',$k,URI::Escape::uri_escape($v));
         }
      } else {
         $attr_string .= ';' if $attr_string;
         $attr_string .= join('=',$k,URI::Escape::uri_escape($attr->{$k}));
      }
   }
   
   return($attr_string);
}


1;

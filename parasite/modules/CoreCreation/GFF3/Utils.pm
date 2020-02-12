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

   substitue( 'foo.gff3', 'source', 'GenBank', 'WormBase_imported' );

=head1 Dependencies

=over

=item Genometools (gt) excuteable in your $PATH

=back

=head1 Functions
   
=cut

use strict;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(validate substitute);

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


=head2 line_by_line

Iterates through GFF3 file, line by line, optionally making changes.

A code reference must be passed.  The code is called for every
B<non-directive> line, and the line is passed; the code should return
the (possibly modified) line.  The code can return C<undef> or an
empty string to indicate that the line should be deleted.

The (possibly modified) GFF3 lines are added to a list which is returned.
Note that lines have been chomp'ed.

=head3 Arguments

Name of GFF3 file, code reference.

=head3 Return values.

A flag indicating success (true) or failure (false) of the validation,
followed by a (hopefully empty) list of validation error messages.

=cut
sub line_by_line {
   my $this = (caller(0))[3];
   my $gff3_file = shift() || die "$this must be passed a GFF3 file name";
   my $code = shift();
   die "$this must be passed a function reference" unless $code && 'CODE' eq ref($code);
   my @new_gff3 = ();
   open(GFF3, $gff3_file) || die die "can't read $gff3_file: $!";
   foreach my $line (<GFF3>) {
      chomp($line);
      if($line =~ m/^##[^#]/) {
         push(@new_gff3,$line);
      } else {
         my $new = $code->($line);
         $new && push(@new_gff3,$new);
      }
   }
   close(GFF3) || die "error whilst reading $gff3_file: $!";

   return(\@new_gff3);
}

1;

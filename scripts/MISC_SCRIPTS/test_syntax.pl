#!/software/bin/perl

my $directory = shift;
die "no dir to test\n" unless $directory;

opendir(DIR,$directory);
my @files = readdir( DIR );
my @failures;

FILE:
foreach $file ( @files ) {
  next if( $file =~ /\#|~/ or -d $file );
  print "$directory/$file\n";
  open (TEST, "head -1 $directory/$file |" ) or die "$!\n";
  while (<TEST>) {
    if( /(perl\S+)/ ) {
      my $return = system ("$1 -c $directory/$file");
      push(@failures,$file) unless $return == 0;
      close TEST;
      next FILE;
    }
  }
}


print "\n\n! ! FAILURES  ! !The following failed perl-c \n", join("\n",@failures),"\n";


__END__

=pod


=head2 NAME - test_syntax.pl

=head1 USAGE
  
=over 4
 
=item test_syntax.pl /scripts
  
=back
  
This script:
  
I<test_syntax.pl MANDATORY arguments:>

B<NONE>

I<test_syntax.pl  Overview:>

This script will run perl -c on every script in the directory given.  The version of perl used will be that specified on the #! line of the script.  Default directory is /wormsrv2/scripts.  A list of those scripts that fail is reported at the end.


I<test_syntax.pl  OPTIONAL arguments:>
directory to test

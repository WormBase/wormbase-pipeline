#!/usr/bin/perl5.6.1 -w

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

use Wormbase;
use File::Basename;
my @files2split = @ARGV;

my %clone2centre = &FetchData('clone2centre') ;

die "no data in clone2centre hash\n" unless %clone2centre;

my $out;
foreach my $file ( @files2split ) {
  ($base, $filename, $ext) = &fileparse($file);
  open ( $blast_file, "<$file" ) or die "cant open input file $file\t$!\n";

  open ($stl_out,">STL_blastx.ace") or die "cant open STL output\t$!\n";
  open ($cam_out,">CAM_blastx.ace") or die "cant open CAM output\t$!\n";

  #Sequence : "cTel33B"
  $out = $cam_out;
  while (<$blast_file>) {
    if( /Sequence \: \"(.*)\"/) {
      my $clone = $1; 
      if ($clone) {
	if( $clone2centre{$clone} eq "HX") {
	  $out = $cam_out;
	}
	elsif( $clone2centre{$clone} eq "RW") {
	  $out = $stl_out;
	}
	else {
	  print "NO CENTRE FOR $clone\n";
	  $out = STDOUT;
	}
      }
    }
    print $out $_;
  }
}

close $out;
close $cam_out,
close $stl_out;

exit(0);


__END__

=pod

=item This will split the blastx data in to CAM and STL based on clone name.  Clone to centre assignments come from the COMMON_DATA file.

=cut

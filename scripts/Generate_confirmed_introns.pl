#!/usr/local/bin/perl5.6.1 -w
                                 
use lib "/wormsrv2/scripts/";  
use Wormbase;

my $ver = &get_wormbase_version;
my $done;
while(! defined($done) ){
  print "which chromosome (eg II, V ) . .(q to exit) ";
  my $csome = <STDIN>;
  chomp $csome;
  exit(0) if "$csome" eq "q";
  
my $file = "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/GFF/CHROMOSOME_$csome.check_intron_cam.gff";

$done = open (GFF, "<$file") or warn "can do that . . either the file is missing or you have entered an invalid chromosome (try I II III IV V or X )\nI'm looking for $file\n\n\n\n\n";
}
my %uniquify;
my $count;
while (<GFF>)
  {
    my @data = split;
    $uniquify{"$data[8]"}++;
    $count++;
  }

my $clone_count;
foreach (sort { $uniquify{$b} <=> $uniquify{$a} } keys %uniquify)
  {
    if( /Clone:(\w+)/ ) {
      print "$1 has $uniquify{$_} introns\n";
      $clone_count++;
    }
  }
print "get busy - you have $count introns to accout for in $clone_count clones\n";
close GFF;
exit(0);



__END__

=pod

=head2 NAME - Generate_confirmed_introns.pl

=head1 USAGE

=over 4

=item Generate_confirmed_introns.pl  [-options]

=back

This script:

Generates a list of clones that have confirmed introns associated with them
The list is ordered by the number of introns / clone.
Hopefully this will mean we can get through larger chunks at time! 

The output goes to STOUT so redirect it in to a file of your choice.

script_template.pl MANDATORY arguments:

=over 4

=item none ( prompts for chromosome )

=back


=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut


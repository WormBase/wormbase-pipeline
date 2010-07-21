#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  overload_rnai.pl
#
#        USAGE:  ./overload_rnai.pl 
#
#  DESCRIPTION:  adds additional information to RNAi GFF lines
#
#       AUTHOR:  $Author: mh6 $
#      VERSION:  $Revision: 1.2 $
#      CREATED:  06/07/10 10:40:04 BST
#===============================================================================

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;

my ($debug,$test,$species,$store,$file,$dontOverwrite);
GetOptions(
   'debug=s'   => \$debug,
   'test'      => \$test,
   'species:s' => \$species,
   'store:s'   => \$store,
   'file:s'    => \$file,
   'dontoverwrite' => \$dontOverwrite,
)||die(@!);


my $wormbase;
if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# fill some reference hashes
my ($r2lab,$r2hist) = &get_rnai2lab();

my @chromosomes = $wormbase->get_chromosome_names(-prefix => 1, -mito => 1);
@chromosomes = qw(contigs) if ($wormbase->assembly_type eq 'contig' || $file);

foreach my $chromosome(@chromosomes){
  my $fh = $file ? (new IO::File $file) : $wormbase->open_GFF_file($chromosome,undef,$log);
  my $outfileName = ($file || $wormbase->GFF_file_name($chromosome,undef));
  my $outf = new IO::File $outfileName.'_rnai','w';

  while (<$fh>){
	unless(/RNAi_(primary|secondary)\s+RNAi_reagent/){
		print $_;
		next;
	}
	chomp;
	print $outf $_;
	my ($rnaid) = /(WBRNAi\d+)/;
	print $outf " ; Laboratory \"$$r2lab{$rnaid}\"" if $$r2lab{$rnaid};
	print $outf " ; History_name \"$$r2hist{$rnaid}\"" if $$r2hist{$rnaid};
	print $outf "\n";
 }
 $wormbase->run_command("mv -f ${outfileName}_rnai $outfileName", $log) unless $dontOverwrite;
}

# Close log files and exit
$log->mail();
exit(0);


# populates two hashes with RNAiID -> Laboratory and RNAiID -> history name
sub get_rnai2lab {
	use Ace;

	my %rnai2lab;
	my %rnai2history;
	my $db = Ace->connect(-path => $wormbase->autoace);
	my $cursor = $db->fetch_many(RNAi => '*');
	while (my $rnai = $cursor->next){
		$rnai2lab{"$rnai"}="${\$rnai->Laboratory}" if $rnai->Laboratory;
		$rnai2history{"$rnai"}="${\$rnai->History_name}" if $rnai->History_name;
	}
	$db->close;
        return \%rnai2lab,\%rnai2history;
}

=pod

=head1 NAME - overload_rnai.pl

=head1 USAGE

=over 4 

=item overload_rnai.pl [-debug USER -test -species WORM -store STORE -file GFF_FILE -dontoverwrite]

=back

this script adds Laboratory and History names to RNAi GFF lines

=head2 Optional arguments

=over 4

=item -debug <username>

send the log mails only to <username>

=item -test

use the test databases

=item -species <species name>

set the species used

=item -store <storable name>

use a storable as wormbase object

=item -file <gff file name>

use a specific GFF file instead of the CHROMOSOME GFFs

=item -dontoverwrite

leave the temp file in place and don't overwrite the original

=back 

=cut

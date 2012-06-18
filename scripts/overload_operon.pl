#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  overload_operon.pl
#
#        USAGE:  ./overload_operon.pl
#
#  DESCRIPTION:  adds additional information to Operon GFF lines
#
#       AUTHOR:  $Author: klh $
#      VERSION:  $Revision: 1.2 $
#      CREATED:  21/05/12 10:40:04 BST
#===============================================================================

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;

my ($debug,$test,$species,$store,$file,$gff_dir);

GetOptions(
	   'debug=s'       => \$debug,
	   'test'          => \$test,
	   'species:s'     => \$species,
	   'store:s'       => \$store,
	   'file:s'        => \$file,
	   'gffdir:s'      => \$gff_dir,
	  );


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
my ($operon_genes) = &get_operon_data();

my @gff_files;
$gff_dir = $wormbase->chromosomes if not defined $gff_dir;

if (defined($file)){
  push(@gff_files,$file);
}
else {
  if($wormbase->assembly_type eq 'contig') {
    push(@gff_files,lc($wormbase->species));
  } else {
    @gff_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
  }
}

foreach my $fileprefix (@gff_files) {
  if ($debug) {print "$fileprefix\n";}
  my $in_file = "$gff_dir/${fileprefix}.gff";
  my $out_file = "$gff_dir/${fileprefix}.operon.gff";

  if (not -e $in_file) {
    $log->log_and_die("GFF file $in_file not found\n");
  }
  if (-z $in_file) { 
    $log->log_and_die("GFF file $in_file zero length\n");
  }

  open(my $in_fh, $in_file) or $log->log_and_die("Could not open $in_file for reading\n");
  open(my $out_fh, ">$out_file") or $log->log_and_die("Could not open $out_file for writing\n");
  
  while (<$in_fh>){
    unless(/operon\s+operon/){
      print $out_fh $_;
      next;
    }
    chomp;
    print $out_fh $_;
    my ($operond) = /(CEOP\d+|CEOP\w+\d+)/;
    if (/Gene/) {
      print "$_";
      $log->log_and_die("\nIt appears that you have already overloaded the CEOP lines for $species\n\n");
    }
    foreach my $genes (@{$$operon_genes{$operond}}) { 
      print $out_fh " ; Gene \"$genes\"" if (defined $$operon_genes{$operond});
    }
    print $out_fh "\n";
  }
  if ($debug) {print "Finished $out_file\n";}
  close($out_fh);
  close($in_fh);
  
  $wormbase->run_command("mv -f $out_file $in_file", $log);
}

# Close log files and exit
$log->mail();
exit(0);


# populates Operon::Gene hash
sub get_operon_data {
	use Ace;

	my %operon_genes;
	my $db = Ace->connect(-path => $wormbase->autoace);
	my $cursor = $db->fetch_many(Operon => '*');
	while (my $operon = $cursor->next){
	  my @tmp_genes = $operon->Contains_gene;
	  my @tmp_genes2 = map{"$_"} @tmp_genes;
	  if ($debug) {print "$operon\n" unless (defined $operon->Contains_gene->name);}
	  push @{$operon_genes{"$operon"}}, @tmp_genes2 if (defined $operon->Contains_gene);
	}
	$db->close;
        return \%operon_genes;
      }

=pod

=head1 NAME - overload_operon.pl

=head1 USAGE

=over 4 

=item overload_operon.pl [-debug USER -test -species WORM -store STORE -file GFF_FILE -dontoverwrite]

=back

this script adds Genes to Operon GFF lines

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

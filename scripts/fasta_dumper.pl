#!/software/bin/perl -w
#
# transcript_builder.pl
# 
# by Anthony Rogers and Gary Williams
#
# Script to make ?Transcript objects
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2013-01-07 17:39:26 $
use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Wormbase;
use File::Path;
use Storable;

my ($debug, $store, $verbose, $database, $test, $wormbase, $species, $class, $method, $out_file, $target_dir);

GetOptions ( "debug:s"      => \$debug,
	     "verbose"      => \$verbose,
	     "class:s"      => \$class,     #specify an acedb class
	     "method:s"     => \$method,    #specify an acedb method
	     "database:s"   => \$database,  #specify a database to query
	     "output:s"     => \$out_file,  #specify an output file name
	     "outdir:s"     => \$target_dir, #specify an output location or else the file will be written to BUILD/species/SEQUENCES/
	     "test"         => \$test,
	     "store:s"      => \$store,
	     "species:s"    => \$species,
	   ) ;



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species,
                           );
}
# Create log file
my $log = Log_files->make_build_log($wormbase);

# Define database to be queried
my $db;
if (not defined $database or $database eq "autoace") {
  $db = $wormbase->autoace;
}else{$db = $database}

my $runtime = $wormbase->runtime;
my $tace = $wormbase->tace;


unless (defined $class) {$log->write_to("!!!!Warning - No class specified so defaulting to Transcript\n");
			 $class = "Transcript";}
unless (defined $method) {$log->write_to("!!!!Warning - No method specified so defaulting to Coding_transcript\n");
			  $method = "Coding_transcript";}

unless (defined $target_dir) {$target_dir = $wormbase->sequences;}
unless (defined $out_file){$out_file     = "coding_transcripts.dna";}

my $out = "$target_dir/$out_file";


############
# Main body
############

$log->write_to("\n$runtime: Started making $class($method) fasta file:\n$out\n");

if (-e "$out") {
$log->write_to("\nOld $out exists in this location, removing.......\n");
$wormbase->run_command("rm -f $out", $log);
}


my $gspecies = $wormbase->full_name('-g_species' => 1);
my $full_name = $wormbase->full_name;


open FAOUT,">$out" ||die($!);

my $object;
my $connection = Ace->connect(-path => "$db") || die (Ace->error);
my $object_it = $connection->fetch_many(-query => "Find $class; Species=\"${full_name}\"; Method=\"$method\";");
while(my $object=$object_it->next){
  if ($verbose){print $object->name."\n";}
  my $dna =$object->asDNA();
  print FAOUT "$dna";
}
  close FAOUT;

if (-e "${out}.gz") {
$log->write_to("\nOld zipped data exists ${out}.gz exists in this location, removing.......\n");
$wormbase->run_command("rm -f ${out}.gz", $log);
}


  $wormbase->run_command("gzip -9 -f $out", $log);
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making $class($method) fasta file\n\nFile location:$out\n\n");

$log->mail();
exit(0);


__END__

=pod

=head2 NAME - fasta_dumper.pl

=head1 USAGE

=over 4

=item fasta_dumper.pl  [-options]

=back

This script "does exactly what it says on the tin". ie it dumps fasta sequences for a given method class combination.

=back

=head2 fasta_dumper.pl arguments:

=over 4

=item * verbose - terminal output

=item * database:s      - either use autoace if used in build process or give the full database path. Basically retrieves paired read info for ESTs from that database.


=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@ebi.ac.uk)

=back

=cut

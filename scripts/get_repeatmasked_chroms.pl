#!/usr/local/bin/perl5.8.0 -w
#
# agp2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# modified for reading in .agp files for worm ensembl


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use lib '/software/worm/lib/bioperl-live';
use lib '/software/worm/ensembl/ensembl/modules';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeq;
use Bio::SeqIO;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $species, $wormbase,$database, $softmask);

my $agp;
my $out_dir;


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "output:s"   => \$out_dir,
	    "database:s" => \$database,
	    "species:s"  => \$species,
            "softmask"   => \$softmask,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}
# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
unless ($database =~ /$species/) {
  $log->log_and_die("are you sure you have the right species / database combo\n Species :$species\nDatabase : $database\n");
}


$out_dir = $wormbase->chromosomes unless $out_dir;
die "cant write to $out_dir\t$!\n" unless (-w $out_dir );

$log->write_to("Connecting to worm_dna\n");

my $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-host'   => 'ia64d',
        '-user'   => 'wormro',
        '-dbname' => $database
        )
    or die "Can't connect to Database $database";

$log->write_to("Building chromosomes\n");

my $file_suffix = $softmask ? "softmasked.dna" : "masked.dna";
my $chr_assembly = ($wormbase->assembly_type eq 'contig') ? 0 : 1;

my ($seqio, $filename);

if (not $chr_assembly) {
  $filename = "$out_dir/${species}_${file_suffix}";
  $seqio = Bio::SeqIO->new(-format => 'fasta',
                           -file => ">$filename");
}


foreach my $seq ( @{$dbobj->get_SliceAdaptor->fetch_all('toplevel')}) {
  my $name = $seq->seq_region_name();

  my ($outfile, $outfh);
  if ($chr_assembly) {
    $filename = "$out_dir/${name}_${file_suffix}";
    $seqio = Bio::SeqIO->new(-format => 'fasta',
                             -file => ">$filename");
  }

  &print_seq($seqio, $seq);

  if ($chr_assembly) {
    $seqio->close();
    $wormbase->run_command("gzip -9 $filename", $log);
  }
}

if (not $chr_assembly){
  $seqio->close();
  $wormbase->run_command("gzip -9 $filename", $log);
}

$log->write_to("Done\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub print_seq {
  my ($seqio,$ens_slice) = @_;
 
  my $seq_name = $ens_slice->seq_region_name;

  $log->write_to("\twriting chromosome $seq_name\n");

  $ens_slice = $softmask
      ? $ens_slice->get_repeatmasked_seq(undef, 1) 
      : $ens_slice->get_repeatmasked_seq;

  my $seqobj = Bio::PrimarySeq->new(-seq => $ens_slice->seq,
                                    -id => $seq_name,
                                    -desc => "1 " . $ens_slice->seq_region_length);

  $seqio->write_seq($seqobj);
}

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head1 NAME

agp2ensembl.pl

=head1 SYNOPSIS

get_repeatmasked_chroms.pl -agp ~wormpipe/Elegans/WSXXX.agp

=head1 DESCRIPTION

extracts RepeatMasked sequence from DNA database using specified AGP file

=head1 OPTIONS

    -agp     agp file

=head1 CONTACT

ar2@sanger.ac.uk

=cut

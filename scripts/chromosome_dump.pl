#!/software/bin/perl -w
#!/usr/local/bin/perl5.8.0 -w
#
# chromosome_dump.pl 
#
# by Keith Bradnam
#
# A script for dumping dna and/or gff files for chromosome objects in autoace
# see pod for more details
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-25 21:53:22 $


use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::PrimarySeq;

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Storable;


######################################################
# Script variables and command-line options          #
######################################################

our ($help, $debug, $dna, $gff, $zipdna, $zipgff, $composition, $database, $dump_dir, $genome_seq_file, $test);
my $store;

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "dna"         => \$dna,
	    "gff"         => \$gff,
	    "zipdna"      => \$zipdna,
	    "zipgff"      => \$zipgff,
	    "composition" => \$composition,
	    "database=s"  => \$database,
	    "dump_dir=s"  => \$dump_dir,
            "genomeseq=s" => \$genome_seq_file,
	    "test"        => \$test,
	    "store:s"     => \$store
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

our $tace   = $wormbase->tace;
our $giface = $wormbase->giface;

##########################################################

# Display help if required
&usage("Help") if ($help);

# Sanity checks
if(!$database && $dump_dir){
  die "You have specified a destination directory to dump to (-dump_dir flag) but not\na source database (-database flag).\n";
}
if(!$dna && !$composition && !$zipdna){
  die "No major option (-dna, -composition, -zipdna) has been specified, try again\n";
}

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wormbase->basedir;
$database = $wormbase->autoace     unless ($database);
$dump_dir = $wormbase->chromosomes unless ($dump_dir);
$genome_seq_file = $wormbase->genome_seq unless ($genome_seq_file);

#####################################################
# Main subroutines
#####################################################

# establish log file.
my $log = Log_files->make_build_log($wormbase);

&dump_dna    if ($dna);
&composition if ($composition );
&zip_files   if ($zipdna);

# say goodnight Barry
$log->mail;

exit(0);




#############################################################################
# Subroutines
#############################################################################


#########################
# dump dna files
#########################

sub dump_dna {

  my $command;

  # command generation
  foreach my $c ($wormbase->get_chromosome_names(-mito => 1,-prefix=> 1)) {
    $command.="find sequence $c\ndna -f $dump_dir/$c.dna\n";
  }
  $command.='quit';

  &execute_ace_command($command,$tace,$database);

  $log->write_to("Removing blank first lines\n");
  foreach ($wormbase->get_chromosome_names(-mito => 1,-prefix=> 1)) {
    $wormbase->remove_blank_lines("$dump_dir/$_.dna", 'no_log');
  }


  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file => ">$genome_seq_file");
  foreach my $chr ($wormbase->get_chromosome_names(-mito => 1,-prefix=> 1)) {
    my $chr_file = "$dump_dir/$chr.dna";
    my $sequence = &read_chromosome($chr_file);
    my $seq = Bio::PrimarySeq->new(-id => $chr,
                                   -seq => uc($sequence));
    $seqio->write_seq($seq);
  }

  $log->write_to("Finished dumping DNA\n\n");
}


###################################
# produce dna composition files
###################################

sub composition {

  $log->write_to("Generating composition.all\n");	
  
  my @chroms = $wormbase->get_chromosome_names(-prefix => 1, -mito => 1);

  my $total;
  my $ns;
  my $gaps;
  my $as;
  my $cs;
  my $gs;
  my $ts;

  $log->write_to("Generating totals file\n");

  foreach my $chr (@chroms) {

    my $file = "$dump_dir/${chr}.dna";

    #print "read $file\n";
    my $seq = &read_chromosome($file);

    my ($a, $c, $g, $t, $gap, $n, $length_dna) = &get_composition($seq);
    $as += $a;
    $cs += $c;
    $gs += $g;
    $ts += $t;
    $total += $length_dna;
    $ns += $n;
    $gaps += $gap;
  }

  open(OUT, ">$dump_dir/composition.all") or $log->log_and_die("Couldn't open composition.all\n");
  print OUT " $total total\n";
  print OUT " a $as\n";
  print OUT " c $cs\n";
  print OUT " g $gs\n";
  print OUT " t $ts\n";
  print OUT " - $gaps\n";
  print OUT " n $ns\n";
  close(OUT);

  my $final_total = $total - $gaps;

  $wormbase->run_command("echo $total $final_total > totals",$log);
  $wormbase->release_composition($log) if ($wormbase->species eq 'elegans');
}

##########################################
# get the composition of a sequence

sub get_composition {
  my ($dna) = @_;

  my ($a, $c, $g, $t, $gap, $n, $length_dna);

  $a = $dna =~ tr/[aA]/A/;
  $c = $dna =~ tr/[cC]/C/;
  $g = $dna =~ tr/[gG]/G/;
  $t = $dna =~ tr/[tT]/T/;
  $gap = $dna =~ tr/\-/-/;

  # the Ns are whatever is not ACGT-
  $length_dna = $n = length $dna;
  $n -= $a;
  $n -= $c;
  $n -= $g;
  $n -= $t;
  $n -= $gap;

  return ($a, $c, $g, $t, $gap, $n, $length_dna);

}

##########################################
# read chromosome from file
#  my $chrom_seq = &read_chromosome($chromosome);

sub read_chromosome {
  my ($chr_file) = @_;

  my $seq = "";
  open (my $fh, $chr_file) 
      or $log->log_and_die("Can't open the dna file $chr_file : $!\n");
  while(<$fh>) {
    /^\>/ and next;
    /^(\S+)$/ and $seq .= $1;
  }

  return $seq
}

##########################
# zip up files
###########################

sub zip_files {
  my @chromosomes = $wormbase->get_chromosome_names(-mito => 1, -prefix => 1);

  foreach my $chr (@chromosomes){
    my $dna_file = "$dump_dir/${chr}.dna";

    if ($zipdna){
      if (-e $dna_file.".gz" ) {
	$log->write_to("\n ${dna_file}.gz exists\n");
      }
      elsif (-e $dna_file) {
	$log->write_to("\n Compressing $dna_file\n");
	system ("/bin/gzip -f $dna_file");
      }
      else {
	$log->write_to("\n ERROR: Couldn't find any dna chromosome files in $dump_dir\n");
      }
    }
  }
}


#####################################################
# execute ace command: tace or giface
#####################################################

sub execute_ace_command {
  my ($command,$exec,$dir)=@_;
  open (WRITEDB,"| $exec $dir ") or do {
    $log->log_and_die("could not find $exec\n"); # throws only errors if it can't find the command
    die();
  };
  print WRITEDB $command;
  close (WRITEDB) or do {
	  $log->log_and_die("execute_ace_command failed\n"); # throws only errors if the command failed
	  die($!);
  }
}

######################################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

__END__

=pod

=head1 NAME - chromosome_dump.pl

chromosome_dump.pl dumps chromosome fasta files from autoace,
calculates composition stats, and compresses output fules. 

=head1 OPTIONS

=over 4

=item -dna

Dump dna files from specified database (see -p flag), dumps one file for each of 
the six nuclear chromosomes, plus one file for the mitochondrial chromosome.

=item -composition (optional)

Calculates composition statistics for any dna files that are generated.
Use in combination with -d option (see above).  Excludes mitochondrial 
chromosome.

=item -database <database>

Specify database that you wish to dump dna/gff files from.  If -database is not specified
the script will dump from $basedir/autoace by default

=item -dump_dir <destination directory for dump files>

Specify destination of the dna and/or gff dump files generated from the -dna or -gff options.
If -dump_dir is not specified, dump files will be written to $basedir/autoace/CHROMOSOMES by default

=item -zipdna (optional)

Compresses any dna files using gzip (will remove any existing files first).  The -composition
option will be run before this stage if -composition is specified.

=item -help

Show these help files.

=item -debug <user>

Only email log file to specified user

=item -test

Run in test mode and use test environment in ~wormpub/TEST_BUILD

=item -quicktest

Same as -test but only runs analysis for one chromosome (III)

=back

=cut

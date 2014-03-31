#!/usr/local/bin/perl5.8.0 -w
#
# landmark_gene2gff.pl
#
# by Keith Bradnam
#
# script for creating extra GFF lines to indicate those genes that are landmark genes
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2014-03-31 22:22:55 $
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my $database;     # choose another database (defaults to autoace)
my $debug;        # debug mode, output only emailed to one person
my $help;         # help mode, show perldoc
my $test;         # use test mode in ~wormpub
my %landmarks;    # hash containing gene IDs as keys and public_name field as value
my $store;                                       # to specify storable file
my $verbose;
my $gff3;

GetOptions(
	   "help"       => \$help,
	   "debug=s"    => \$debug,
	   "test"       => \$test,
	   "verbose"    => \$verbose,
	   "database=s" => \$database,
	   'store=s'    => \$store,
           'gff3'       => \$gff3,
	   );

############################
# recreate configuration   #
############################
my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test, ) }

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


###########################################
# Variables Part II (depending on $wormbase)    #
###########################################

# database/file paths and locations
my $basedir = $wormbase->basedir;

# set default database if -database not used
$database = $wormbase->autoace if !$database;

###############################
#
#  main part of script
#
###############################

# get list of landmark genes from database, and add to hash
&get_landmark_genes;

# now loop through GFF files to look for existing gene spans
my @chromosomes = $wormbase->get_chromosome_names(-mito => 0, -prefix => 1);

foreach my $chromosome (@chromosomes) {

  $log->write_to("Processing chromosome $chromosome\n");
  
  my $infile = ($gff3) ? $wormbase->GFF3_file_name($chromosome) : $wormbase->GFF_file_name($chromosome);
  my $outfile_base = $wormbase->gff_splits . "/${chromosome}_landmarks";
  my $outfile = ($gff3) ? "${outfile_base}.gff3" : "${outfile_base}.gff";

  # open input/output streams
  open( my $out_fh, ">$outfile" ) or die "Could not open $outfile for writing\n";
  open( my $in_fh, $infile) or die "Could not open $infile for reading\n";
  
  while (<$in_fh>) {
    /^\#/ and next;

    my @data = split(/\t+/, $_);
    next unless $data[2] eq 'gene';
    
    my $gene;
    if ($gff3) {
      ($gene) = $data[8] =~ /Name=Gene:(WBGene\d+)/;
    } else {
      ($gene) = $data[8] =~ /Gene \"(WBGene\d+)\"/;
    }
    
    # check gene from GFF file with genes in hash to see if it is a landmark gene, write to output if so
    if (defined $gene and  $landmarks{$gene} ) {
      my $lm = $landmarks{$gene};
      print $out_fh "$data[0]\tlandmark\tgene\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t";
      if ($gff3) {
        print $out_fh "ID=landmark:$lm;locus=$lm\n";
      } else {
        print $out_fh "Locus \"$lm\"\n";
      }
    }
  }
  close($in_fh);
  close($out_fh) or die "Could not cleanly close output file\n";
}

##################
# Check the files
##################

my $gff_dir = $wormbase->gff_splits;
foreach my $chromosome (@chromosomes) {
  if ($gff3) {
    $wormbase->check_file("$gff_dir/${chromosome}_landmarks.gff3", $log,
                          minsize => 900,
                          maxsize => 2000,
                          lines => ["^${chromosome}\\s+landmark\\s+gene\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+]\\s+\\S+\\s+.*Locus=\\S+"],
        );
  } else {
    $wormbase->check_file("$gff_dir/${chromosome}_landmarks.gff", $log,
                          minsize => 900,
                          maxsize => 2000,
                          lines => ["^${chromosome}\\s+landmark\\s+gene\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+]\\s+\\S+\\s+Locus\\s+\\S+"],
        );
  }
}


# tidy up and exit
$log->mail();
exit(0);


##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

##############################
# get list of landmark genes
##############################

sub get_landmark_genes {

    $log->write_to("Getting list of landmark genes from $database\n");

    my $tace = $wormbase->tace;
    my $db   = Ace->connect(
        -path    => $database,
        -program => $tace
      )
      or $log->log_and_die("Connection failure: ". Ace->error);

    # only want landmark genes which will be in GFF files, i.e. those with Sequence name
    # this is virtually all of them
    my @landmarks = $db->fetch( -query => "Find Gene WHERE Landmark_gene AND Sequence_name" );

    # build hash
    foreach my $gene (@landmarks) {
        my $public_name = $gene->Public_name->name;
        $landmarks{ $gene->name } = $public_name;
    }
    $db->close;

}

__END__
                                                                                                   
=pod
                                                                                                   
=head2 NAME - landmark_genes2gff.pl
                                                                                                   
=head1 USAGE
                                                                                                   
=over 4
                                                                                                   
=item landmark_genes2gff.pl  [-debug -help -database=s -test]
                                                                                                   
=back

This script takes the list of 100 or so 'landmark' genes (those gene objects with a Landmark_gene tag)
and if that gene appears in the GFF files, then the script will write extra GFF output, i.e.

From this:

CHROMOSOME_II   gene    gene    23347   24428   .       +       .       Gene "WBGene00005017"

To this:

CHROMOSOME_II framework gene    23347   24428   .       -       .       Locus bli-2


The final field of the GFF file will contain the 'Public_name' field for the gene object.  This
should always be a CGC style name.

Extra GFF output files will be written to $database/GFF_SPLITS/WSXXX and then need to be appended
to the main GFF files at the end of the build.  This process is for CSHL who need this information.

                                                                                                   
landmark_genes2gff.pl MANDATORY arguments:
                                                                                                   
=over 4

=item none
                                                                                                   
=back

map_Alleles.pl  OPTIONAL arguments:
                                                                                                   
=over 4
                                                                                                   
=item -help
                                                                                                   
this stuff
                                                                                                   
=back
                                                                                                   
=over 4
                                                                                                   
=item -debug <user>

log output goes to <user>
                                                                                                   
=back
                                                                                                   
=over 4
                                                                                                   
=item -database
                                                                                                   
specify a database path (defaults to /wormsrv2/autoace)
                                                                                                   
=back

=over 4
                                                                                                   
=item -test
                                                                                                   
use test build environment in ~wormpub
                                                                                                   
=back
                                                                                                   
=head1 AUTHOR
 
=over 4
 
=item Keith Bradnam (krb@sanger.ac.uk)
 
=back

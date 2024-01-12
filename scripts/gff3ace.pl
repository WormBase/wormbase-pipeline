#!/software/bin/perl -w
#
# Small script to convert GFF3 gene predictions to ace
# It remaps the genomic locations from WS160 to the current version
# To not remap, set the release to be 0.
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2012-06-22 08:56:52 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
use Coords_converter;
use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $output, $source, $release, $method, $species, $addend);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
	    "output:s"   => \$output,
	    "source:s"   => \$source,
	    "release:i"  => \$release, # version to map from, set to zero to not remap at all
	    "method:s"   => \$method, # method to specify in ace output
	    "species:s"  => \$species,
	    "addend"     => \$addend, # use this if the CDS in the GFF does not include the STOP codon (acedb models should include it)
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


#################################
# check input arguments
#################################

die "-method name must be set to give the name of the Method in the ace file\n" unless $method;

$species = $wormbase->full_name;

# the default remapping action is not to remap from a past release
$release = 0 unless $release;

#################################

my $database = $wormbase->autoace;
# my $coords = Coords_converter->invoke($database, 0, $wormbase);
my %clone_child;
my %gene_info;
my %sources;			# keep a note of the sources used

# read in the mapping data
my $assembly_mapper;
if ($release != 0) {
  $assembly_mapper = Remap_Sequence_Change->new($release, $wormbase->get_wormbase_version, $wormbase->species, $wormbase->genome_diffs);
}

# suck the data in
open (IN, "<$input") || die "Can't open $input\n";
open (OUT, ">$output") || die "Can't open $output\n";

while (my $line = <IN>) {
  
  if ($line =~ /^#/) {next;}

  my @fields = split /\s+/, $line;

  # ignore lines that are not from the source we want
  $sources{$fields[1]} = 1;	# keep a note of the sources
  if (defined $source && $fields[1] ne $source) {next;}

  # ignore everything except the CDS, exon and mRNA/transcript lines
  if ($fields[2] ne "CDS" && $fields[2] !~ /exon/ && $fields[2] ne "mRNA" && $fields[2] ne "transcript") {next;}

  my ($clone, $start, $end);

  if ($fields[2] eq "mRNA" || $fields[2] eq "transcript") {
    my @other = split /;/, $fields[8];
    my @id = map { /ID=(\S+)/ ? $1 : ()} @other; 

    #print "$id[0] $fields[0] $fields[3] $fields[4] $fields[6]\n";
    my $id = $id[0];

    # for converting the genBlastG files from Jack Chen's lab we would like to populate the Remark field
    if ($method eq "genBlastG") {
      my @name = map { /Name=(\S+)/ ? $1 : ()} @other; 
      my @pid = map { /PID=(\S+)/ ? $1 : ()} @other; 
      my @coverage = map { /Coverage=(\S+)/ ? $1 : ()} @other; 
      if (@name) {
	$gene_info{$id}{'remark'} = "Elegans homolog: $name[0], homology % ID: $pid[0], homology coverage: $coverage[0]";
      } else {
	$gene_info{$id}{'remark'} = "curated gene used";
      }
    }

    # remap to current genome positions
    if ($release != 0) {
      if ($fields[0] !~ /^CHROMOSOME_/) {$fields[0] = "CHROMOSOME_$fields[0]"};
      my ($indel, $change);	# not used
      ($fields[3], $fields[4], $fields[6], $indel, $change) = $assembly_mapper->remap_gff($fields[0], $fields[3], $fields[4], $fields[6]);
    }

    # get the clone coords
    ($clone, $start, $end) =  ($fields[0],$fields[3],$fields[4]);#$coords->LocateSpan($fields[0], $fields[3], $fields[4]);

    # adjust the end to include the STOP codon if it has been copied over from elegans
    if ($addend) {
      if ($fields[6] eq '+') {$end += 3;} else {$start -= 3;}
    }

    # store the CDS_child data
    $clone_child{$clone}{$id}{'start'} = $start;
    $clone_child{$clone}{$id}{'end'} = $end;
    $clone_child{$clone}{$id}{'sense'} = $fields[6];

    # store some useful information for the gene
    $gene_info{$id}{'clone'} = $clone;
    $gene_info{$id}{'start'} = $start;
    $gene_info{$id}{'end'} = $end;
    $gene_info{$id}{'sense'} = $fields[6];
    if ($fields[6] eq '+') {
      $gene_info{$id}{'chrom_start'} = $fields[3];
    } else {			# reverse sense - we work from the end of the gene when counting exons
      $gene_info{$id}{'chrom_start'} = $fields[4];
    }
    

  } else {
    my @other = split /;/, $fields[8];
    my @id = map { /Parent=(\S+)/ ? $1 : ()} @other; 

    #print "$id[0] $fields[0] $fields[3] $fields[4] $fields[6]\n";
    my $id = $id[0];

    # remap to current genome positions
    if ($release != 0) {
      if ($fields[0] !~ /^CHROMOSOME_/) {$fields[0] = "CHROMOSOME_$fields[0]"};
      my ($indel, $change);	# not used
      ($fields[3], $fields[4], $fields[6], $indel, $change) = $assembly_mapper->remap_gff($fields[0], $fields[3], $fields[4], $fields[6]);
    }

    # get start/end of exon relative to the start of the gene
    my $chrom_start = $gene_info{$id}{'chrom_start'};
    if (! defined $chrom_start) {die "The chrom_start for gene $id is not defined\n";}
    if ($gene_info{$id}{'sense'} eq '+') {
      $start = $fields[3] - $chrom_start + 1;
      $end = $fields[4] - $chrom_start + 1;

      # store the exon
      push @{$gene_info{$id}{'starts'}}, $start;
      push @{$gene_info{$id}{'ends'}}, $end;
    } else {			# reverse sense
      $start = $chrom_start - $fields[3] + 1;
      $end = $chrom_start - $fields[4] + 1;

      # store the exon
      push @{$gene_info{$id}{'starts'}}, $end;
      push @{$gene_info{$id}{'ends'}}, $start;
    }


  }
}

# write out the exons of the gene
foreach my $id (keys %gene_info) {
  if ($#{$gene_info{$id}{'starts'}} == -1) {die "There is no exon information avalable for gene $id\n"}

  # sort the start and end exons into ascending order as there is no guarantee that they are in any order
  @{$gene_info{$id}{'starts'}} = sort {$a <=> $b} @{$gene_info{$id}{'starts'}};
  @{$gene_info{$id}{'ends'}} = sort {$a <=> $b} @{$gene_info{$id}{'ends'}};

  my $clone = $gene_info{$id}{'clone'};
  my $remark = $gene_info{$id}{'remark'};
  print OUT "\nCDS : \"$id\"\n";
  print OUT "Sequence \"$clone\"\n";
  print OUT "Species \"$species\"\n";
  print OUT "CDS\n";
  print OUT "Method \"$method\"\n";
  print OUT "Remark \"$remark\"\n";
  for (my $i = 0; $i < @{$gene_info{$id}{'starts'}}; $i++) {
    my $start = $gene_info{$id}{'starts'}[$i];
    my $end = $gene_info{$id}{'ends'}[$i];
    if ($addend && $i ==  $#{$gene_info{$id}{'ends'}}) {$end += 3;}
    print OUT "Source_exons $start $end\n";
  }
}

# now write out the CDS_child data
foreach my $clone (keys %clone_child) {
  print OUT "\nSequence : \"$clone\"\n";
  foreach my $id (keys %{$clone_child{$clone}}){
    my $start = $clone_child{$clone}{$id}{'start'};
    my $end = $clone_child{$clone}{$id}{'end'};
    my $sense = $clone_child{$clone}{$id}{'sense'};
    if ($sense eq '+') {
      print OUT "CDS_child \"$id\" $start $end\n";
    } else {
      print OUT "CDS_child \"$id\" $end $start\n";
    }
  }
}


close(OUT);
close(IN);


# report the sources used
print "\n\nSources read in:\n";
foreach my $s (keys %sources) {
  print "$s";
  if (! defined $source || $s eq $source) {print "\t- processed to ace";}
  print "\n";
}


$log->mail();
print "Finished.\n" if ($verbose);
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




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - gff3ace.pl

=head1 USAGE

=over 4

=item gff3ace.pl  [-options]

=back

This script reads in a GF3 file of gene predictions and writes out an ACE file of the CDS regions.

For example:

perl gff3ace.pl -in AUGUSTUS_validated.gff3.2 -out augustus.ace -source AUGUSTUS_cat3_ver1 -method augustus

perl gff3ace.pl -in JIGSAW.gff3  -out jigsaw.ace -source jigsaw-3.2.8-lin_phaseII_cat4_ver2 -method jigsaw



script_template.pl MANDATORY arguments:

=over 4

=item -input input file of gene predictions in GFF3 format

=back

=over 4

=item -output output ACE file of CDS objects

=back

=over 4

=item -method name of Method class of these predictions

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -release version_number to remap the genomic positions from. It is assumed by default that there is no remapping. 

=back

=over 4

=item -source source GFF3 field name to use. By default this script will read all source lines in the input file. specifying an explicit source name will restrict the processing to just those lines with the specified source.

=back

=over 4

=item -species species_name. By default, this script will write ACE data specifying the species as 'elegans'. This specifies a different species.

=back

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut

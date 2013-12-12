#!/usr/bin/env perl
# 
# make_pseudo_map_positions.pl
#
# Script to identify genes which can have their Interpolated_map_position tag promoted to a Map position
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-12-12 18:58:42 $

use strict;
use warnings;
use lib $ENV{CVS_DIR};
use Wormbase;
use Ace;
use Getopt::Long;


###################################################
# command line options                            # 
###################################################

my $maintainers = "All";
my ($help, $database, $verbose, $test, $noload,$store,$debug,$acefile);

GetOptions (	
  "help"        => \$help,
  "database=s"  => \$database,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "noload"      => \$noload,
  'store=s'     => \$store,
  'debug=s'     => \$debug,
  'acefile=s'   => \$acefile,
    );

# Display help if required
&usage("Help") if $help;

###################################################
# Miscellaneous important variables               # 
###################################################

# recreate configuration
my $wb;
if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n")
} else { 
  $wb = Wormbase->new( 
    -debug => $debug, 
    -test => $test, ); 
}

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wb->basedir;

# create log file, open output file handles
my $log = Log_files->make_build_log($wb);

$database = $wb->autoace if (!$database);           # specify autoace as the default database
$acefile  = $wb->acefiles."/pseudo_map_positions.ace" if not defined $acefile;


###############################################################################################
# Look for CGC-named Gene objects without a 'Map' tag or mapping_data but which have an allele
# and a sequence connection and an 'Interpolated_map_position' tag.  This is for creating 
# inferred multipoint data for Jonathan Hodgkin to approve
###############################################################################################

$log->write_to("Genes with newly promoted Map positions:\n\n");

# open a connection to database
my $db = Ace->connect(-path  => $database ) ||  $log->log_and_die("Connection failure: ",Ace->error);

# build query to find candidate genes
my $query  = "find Gene * WHERE !Mapping_data & Allele & CGC_name & Sequence_name & Interpolated_map_position & Species =\"*elegans\"";
push(my @candidate_genes, $db->find($query) );

my @to_promote;
foreach my $gene (@candidate_genes){

  # get basic details from gene
  my $map      = $gene->Interpolated_map_position(1);
  my $position = $gene->Interpolated_map_position(2);
  my $cgc_name = $gene->CGC_name;
  my @alleles  = $gene->Allele(1);

  my @supporting_alleles;
  foreach my $allele (@alleles){

    # only want to consider alleles which are not Transposon insertions or SNPs
    next if(defined($allele->Transposon_insertion));

    my $method = $allele->Method;
    next if not grep { $method eq $_ } ('Allele',
                                        'Deletion_allele',
                                        'Deletion_and_insertion_allele',
                                        'Insertion_allele',
                                        'Substitution_allele',
                                        'KO_consortium_allele',
                                        'NBP_knockout_allele',
                                        'NemaGENETAG_consortium_allele');

    push @supporting_alleles, $allele;
  }
  if (@supporting_alleles) {
    push @to_promote, {
      gene => $gene,
      cgc  => $cgc_name,
      map  => $map,
      pos  => $position,
      alleles => \@supporting_alleles,
    };
  }
}


# open main output file to be loaded into autoace (and subsequently geneace)
open(OUT, ">$acefile") || die $!;
print OUT "// this file contains details of Genes with interpolated map positions\n";
print OUT "// that can be 'upgraded' to a (pseudo) Map position.  We do this only\n";
print OUT "// for C. elegans genes that have a CGC-name, allele connection (but not a\n";
print OUT "// transposon insertion), and a sequence connection\n\n";

print OUT "// This is done at the request of Jonathan Hodgkin so that more genes\n";
print OUT "// can be placed on the genetic map.  A later build script will have to\n";
print OUT "// make associated Multi_pt_data objects to accompany these genes (all\n";
print OUT "// genes with map positions need to have mapping data).\n\n";

print OUT "// This file should also be uploaded back into geneace...normally this\n";
print OUT "// happen automatically by a later build script\n\n";

foreach my $obj (sort { $a->{map} cmp $b->{map} or $a->{pos} <=> $b->{pos} } @to_promote) {
  my @allele_names = map { $_->Public_name } @{$obj->{alleles}};
  $log->write_to(sprintf("%-10s\t%-6s\t'promoted' map position: %3s\t%10.5f\talleles: %s\n", 
                         $obj->{gene}, 
                         $obj->{cgc},
                         $obj->{map}, 
                         $obj->{pos}, 
                         join(",", sort @allele_names)));
  # first remove existing Interpolated_map_position before adding new data
  print OUT "\nGene : \"$obj->{gene}\" \/\/ $obj->{cgc}\n";
  print OUT "-D Interpolated_map_position\n\n";
  
  print OUT "Gene : \"$obj->{gene}\"\n";
  print OUT "Map \"$obj->{map}\" Position $obj->{pos}\n";
  print OUT "Pseudo_map_position\n";
  print OUT "Remark \"Map position created from combination of previous interpolated map position (based on known location of sequence) and allele information.  Therefore this is not a genetic map position based on recombination frequencies or genetic experiments.  This was done on advice of the CGC.\" CGC_data_submission\n\n\n";
}
close(OUT);
$db->close;

$log->write_to("\nTotal: ". scalar(@to_promote) . " genes to become inferred genetic marker(s)\n\n");

# we do not wish to throw an error if there are large differences in
# the number of objects loaded from this file between on Build and the
# next

unless ($noload) {
  my $accept_large_differences = 1; 
  $wb->load_to_database($wb->autoace, $acefile, 'pseudo_map_postn', $log, 0, $accept_large_differences);
}

$log->mail();
$log->mail("genenames\@wormbase.org", "list of promoted map positions") unless $wb->test;

exit(0);


sub usage {
  my $error = shift;
  if ($error == 0) {
    # Normal help menu
    exec ('perldoc',$0);
  }
}



#--------------------------------------------------------------------------------------------------------------------




__END__

                                                                                                       
=pod

=head2   NAME - make_pseudo_map_positions.pl

=head1 USAGE

=over 4

=item make_pseudo_map_positions.pl -[options]


=back

=head1 DESCRIPTION


In every build there are usually some Gene objects with 'Interpolated_map_positions' which can be
'promoted' into having 'pseudo' map positions.  I.e. they will now use a Map tag and appear on 
the genetic map as real markers.  Jonathan Hodgkin gets us to do this to add more markers to the
map.

The critera for genes to be promoted is that they have to be attached to a sequence, have a CGC-name,
and have a (non-Transposon insertion) allele attached.  Finally, we only do this for C. elegans 
genes (as we only have a C. elegans genetic map).

This script makes an acefile that gets loaded into the build which will remove existing 
Interpolated_map_position tags for some genes and transfer the map value and position to the Map tag
instead.  We also flag these genes with the 'Pseudo_map_position' tag and leave a suitable remark.

The script will make an acefile at $wb->autoace/acefiles/pseudo_map_positions.ace
This acefile will be loaded to autoace by default,  but it also needs to be loaded back 
into geneace.  Jonathan receives a separate email listing the genes which have been changed.


=back

=head1 MANDATORY arguments: <none>

=back

=head1 OPTIONAL arguments: -help, -database, -verbose, -test, -store


=over 4

=item -help

Displays this help

=item -database

Specify a path to a valid acedb database to query (defaults to /wormsrv2/autoace)

=item -verbose

Display names of promoted genes as the script finds them

=item -test

Use the test build environment

=item -store

Use a stored configuration

=back

=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)


but based on original code in another script by Chao-kung Chen.

=back

=cut

/#!/usr/bin/perl5.8.0 -w
# creating WBGene objects from gene data
# mt3 June 2007

use strict;

# declare variables
my $WBGene;
my $genename;
my $person_id;
my $genenamefile = "processgene.txt";
my $output = "processgene.ace";
my @acefile_values;
my $i;

# read in input file
open ( INFILE, "<$genenamefile" ) or die ( "Couldn't open file
$genenamefile: $!\n" );

# load $genenamefile into array line by line
@acefile_values = <INFILE>;
close INFILE;

# open output file
open (OUTFILE,  ">$output");

for $i ( 0 .. $#acefile_values ) {

# get first line of acefile values array
#print $acefile_values[ $i ];
# split first line of array into scalars
( $WBGene, $genename, $person_id ) = split( / / , $acefile_values[ $i ]);
# print these values out.
 print ( "WBGene $WBGene, Gene name $genename, Person $person_id" );

# check format of input file (WBGeneID, genename, person_evidence)
$WBGene =~ ( m/WBGene\w+/ ) or die;
$genename =~ ( m/\w+\-\d+/ ) or die;
$person_id =~ ( m/WBPerson\d+/ ) or die;
chomp $person_id;

# Add WBGeneID to file
print ( OUTFILE "Gene : $WBGene\n" );

# Add Version = 2
print ( OUTFILE "Version 2\n" );

# Add History - includes genename
print ( OUTFILE "History Version_change 2 now WBPerson2970 CGC_name $genename\n");

# Add CGC_name and person_evidence
print ( OUTFILE "CGC_name $genename Person_evidence $person_id\n");

# Add Public_name
print (OUTFILE "Public_name $genename\n");

# Add Gene_class
$genename =~ ( m/(\w+)\-\d+/ ) ;
my $gene_class = $1;
print (OUTFILE "Gene_class $gene_class\n\n");
}
close ( OUTFILE );
print ( "\nFinished.\n" );

=pod

=head2 NAME - processgene.pl

=head1 USAGE

=over 4

=item processgene.pl 

=back

=head1 DESCRIPTION

A script designed to create an .ace file of the tags needed in Gene objects 
when CGC_names are assigned.
The script requires a file called processgene.txt as input. The file must
be in the format of GeneID CGC_name WBPerson e.g. 

WBGene00021073 nol-1 WBPerson733

The script must be run from the directory where processgene.txt resides and takes processgene.txt 
as input and writes the output into processgene.ace

e.g 

> perl ~wormpub/wormbase/scripts/GENEACE/processgene.pl

The output processgene.txt would look like this:

 Gene : WBGene00021073
 Version 2
 History Version_change 2 now WBPerson2970 CGC_name nol-1
 CGC_name nol-1 Person_evidence WBPerson733
 Public_name nol-1
 Gene_class nol

The script assumes the Version to be 2 (based on the most common update), and
this MUST be checked before loading into Geneace.

=head1 AUTHOR Mary Ann Tuli (mt3@sanger.ac.uk)

=back

=cut

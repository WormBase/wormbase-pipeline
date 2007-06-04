#!/usr/bin/perl5.8.0 -w
# creating WBGene objects from gene data
# mt3 June 2007

use strict;

# declare variables
my $WBGene;
my $genename;
my $person_id;
my $genenamefile = "newgenes.txt";
my $output = "newgene.ace";
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
print ( OUTFILE "History now WBPerson0002970 CGC_name $genename\n");

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

#!/usr/local/bin/perl5.8.0 -w
#
# update_web_gene_names.pl
#
# completely rewritten by Keith Bradnam from list_loci_designations
#
# Last updated by: $Author: wormpub $     
# Last updated on: $Date: 2006-02-13 14:57:25 $      
#
# This script should be run under a cron job and simply update the webpages that show
# current gene names and sequence connections.  Gets info from geneace.  

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Carp;
use Getopt::Long;



##############################
# Command line options       #
##############################

my $weekly; # for weekly cronjob which will interrogate current_DB
my $daily;  # daily updates which will interrogate /nfs/disk100/wormpub/DATABASES/geneace
 
GetOptions ("weekly" => \$weekly,
            "daily"  => \$daily);


##############################
# Script variables (run)     #
##############################

my $tace  = &tace;
my $www = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/LOCI"; # where output will be going

my $rundate = &rundate;
my $log = "/tmp/update_web_gene_names";
my $database; 



# Set up log file

open(LOG,">$log") || carp "Couldn't open tmp log file\n";
print LOG &runtime, " :Started running update_web_gene_names.pl\n\n";

die "Can't run both -weekly and -daily at the same time!\n" if ($weekly && $daily);

# make the a-z lists based on CGC_name using current_DB
if($weekly){
  print LOG "Creating loci pages based on current_DB\n";
  $database = "/nfs/disk100/wormpub/DATABASES/current_DB";
  &create_currentDB_loci_pages;
}

# make lists of gene2molecular_name and molecular_name2gene
if($daily){
  print LOG "Making daily update lists\n";
  $database = "/nfs/disk100/wormpub/DATABASES/geneace";
  &make_gene_lists;
}


###################################################
# Tidy up - close things, mail log, run webpublish
###################################################

# now update pages using webpublish
chdir($www) || print LOG "Couldn't run chdir\n";


if($weekly){
  system("/usr/local/bin/webpublish -f -q *.shtml") && print LOG "Couldn't run webpublish on html files\n";
}
system("/usr/local/bin/webpublish -f -q *.txt") && print LOG "Couldn't run webpublish on text file\n";

print LOG &runtime, " : Finished running script\n";

&mail_maintainer("update_web_gene_names.pl","mt3\@sanger.ac.uk","$log");

close(LOG);
exit(0);




#######################################################################################
#
#                    T  H  E    S  U  B  R  O  U  T  I  N  E  S
#
#######################################################################################


sub create_currentDB_loci_pages{

  # query against current_DB
  my $db = Ace->connect(-path  => "$database",
                        -program =>$tace) || do { print "Connection failure: ",Ace->error; croak();};


  # open text file which will contain all genes
  open (TEXT, ">$www/loci_all.txt") || croak "Couldn't open text file for writing to\n";
  print TEXT "Gene name, WormBase gene ID, Gene version, CDS name (Wormpep ID), RNA gene name, pseudogene name, other names, cgc approved?\n";

  foreach my $letter ("a".."z"){
    # Get all Loci
    my @gene_names = $db->fetch(-query=>"Find Gene_name $letter\* WHERE Public_name_for");
    
    # loop through each file (one for each letter a-z)
    open (HTML, ">$www/loci_designations_${letter}.shtml") || croak "Couldn't open file for writing to\n";
    
    my $line = 0;

    # cycle through each locus in database
    foreach my $gene_name (@gene_names){


      # skip gene names that are just sequence names
      next unless ($gene_name->CGC_name_for || $gene_name->Other_name_for);

      # need to get Gene object via CGC_name_for or Other_name_for tags
      my $gene;
      my $public_name;

      if($gene_name->CGC_name_for){
	$gene = $gene_name->CGC_name_for;
	$public_name = $gene->CGC_name;
      }
      else{
	$gene = $gene_name->Other_name_for;
	$public_name = $gene->Other_name;
      }
      
      # ignore non C. elegans genes for now
      my $species = $gene->Species;
      next unless ($species eq "Caenorhabditis elegans");

      # ignore dead genes
      next if($gene->Status->name ne "Live");

      # Set alternating colours for each row of (HTML) output 
      if (($line % 2) == 0) { 
	print HTML "<TR BGCOLOR=\"lightblue\">\n";
      }
      else {
	print HTML "<TR BGCOLOR=\"white\">\n";
      }
      
      # Column 1 - ?Gene name
      print HTML "<TD align=center><A HREF=\"http://www.wormbase.org/db/gene/gene?name=${public_name}\">${public_name}</a></TD>";
      print TEXT "$public_name,";
      
      
      # Column 2 - Gene ID
      print HTML "<TD align=center><A HREF=\"http://www.wormbase.org/db/gene/gene?name=${gene};class=Gene\">${gene}</a></TD>";
      print TEXT "$gene,";
      
      # Column 3 - Gene Version
      my $version = $gene->Version;
      print HTML "<TD align=center>$version</TD>";
      print TEXT "$version,";
      
      # Column 4 - ?CDS connections
      if(defined($gene->at('Molecular_info.Corresponding_CDS'))){
	my @CDSs = $gene->Corresponding_CDS;
	print HTML "<TD>";
	foreach my $cds (@CDSs){
	  # also get wormpep identifier for each protein
	  my $protein = $cds->Corresponding_protein;
	  print HTML "<A HREF=\"http://www.wormbase.org/db/gene/gene?name=${cds};class=CDS\">${cds}</a> ";
	  print HTML "(<A HREF=\"http://www.wormbase.org/db/seq/protein?name=${protein};class=Protein\">${protein}</a>) ";
	  print TEXT "$cds ($protein) ";
	}
	print TEXT ",,,";
	print HTML "</TD><TD>&nbsp</TD><TD>&nbsp</TD>";
      }
      
      
      # Column 5 - ?Transcript connections
      elsif(defined($gene->at('Molecular_info.Corresponding_transcript'))){
	print HTML "<TD>&nbsp</TD>";
	my @transcripts = $gene->Corresponding_transcript;
	print HTML "<TD>";
	print TEXT ",";
	foreach my $i (@transcripts){
	  print HTML "<A HREF=\"http://www.wormbase.org/db/seq/sequence?name=${i}\">${i}</a> ";
	  print TEXT "$i ";
	}
	print TEXT ",,";
	print HTML "</TD><TD>&nbsp</TD>";
      }
      
      # Column 6 - ?Pseudogene connections
      elsif(defined($gene->at('Molecular_info.Corresponding_pseudogene'))){
	my @pseudogenes = $gene->Corresponding_pseudogene;
	print HTML "<TD>&nbsp</TD><TD>&nbsp</TD><TD>";
	print TEXT ",,";
	foreach my $i (@pseudogenes){
	  print HTML "<A HREF=\"http://www.wormbase.org/db/seq/sequence?name=${i}\">${i}</a> ";
	  print TEXT "$i ";
	}
	print HTML "</TD>";
	print TEXT ",";
      }
      
      # Blank columns if no ?Sequence, ?Transcript, or ?Pseudogene
      else{
	print HTML "<TD>&nbsp</TD><TD>&nbsp</TD><TD>&nbsp</TD>";
	print TEXT ",,,";
      }
      
      
      # Column 7 - Other names for ?Gene
      if(defined($gene->at('Identity.Name.Other_name'))){
	my @other_names = $gene->Other_name;
	print HTML "<TD>";
	foreach my $i (@other_names){
	  print HTML "${i} ";
	  print TEXT "$i ";
	}
	print HTML "</TD>";
      }
      else{
	print HTML "<TD>&nbsp</TD>";
      }
      print TEXT ",";
      
      
      # Column 8 CGC approved?
      if(defined($gene->at('Identity.Name.CGC_name'))){
	print HTML "<TD align=center>approved</TD>\n";
	print TEXT "approved"
	  }
      else{
	print HTML"<TD>&nbsp<TD>\n";
      }
      
      $line++;
      print HTML "</TR>\n";
      print TEXT "\n";
      $gene->DESTROY();
    }
    close(HTML);    
  }
  close(TEXT);

  $db->close;


}


############################################################

sub make_gene_lists{

  my %molecular_name2gene;
  my %gene2molecular_name;
  my %transposon_genes;
  
  # connect to AceDB using TableMaker, 
  my $command="Table-maker -p /nfs/disk100/wormpub/DATABASES/geneace/wquery/gene2molecular_name.def\nquit\n";
  open (TACE, "echo '$command' | $tace /nfs/disk100/wormpub/DATABASES/geneace |") || print LOG "ERROR: Can't open tace connection to /nfs/disk100/wormpub/DATABASES/geneace\n";
  while (<TACE>) {
    chomp;
    # skip any acedb banner text (table maker output has all fields surrounded by "")
    next if ($_ !~ m/^\"/);
    # skip acedb prompts
    next if (/acedb/);
    # skip empty fields
    next if ($_ eq "");
                                                                                           
    # get rid of quote marks
    s/\"//g;
                                                                                           
    # split the line into various fields
    my ($gene,$sequence_name,$cgc_name) = split(/\t/, $_) ;

    # populate hashes, appending CGC name if present
    if(defined($cgc_name)){
      $molecular_name2gene{$sequence_name} = "$gene $cgc_name";
      $gene2molecular_name{$gene} = "$sequence_name $cgc_name";
    }
    else{
      $molecular_name2gene{$sequence_name} = $gene;
      $gene2molecular_name{$gene} = $sequence_name;
    }
  }
  close TACE;

  # now fire off second query to get dead genes which were made into Transposons
  # this is to help Darin
  $command = "Table-maker -p /nfs/disk100/wormpub/DATABASES/geneace/wquery/genes_made_into_transposons.def\nquit\n";
  open (TACE, "echo '$command' | $tace /nfs/disk100/wormpub/DATABASES/geneace |") || print LOG "ERROR: Can't open tace connection to /nfs/disk100/wormpub/DATABASES//geneace\n";
  while (<TACE>) {
    chomp;
    # skip any acedb banner text (table maker output has all fields surrounded by "")
    next if ($_ !~ m/^\"/);
    # skip acedb prompts
    next if (/acedb/);
    # skip empty fields
    next if ($_ eq "");
                                                                                           
    # get rid of quote marks
    s/\"//g;
                                                                                           
    # split the line into various fields
    my ($gene,$public_name) = split(/\t/, $_) ;
    $transposon_genes{$gene} = $public_name;
  }

  # set up various output files (first two are reverse of each other)

  open (GENE2MOL, ">$www/genes2molecular_names.txt") || die "ERROR: Couldn't open genes2molecular_names.txt  $!\n";
  foreach my $key (sort keys %gene2molecular_name){
    print GENE2MOL "$key\t$gene2molecular_name{$key}\n";	      
  }
  close(GENE2MOL);


  open (MOL2GENE, ">$www/molecular_names2genes.txt") || die "ERROR: Couldn't open molecular_names2genes.txt $!\n";
  foreach my $key (sort keys %molecular_name2gene){
    print MOL2GENE "$key\t$molecular_name2gene{$key}\n";
  }
  close(MOL2GENE);

  open (TRANSPOSONS, ">$www/transposon_genes.txt") || die "ERROR: Couldn't open transposon_genes.txt $!\n";
  foreach my $key (sort keys %transposon_genes){
    print TRANSPOSONS "$key\t$transposon_genes{$key}\n";
  }
  close(TRANSPOSONS);

}



__END__

=pod

=head1 NAME - update_web_gene_names.pl

=back

=head1 USAGE

=over 4

=item update_web_gene_names.pl 

Simply takes the latest set of gene names in geneace and writes to the development web site
a set of HTML pages (one for each letter of the alphabet) containing all gene names starting
with that letter.  Makes these tables hyperlinked to WormBase and also includes other names
and CDS (with wormpep)/transcript/pseudogene connections.  Also includes Gene IDs and 
version numbers as extra columns.

Additionally makes two other files which is just gene 2 molecular name and vice versa.

When script finishes it copies across to the live web site.  This script should normally be
run every night on a cron job for the genes2molecular_names.txt file and weekly for the a-z pages.

=back

=head2 camcheck.pl MANDATORY arguments:

=over 4

=item none

=back

=head2 camcheck.pl OPTIONAL arguments:

=over 4

=item -daily (queries geneace: makes gene2molecular_name and molecular_name2gene files)

=item -weekly (queries current_DB, makes a-z lists)

=back


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

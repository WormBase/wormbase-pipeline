#!/usr/local/bin/perl5.8.0 -w
#
# update_web_gene_names.pl
#
# completely rewritten by Keith Bradnam from list_loci_designations
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-07-14 09:52:34 $      
#
# This script should be run under a cron job and simply update the webpages that show
# current gene names and sequence connections.  Gets info from geneace.  


#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Carp;

##############################
# Script variables (run)     #
##############################

my $tace  = &tace;
my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans/LOCI"; # where output will be going

my $rundate = &rundate;
my $log = "/tmp/update_web_gene_names";

open(LOG,">$log") || carp "Couldn't open tmp log file\n";
print LOG "Running update_web_gene_names.pl on $rundate\n\n";


# make the a-z lists based on CGC_name
 print LOG "Creating loci pages\n";
&create_loci_pages;

# now make lists of gene2molecular_name and molecular_name2gene
print LOG "Making gene lists\n";
&make_gene_lists;

###################################################
# Tidy up - close things, mail log, run webpublish
###################################################


# now update pages using webpublish
chdir($www) || print LOG "Couldn't run chdir\n";

system("/usr/local/bin/webpublish -f -q *.shtml") && print LOG "Couldn't run webpublish on html files\n";
system("/usr/local/bin/webpublish -f -q *.txt") && print LOG "Couldn't run webpublish on text file\n";
system("/usr/local/bin/webpublish -f -q *.ace") && print LOG "Couldn't run webpublish on ace file\n";
#
&mail_maintainer("update_web_gene_names.pl","krb\@sanger.ac.uk","$log");

close(LOG);
exit(0);




#######################################################################################
#
#                    T  H  E    S  U  B  R  O  U  T  I  N  E  S
#
#######################################################################################


sub create_loci_pages{

  # query against active geneace database
  my $db = Ace->connect(-path  => "/wormsrv1/geneace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; croak();};

  # open text file which will contain all genes
  open (TEXT, ">$www/loci_all.txt") || croak "Couldn't open text file for writing to\n";
  print TEXT "CGC gene name, WormBase gene ID, CDS name, transcript name, pseudogene name, other names, cgc approved?\n";

  foreach my $letter ("a".."z"){
    # Get all Loci
    my @genes = $db->fetch(-query=>"Find Gene_name $letter\* WHERE CGC_name_for; Follow CGC_name_for; \(Species = \"Caenorhabditis elegans\"\)");
    
    # loop through each file (one for each letter a-z)
    open (HTML, ">$www/loci_designations_${letter}.shtml") || croak "Couldn't open file for writing to\n";
    # Text file for simpler handling
    
    my $line = 0;

    # cycle through each locus in database
    foreach my $gene (@genes){
      
      my $cgc_name = $gene->CGC_name;
      
      # Set alternating colours for each row of (HTML) output 
      if (($line % 2) == 0) { 
	print HTML "<TR BGCOLOR=\"lightblue\">\n";
      }
      else {
	print HTML "<TR BGCOLOR=\"white\">\n";
      }
      
      # Column 1 - ?Gene name
      print HTML "<TD><A HREF=\"http://www.wormbase.org/db/gene/gene?name=${cgc_name}\">${cgc_name}</a></TD>";
      print TEXT "$cgc_name,";
      
      
      # Column 2 - Gene ID
      print HTML "<TD>$gene</TD>";
      print TEXT "$gene,";
      
      
      # Column 3 - ?CDS connections
      if(defined($gene->at('Molecular_info.CDS'))){
	my @CDSs = $gene->CDS;
	print HTML "<TD>";
	foreach my $i (@CDSs){
	  print HTML "<A HREF=\"http://www.wormbase.org/db/seq/sequence?name=${i}\">${i}</a> ";
	  print TEXT "$i ";
	}
	print TEXT ",,,";
	print HTML "</TD><TD>&nbsp</TD><TD>&nbsp</TD>";
      }
      
      
      # Column 4 -  ?Transcript connections
      elsif(defined($gene->at('Molecular_info.Transcript'))){
	print HTML "<TD>&nbsp</TD>";
	my @transcripts = $gene->Transcript;
	print HTML "<TD>";
	print TEXT ",";
	foreach my $i (@transcripts){
	  print HTML "<A HREF=\"http://www.wormbase.org/db/seq/sequence?name=${i}\">${i}</a> ";
	  print TEXT "$i ";
	}
	print TEXT ",,";
	print HTML "</TD><TD>&nbsp</TD>";
      }
      
      # Column 5 - ?Pseudogene connections
      elsif(defined($gene->at('Molecular_info.Pseudogene'))){
	my @pseudogenes = $gene->Pseudogene;
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
      
      
      # Column 6 - Other names for ?Gene
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
      
      
      # Column 7 CGC approved?
      if(defined($gene->at('Identity.Name.CGC_name'))){
	print HTML "<TD>CGC approved</TD>\n";
	print TEXT "CGC approved"
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
  my %gene2public_name;

  # connect to AceDB using TableMaker, 
  my $command="Table-maker -p /wormsrv1/geneace/wquery/gene2molecular_name.def\nquit\n";
  print LOG "Using command $command\n";
  print LOG "tace is set to $tace\n";
  open (TACE, "echo '$command' | $tace /wormsrv1/geneace |") || print LOG "ERROR: Can't open tace connection to /wormsrv1/geneace\n";
  while (<TACE>) {
    chomp;
    next if (/acedb/);
    next if ($_ eq "");
    last if (/\/\//);
                                                                                           
    # get rid of quote marks
    s/\"//g;
                                                                                           
    # split the line into various fields
    my ($gene,$cds,$transcript,$pseudogene,$public_name) = split(/\t/, $_) ;

    $gene2public_name{$gene} = "$public_name";
    print LOG "$public_name\n";
    if($cds){
      $gene2molecular_name{$gene} .= "CDS_$cds ";
      $molecular_name2gene{$cds} = $gene;
    }
    if($transcript){
      $gene2molecular_name{$gene} .= "Transcript_$transcript ";
      $molecular_name2gene{$transcript} = $gene;
    }
    if($pseudogene){
      $gene2molecular_name{$gene} .= "Pseudogene_$pseudogene ";
      $molecular_name2gene{$pseudogene} = $gene;
    }
  }
  close TACE;
  



  # set up three output files, two are reverse of each other then one ace file so Darin or Paul
  # can upload to camace/stlace if necessary
  open (ACE, ">$www/genes2molecular_names.ace") || print LOG "Couldn't open genes2molecular_names.ace\n";
  open (GENE2MOL, ">$www/genes2molecular_names.txt") || print LOG "Couldn't open genes2molecular_names.txt\n";

  foreach my $key (sort keys %gene2molecular_name){
    print ACE "Gene : $key\n";

    # need to check if there is more than one worm_gene per Gene object
    my @names = split(/ /,$gene2molecular_name{$key});
    foreach my $name (@names){
      # treat differently depending on class
      if($name =~ m/^CDS_/){
	$name =~ s/CDS_//;
	print ACE "CDS $name\n";
      }
      elsif($name =~ m/^Transcript_/){
	$name =~ s/Transcript_//;
	print ACE "Transcript $name\n";
      }
      elsif($name =~ m/^Pseudogene_/){
	$name =~ s/Pseudogene_//;
	print ACE "Pseudogene $name\n";
      }    
      print GENE2MOL "$key\t$name\n";	      
    }
    # one Public_name field always per Gene
    print ACE "Public_name $gene2public_name{$key}\n\n";
  }

  open (MOL2GENE, ">$www/molecular_names2genes.txt") || croak "Couldn't open molecular_names2genes.txt\n";
  foreach my $key (sort keys %molecular_name2gene){
    print MOL2GENE "$key\t$molecular_name2gene{$key}\n";
  }

  close(GENE2MOL);
  close(MOL2GENE);
  close(ACE);
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
and sequence/transcript/pseudogene connections.  Now includes Gene IDs as an extra column
but this means that Loci are not sorted alphabetically at the moment.

When script finishes it copies across to the live web site.  This script should normally be
run every night on a cron job.

=back

=head2 camcheck.pl MANDATORY arguments:

=over 4

=item none

=back

=head2 camcheck.pl OPTIONAL arguments:

=over 4

=item none

=back


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

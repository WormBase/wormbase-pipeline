#!/usr/local/bin/perl5.8.0 -w
#
# update_web_gene_names.pl
#
# completely rewritten by Keith Bradnam from list_loci_designations
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-04-14 15:29:49 $      
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

# query against active geneace database
my $db = Ace->connect(-path  => "/wormsrv1/geneace",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; croak();};

# Get all Loci
my @genes = $db->fetch(-query=>"Find Gene_name WHERE CGC_name_for; Follow CGC_name_For; \(Species = \"Caenorhabditis elegans\"\)");


# Need to open first file in series to make this all work
open (HTML, ">$www/loci_designations_a.shtml") || croak "Couldn't open file for writing to\n";
# Text file for simpler handling
open (TEXT, ">$www/loci_all.txt") || croak "Couldn't open text file for writing to\n";
print TEXT "Locus, sequence name, transcript name, pseudogene name, other names, cgc approved?\n";
my $prev_initial = "a";


my $line = 0;

# cycle through each locus in database
foreach my $gene (@genes){
  
  my $cgc_name = $gene->CGC_name;
  print "$gene -> $cgc_name\n";
  my $initial = substr($cgc_name,0,1);
  $initial = lc($initial);
  if ($initial ne $prev_initial){
    close(HTML);
    open (HTML, ">$www/loci_designations_$initial.shtml") || croak "Couldn't open file for writing to\n";
    $prev_initial = $initial;
  }

  # Set alternating colours for each row of (HTML) output 
  if (($line % 2) == 0) { 
    print HTML "<TR BGCOLOR=\"lightblue\">\n";
  }
  else {
    print HTML "<TR BGCOLOR=\"white\">\n";
  }
  
  # Column 1 - ?Locus name
  print HTML "<TD><A HREF=\"http://www.wormbase.org/db/gene/gene?name=${cgc_name}\">${cgc_name}</a></TD>";
  print TEXT "$gene,";


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


  # Column 6 - Other names for ?Locus
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


###################################################
# Tidy up - close things, mail log, run webpublish
###################################################

close(HTML);  

$db->close;

my $rundate = &rundate;
my $log = "/tmp/update_web_gene_names";

open(LOG,">$log") || carp "Couldn't open tmp log file\n";

print LOG "Running update_web_gene_names.pl on $rundate\n\n";

# now update pages using webpublish
chdir($www) || print LOG "Couldn't run chdir\n";

#system("/usr/local/bin/webpublish -f -q *.shtml") && print LOG "Couldn't run webpublish on html files\n";
#system("/usr/local/bin/webpublish -f -q *.txt") && print LOG "Couldn't run webpublish on text file\n";

&mail_maintainer("update_web_gene_names.pl","krb\@sanger.ac.uk","$log");

close(LOG);
exit(0);



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

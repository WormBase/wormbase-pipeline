#!/usr/local/bin/perl5.8.0 -w
#
# fetch_EMBL_seqs_for_blatting.pl
#
# by Keith Bradnam
# 
# Attempt to unify all of the diverse scripts to fetch ESTs, OSTs, mRNAs etc. used by blat 
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2003-09-10 18:11:04 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;
use Data::Dumper;

##############################
# command-line options       #
##############################

my $help;               # Help/Usage page
my $ace;                # Dump acefile
my $debug;              # for debug mode
my $verbose;            # turn on extra output
my $blastdb;            # make blast database using pressdb?
my ($est, $mrna, $ost, $nematode, $embl);


GetOptions (
	    "ace"      => \$ace,
	    "est"      => \$est,
	    "mrna"     => \$mrna,
	    "ost"      => \$ost,
	    "nematode" => \$nematode,
	    "embl"     => \$embl,
	    "debug=s"  => \$debug,
	    "verbose"  => \$verbose,
	    "blastdb"  => \$blastdb,
            "help"     => \$help
	    );

# Help pod if needed
&usage(0) if ($help);


##############################
# Script variables (run)     #
##############################

my $acc;      # EMBL accession
my $id;       # EMBL ID   
my $sv;       # EMBL sequence version
my $def;      # EMBL description
my $protid;   # EMBL Protein_ID
my $protver;  # EMBL Protein_ID version
my $org;      # EMBL species

my $ost_seq;  # tag for OST/EST split
my %EST_name; # EST accession => name
my %EST_dir;  # EST accession => orientation [5|3]

my $dir = "/nfs/disk100/wormpub/analysis/ESTs"; # path for output files
my $getz = "/usr/local/pubseq/bin/getzc"; # getz binary


&make_ests          if ($est || $ost);
&make_mrnas         if ($mrna);
&make_embl_cds      if ($embl);
&make_nematode_ests if ($nematode);

exit(0);


#################################################################
#
# C. elegans ESTs and OSTs
#
#################################################################

sub make_ests{

  # read accession->yk name hash for EST names 
  open (FH, "</wormsrv2/autoace/BLAT/EST.dat") or die "EST.dat : $!\n";
  undef $/;
  my $data = <FH>;
  eval $data;
  die if $@;
  $/ = "\n";
  close(FH);

  open (OUT_EST, ">$dir/elegans_ESTs")     if ($est);
  open (OUT_OST, ">$dir/elegans_OSTs")     if ($ost);
  open (OUT_ACE, ">$dir/elegans_ESTs.ace") if ($ace);
  open (SEQUENCES, "$getz -sf fasta -f \"id acc des seq sv\" \'([embl-org:caenorhabditis elegans] \& [embl-div:est])\' |") ;

  while (<SEQUENCES>) {

    unless (/^AC\s+/ || /^DE\s+/ || /^>/ || /^ID\s+/ || /^SV\s+/) {
      print OUT_EST  if ($ost_seq == 0);
      print OUT_ACE  if (($ost_seq == 0) && ($ace));
      print OUT_OST  if ($ost_seq == 1);      
    }
    
    ($acc = 1) if (/^AC\s+(\S+);/);
    ($id = $1) if (/^ID\s+(\S+)/);
    ($sv = $1) if (/^SV\s+\S+\.(\d+)/);

    if (/^DE\s+(.+)/)   {
      $def = $def." ".$1;
      
      if ($def =~ /^ OST/) {
	$ost_seq = 1;
      }
      else {
	$ost_seq = 0;
      }
    }
    if (/^>/) {
      if ($ost_seq == 0) {
	print OUT_EST ">$acc $id $def\n";
	
	if ($ace) {
	  if (exists $EST_name{$acc}) {
	    print OUT_ACE "\nSequence : \"$EST_name{$acc}\"\nDatabase EMBL $id $acc $sv\n" if ($ace);
	  }
	  else {	
	    print OUT_ACE "\nSequence : \"$acc\"\nDatabase EMBL $id $acc $sv\n" if ($ace);
	  }
	  print OUT_ACE   "Species \"Caenorhabditis elegans\"\n" if ($ace);
	  print OUT_ACE   "Title \"$def\"\nMethod EST_elegans\n" if ($ace);
	  
	  if (exists $EST_name{$acc}) {
	    print OUT_ACE "\nDNA \"$EST_name{$acc}\"\n" if ($ace);
	  }
	  else {	
	    print OUT_ACE   "\nDNA \"$acc\"\n" if ($ace);
	  }
	}
	
	$def = "";
	$id  = "";
	$sv  = "";
	$acc = "";
      } 
      else {
	print OUT_OST ">$acc $id $def\n";
	$def = "";
	$id  = "";
	$sv  = "";
	$acc = "";
      }
    }
  }
  close (SEQUENCES);
  close (OUT_EST) if ($est);
  close (OUT_OST) if ($ost);
  close (OUT_ACE) if ($ace);

  # pressdb fasta database
  system ("/usr/local/pubseq/bin/pressdb $dir/elegans_ESTs > /dev/null") if ($blastdb && $est);
  system ("/usr/local/pubseq/bin/pressdb $dir/elegans_OSTs > /dev/null") if ($blastdb && $ost);

}



#################################################################
#
# C. elegans mRNAs
#
#################################################################

sub make_mrnas{

#  my $ftp_file   = "/nfs/disk69/ftp/pub/wormbase/sequences/ESTS/C.elegans_nematode_mRNAs";
    
  # open filehandles for output files 
  open (OUT_MRNA, ">$dir/elegans_mRNAs");
  open (OUT_ACE,  ">$dir/elegans_mRNAs.ace") if ($ace);

  open (SEQUENCES, "$getz -f \"acc\" \'([emblnew-org:caenorhabditis elegans*] \& [emblnew-mol:mrna]  \! [emblnew-div:EST]) | ([embl-org:caenorhabditis elegans*] \& [embl-mol:mrna] \! [embl-div:EST])\' |") ;

  while (<SEQUENCES>) {
    chomp;
    ($acc) = (/^AC\s+(\S+)\;/);
    print "Parsing EMBL accession: '$acc'\n" if ($verbose);
    next if ($acc eq "");
    $def = "";
    # pfetch each sequence individually
    open (LOOK, "/usr/local/pubseq/bin/pfetch -F $acc |");
    while (<LOOK>) {
      print if ($verbose);
      
      if (/^\s/) {
	s/\s+//g;
	s/\d+//g;
	print OUT_MRNA "$_\n";
	print OUT_ACE  "$_\n" if ($ace);
      }
      
      if (/^ID\s+(\S+)/)                         {$id  = $1;}
      if (/^SV\s+\S+\.(\d+)/)                    {$sv  = $1;}
      if (/^DE\s+(.+)/)                          {$def = $def." ".$1;}
      if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {$protid=$1; $protver=$2;}
      if (/^SQ/) {
	print OUT_MRNA ">$acc $id $def\n";
	if ($ace){
	  print OUT_ACE   "\nSequence : \"$acc\"\nDatabase EMBL $id $acc $sv\n";
	  print OUT_ACE   "Protein_id $acc $protid $protver\n";
	  print OUT_ACE   "Species \"Caenorhabditis elegans\"\n";
	  print OUT_ACE   "Title \"$def\"\nMethod NDB\n";
	  print OUT_ACE   "\nDNA \"$acc\"\n";
	}
	# reset vars
	$def = ""; $id = ""; $acc = ""; $sv = ""; $protid = ""; $protver ="";
      }
    }
  }
  # close filehandles
  close (SEQUENCES);  
  close(OUT_MRNA);
  close(OUT_ACE) if ($ace);
  

  # pressdb fasta database
  system ("/usr/local/pubseq/bin/pressdb $dir/elegans_mRNAs > /dev/null") if ($blastdb);
  
# clean FTP copy and transfer new ones
#  system ("/bin/rm -f ${ftp_file}.gz");
#  system ("cp $fasta_file $ftp_file");
#  system ("/bin/gzip $ftp_file");


}


#################################################################
#
# Other nematode species ESTs
#
#################################################################

sub make_nematode_ests{

#  my $ftp_file   = "/nfs/disk69/ftp/pub/wormbase/sequences/ESTS/non_C.elegans_nematode_ESTs";
  # remove old files
#  system ("/bin/rm -f $fasta_file");
#  system ("/bin/rm -f $ace_file");
  
  # open filehandles for output files 
  open (OUT_NEM, ">$dir/other_nematode_ESTs");
  open (OUT_ACE, ">$dir/other_nematode_ESTs.ace") if ($ace);
  
  my @species_list = (
		      'Ancylostoma caninum',
		      'Ancylostoma ceylanicum',
		      'Angiostrongylus cantonensis',
		      'Ascaris lumbricoides',
		      'Ascaris suum',
		      'Brugia malayi',
		      'Brugia pahangi',
		      'Caenorhabditis briggsae',
		      'Caenorhabditis marpasi',
		      'Caenorhabditis vulgarensis',
		      'Dictyocaulus viviparus',
		      'Dirofilaria immitis',
		      'Globodera mexicana',
		      'Globodera pallida',
		      'Globodera rostochiensis',
		      'Haemonchus contortus',
		      'Heterodera glycines',
		      'Heterodera schachtii',
		      'Litomosoides sigmodontis',
		      'Loa loa',
		      'Meloidogyne arenaria',
		      'Meloidogyne incognita',
		      'Meloidogyne hapla',
		      'Meloidogyne javanica',
		      'Necator americanus',
		      'Nippostrongylus brasiliensis',
		      'Oesophagostomum dentatum',
		      'Onchocerca ochengi',
		      'Onchocerca volvulus',
		      'Ostertagia ostertagi',
		      'Parastrongyloides trichosuri',
		      'Pratylenchus penetrans',	
		      'Pristionchus pacificus',
		      'Steinernema feltiae',
		      'Strongyloides ratti',
		      'Strongyloides stercoralis',
		      'Teladorsagia circumcincta',
		      'Toxocara canis',
		      'Trichinella spiralis',
		      'Trichostrongylus vitrinus',
		      'Trichuris muris',
		      'Wuchereria bancrofti',
		      'Zeldia punctata'
		      );

  foreach my $species (@species_list) {

    my $species_count = 0;
    $def = "";
    open (SEQUENCES, "$getz -e \'([embl-org:$species*]) \& ([embl-key:est])' |") ;
    while (<SEQUENCES>) {
      chomp;
      
      if (/^\s/) {
	s/\s+//g;
	s/\d+//g;
	print OUT_NEM "$_\n";
	print OUT_ACE "$_\n" if ($ace);
      }
      
      if (/^ID\s+(\S+)/)                         {$id  = $1;}
      if (/^SV\s+(\S+)\.(\d+)/)                  {$acc = $1; $sv = $2;}
      if (/^DE\s+(.+)/)                          {$def = $def." ".$1;}
      if (/^SQ/) {
	print OUT_NEM ">$acc $id $def\n";
	if ($ace){
	  print OUT_ACE "\nSequence : \"$acc\"\nDatabase EMBL $id $acc $sv\n";
	  print OUT_ACE "Species \"$species\"\n";
	  print OUT_ACE "Title \"$def\"\nMethod EST_nematode\n";
	  print OUT_ACE "\nDNA \"$acc\"\n";
	}
	$species_count++;
	# reset vars
	$def = ""; $id = ""; $acc = ""; $sv = ""; $protid = ""; $protver ="";
      }
    }    
    
    close SEQUENCES;

    print "// Parsed $species_count entries for $species\n" if ($verbose);

  }
    

  # close filehandles
  close(OUT_NEM);
  close(OUT_ACE) if ($ace);


  # pressdb fasta database
  system ("/usr/local/pubseq/bin/pressdb $dir/other_nematode_ESTs > /dev/null") if ($blastdb);

# clean FTP copy and transfer new ones
#system ("/bin/rm -f ${ftp_file}.gz");
#system ("cp $fasta_file $ftp_file");
#system ("/bin/gzip $ftp_file");

}

#################################################################
#
# Non-wormbase C. elegans CDS in EMBL
#
#################################################################

sub make_embl_cds{

  open (OUT_EMBL, ">$dir/elegans_embl_cds");
  open (OUT_ACE, ">$dir/elegans_embl_cds.ace") if ($ace);

  open (SEQUENCES, "$getz -e \'(([emblnew-Division:inv] & [emblnew-Molecule:genomic dna*] & [emblnew-org:Caenorhabditis elegans*] ! [emblnew-Keywords:HTG*]) & ([emblnew-FtKey:cds] > emblnew)) | (([emblrelease-Division:inv] & [emblrelease-Molecule:genomic dna*] & [emblrelease-org:Caenorhabditis elegans*] ! [emblrelease-Keywords:HTG*]) & ([emblrelease-FtKey:cds] >emblrelease))\' |");


  $def = "";
  while (<SEQUENCES>) {
    chomp;
    if (/^\s/) {                   # matches only the sequence lines
      s/\s+//g;                  # remove white space
      s/\d+//g;                  # remove numerals
      print OUT_ACE "$_\n" if ($ace);   # print to output
    }
    
    if (/^ID\s+(\S+)/)                         {$id  = $1;}
    if (/^SV\s+(\S+)\.(\d+)/)                  {$acc = $1; $sv  = $2;}
    if (/^DE\s+(.+)/)                          {$def = $def." ".$1;}
    if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {$protid=$1; $protver=$2;}
    if (/^SQ/) {
      print "// Parsed accession $acc\n" if ($verbose);
      
      open (CDS, "CDS.pl $acc |");   # Ugly. this needs rationalising dl
      while  (<CDS>) {
	print OUT_EMBL $_;
      }
      close(CDS);
      
      # print out the acefile version
      if ($ace){
	print OUT_ACE   "\nSequence : \"$acc\"\nDatabase EMBL $id $acc $sv\n";
	print OUT_ACE   "Species \"Caenorhabditis elegans\"\n";
	print OUT_ACE   "Title \"$def\"\nMethod NDB\n";
	print OUT_ACE  "\nDNA \"$acc\"\n";
      }
      # reset vars
      $def = ""; $id = ""; $acc = ""; $sv = ""; $protid = ""; $protver ="";
    }
  }
  
  close(SEQUENCES);
  close(OUT_EMBL);
  close(OUT_ACE) if ($ace);

  # pressdb fasta database
  system ("/usr/local/pubseq/bin/pressdb $dir/elegans_embl_cds > /dev/null") if ($blastdb);

}


#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub usage {
  my $error = shift;
  if ($error == 0) {
    # Normal help menu
    exec ('perldoc',$0);
  }
}
###############################################################################






__END__

=pod

=head2   NAME - make_C.elegans_nematode_ESTs

=head1 USAGE

=over 4

=item make_C.elegans_nematode_ESTs

=back

make_C.elegans_nematode_ESTs runs an SRS query on the EMBL database
to return all C.elegans sequences from the EST division. This set
is then split using a simple perl regex into the OST set (Marc Vidal's
Orfeome project 1.1) and everything else. These files are used for the
BLAT analysis as part of each WormBase build.

=back

=cut


=pod

=head2   NAME - make_C.elegans_nematode_mRNAs

=head1 USAGE

=over 4

=item make_C.elegans_nematode_mRNAs [-options]

=back

make_C.elegans_nematode_mRNAs is a wrapper around a getz
call to make a fasta file for all C.elegans mRNAs in the
EMBL database. The script produces two files:

fasta file for BLAT similarity mapping

ace file for loading data into ACEDB

make_C.elegans_nematode_mRNAs MANDATORY arguments:

=over 4

=item none, 

=back

make_C.elegans_nematode_mRNAs OPTIONAL arguments:

=over 4

=item -debug, verbose report (e-mails dl1 only)

=item -help, this help page

=back

make_C.elegans_nematode_mRNAs DEFAULT behaviour:

make_C.elegans_nematode_mRNAs will write the two 
data files into /nfs/disk100/wormpub/analysis/ESTs.

=head1 REQUIREMENTS

=back

=item none

=back

This script is called by the following scripts:

=over 4

=item none

=back

=head1 AUTHOR

=back

=item Dan Lawson (dl1@sanger.ac.uk)

=back

=cut

=pod

=head2   NAME - make_non_C.elegans_nematode_ESTs

=head1 USAGE

=over 4

=item make_non_C.elegans_nematode_ESTs [-options]

=back

make_non_C.elegans_nematode_ESTs is a wrapper around a getz
call to make a fasta file for all non-elegans nematode ESTs in 
the EMBL database. The script produces two files:

fasta file for BLAT similarity mapping

ace file for loading data into ACEDB

make_non_C.elegans_nematode_ESTs MANDATORY arguments:

=over 4

=item none, 

=back

make_non_C.elegans_nematode_ESTs OPTIONAL arguments:

=over 4

=item -debug, verbose report (e-mails dl1 only)

=item -help, this help page

=back

make_non_C.elegans_nematode_ESTs DEFAULT behaviour:

make_non_C.elegans_nematode_ESTs will write the two 
data files into /nfs/disk100/wormpub/analysis/ESTs.

=head1 REQUIREMENTS

=back

=item none

=back

This script is called by the following scripts:

=over 4

=item none

=back

=head1 AUTHOR

=back

=item Dan Lawson (dl1@sanger.ac.uk)

=back

=cut




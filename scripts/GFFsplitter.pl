#!/usr/local/bin/perl5.8.0 -w
#
# GFFsplitter.pl
# 
# by Dan Lawson
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-03-25 11:12:13 $
#
# Usage GFFsplitter.pl [-options]


#################################################################################
# variables                                                                     #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Ace;


##################################################
# Script variables and command-line options      #
##################################################
my $maintainers = "All";
my $WS_version = &get_wormbase_version_name;
our $lockdir = "/wormsrv2/autoace/logs/";


my $help;      # Help/Usage page
my $archive;   # archive GFF_splits directory into a WSxx directory
my $debug;     # debug
my $verbose;   # verbose mode
our $log;

GetOptions (
	    "help"      => \$help,
	    "archive"   => \$archive,
	    "debug:s"   => \$debug
	    );

# help 
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;


##############################
# Paths etc                  #
##############################

my $datadir = "/wormsrv2/autoace/GFF_SPLITS";

# create GFF_SPLITS subdirectory if it doesn't already exist
if (! -e "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS"){
  system("mkdir /wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS") && die "Couldn't create directory\n";
}

# check to see if full chromosome gff dump files exist
if (! -e "/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_I.gff"){
  print "chromosome_dump.pl may have failed/n";
}



##########################################################
# Archive the GFF_splits directory into a WSxx directory
##########################################################

# runs only if -a is specified
if($archive){
  print "Renaming $datadir/GFF_SPLITS to $datadir/$WS_version\n";
  system("mv $datadir/GFF_SPLITS ${datadir}/${WS_version}") && die "Couldn't rename directory\n";
  exit(0);
}

#################################################################################
# Main Loop                                                                     #
#################################################################################


# prepare array of file names and sort names
our @files = (
	      'CHROMOSOME_I',
	      'CHROMOSOME_II',
	      'CHROMOSOME_III',
	      'CHROMOSOME_IV',
	      'CHROMOSOME_V',
	      'CHROMOSOME_X',
	      );

our @gff_files = sort @files; 
undef @files; 

our @gff_classes;

# prepare array of data types to dump
READARRAY:   while (<DATA>) {
    chomp $_;
    last READARRAY if $_ =~ /END/;
    push (@gff_classes,$_);
}


############################################################
# loop through each GFF file                               #
############################################################

my %GFF;
my $file;

foreach $file (@gff_files) {
  undef(%GFF);
  next if ($file eq "");
    
  print LOG "# File $file\n";
  print     "\n# File $file\n" if ($verbose);
  
  my $line_count = 0;
  
  #########################################################
  # open the gff file                                     #
  #########################################################
  
  my $running_total = 0;
  my ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name);
  my @header = "";
  
  open (GFF, "</wormsrv2/autoace/CHROMOSOMES/$file.gff")  || die "Cannot open /wormsrv2/autoace/CHROMOSOMES/$file.gff\n";
  while (<GFF>) {
    chomp;
    $line_count++;
    
    print "." if ((($line_count % 5000) == 0) && $verbose);
    
    #skip header lines of file
    if (/^\#/) {push (@header,$_); next;}
    
    ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name) = split /\t/;
    
    # Clone path
    if    ( ($source eq "Genomic_canonical")    && ($feature eq "region"))      {push (@{$GFF{$file}{clone_path}},$_);}
    # Genes (CDSs)  
    elsif ((($source eq "curated")              && ($feature eq "CDS")))        {push (@{$GFF{$file}{genes}},$_);}
    # Pseudogenes
    elsif ( ($source eq "Pseudogene")           && ($feature eq "Pseudogene"))  {push (@{$GFF{$file}{pseudogenes}},$_);}
    # RNA genes 
    elsif ($feature =~ m/_primary_transcript/)                                  {push (@{$GFF{$file}{rna}},$_);}

    # coding_exon (used to be CDS)
    elsif (($source eq "curated")               && ($feature eq "coding_exon")) {push (@{$GFF{$file}{coding_exon}},$_);}
    # Exon      
    elsif (($source eq "curated")               && ($feature eq "exon"))        {push (@{$GFF{$file}{exon}},$_);}
    elsif (($source eq "Pseudogene")            && ($feature eq "exon"))        {push (@{$GFF{$file}{exon_pseudogene}},$_);}
    elsif (($source eq "Non_coding_transcript") && ($feature eq "exon"))        {push (@{$GFF{$file}{exon_noncoding}},$_);}
    elsif (($source eq "tRNAscan-SE-1.23")      && ($feature eq "exon"))        {push (@{$GFF{$file}{exon_tRNA}},$_);}


    # Intron    
    elsif (($source eq "curated")               && ($feature eq "intron"))      {push (@{$GFF{$file}{intron}},$_);}
    # all other introns
    elsif ($feature eq "intron")                                                {push (@{$GFF{$file}{intron_all}},$_);}

    # Repeats
    elsif ($source eq "RepeatMasker")                                              {push (@{$GFF{$file}{repeats}},$_);}

    # Assembly tags
    elsif ($source eq "assembly_tag")                                           {push (@{$GFF{$file}{assembly_tags}},$_);}
    # SL1/SL2 acceptor sites
    elsif ($feature =~ m/_acceptor_site/)                                       {push (@{$GFF{$file}{ts_site}},$_);}
    # Oligo mapping
    elsif ($feature eq "OLIGO")                                                 {push (@{$GFF{$file}{oligos}},$_);}
    # RNAi
    elsif ($feature eq "RNAi_reagent")                         {push (@{$GFF{$file}{RNAi}},$_);}
    # PCR_products
    elsif ($feature eq "PCR_product")                 {push (@{$GFF{$file}{PCR_products}},$_);}
    # Alleles
    elsif (($source eq "Allele") || ($source eq "Mos_insertion_allele")) {push (@{$GFF{$file}{allele}},$_);}
    # operons
    elsif ($source eq "operon")                       {push (@{$GFF{$file}{operon}},$_);}
    # Oligo_set
    elsif (($source eq "Oligo_set") && ($feature eq "reagent")) {push (@{$GFF{$file}{Oligo_set}},$_);}
    # Clone ends
    elsif ((/Clone_left_end/) || (/Clone_right_end/)) {push (@{$GFF{$file}{clone_ends}},$_);}
    # cDNA for RNAi
    elsif (/cDNA_for_RNAi/)                           {push (@{$GFF{$file}{cDNA_for_RNAi}},$_);}
    # BLAT_EST
    elsif (/BLAT_EST_BEST/)                           {push (@{$GFF{$file}{BLAT_EST_BEST}},$_);
						       push (@{$GFF{$file}{BLAT_TRANSCRIPT_BEST}},$_);}
    elsif (/BLAT_EST_OTHER/)                          {push (@{$GFF{$file}{BLAT_EST_OTHER}},$_);}
    # BLAT_OST
    elsif (/BLAT_OST_BEST/)                           {push (@{$GFF{$file}{BLAT_OST_BEST}},$_);
						       push (@{$GFF{$file}{BLAT_TRANSCRIPT_BEST}},$_);}
    elsif (/BLAT_OST_OTHER/)                          {push (@{$GFF{$file}{BLAT_OST_OTHER}},$_);}
    # BLAT_mRNA
    elsif (/BLAT_mRNA_BEST/)                          {push (@{$GFF{$file}{BLAT_mRNA_BEST}},$_);}
    elsif (/BLAT_mRNA_OTHER/)                         {push (@{$GFF{$file}{BLAT_mRNA_OTHER}},$_);}
    # BLAT_EMBL
    elsif (/BLAT_EMBL_BEST/)                          {push (@{$GFF{$file}{BLAT_EMBL_BEST}},$_);}
    elsif (/BLAT_EMBL_OTHER/)                         {push (@{$GFF{$file}{BLAT_EMBL_OTHER}},$_);}
    # BLAT_NEMATODE
    elsif (/BLAT_NEMATODE/)                           {push (@{$GFF{$file}{BLAT_NEMATODE}},$_);}
    # BLAT_TC1_BEST
    elsif (/BLAT_TC1_BEST/)                           {push (@{$GFF{$file}{BLAT_TC1_BEST}},$_);}
    # Expr_profile
    elsif (/Expr_profile/)                            {push (@{$GFF{$file}{Expr_profile}},$_);}
    # UTR         
    elsif (/UTR/)                                     {push (@{$GFF{$file}{UTR}},$_);}
    # Protein similarities
    elsif (/wublastx/)                                {push (@{$GFF{$file}{BLASTX}},$_);}
    # C. briggsae
    elsif (/:waba/)                                   {push (@{$GFF{$file}{WABA_BRIGGSAE}},$_);}

    # REST OF LINES
    else                                              {push (@{$GFF{$file}{rest}},$_); $running_total--;}
    
    # increment no of lines sorted to bin
    $running_total++; 
    
  }
  close(GFF);
  ########################################
  # write some output                    #
  ########################################
  
  foreach my $tag (@gff_classes) {
    print "# $file $tag\n" if ($debug);
    
    open (OUT, ">$datadir/GFF_SPLITS/$file.$tag.gff") or die "Can't open file\n";
    
    # gff header lines
    foreach my $line (@header) {
      next if ($line eq "");
      print OUT "$line\n";
    }
    
    # gff split data lines
    foreach my $line (@{$GFF{$file}{$tag}}) {
      next if ($line eq "");
      print OUT "$line\n";
    }
    close OUT;
  }
  
  
  #########################################
  # Alter clone file to include accession #
  #########################################
  
  # GFF clone_path with EMBL accessions and sequence versions
  my $input_file = "$datadir/GFF_SPLITS/$file.clone_path.gff";
  my $output_file = "$datadir/GFF_SPLITS/$file.clone_acc.gff";
  &GFF_clones_with_accessions("$input_file", "$output_file");
  
  
  # GFF genes with wormpep CE accessions
  # Shouldn't do this unless Wormpep has been made else no Corresponding_protein tags in database
  if(-e "$lockdir/D1:Build_wormpep_final"){
    $input_file = "$datadir/GFF_SPLITS/$file.genes.gff";
    $output_file = "$datadir/GFF_SPLITS/$file.genes_acc.gff";
    &GFF_genes_with_accessions("$input_file", "$output_file");
    system ("mv -f $output_file $input_file");
  }
  
  # GFF UTRs with CDS names
  # Shouldn't attempt to do this if UTR data has not been generated
  if(-e "$lockdir/B10:Generate_UTR_data" ) {
    my $utr_file = "$datadir/GFF_SPLITS/$file.UTR.gff";
    my $utr_cds_file = "$datadir/GFF_SPLITS/$file.UTR_CDS.gff";
    &GFF_with_UTR("$utr_file","$utr_cds_file");
  }
  
}

# Tidy up
close (LOG);

&mail_maintainer("GFFsplitter.pl finished",$maintainers,$log);

exit(0);






###############################
# subroutines                 #
###############################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}

##########################################


sub file_size {
    my $file = shift;
    my $file_length;
    open (FILE_SIZE, "/bin/cat $file | wc -l |");
    while (<FILE_SIZE>) {
	chomp;
	s/\s//g;
	$file_length = $_;
    }
    close FILE_SIZE;
    print LOG "Counting No. of lines in file '$file' : $file_length\n" if ($debug);
    return $file_length;
}


##########################################
sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  elsif ($error eq "Debug") {
    # No debug person named
    print "You haven't supplied your name\nI won't run in debug mode
         until I know who you are\n";
    exit (0);
  }
}

################################################
#
# Post-processing GFF routines
#
################################################

sub GFF_with_UTR{
  my $utr = shift;
  my $utr_cds = shift;
  open( UTR, "<$utr" );
  open( UTR_CDS, ">$utr_cds");
  while (<UTR>) {
    if (/^\#/){
      print UTR_CDS $_;
    }
    else{
      (/_UTR:(\S+)\"/);
      print UTR_CDS "$`" . "_UTR:$1\" CDS=\"$1\"\n";
    }
  }
  close UTR;
  close UTR_CDS;
  
  system ("mv -f $utr_cds $utr");
}

########################################

sub GFF_clones_with_accessions{
  my $infile   = shift;
  my $outfile = shift;
  my $wormdb = "/wormsrv2/autoace";

  my $db = Ace->connect(-path=>$wormdb) || do { print "Connection failure to $wormdb: ",Ace->error; die();};

  open (GFF, "<$infile")  || die "Can't open GFF file: $infile\n\n";
  open (OUT, ">$outfile") || die "Can't write to $outfile\n";
  while (<GFF>) {

    next if (/^\#/);
    chomp;

    my @gff = split (/\t/,$_);
    my ($gen_can) = $gff[8] =~ /Sequence \"(\S+)\"/; 
    my $obj = $db->fetch(Sequence=>$gen_can);
    if (!defined ($obj)) {
      print "Could not fetch sequence '$gen_can'\n";
      next;
    }

    # Grab accession and sequence version info from the database
    my $acc;
    my $sv;

    if($obj->at('DB_info.Database.EMBL.NDB_AC')){
      ($acc) = $obj->at('DB_info.Database.EMBL.NDB_AC');
      ($sv) = $obj->at('DB_info.Database.EMBL.NDB_SV');
    }
    elsif($obj->at('DB_info.Database.GenBank.NDB_AC')){
      ($acc) = $obj->at('DB_info.Database.GenBank.NDB_AC');
      ($sv) = $obj->at('DB_info.Database.GenBank.NDB_SV');
    }

    # now just want the numerical suffix of sequence version field
    $sv =~ s/.*\.//;

     print OUT "$_ acc=$acc ver=$sv\n";
    
    $obj->DESTROY();
    
  }
  close(GFF);
  close(OUT);
  $db->close;
 
}

########################################

sub GFF_genes_with_accessions{
  my $infile  = shift;
  my $outfile = shift;
  my $wormdb = "/wormsrv2/autoace";

  my $db = Ace->connect(-path=>$wormdb) || do { print "Connection failure to $wormdb: ",Ace->error; die();};

  open (GFF, "<$infile") || die "Can't open GFF file: $infile\n\n";
  open (OUT, ">$outfile") || die "Can't write to $outfile\n";
  while (<GFF>) {

    next if (/^\#/);    
    chomp;

    my @gff = split (/\t/,$_);
    next if ($gff[1] ne "curated");
    (my $gen_can) = $gff[8] =~ /CDS \"(\S+)\"/; 
    my $obj = $db->fetch(CDS=>$gen_can);
    if (!defined ($obj)) {
      print "Could not fetch CDS '$gen_can'\n";
      next;
    }

    my $acc = $obj->Corresponding_protein(1);
    $acc =~ s/WP\://g;

    if ($acc ne "") {
      $acc =~ s/WP\://g;
      print OUT "$_ wp_acc=$acc\n";
    }
    else {  # could be a tRNA
      print OUT "$_\n";
    }
    undef ($acc); 

    $obj->DESTROY();

  }
  close (GFF);
  close (OUT);
  $db->close;

 
}

__DATA__
clone_path
genes
pseudogenes
rna
coding_exon
exon
exon_tRNA
exon_pseudogene
exon_noncoding
intron
intron_all
repeats
assembly_tags
ts_site
oligos
RNAi
allele
clone_ends
PCR_products
cDNA_for_RNAi
BLAT_EST_BEST
BLAT_EST_OTHER
BLAT_OST_BEST
BLAT_OST_OTHER
BLAT_TRANSCRIPT_BEST
BLAT_mRNA_BEST
BLAT_mRNA_OTHER
BLAT_EMBL_BEST
BLAT_EMBL_OTHER
BLAT_NEMATODE
BLAT_TC1_BEST
Expr_profile
UTR
BLASTX
WABA_BRIGGSAE
operon
Oligo_set
rest
__END__



=pod

=head2 NAME - GFFsplitter.pl

=back 

=head1 USAGE

=over 4

=item GFFsplitter.pl <options>

=back

This script splits the large GFF files produced during the build process into
smaller files based on a named set of database classes to be split into.
Output written to /wormsrv2/autoace/GFF_SPLITS/WSxx

=over 4

=item MANDATORY arguments: 

None.

=back

=over 4

=item OPTIONAL arguments: -help, this help page.

= item -debug <user>, only email report/logs to <user>

= item -archive, archives (gzips) older versions of GFF_SPLITS directory

= item -verbose, turn on extra output to screen to help track progress


=back


=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut

#!/usr/local/bin/perl5.8.0 -w
#
# GFFsplitter.pl
# 
# by Dan Lawson
#
# Last updated on: $Date: 2006-02-17 11:32:47 $
#
# Usage GFFsplitter.pl [-options]


#################################################################################
# variables                                                                     #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use File::Path;
use Getopt::Long;
use Ace;
use Carp;
use Log_files;
use Storable;

##################################################
# Script variables and command-line options      #
##################################################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $archive;   # archive GFF_splits directory into a WSxx directory
my $chrom;     # single chromosome mode
my $splitdir;
my $gffdir;

GetOptions (
	    "help"      => \$help,
	    "archive"   => \$archive,
	    "debug:s"   => \$debug,
	    "chrom:s"   => \$chrom,
	    "gffdir:s"  => \$gffdir,
	    "splitdir:s"=> \$splitdir,
	    "test"      => \$test,
	    "verbose"   => \$verbose,
	    "store:s"     => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
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

##############################
# Paths etc                  #
##############################

my $database        = $wormbase->autoace;
$splitdir        = $wormbase->gff_splits unless $splitdir;
$gffdir          = $wormbase->chromosomes unless $gffdir;
my $lockdir      = $wormbase->logs;

my $WS_version = $wormbase->get_wormbase_version_name;

$log->write_to("splitting GFF from $gffdir to $gffdir\n");

my @files;
# prepare array of file names and sort names


if (defined($chrom)){
    unless (grep { $chrom eq $_ } ('I','II','III','IV','V','X','MtDNA')) {
	die "ERROR: $chrom is an incorrect chromosome number, please use I, II, III etc.\n";
    }
    @files = (
		  "CHROMOSOME_${chrom}"
		  );
}
else {
    @files = (
		  'CHROMOSOME_I',
		  'CHROMOSOME_II',
		  'CHROMOSOME_III',
		  'CHROMOSOME_IV',
		  'CHROMOSOME_V',
		  'CHROMOSOME_X',
		  'CHROMOSOME_MtDNA'
		  );
}

our @gff_files = sort @files; 
undef @files; 

our @gff_classes;


# create GFF_SPLITS subdirectory if it doesn't already exist
if (! -e $splitdir){
  mkpath($splitdir) or die "Couldn't create directory $splitdir\n";
}

##########################################################
# Archive the GFF_splits directory into a WSxx directory
##########################################################

# runs only if -a is specified
if($archive){
  my $archive_split = $wormbase->basedir . "$WS_version/GFF_SPLITS";
  # create archive GFF_SPLITS subdirectory if it doesn't already exist
  if (! -e $archive_split){
    mkpath($archive_split) or die "Couldn't create directory $archive_split\n";
  }
  print "Renaming $splitdir to $archive_split\n" if ($verbose);
  system("mv $splitdir $archive_split") && die "Couldn't rename directory\n";
  exit(0);
}

# check to see if full chromosome gff dump files exist
foreach my $file (@gff_files) {
    unless (-e "$gffdir/$file.gff") {
      die "$gffdir/$file.gff is missing\n";
    }
    if (-e -z "$gffdir/$file.gff") {
      die "$gffdir/$file.gff is zero length\n";
    }
}


#################################################################################
# Main Loop                                                                     #
#################################################################################

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
    
  $log->write_to("# File $file\n");
  print     "\n# File $file\n" if ($verbose);
  
  my $line_count = 0;
  
  #########################################################
  # open the gff file                                     #
  #########################################################
  
  my $running_total = 0;
  my ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name);
  my @header = "";
  
  open (GFF, "<$gffdir/$file.gff")  || die "Cannot open $gffdir/$file.gff\n";
  while (<GFF>) {
    chomp;
    $line_count++;
    
    print "." if ((($line_count % 5000) == 0) && $verbose);
    
    #skip header lines of file
    if (/^\#/) {push (@header,$_); next;}
    
    ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name) = split /\t/;
    
    # Clone path
    if    ( ($source eq "Genomic_canonical")    && ($feature eq "region"))      {push (@{$GFF{$file}{clone_path}},$_);}
    # Genes (WBGene)
    elsif ($source eq "gene")                                                   {push (@{$GFF{$file}{WBGene}},$_);}
    # Genes (CDSs)  
    elsif ((($source eq "curated")              && ($feature eq "CDS")))        {push (@{$GFF{$file}{CDS}},$_);
									         push (@{$GFF{$file}{worm_genes}},$_);}
    # Pseudogenes
    elsif ( ($source eq "Pseudogene")           && ($feature eq "Pseudogene"))  {push (@{$GFF{$file}{pseudogenes}},$_);
									         push (@{$GFF{$file}{worm_genes}},$_);}
    # RNA genes 
    elsif ($feature =~ m/_primary_transcript/)                                  {push (@{$GFF{$file}{rna}},$_);
									         push (@{$GFF{$file}{worm_genes}},$_);}
    # Transposon
    elsif( ($source eq "Transposon")         or ($source eq "Transposon_CDS") ) {push (@{$GFF{$file}{transposon}},$_);}

    # Coding_transcripts
    elsif ($source eq "Coding_transcript")                                      {push (@{$GFF{$file}{Coding_transcript}},$_);}
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

    # Genefinder predictions
    elsif ($source eq "Genefinder")                                             {push (@{$GFF{$file}{Genefinder}},$_);}
    # History predictions
    elsif ($source eq "history")                                                {push (@{$GFF{$file}{history}},$_);}
    # Twinscan predictions
    elsif ($source eq "twinscan")                                               {push (@{$GFF{$file}{twinscan}},$_);}
    # Repeats
    elsif ( ($source eq "RepeatMasker") || ($source eq "inverted") || ($source eq "tandem"))  {push (@{$GFF{$file}{repeats}},$_);}

    # Assembly tags
    elsif ($source eq "assembly_tag")                                           {push (@{$GFF{$file}{assembly_tags}},$_);}
    
    ## Features ##

    # SL1/SL2 acceptor sites
    elsif ( ($source eq "SL1") || ($source eq "SL2") )                          {push (@{$GFF{$file}{TSL_site}},$_);}

    # polyA signal and polyA sites
    elsif ( ($source eq "polyA_signal_sequence") || ($source eq "polyA_site") ) {push (@{$GFF{$file}{polyA}},$_);}


    # Oligo mapping
    elsif ($feature eq "oligo")                                                 {push (@{$GFF{$file}{oligos}},$_);}
    # RNAi
    elsif ($source eq "RNAi_primary")                                           {push (@{$GFF{$file}{RNAi_primary}},$_);}
    elsif ($source eq "RNAi_secondary")                                         {push (@{$GFF{$file}{RNAi_secondary}},$_);}

    # PCR_products
    elsif ($feature eq "PCR_product")                                           {push (@{$GFF{$file}{PCR_products}},$_);}
    # Alleles
    elsif (($source eq "Allele") || ($source eq "Mos_insertion_allele"))        {push (@{$GFF{$file}{allele}},$_);}
    # operons
    elsif ($source eq "operon")                                                 {push (@{$GFF{$file}{operon}},$_);}
    # Oligo_set
    elsif (($source eq "Oligo_set") && ($feature eq "reagent"))                 {push (@{$GFF{$file}{Oligo_set}},$_);}
    # Clone ends
    elsif ((/Clone_left_end/) || (/Clone_right_end/))                           {push (@{$GFF{$file}{clone_ends}},$_);}
    # cDNA for RNAi
    elsif (/cDNA_for_RNAi/)                                                     {push (@{$GFF{$file}{cDNA_for_RNAi}},$_);}
    # BLAT_EST
    elsif (/BLAT_EST_BEST/)                                                     {push (@{$GFF{$file}{BLAT_EST_BEST}},$_);
						                                 push (@{$GFF{$file}{BLAT_TRANSCRIPT_BEST}},$_);}
    elsif (/BLAT_EST_OTHER/)                                                    {push (@{$GFF{$file}{BLAT_EST_OTHER}},$_);}
    # BLAT_OST
    elsif (/BLAT_OST_BEST/)                                                     {push (@{$GFF{$file}{BLAT_OST_BEST}},$_);
						                                 push (@{$GFF{$file}{BLAT_TRANSCRIPT_BEST}},$_);}
    elsif (/BLAT_OST_OTHER/)                                                    {push (@{$GFF{$file}{BLAT_OST_OTHER}},$_);}
    # BLAT_mRNA
    elsif (/BLAT_mRNA_BEST/)                                                    {push (@{$GFF{$file}{BLAT_mRNA_BEST}},$_);
						                                 push (@{$GFF{$file}{BLAT_TRANSCRIPT_BEST}},$_);}
    elsif (/BLAT_mRNA_OTHER/)                                                   {push (@{$GFF{$file}{BLAT_mRNA_OTHER}},$_);}
    # BLAT_EMBL
    elsif (/BLAT_EMBL_BEST/)                                                    {push (@{$GFF{$file}{BLAT_EMBL_BEST}},$_);}
    elsif (/BLAT_EMBL_OTHER/)                                                   {push (@{$GFF{$file}{BLAT_EMBL_OTHER}},$_);}
    # BLAT_NEMATODE
    elsif (/BLAT_NEMATODE/)                                                     {push (@{$GFF{$file}{BLAT_NEMATODE}},$_);}
    # BLAT_NEMBASE
    elsif (/BLAT_NEMBASE/)                                                      {push (@{$GFF{$file}{BLAT_NEMBASE}},$_);}
    # BLAT_TC1_BEST
    elsif (/BLAT_TC1_BEST/)                                                     {push (@{$GFF{$file}{BLAT_TC1_BEST}},$_);}
    elsif (/BLAT_TC1_OTHER/)                                                    {push (@{$GFF{$file}{BLAT_TC1_OTHER}},$_);}
    # BLAT_ncRNA_BEST
    elsif (/BLAT_ncRNA_BEST/)                                                   {push (@{$GFF{$file}{BLAT_ncRNA_BEST}},$_);}
    elsif (/BLAT_ncRNA_OTHER/)                                                  {push (@{$GFF{$file}{BLAT_ncRNA_OTHER}},$_);}
    # BLAT_WASHU
    elsif (/BLAT_WASHU/)                                                        {push (@{$GFF{$file}{BLAT_WASHU}},$_);}
    # Expr_profile
    elsif (/Expr_profile/)                                                      {push (@{$GFF{$file}{Expr_profile}},$_);}
    # Protein similarities
    elsif (/wublastx/)                                                          {push (@{$GFF{$file}{BLASTX}},$_);}
    # C. briggsae
    elsif ($source =~ /waba/)                                                   {push (@{$GFF{$file}{WABA_BRIGGSAE}},$_);}

    # TEC_RED similarities
    elsif ($source eq "TEC_RED")                                                {push (@{$GFF{$file}{TEC_RED}},$_);}
    # SAGE transcript
    elsif ($source eq "SAGE_transcript")                                        {push (@{$GFF{$file}{SAGE}},$_);}



    # REST OF LINES
    else                                                                        {push (@{$GFF{$file}{rest}},$_); $running_total--;}
    
    # increment no of lines sorted to bin
    $running_total++; 
    
  }
  close(GFF);
  ########################################
  # write some output                    #
  ########################################
  
  foreach my $tag (@gff_classes) {
    print "# $file $tag\n" if ($verbose);
    open (OUT, ">$splitdir/$file.$tag.gff") or die "Can't open file\n";
    
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
  my $input_file = "$splitdir/$file.clone_path.gff";
  my $output_file = "$splitdir/$file.clone_acc.gff";
  &GFF_clones_with_accessions("$input_file", "$output_file");
  
  
  # GFF genes with wormpep CE accessions
  # Shouldn't do this unless Wormpep has been made else no Corresponding_protein tags in database
  if(-e "$lockdir/D1:Build_wormpep_final"){
    $input_file = "$splitdir/$file.CDS.gff";
    $output_file = "$splitdir/$file.CDS_acc.gff";
    &GFF_CDS_with_accessions("$input_file", "$output_file");
    system ("mv -f $output_file $input_file");
  }
  
}

# Tidy up
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


###############################
# subroutines                 #
###############################

##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  elsif ($error eq "No GFF file") {
      # No GFF file to work from
      print "One (or more) GFF files are absent from $gffdir\n\n";
      exit(0);
  }
  elsif ($error eq "Zero length GFF file") {
      # Zero length GFF file
      print "One (or more) GFF files are zero length. The GFF dump may not have worked\n\n";
      exit(0);
  }
  elsif ($error eq "Debug") {
    # No debug person named
    print "You haven't supplied your name\nI won't run in debug mode until I know who you are\n\n";
    exit (0);
  }
}

################################################
#
# Post-processing GFF routines
#
################################################

########################################

sub GFF_clones_with_accessions{
  my $infile   = shift;
  my $outfile = shift;

  my $db = Ace->connect(-path=>$database) || do { print "Connection failure to $database: ",Ace->error; die();};

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

sub GFF_CDS_with_accessions{
  my $infile  = shift;
  my $outfile = shift;

  my $db = Ace->connect(-path=>$database) || do { print "Connection failure to $database: ",Ace->error; die();};

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
CDS
WBGene
pseudogenes
transposon
rna
worm_genes
Coding_transcript
coding_exon
exon
exon_tRNA
exon_pseudogene
exon_noncoding
intron
intron_all
Genefinder
history
repeats
assembly_tags
TSL_site
polyA
oligos
RNAi_primary
RNAi_secondary
TEC_RED
SAGE
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
BLAT_NEMBASE
BLAT_TC1_BEST
BLAT_TC1_OTHER
BLAT_ncRNA_BEST
BLAT_ncRNA_OTHER
BLAT_WASHU
Expr_profile
BLASTX
WABA_BRIGGSAE
operon
Oligo_set
rest
twinscan
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
Output written to ~wormpub/BUILD/autoace/GFF_SPLITS/GFF_SPLITS

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

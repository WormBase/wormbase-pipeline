#!/usr/local/bin/perl5.6.1 -w
#
# GFFsplitter.pl
# 
# by Dan Lawson
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-01-17 15:28:41 $
#
# Usage GFFsplitter.pl [-options]


#################################################################################
# variables                                                                     #
#################################################################################

use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Ace;

$|=1;


##################################################
# Script variables and command-line options      #
##################################################
my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;
my $WS_version = &get_wormbase_version_name;
our $lockdir = "/wormsrv2/autoace/logs";


my $help;      # Help/Usage page
my $archive;   # archive GFF_splits directory into a WSxx directory
my $debug;     # debug

GetOptions (
	    "help"      => \$help,
	    "archive"   => \$archive,
	    "debug:s"   => \$debug
	    );

# help 
&usage("Help") if ($help);

# no debug name
print "DEBUG = \"$debug\"\n\n" if $debug;
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");



##############################
# Paths etc                  #
##############################

my $datadir = "/wormsrv2/autoace/GFF_SPLITS";

# create GFF_SPLITS subdirectory if it doesn't already exist
if (! -e "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS"){
  system("mkdir /wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS") && die "Couldn't create directory\n";
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

 ########################################
 # Open logfile                         #
 ########################################

my $logfile = "/wormsrv2/logs/GFFsplitter.$rundate.$$";
open (LOG,">$logfile");
LOG->autoflush();

print LOG "# GFFsplitter\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";

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

our $debug = 1;
my %GFF;
my $file;

foreach $file (@gff_files) {
    
    next if ($file eq "");
    
    print LOG "# File $file\n";
    print     "# File $file\n" if ($debug);

    my $line_count = 0;

    #########################################################
    # open the gff file                                     #
    #########################################################
    
    my $running_total = 0;
    my ($name,$feature,$method,$start,$stop,$score,$strand,$other);
    my @header = "";

    open (GFF, "</wormsrv2/autoace/CHROMOSOMES/$file.gff");
    while (<GFF>) {
	chomp;
	$line_count++;
	
	print "." if (($line_count % 5000) == 0);

	#skip header lines of file
	if (/^\#/) {push (@header,$_); next;}

	($name,$method,$feature,$start,$stop,$score,$strand,$other,$name) = split /\t/;

	# Clone path
	if (($method eq "Genomic_canonical") 
	    && ($feature eq "Sequence"))                  {push ( @{$GFF{$file}{clone_path}},$_);}
	# Genes     
	elsif ( (($method eq "curated") && ($feature eq "Sequence")) ||
		(($method eq "provisional") && ($feature eq "Sequence")) )
		{push (@{$GFF{$file}{genes}},$_);}
	elsif ( ($method eq "Pseudogene")   && ($feature eq "Sequence"))              
		{push (@{$GFF{$file}{pseudogenes}},$_);}
	# RNA genes 
	elsif ((($method eq "RNA")    || ($method eq "tRNAscan-SE-1.11") ||
		($method eq "rRNA")   || ($method eq "scRNA") ||
		($method eq "snRNA")  || ($method eq "snoRNA") || 
		($method eq "miRNA")) && ($feature eq "Transcript"))               
		{push (@{$GFF{$file}{rna}},$_);}
	# CDS Exon  
	elsif ( (($method eq "curated") || ($method eq "provisional")) 
		&& ($feature eq "CDS"))                   {push (@{$GFF{$file}{CDS_exon}},$_);}
	# Exon      
	elsif ((($method eq "curated") 
		|| ($method eq "provisional")) 
	       && ($feature eq "exon"))                   {push (@{$GFF{$file}{exon}},$_);}
       	elsif (($method eq "tRNAscan-SE-1.11") 
	       && ($feature eq "exon"))                   {push (@{$GFF{$file}{exon_tRNA}},$_);}
	elsif (($method eq "Pseudogene") 
	       && ($feature eq "exon"))                   {push (@{$GFF{$file}{exon_pseudogene}},$_);}
	# Intron    
	elsif ((($method eq "curated") 
		|| ($method eq "provisional")) 
	       && ($feature eq "intron"))                 {push (@{$GFF{$file}{intron}},$_);}
	elsif (($method eq "tRNAscan-SE-1.11") 
	       && ($feature eq "intron"))                 {push (@{$GFF{$file}{intron_tRNA}},$_);}
	elsif (($method eq "Pseudogene") 
	       && ($feature eq "intron"))                 {push (@{$GFF{$file}{intron_pseudogene}},$_);}
	# Intron confirmed by ESTs
	elsif (/Confirmed_by_EST/)                        {push (@{$GFF{$file}{intron_confirmed_CDS}},$_);}
	elsif (/Confirmed_by_cDNA/)                       {push (@{$GFF{$file}{intron_confirmed_CDS}},$_);}
	elsif (/Confirmed_in_UTR/)                        {push (@{$GFF{$file}{intron_confirmed_UTR}},$_);}
        # Repeats
	elsif (($method eq "tandem") 
            || ($method eq "inverted") 
            || ($method eq "hmmfs.3") 
            || ($method eq "scan") 
            || ($feature eq "repeat_region"))             {push (@{$GFF{$file}{repeats}},$_);}
   	# TC1 insertions
	elsif ($method eq "BLASTN_TC1")                   {push (@{$GFF{$file}{tc_insertions}},$_);}
   	# Protein similarities
	elsif ($method eq "BLASTX")                       {push (@{$GFF{$file}{BLASTX}},$_);}
   	# C.briggsae similarities
	elsif (($method eq "WABA_coding") 
	       || ($method eq "WABA_strong")
               || ($method eq "WABA_weak"))               {push (@{$GFF{$file}{WABA_BRIGGSAE}},$_);}
   	# Assembly tags
	elsif ($method eq "assembly_tag")                 {push (@{$GFF{$file}{assembly_tags}},$_);}
	# TS site
	elsif ((/misc_feature/) && 
               (/trans-splice site/))                     {push (@{$GFF{$file}{ts_site}},$_);}
   	# Oligo mapping
	elsif ($feature eq "OLIGO")                       {push (@{$GFF{$file}{oligos}},$_);}
   	# RNAi
	elsif ($method eq "RNAi")                         {push (@{$GFF{$file}{RNAi}},$_);}
   	# GENEPAIR
	elsif ($method eq "GenePair_STS")                 {push (@{$GFF{$file}{genepair}},$_);}
   	# Alleles
	elsif (/Allele/)                                  {push (@{$GFF{$file}{allele}},$_);}
   	# Clone ends
	elsif ((/Clone_left_end/)  
               || (/Clone_right_end/))                    {push (@{$GFF{$file}{clone_ends}},$_);}
   	# PCR Products
	elsif (/PCR_product/)                             {push (@{$GFF{$file}{PCR_products}},$_);}
   	# cDNA for RNAi
	elsif (/cDNA_for_RNAi/)                           {push (@{$GFF{$file}{cDNA_for_RNAi}},$_);}
   	# BLAT_EST
	elsif (/BLAT_EST_BEST/)                           {push (@{$GFF{$file}{BLAT_EST_BEST}},$_);}
	elsif (/BLAT_EST_OTHER/)                          {push (@{$GFF{$file}{BLAT_EST_OTHER}},$_);}
   	# BLAT_mRNA
	elsif (/BLAT_mRNA_BEST/)                          {push (@{$GFF{$file}{BLAT_mRNA_BEST}},$_);}
	elsif (/BLAT_mRNA_OTHER/)                         {push (@{$GFF{$file}{BLAT_mRNA_OTHER}},$_);}
   	# BLAT_EMBL
	elsif (/BLAT_EMBL_BEST/)                          {push (@{$GFF{$file}{BLAT_EMBL_BEST}},$_);}
	elsif (/BLAT_EMBL_OTHER/)                         {push (@{$GFF{$file}{BLAT_EMBL_OTHER}},$_);}
   	# BLATX_NEMATODE
	elsif (/BLATX_NEMATODE/)                          {push (@{$GFF{$file}{BLATX_NEMATODE}},$_);}
   	# Expr_profile
	elsif (/Expr_profile/)                            {push (@{$GFF{$file}{Expr_profile}},$_);}
   	# UTR         
	elsif (/UTR/)                                     {push (@{$GFF{$file}{UTR}},$_);}
	# REST OF LINES
	else                                              {push (@{$GFF{$file}{rest}},$_); $running_total--;}
	
	# increment no of lines sorted to bin
	$running_total++; 
	
    }

    ########################################
    # write some output                    #
    ########################################

    foreach my $tag (@gff_classes) {
	print "# $file $tag\n";

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
    unless(-e "$lockdir/D1:Build_wormpep_final"){
      system ("GFF_with_wormpep_accessions.pl $datadir/GFF_SPLITS/$file.genes.gff > $datadir/GFF_SPLITS/$file.genes_acc.gff");
      $input_file = "$datadir/GFF_SPLITS/$file.genes.gff";
      $output_file = "$datadir/GFF_SPLITS/$file.genes_acc.gff";
      &GFF_genes_with_accessions("$input_file", "$output_file");
      system ("mv -f $output_file $input_file");
    }
    
    # GFF UTRs with CDS names
    # Shouldn't attempt to do this if UTR data has not been generated
    unless( -e "$lockdir/B10:Generate_UTR_data" ) {
      my $utr_file = "$datadir/GFF_SPLITS/$file.UTR.gff";
      my $utr_cds_file = "$datadir/GFF_SPLITS/$file.UTR_CDS.gff";
      &GFF_with_UTR("$utr_file","$utr_cds_file");
  }
    
}
close LOG;


###############################
# Mail log to curator         #
###############################

open (OUTLOG,  "|/usr/bin/mailx -s \"WormBase Report: GFFsplitter\" $maintainers ");
open (READLOG, "<$logfile");
while (<READLOG>) {
    print OUTLOG "$_";
}
close READLOG;
close OUTLOG;

###############################
# hasta luego                 #
###############################

exit(0);

###############################
# subroutines                 #
###############################

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
  open( UTR, "$utr" );
  open( UTR_CDS, "$utr_cds");
  while (<UTR>) {
    
    print UTR_CDS $_ if (/^\#/);
    (/_UTR:(\S+)\"/);
    print UTR_CDS "$`" . "_UTR:$1\" CDS=\"$1\"\n";
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
  
    my @acc = $obj->DB_info->row();
    
    print OUT "$_ acc=$acc[3] ver=$acc[4]\n";
    
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
    (my $gen_can) = $gff[8] =~ /Sequence \"(\S+)\"/; 

    my $obj = $db->fetch(Sequence=>$gen_can);
    if (!defined ($obj)) {
      print "Could not fetch sequence '$gen_can'\n";
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
CDS_exon
exon
exon_tRNA
exon_pseudogene
intron
intron_tRNA
intron_pseudogene
intron_confirmed_CDS
intron_confirmed_UTR
repeats
tc_insertions
BLASTX
WABA_BRIGGSAE
assembly_tags
ts_site
oligos
RNAi
genepair
alleles
clone_ends
PCR_products
cDNA_for_RNAi
BLAT_EST_BEST
BLAT_EST_OTHER
BLAT_mRNA_BEST
BLAT_mRNA_OTHER
BLAT_EMBL_BEST
BLAT_EMBL_OTHER
BLATX_NEMATODE
Expr_profile
UTR
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



=back


=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut

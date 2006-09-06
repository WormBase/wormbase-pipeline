#!/usr/local/bin/perl5.8.0 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase & wormpep FTP sites
# 
# Originally written by Dan Lawson
#
# Last updated by: $Author: wormpub $
# Last updated on: $Date: 2006-09-06 13:14:21 $
#
# see pod documentation (i.e. 'perldoc make_FTP_sites.pl') for more information.
#
##########################################################################################


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use IO::Handle;


#################################################################################
# Command-line options and variables                                            #
#################################################################################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $release; # only copy across release files
my $chroms;  # only copy across chromosome files
my $misc;    # only copy misc files
my $wormpep; # only copy wormpep files
my $genes;   # only copy confirmed genes
my $cDNA;    # only copy cDNA2orf file
my $geneIDs; # only copy file of gene IDs
my $pcr;     # only copy file of PCR products
my $homols;  # only copy best blast hits 
my $all;     # copy everything across

GetOptions ("help"     => \$help,
	    "debug=s"  => \$debug,
	    "test"     => \$test,
	    "verbose"  => \$verbose,
	    "store:s"    => \$store,
	    "release"  => \$release,
	    "chroms"   => \$chroms,
	    "misc"     => \$misc,
	    "wormpep"  => \$wormpep,
	    "genes"    => \$genes,
	    "cDNAlist" => \$cDNA,
	    "geneIDs"  => \$geneIDs,
	    "pcr"      => \$pcr,
	    "homols"   => \$homols,
	    "all"      => \$all);


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


# using -all option?
($release=$chroms=$misc=$wormpep=$genes=$cDNA=$geneIDs=$pcr=$homols = 1 ) if ($all);

my $base_dir = $wormbase->basedir;    # The BUILD directory
my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $targetdir = "/nfs/disk69/ftp/pub/wormbase";  # default directory, can be overidden

my $WS              = $wormbase->get_wormbase_version();      # e.g.   132
my $WS_name         = $wormbase->get_wormbase_version_name(); # e.g. WS132
my $maintainers     = "All";
my $runtime;



#################################################################################
# Main                                                                          #
#################################################################################


&copy_release_files if ($release);    # make a new directory for the WS release and copy across release files

&copy_chromosome_files if ($chroms);  # make a new /CHROMOSOMES directory for the DNA, GFF, and agp files and copy files across

&copy_misc_files if ($misc);          # copy across models.wrm and other misc. files, e.g. wormRNA

&copy_wormpep_files if ($wormpep);    # copied from ~wormpub/WORMPEP

&extract_confirmed_genes if ($genes); # make file of confirmed genes from autoace and copy across

&make_cDNA2ORF_list if ($cDNA);       # make file of cDNA -> ORF connections and add to FTP site

&make_geneID_list if ($geneIDs);      # make file of WBGene IDs -> CGC name & Sequence name and add to FTP site

&make_pcr_list if ($pcr);             # make file of PCR products -> WBGene IDs, CDS, CGC name

&copy_homol_data if ($homols);        # copies best blast hits files across



################################
#
# Tidy up and exit
#
################################


# warn about errors in subject line if there were any
my $errors = $log->report_errors;
$log->mail;
print "Finished.\n" if ($verbose);
exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################


##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying release files\n");

  $wormbase->run_command("mkdir $targetdir/$WS_name", $log) unless -e "$targetdir/$WS_name";

  my $filename;

  opendir (RELEASE,"$ace_dir/release") or croak ("Could not open directory $ace_dir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    $wormbase->run_command("scp $ace_dir/release/$filename $targetdir/$WS_name/$filename", $log);

    my $O_SIZE = (stat("$ace_dir/release/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$WS_name/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir RELEASE;
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished\n\n");
  
}

##################################################
# copy the DNA, GFF, and agp files across
##################################################

sub copy_chromosome_files{

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying chromosome files\n");
  my $filename;
  $wormbase->run_command("mkdir $targetdir/$WS_name/CHROMOSOMES", $log) unless -e "$targetdir/$WS_name/CHROMOSOMES";

  opendir (DNAGFF,"$ace_dir/CHROMOSOMES") or croak ("Could not open directory $ace_dir/CHROMOSOMES");
  while (defined($filename = readdir(DNAGFF))) {
    if (($filename eq ".")||($filename eq "..")||($filename eq "SUPPLEMENTARY_GFF")) { next;}
    $wormbase->run_command("scp $ace_dir/CHROMOSOMES/$filename $targetdir/$WS_name/CHROMOSOMES/$filename", $log);
    my $O_SIZE = (stat("$ace_dir/CHROMOSOMES/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$WS_name/CHROMOSOMES/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir DNAGFF;

  $wormbase->run_command("mkdir $targetdir/$WS_name/CHROMOSOMES/SUPPLEMENTARY_GFF", $log) unless -e "$targetdir/$WS_name/CHROMOSOMES/SUPPLEMENTARY_GFF";
    opendir (DNAGFFSUP,"$ace_dir/CHROMOSOMES/SUPPLEMENTARY_GFF") or croak ("Could not open directory $ace_dir/CHROMOSOMES/SUPPLEMENTARY_GFF");
  while (defined($filename = readdir(DNAGFFSUP))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    $wormbase->run_command("scp $ace_dir/CHROMOSOMES/SUPPLEMENTARY_GFF/$filename $targetdir/$WS_name/CHROMOSOMES/SUPPLEMENTARY_GFF/$filename", $log);
    my $O_SIZE = (stat("$ace_dir/CHROMOSOMES/SUPPLEMENTARY_GFF/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$WS_name/CHROMOSOMES/SUPPLEMENTARY_GFF/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy SUPPLEMENTARY_GFF/$filename\n";
    } 
  }
  closedir DNAGFFSUP;

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

############################################
# copy across models.wrm and misc files
#############################################

sub copy_misc_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying misc files\n");

  # Copy across the models.wrm file
  $wormbase->run_command("scp $ace_dir/wspec/models.wrm $targetdir/$WS_name/models.wrm.$WS_name", $log);

  # copy some miscellaneous files across
  my $old_release = $WS -1;
  $wormbase->run_command("scp ".	$wormbase->compare."/WS$old_release-$WS_name.dbcomp $targetdir/$WS_name/", $log);
  $wormbase->run_command("scp $base_dir/autoace_config/INSTALL $targetdir/$WS_name/", $log);

  # tar, zip, and copy WormRNA files across from BUILD/WORMRNA
  my $dest = "$base_dir/WORMRNA/wormrna$WS";
  chdir "$dest" or croak "Couldn't cd $dest\n";
  $wormbase->run_command("/bin/tar -cf $targetdir/$WS_name/wormrna$WS.tar README wormrna$WS.rna", $log);
  $wormbase->run_command("/bin/gzip $targetdir/$WS_name/wormrna$WS.tar", $log);

  # zip and copy the microarray oligo mapping files.
  chdir "$ace_dir";
  $wormbase->run_command("/bin/gzip *oligo_mapping", $log);
  $wormbase->run_command("cp *oligo_mapping.gz $targetdir/$WS_name/", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");

}


############################################
# copy across wormpep files
#############################################

sub copy_wormpep_files{

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying wormpep files\n");

  my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";
  my $wp_source_dir = $wormbase->wormpep;
  my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
  my $wp_ftp_dir = "$wormpep_ftp_root/wormpep$WS";
  mkdir $wp_ftp_dir unless -e $wp_ftp_dir;

 
  my @wormpep_files = $wormbase->wormpep_files;
  
  foreach my $file ( @wormpep_files ){
  # copy the wormpep release files across
    $wormbase->run_command("scp $wp_source_dir/$file$WS $wp_ftp_dir/$file$WS", $log);
    &CheckSize("$wp_source_dir/$file$WS","$wp_ftp_dir/$file$WS");
  }

  # tar up the latest wormpep release and copy across
  my $tgz_file = "$wp_source_dir/wormpep$WS.tar";
  my $command = "/bin/tar -c -h -P \"$base_dir/WORMPEP/\" -f $tgz_file";

  # grab list of wormpep file names from subroutine
  foreach my $file (@wormpep_files){
      $command .= " $wp_source_dir/$file$WS";
  }
  $wormbase->run_command("$command", $log);
  $wormbase->run_command("/bin/gzip $tgz_file", $log);
  $tgz_file .= ".gz";
  $wormbase->run_command("mv $tgz_file $targetdir/$WS_name", $log);


  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting confirmed genes\n");

  my $db = Ace->connect(-path  => "$ace_dir/");
  my $query = "Find elegans_CDS; Confirmed";
  my @confirmed_genes   = $db->fetch(-query=>$query);

  open(OUT,">${targetdir}/$WS_name/confirmed_genes.$WS_name") || croak "Couldn't write to ${targetdir}/$WS_name/confirmed_genes.$WS_name\n";

  foreach my $seq (@confirmed_genes){
    my $dna = $seq->asDNA();
    print OUT "$dna";
  }

  close(OUT);
  $wormbase->run_command("/bin/gzip ${targetdir}/$WS_name/confirmed_genes.$WS_name", $log);

  $db->close;

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished extracting\n\n");
  return(0);
}

################################################################################
# make list of cDNA -> orf connections
################################################################################

sub make_cDNA2ORF_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making cDNA2ORF files\n");
  # simple routine to just get cDNA est names and their correct ORFs and make an FTP site file
  # two columns, second column supports multiple ORF names

  my $tace = $wormbase->tace;
  my $command=<<EOF;
Table-maker -p "$ace_dir/wquery/cDNA2CDS.def" quit
EOF

  my $dir = "$ace_dir";

  my %cDNA2orf;
  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    chomp;
    if (m/^\"/){
      s/\"//g;
      m/(.+)\s+(.+)/;     
      $cDNA2orf{$1} .= "$2 ";
    }
  }
  close(TACE);
  
  # output to ftp site

  
  my $out = "$targetdir/$WS_name/cDNA2orf.$WS_name";
  open(OUT, ">$out") || croak "Couldn't open $out\n";

  foreach my $key (keys %cDNA2orf){
    print OUT "$key,$cDNA2orf{$key}\n";
  }
  close(OUT);

  $wormbase->run_command("/bin/gzip $out", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making files\n\n");  
  
}

################################################################################
# make list of WBGene IDs to CGC name and Sequence name
################################################################################

sub make_geneID_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making Gene ID list\n");
  # For each 'live' Gene object, extract 'CGC_name' and 'Sequence_name' fields (if present)

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/gene2cgc_name_and_sequence_name.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$targetdir/$WS_name/geneIDs.$WS_name";

  open (OUT, ">$out") || croak "Couldn't open $out\n";
  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
      if (m/^\"/){
	  s/\"//g;
	  tr/\t/,/;
	  print OUT "$_";
      }
  }
  close(TACE);
  close(OUT);

  $wormbase->run_command("/bin/gzip $out", $log);
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making list\n\n");
}


################################################################################
# make list of PCR_product connections to CDS and Gene ID plus CGC name
################################################################################

sub make_pcr_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making PCR product 2 gene list list\n");

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/pcr_product2gene.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$targetdir/$WS_name/pcr_product2gene.$WS_name";

  # hashes needed because one pcr product may hit two or more genes
  my %pcr2gene;
  my %gene2sequence;
  my %gene2cgc;

  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    if (m/^\"/){
      s/\"//g;
      chomp;
      # split the line into various fields
      my ($pcr,$gene,$cgc,$sequence) = split(/\t/, $_) ;

      # fill hashes and arrays, allow for multiple genes per PCR product
      $pcr2gene{$pcr}      .= "$gene,";
      ($gene2cgc{$gene} = $cgc) if($cgc);
      $gene2sequence{$gene} = "$sequence";
    }
  }
  close(TACE);



# Now write output, cycling through list of pcr products in %pcr2gene

  open (OUT, ">$out") || croak "Couldn't open $out\n";

  foreach my $pcr (keys %pcr2gene){
    my @genes = split(/,/,$pcr2gene{$pcr});
    my $counter =0;
    print OUT "$pcr";
    foreach my $gene (@genes){

      # remove next element if it is the same gene ID and start loop again
      if (defined($genes[$counter+1]) && $genes[$counter+1] eq $genes[$counter]){
	splice(@genes, $counter+1,1);
	redo;
      }
      $counter++;
      print OUT "\t$gene";
      print OUT "($gene2cgc{$gene})" if (exists($gene2cgc{$gene}));
      
      # now print sequence name
      print OUT ",$gene2sequence{$gene}";
      
    }
    print OUT "\n";
  }

  close(OUT);

  $wormbase->run_command("/bin/gzip $out", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making list\n\n");
}

########################

sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    $log->write_to("\tERROR: $first SRC: $F_SIZE TGT: $S_SIZE - not same size, please check\n");
  } 
}

############################################################
# copy best blast hits file to ftp site
############################################################

sub copy_homol_data{


  my $blat_dir  = "$ace_dir/BLAT";
  my $blast_dir = $wormbase->acefiles;

  # does this to tidy up???? Not sure why these lines are here, krb
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_pep_blastp.ace", $log);
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_brigpep_blastp.ace", $log);
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_dna_blastx.ace", $log);
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_pep_motif_info.ace", $log);
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_brigpep_motif_info.ace", $log);

  # compress best blast hits files and then copy to FTP site
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_pep_best_blastp_hits", $log);
  $wormbase->run_command("/bin/gzip -f $blast_dir/worm_brigpep_best_blastp_hits", $log);
  $wormbase->run_command("scp $blast_dir/worm_pep_best_blastp_hits.gz      $targetdir/$WS_name/best_blastp_hits.$WS_name.gz", $log);
  $wormbase->run_command("scp $blast_dir/worm_brigpep_best_blastp_hits.gz  $targetdir/$WS_name/best_blastp_hits_brigpep.$WS_name.gz", $log);


}

##################################################################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##################################################################################





__END__

=pod

=head1 NAME - make_FTP_sites.pl

=back


=head1 USAGE

=over 4

=item make_FTP_sites.pl

=back


This script does :

 [01] - make a new directory for the WS release
 [02] - copy the WS release files to the target directory
 [03] - make a new directory for the chromosome DNA/GFF/AGP files
 [04] - copy the chromosome DNA/GFF/AGP files to the target directory
 [05] - copy the models.wrm file across (also present in the database.*.tar.gz files)
 [06] - copy the relevant dbcomp file across
 [07] - copy across latest wormpep release
 [08] - make wormpep FTP site
 [09] - copy WormRNA release across
 [10] - extract confirmed genes from autoace and make a file on FTP site
 [11] - delete the old symbolic link and make the new one
 [12] - delete the old WS release version directory
 [13] - makes a file of cDNA2orf connections
 [14] - makes a file of all gene IDs with CGC names and Sequence names (where present)
 [15] - exit gracefully


=over 4

=item MANDATORY arguments:

none

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)


=head1 AUTHOR - Originally written by Dan Lawson (dl1@sanger.ac.uk) with extensive rewrite
and upgrade by Keith Bradnam (krb@sanger.ac.uk)

=cut



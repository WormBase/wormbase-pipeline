#!/usr/local/bin/perl5.6.1 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase (public), wormbase (private), wormpep sites
# 
# Originally written by Dan Lawson
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-01-23 18:15:22 $
#
# see pod documentation (i.e. 'perldoc make_FTP_sites.pl') for more information.
#
##########################################################################################


use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use Getopt::Long;
use Ace;
use IO::Handle;
$|=1;


#################################################################################
# Command-line options and variables                                            #
#################################################################################

our $sourcedir = "/wormsrv2/autoace";
our $targetdir = "/nfs/disk69/ftp/pub/wormbase";  # default directory, can be overidden
our $runtime;
our $log;

our $release            = &get_wormbase_version_name(); # e.g. WS89
our $release_number     = &get_wormbase_version();      # e.g.   89
our $old_release        = $release_number - 1;
our $wormrna_release    = $release_number;
our $wormpep            = $release_number;
our ($help,$debug);
my $maintainers = "All";

GetOptions (
	    "help"     => \$help,
	    "debug=s"   => \$debug
	   );

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


#################################################################################
# Main                                                                          #
#################################################################################


&create_log_file;

&copy_release_files;    # make a new directory for the WS release and copy across release files

&copy_chromosome_files; # make a new /CHROMOSOMES directory for the DNA, GFF, and agp files and copy files across

&copy_misc_files;       # copy across models.wrm and other misc. files, e.g. wormRNA

&copy_wormpep_files;    # copied from ~wormpub/WORMPEP

&extract_confirmed_genes; # make file of confirmed genes from autoace and copy across

&copy_homol_data;        # copy blat and blast data to private ftp site for St. Louis

$runtime = &runtime;
print LOG "\n\n\nmake_FTP_sites.pl finished at $runtime\n";

close(LOG);
&mail_maintainer("WormBase Report: make_FTP_sites.pl","$maintainers",$log);


exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################

########################################
# Open logfile                         #
########################################

sub create_log_file{
  my $rundate = `date +%y%m%d`; chomp $rundate;
  $runtime = &runtime;
  $log="/wormsrv2/logs/make_FTP_sites.pl.$rundate.$$";
  open (LOG,">$log");
  LOG->autoflush();

  print LOG "make_FTP_sites.pl started at $runtime\n\n";
}

##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{
  $runtime = &runtime;
  print LOG "$runtime: Starting to copy release directory to FTP site\n\n";

  system ("mkdir $targetdir/$release") unless -e "$targetdir/$release";

  my $filename;

  opendir (RELEASE,"$sourcedir/release") or die ("Could not open directory $sourcedir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    system ("cp $sourcedir/release/$filename $targetdir/$release/$filename") && die "ERROR: can't copy source directory\n";
    my $O_SIZE = (stat("$sourcedir/release/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$release/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      print LOG "\tError - file $filename not transferred regularly - please check\n";
      die();
    } 
    else {
      print LOG "\tCopied filename: $filename SRC: $O_SIZE TGT: $N_SIZE\n";
    }
  }
  closedir RELEASE;

  $runtime = &runtime;
  print LOG "$runtime: Finished copying the release directory to FTP site\n";
  
}

##################################################
# copy the DNA, GFF, and agp files across
##################################################

sub copy_chromosome_files{

  $runtime = &runtime;
  print LOG "Starting to copy CHROMOSOMES directory : '$targetdir/$release/CHROMOSOMES'\n";

  my $filename;
  system ("mkdir $targetdir/$release/CHROMOSOMES") unless -e "$targetdir/$release/CHROMOSOMES";

  opendir (DNAGFF,"$sourcedir/CHROMOSOMES") or die ("Could not open directory $sourcedir/CHROMOSOMES");
  while (defined($filename = readdir(DNAGFF))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    system ("cp $sourcedir/CHROMOSOMES/$filename $targetdir/$release/CHROMOSOMES/$filename") && die "ERROR: Can't copy chromosome file\n";
    my $O_SIZE = (stat("$sourcedir/CHROMOSOMES/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$release/CHROMOSOMES/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      print LOG "Error - file $filename not transferred regularly - please check\n";
      die();
    } 
    else {
      print LOG "\tCopied filename: $filename SRC: $O_SIZE TGT: $N_SIZE\n";
    }
  }
  closedir DNAGFF;
  
  $runtime = &runtime;
  print LOG "$runtime: Finished copying the CHROMOSOMES directory\n";

}

############################################
# copy across models.wrm and misc files
#############################################

sub copy_misc_files{

  $runtime = &runtime;
  print LOG "$runtime: Started copying misc files\n";


  # Copy across the models.wrm file
  system("cp $sourcedir/wspec/models.wrm $targetdir/$release/models.wrm.$release") && die "ERROR: can't copy models.wrm\n";

  # copy some miscellaneous files across
  system("cp /wormsrv2/autoace/COMPARE/WS$old_release-$release.dbcomp $targetdir/$release/") && die "ERROR: can't copy dbcomp file\n";
  
  system("cp /wormsrv2/autoace_config/INSTALL $targetdir/$release/") && die "ERROR: can't copy INSTALL script\n";

  # tar, zip, and copy WormRNA files across from wormsrv2/WORMRNA
  my $dest = "/wormsrv2/WORMRNA/wormrna${wormrna_release}";
  chdir "$dest" or die "Couldn't cd $dest\n";
  system("/bin/tar -cf $targetdir/$release/wormrna${wormrna_release}.tar README wormrna${wormrna_release}.rna") && die "ERROR: can't create wormrna tar file\n";
  system("/bin/gzip $targetdir/$release/wormrna${wormrna_release}.tar") && die "ERROR: can't create wormrna gzip file\n";


  $runtime = &runtime;
  print LOG "$runtime: Finished copying misc files\n";

}


############################################
# copy across wormpep files
#############################################

sub copy_wormpep_files{

  $runtime = &runtime;
  print LOG "$runtime: Started copying wormpep files\n";

  # move wormpep release from /wormsrv2 to /disk100/wormpub
  my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";

  print LOG "\tRemoving old files in $wormpub_dir\n";

  unlink("$wormpub_dir/wormpep_current")           || print LOG "ERROR: Cannot delete files in $wormpub_dir\n";
  unlink("$wormpub_dir/wormpep.accession_current") || print LOG "ERROR: Cannot delete files in $wormpub_dir\n";
  unlink("$wormpub_dir/wormpep.dna_current")       || print LOG "ERROR: Cannot delete files in $wormpub_dir\n";
  unlink("$wormpub_dir/wormpep.history_current")   || print LOG "ERROR: Cannot delete files in $wormpub_dir\n";
  unlink("$wormpub_dir/wp.fasta_current")          || print LOG "ERROR: Cannot delete files in $wormpub_dir\n";


  print LOG "\tCopying new files to $wormpub_dir\n";
  my $new_wpdir = "/wormsrv2/WORMPEP/wormpep${wormpep}";
  system ("cp $new_wpdir/wormpep_current $wormpub_dir/wormpep_current")                     && die "ERROR: Cannot copy file to $wormpub_dir\n";
  system ("cp $new_wpdir/wormpep.accession$wormpep $wormpub_dir/wormpep.accession_current") && die "ERROR: Cannot copy file to $wormpub_dir\n";
  system ("cp $new_wpdir/wormpep.dna$wormpep $wormpub_dir/wormpep.dna_current")             && die "ERROR: Cannot copy file to $wormpub_dir\n";
  system ("cp $new_wpdir/wormpep.history$wormpep $wormpub_dir/wormpep.history_current")     && die "ERROR: Cannot copy file to $wormpub_dir\n";
  system ("cp $new_wpdir/wp.fasta$wormpep $wormpub_dir/wp.fasta_current")                   && die "ERROR: Cannot copy file to $wormpub_dir\n";
  system ("/usr/local/pubseq/bin/setdb $wormpub_dir/wormpep_current")                       && die "ERROR: setdb in wormpub failed\n";
  system ("chmod +rw $new_wpdir/*") && die "ERROR: chmod command failed\n";


  # tar up the latest wormpep release and copy across
  print LOG "Creating tar file and copying across to ftp site\n";
  system ("/bin/tar -hcf $wormpub_dir/wormpep${wormpep}.tar $new_wpdir/wormpep${wormpep} $new_wpdir/wormpep.accession${wormpep} $new_wpdir/wormpep.diff${wormpep} $new_wpdir/wormpep.dna${wormpep} $new_wpdir/wormpep.history${wormpep} $new_wpdir/wormpep.table${wormpep} $new_wpdir/wp.fasta${wormpep}") && die "ERROR: Couldn't run tar command\n";
  system ("/bin/gzip $wormpub_dir/wormpep${wormpep}.tar") && die "ERROR: Cannot gzip tar file\n";

  system("mv $wormpub_dir/wormpep${wormpep}.tar.gz $targetdir/$release") && die "ERROR: Couldn't move wormpep gzip tar file from $wormpub_dir to $targetdir/$release\n";


  $runtime = &runtime;
  print LOG "$runtime: Starting to make wormpep ftp directory\n";


  my $wp_source_dir = "/wormsrv2/WORMPEP/wormpep${wormpep}";
  my $wp_ftp_dir = glob("~ftp/pub/databases/wormpep");
  my $wp_old_release = "old_wormpep"."$old_release";

  print LOG "\tMaking $wp_ftp_dir/$wp_old_release directory\n";
  system("mkdir $wp_ftp_dir/$wp_old_release") unless -e "$wp_ftp_dir/$wp_old_release";


  # move the actual files to old_ directory
  print LOG "\tstarting to move files from current wormpep directory to old wormpep directory on FTP site\n";

  system("mv $wp_ftp_dir/wormpep.accession $wp_ftp_dir/$wp_old_release/wormpep.accession$old_release");
  system("mv $wp_ftp_dir/wormpep.diff $wp_ftp_dir/$wp_old_release/wormpep.diff$old_release");
  system("mv $wp_ftp_dir/wormpep.dna $wp_ftp_dir/$wp_old_release/wormpep.dna$old_release"); 
  system("mv $wp_ftp_dir/wormpep.history $wp_ftp_dir/$wp_old_release/wormpep.history$old_release");  
  system("mv $wp_ftp_dir/wormpep.table $wp_ftp_dir/$wp_old_release/wormpep.table$old_release");   
  system("mv $wp_ftp_dir/wormpep${old_release} $wp_ftp_dir/$wp_old_release/");
  system("mv $wp_ftp_dir/wp.fasta $wp_ftp_dir/$wp_old_release/wp.fasta$old_release");
  print LOG "\tfinished moving files to old wormpep directory on FTP site\n";



  # copy the wormpep release files across
  system("cp $wp_source_dir/wormpep.accession$wormpep $wp_ftp_dir/wormpep.accession");
  &CheckSize("$wp_source_dir/wormpep.accession$wormpep","$wp_ftp_dir/wormpep.accession");
  print LOG "\tmove & rename wormpep.accession$wormpep => wormpep.accession\n";
  
  system("cp $wp_source_dir/wormpep.diff$wormpep $wp_ftp_dir/wormpep.diff");
  &CheckSize("$wp_source_dir/wormpep.diff$wormpep","$wp_ftp_dir/wormpep.diff");
  print LOG "\tmove & rename wormpep.diff$wormpep => wormpep.diff\n";
  
  system("cp $wp_source_dir/wormpep.dna$wormpep $wp_ftp_dir/wormpep.dna");
  &CheckSize("$wp_source_dir/wormpep.dna$wormpep","$wp_ftp_dir/wormpep.dna");
  print LOG "\tmove & rename wormpep.dna$wormpep => wormpep.dna\n";
  
  system("cp $wp_source_dir/wormpep.history$wormpep $wp_ftp_dir/wormpep.history");
  &CheckSize("$wp_source_dir/wormpep.history$wormpep","$wp_ftp_dir/wormpep.history");
  print LOG "\tmove & rename wormpep.history$wormpep => wormpep.history\n";
  
  system("cp $wp_source_dir/wormpep.table$wormpep $wp_ftp_dir/wormpep.table");
  &CheckSize("$wp_source_dir/wormpep.table$wormpep","$wp_ftp_dir/wormpep.table");
  print LOG "\tmove & rename wormpep.table$wormpep => wormpep.table\n";
  
  system("cp $wp_source_dir/wormpep$wormpep $wp_ftp_dir/wormpep${wormpep}");
  &CheckSize("$wp_source_dir/wormpep$wormpep","$wp_ftp_dir/wormpep${release}");
  print LOG "\tmove wormpep$wormpep =>  wormpep$wormpep\n";
  
  system("cp $wp_source_dir/wp.fasta$wormpep $wp_ftp_dir/wp.fasta");
  &CheckSize("$wp_source_dir/wp.fasta$wormpep","$wp_ftp_dir/wp.fasta");
  print LOG "\tmove & rename wp.fasta$wormpep => wp.fasta\n\n";
  

  # delete the old symbolic link and make the new one  
  system("cd $wp_ftp_dir; ln -fs wormpep$wormpep wormpep");
  system("cd $wp_ftp_dir; ln -fs $wp_old_release/wormpep$old_release wormpep.prev");
  print LOG "Deleted the old sym_link and created the new one\n\n";


  $runtime = &runtime;
  print LOG "$runtime: Finished copying wormpep files to ~wormpub, and creating wormbase and wormpep ftp directories\n";
}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  $runtime = &runtime;
  print LOG "$runtime: Starting to extract confirmed genes from autoace\n";


  my $db = Ace->connect(-path  => "/wormsrv2/autoace/");
  my $query = "Find Sequence; Confirmed_by";
  my @confirmed_genes   = $db->fetch(-query=>$query);


  open(OUT,">${targetdir}/${release}/confirmed_genes.${release}") || die "Couldn't write to ${targetdir}/${release}/confirmed_genes.${release}\n";

  foreach my $seq (@confirmed_genes){
    my $dna = $seq->asDNA();
    
    my (@type) = $seq->get('Confirmed_by');
    
    if(defined($type[1])){
      $dna =~ s/(>\w+\.\w+)/$1 Confirmed_by_EST_and_cDNA/;    
    }       
    else{
      $dna =~ s/(>\w+\.\w+)/$1 Confirmed_by_$type[0]/;        
    }       
    print OUT "$dna";
  }

  close(OUT);
  print LOG "gzipping confirmed_genes.${release}\n";
  system("/bin/gzip ${targetdir}/${release}/confirmed_genes.${release}") && print LOG "ERROR: couldn't gzip confirmed_genes\n";


  $runtime = &runtime;
  print LOG "$runtime: Finishing extract confirmed genes from autoace\n";
  
  $db->close;
  return(0);
  
}

########################

sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    print LOG "ERROR - file $first not transferred regularly - please check\n";
  } 
  else {
    print LOG "Copied filename: $first SRC: $F_SIZE TGT: $S_SIZE\n";
  }
}

############################################################
# copy BLAT results to the private FTP site for St. Louis
############################################################

sub copy_homol_data{

  $runtime = &runtime;
  print LOG "$runtime: Starting to copy BLAT data to private ftp site\n";

  my $blat_dir  = "/wormsrv2/autoace/BLAT";
  my $blast_dir = "/wormsrv2/wormbase/ensembl_dumps";
  my $private_ftp = "/nfs/privateftp/ftp-wormbase/pub/data/st_louis_homol_data";
  system("rm -f $private_ftp/*gz");
  
  system("cp $blat_dir/stlace.blat.EST.ace                  $private_ftp/${release}_stlace.blat.EST.ace");
  system("cp $blat_dir/stlace.blat.mRNA.ace                 $private_ftp/${release}_stlace.blat.mRNA.ace");
  system("cp $blat_dir/stlace.good_introns.EST.ace          $private_ftp/${release}_stlace.blat.good_introns.EST.ace");
  system("cp $blat_dir/stlace.good_introns.mRNA.ace         $private_ftp/${release}_stlace.blat.good_introns.mRNA.ace");
  
  system("cp $blat_dir/virtual_objects.stlace.BLAT_EST.ace  $private_ftp/${release}_virtual_objects.stlace.BLAT_EST.ace");
  system("cp $blat_dir/virtual_objects.stlace.BLAT_mRNA.ace $private_ftp/${release}_virtual_objects.stlace.BLAT_mRNA.ace");
  system("cp $blat_dir/virtual_objects.stlace.ci.EST.ace    $private_ftp/${release}_virtual_objects.stlace.ci.EST.ace");
  system("cp $blat_dir/virtual_objects.stlace.ci.mRNA.ace   $private_ftp/${release}_virtual_objects.stlace.ci.mRNA.ace");
  
  $runtime = &runtime;
  print LOG "$runtime: Starting to copy BLAST data to private ftp site\n";

  system("cp $blast_dir/blastp_ensembl.ace          $private_ftp/${release}_blastp_data.ace");
  system("cp $blast_dir/blastx_ensembl.ace          $private_ftp/${release}_blastx_data.ace");
  system("cp $blast_dir/ensembl_motif_info.ace      $private_ftp/${release}_protein_motif_data.ace");

  system("/bin/gzip $private_ftp/*ace");

  $runtime = &runtime;
  print LOG "$runtime: Finishing copying BLAST data to private ftp site\n";

}

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


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
 [11] - copies BLAT and BLAST data to private FTP site for St. Louis
 [12] - delete the old symbolic link and make the new one
 [13] - delete the old WS release version directory
 [14] - exit gracefully


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



#!/usr/local/bin/perl5.8.0 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase (public), wormbase (private), wormpep sites
# 
# Originally written by Dan Lawson
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2004-05-07 11:28:54 $
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
use Carp;

#################################################################################
# Command-line options and variables                                            #
#################################################################################

our $sourcedir = "/wormsrv2/autoace";
our $targetdir = "/nfs/disk69/ftp/pub/wormbase";  # default directory, can be overidden
our $log;

our $release            = &get_wormbase_version_name(); # e.g. WS89
our $release_number     = &get_wormbase_version();      # e.g.   89
our $old_release        = $release_number - 1;
our $wormrna_release    = $release_number;
our $wormpep            = $release_number;
our ($help,$debug,$norelease);
my $maintainers = "All";
my $errors = 0;    # tracking system call errors

GetOptions (
	    "help"       => \$help,
	    "debug=s"    => \$debug,
	    "norelease"  => \$norelease
	   );

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ($maintainers = $debug . '\@sanger.ac.uk');
}


#################################################################################
# Main                                                                          #
#################################################################################


&create_log_file;

&copy_release_files unless ($norelease);    # make a new directory for the WS release and copy across release files

&copy_chromosome_files; # make a new /CHROMOSOMES directory for the DNA, GFF, and agp files and copy files across

&copy_misc_files;       # copy across models.wrm and other misc. files, e.g. wormRNA

&copy_wormpep_files;    # copied from ~wormpub/WORMPEP

&extract_confirmed_genes; # make file of confirmed genes from autoace and copy across

&make_yk2ORF_list;       # make file of yk EST -> ORF connections and add to FTP site

#&copy_homol_data;        # copy blat and blast data to private ftp site for St. Louis



################################
#
# Tidy up and exit
#
################################


print LOG "\n\n\n",&runtime," make_FTP_sites.pl finished\n";

close(LOG);

# warn about errors in subject line if there were any
if($errors == 0){
  &mail_maintainer("BUILD REPORT: make_FTP_sites.pl",$maintainers,$log);
}
elsif ($errors ==1){
  &mail_maintainer("BUILD REPORT: make_FTP_sites.pl : $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer("BUILD REPORT: make_FTP_sites.pl : $errors ERRORS!!!",$maintainers,$log);
}

exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################

########################################
# Open logfile                         #
########################################

sub create_log_file{
  my $rundate = `date +%y%m%d`; chomp $rundate;
  $log="/wormsrv2/logs/make_FTP_sites.pl.$rundate.$$";
  open (LOG,">$log");
  LOG->autoflush();

  print LOG &runtime, " make_FTP_sites.pl started\n\n";
}

##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{

  &run_command("mkdir $targetdir/$release") unless -e "$targetdir/$release";

  my $filename;

  opendir (RELEASE,"$sourcedir/release") or croak ("Could not open directory $sourcedir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    &run_command("scp $sourcedir/release/$filename $targetdir/$release/$filename");

    my $O_SIZE = (stat("$sourcedir/release/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$release/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      print LOG "\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n";
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir RELEASE;

  
}

##################################################
# copy the DNA, GFF, and agp files across
##################################################

sub copy_chromosome_files{


  my $filename;
  &run_command("mkdir $targetdir/$release/CHROMOSOMES") unless -e "$targetdir/$release/CHROMOSOMES";

  opendir (DNAGFF,"$sourcedir/CHROMOSOMES") or croak ("Could not open directory $sourcedir/CHROMOSOMES");
  while (defined($filename = readdir(DNAGFF))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    &run_command("scp $sourcedir/CHROMOSOMES/$filename $targetdir/$release/CHROMOSOMES/$filename");
    my $O_SIZE = (stat("$sourcedir/CHROMOSOMES/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$release/CHROMOSOMES/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      print LOG "\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n";
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir DNAGFF;
  

}

############################################
# copy across models.wrm and misc files
#############################################

sub copy_misc_files{


  # Copy across the models.wrm file
  &run_command("scp $sourcedir/wspec/models.wrm $targetdir/$release/models.wrm.$release");

  # copy some miscellaneous files across
  &run_command("scp /wormsrv2/autoace/COMPARE/WS$old_release-$release.dbcomp $targetdir/$release/");
  &run_command("scp /wormsrv2/autoace_config/INSTALL $targetdir/$release/");

  # tar, zip, and copy WormRNA files across from wormsrv2/WORMRNA
  my $dest = "/wormsrv2/WORMRNA/wormrna${wormrna_release}";
  chdir "$dest" or croak "Couldn't cd $dest\n";
  &run_command("/bin/tar -cf $targetdir/$release/wormrna${wormrna_release}.tar README wormrna${wormrna_release}.rna");
  &run_command("/bin/gzip $targetdir/$release/wormrna${wormrna_release}.tar");

  # zip and copy interpolated map across from /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/

  chdir "/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/";
  &run_command("/bin/gzip WS*interpolated_map.txt"); 
  &run_command("scp ${release}_CDSes_interpolated_map.txt.gz $targetdir/$release/gene_interpolated_map_positions.${release}.gz");
  &run_command("scp ${release}_Clones_interpolated_map.txt.gz $targetdir/$release/clone_interpolated_map_positions.${release}.gz");


}


############################################
# copy across wormpep files
#############################################

sub copy_wormpep_files{


  my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";
  my $wp_source_dir = "/wormsrv2/WORMPEP/wormpep${wormpep}";
  my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
  my $wp_ftp_dir = "$wormpep_ftp_root/wormpep${wormpep}";
  mkdir $wp_ftp_dir unless -e $wp_ftp_dir;

  foreach my $file ( &wormpep_files ){
  # copy the wormpep release files across
    &run_command("scp $wp_source_dir/$file$wormpep $wp_ftp_dir/$file");
    &CheckSize("$wp_source_dir/$file$wormpep","$wp_ftp_dir/$file");
  }

  &run_command("cd $wormpep_ftp_root; ln -fs $wp_ftp_dir wormpep_dev");
}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{


  my $db = Ace->connect(-path  => "/wormsrv2/autoace/");
  my $query = "Find elegans_CDS; Confirmed_by";
  my @confirmed_genes   = $db->fetch(-query=>$query);


  open(OUT,">${targetdir}/${release}/confirmed_genes.${release}") || croak "Couldn't write to ${targetdir}/${release}/confirmed_genes.${release}\n";

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
  &run_command("/bin/gzip ${targetdir}/${release}/confirmed_genes.${release}");

  $db->close;
  return(0);
  
}

################################################################################
# make list of yk EST -> orf connections
################################################################################

sub make_yk2ORF_list {

  # simple routine to just get yk est names and their correct ORFs and make an FTP site file
  # two columns, second column supports multiple ORF names

  my $tace = &tace;
  my $command=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/yk2ORF.def" quit
EOF

  my $dir = "/wormsrv2/autoace";
  my %est2orf;
  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    chomp;
    if (m/^\"/){
      s/\"//g;
      m/(.+)\s+(.+)/;     
      $est2orf{$1} .= "$2 ";
    }
  }
  close(TACE);
  
  # output to ftp site

  
  my $out = "$targetdir/$release/yk2orf.$release";
  open(OUT, ">$out") || croak "Couldn't open $out\n";

  foreach my $key (keys %est2orf){
    print OUT "$key,$est2orf{$key}\n";
  }
  close(OUT);

  &run_command("/bin/gzip $out");
  
  
}

########################

sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    print LOG "\tERROR: $first SRC: $F_SIZE TGT: $S_SIZE - not same size, please check\n";
  } 
}

############################################################
# copy BLAT results to the private FTP site for St. Louis
############################################################

sub copy_homol_data{


  my $blat_dir  = "/wormsrv2/autoace/BLAT";
  my $blast_dir = "/wormsrv2/wormbase/ensembl_dumps";
  my $private_ftp = "/nfs/privateftp/ftp-wormbase/pub/data/st_louis_homol_data";
  &run_command("rm -f $private_ftp/*gz");


  &run_command("/bin/gzip -f $blast_dir/worm_pep_blastp.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_blastp.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_dna_blastx.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_pep_motif_info.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_motif_info.ace");

  &run_command("/bin/gzip -f $blast_dir/worm_pep_best_blastp_hits");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_best_blastp_hits");

  &run_command("scp $blast_dir/worm_pep_best_blastp_hits.gz      $targetdir/$release/best_blastp_hits.$release.gz");
  &run_command("scp $blast_dir/worm_brigpep_best_blastp_hits.gz  $targetdir/$release/best_blastp_hits_brigpep.$release.gz");

  &run_command("/bin/gzip $private_ftp/*ace");


}

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##################################################################################
#
# Simple routine which will run commands via system calls but also check the 
# return status of a system call and complain if non-zero, increments error check 
# count, and prints a log file error
#
##################################################################################

sub run_command{
  my $command = shift;
  print LOG &runtime, ": started running $command\n";
  my $status = system($command);
  if($status != 0){
    $errors++;
    print LOG "ERROR: $command failed\n";
  }

  # for optional further testing by calling subroutine
  return($status);
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



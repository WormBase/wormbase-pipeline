#!/usr/local/bin/perl5.8.0 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase & wormpep FTP sites
# 
# Originally written by Dan Lawson
#
# Last updated by: $Author: dl1 $
# Last updated on: $Date: 2004-10-08 15:10:10 $
#
# see pod documentation (i.e. 'perldoc make_FTP_sites.pl') for more information.
#
##########################################################################################


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;
use IO::Handle;
use Carp;

#################################################################################
# Command-line options and variables                                            #
#################################################################################

my $sourcedir = "/wormsrv2/autoace";
my $targetdir = "/nfs/disk69/ftp/pub/wormbase";  # default directory, can be overidden

my $WS              = &get_wormbase_version();      # e.g.   132
my $WS_name         = &get_wormbase_version_name(); # e.g. WS132
my $maintainers     = "All";
my $runtime;

my $help;
my $debug;
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

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ($maintainers = $debug . '\@sanger.ac.uk');
}

# using -all option?
($release=$chroms=$misc=$wormpep=$genes=$cDNA=$geneIDs=$pcr= 1) if ($all);


#################################################################################
# Main                                                                          #
#################################################################################

# open log file
my $log = Log_files->make_build_log();


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

$log->write_to("\n$errors errors found\n");

if($errors == 0){
  $log->mail("$maintainers","BUILD REPORT: make_FTP_sites.pl");
}
elsif ($errors ==1){
  $log->mail("$maintainers","BUILD REPORT: make_FTP_sites.pl : $errors ERROR!");
}
else{
  $log->mail("$maintainers","BUILD REPORT: make_FTP_sites.pl : $errors ERRORS!!!");
}

exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################


##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{
  $runtime = &runtime;
  $log->write_to("$runtime: copying release files\n");

  &run_command("mkdir $targetdir/$WS_name") unless -e "$targetdir/$WS_name";

  my $filename;

  opendir (RELEASE,"$sourcedir/release") or croak ("Could not open directory $sourcedir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    &run_command("scp $sourcedir/release/$filename $targetdir/$WS_name/$filename");

    my $O_SIZE = (stat("$sourcedir/release/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$WS_name/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir RELEASE;
  $runtime = &runtime;
  $log->write_to("$runtime: Finished\n\n");
  
}

##################################################
# copy the DNA, GFF, and agp files across
##################################################

sub copy_chromosome_files{

  $runtime = &runtime;
  $log->write_to("$runtime: copying chromosome files\n");
  my $filename;
  &run_command("mkdir $targetdir/$WS_name/CHROMOSOMES") unless -e "$targetdir/$WS_name/CHROMOSOMES";

  opendir (DNAGFF,"$sourcedir/CHROMOSOMES") or croak ("Could not open directory $sourcedir/CHROMOSOMES");
  while (defined($filename = readdir(DNAGFF))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    &run_command("scp $sourcedir/CHROMOSOMES/$filename $targetdir/$WS_name/CHROMOSOMES/$filename");
    my $O_SIZE = (stat("$sourcedir/CHROMOSOMES/$filename"))[7];
    my $N_SIZE = (stat("$targetdir/$WS_name/CHROMOSOMES/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir DNAGFF;
  $runtime = &runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

############################################
# copy across models.wrm and misc files
#############################################

sub copy_misc_files{
  $runtime = &runtime;
  $log->write_to("$runtime: copying misc files\n");

  # Copy across the models.wrm file
  &run_command("scp $sourcedir/wspec/models.wrm $targetdir/$WS_name/models.wrm.$WS_name");

  # copy some miscellaneous files across
  my $old_release = $WS -1;
  &run_command("scp /wormsrv2/autoace/COMPARE/WS$old_release-$WS_name.dbcomp $targetdir/$WS_name/");
  &run_command("scp /wormsrv2/autoace_config/INSTALL $targetdir/$WS_name/");

  # tar, zip, and copy WormRNA files across from wormsrv2/WORMRNA
  my $dest = "/wormsrv2/WORMRNA/wormrna$WS";
  chdir "$dest" or croak "Couldn't cd $dest\n";
  &run_command("/bin/tar -cf $targetdir/$WS_name/wormrna$WS.tar README wormrna$WS.rna");
  &run_command("/bin/gzip $targetdir/$WS_name/wormrna$WS.tar");

  # zip and copy interpolated map across from /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/

  chdir "/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/";
  &run_command("/bin/gzip WS*interpolated_map.txt"); 
  &run_command("scp ${WS_name}_CDSes_interpolated_map.txt.gz $targetdir/$WS_name/gene_interpolated_map_positions.${WS_name}.gz");
  &run_command("scp ${WS_name}_Clones_interpolated_map.txt.gz $targetdir/$WS_name/clone_interpolated_map_positions.${WS_name}.gz");

  $runtime = &runtime;
  $log->write_to("$runtime: Finished copying\n\n");

}


############################################
# copy across wormpep files
#############################################

sub copy_wormpep_files{

  $runtime = &runtime;
  $log->write_to("$runtime: copying wormpep files\n");

  my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";
  my $wp_source_dir = "/wormsrv2/WORMPEP/wormpep$WS";
  my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
  my $wp_ftp_dir = "$wormpep_ftp_root/wormpep$WS";
  mkdir $wp_ftp_dir unless -e $wp_ftp_dir;

 
  my @wormpep_files = &wormpep_files;
  
  foreach my $file ( @wormpep_files ){
  # copy the wormpep release files across
    &run_command("scp $wp_source_dir/$file$WS $wp_ftp_dir/$file$WS");
    &CheckSize("$wp_source_dir/$file$WS","$wp_ftp_dir/$file$WS");
  }

  # tar up the latest wormpep release and copy across
  my $tgz_file = "$wp_source_dir/wormpep$WS.tar";
  my $command = "/bin/tar -c -h -P \"/wormsrv2/WORMPEP/\" -f $tgz_file";

  # grab list of wormpep file names from subroutine
  foreach my $file (@wormpep_files){
      $command .= " $wp_source_dir/$file$WS";
  }
  &run_command("$command");
  &run_command("/bin/gzip $tgz_file");
  $tgz_file .= ".gz";
  &run_command("mv $tgz_file $targetdir/$WS_name");


  $runtime = &runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  $runtime = &runtime;
  $log->write_to("$runtime: Extracting confirmed genes\n");

  my $db = Ace->connect(-path  => "/wormsrv2/autoace/");
  my $query = "Find elegans_CDS; Confirmed_by";
  my @confirmed_genes   = $db->fetch(-query=>$query);


  open(OUT,">${targetdir}/$WS_name/confirmed_genes.$WS_name") || croak "Couldn't write to ${targetdir}/$WS_name/confirmed_genes.$WS_name\n";

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
  &run_command("/bin/gzip ${targetdir}/$WS_name/confirmed_genes.$WS_name");

  $db->close;

  $runtime = &runtime;
  $log->write_to("$runtime: Finished extracting\n\n");
  return(0);

}

################################################################################
# make list of cDNA -> orf connections
################################################################################

sub make_cDNA2ORF_list {

  $runtime = &runtime;
  $log->write_to("$runtime: making cDNA2ORF files\n");
  # simple routine to just get cDNA est names and their correct ORFs and make an FTP site file
  # two columns, second column supports multiple ORF names

  my $tace = &tace;
  my $command=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/cDNA2CDS.def" quit
EOF

  my $dir = "/wormsrv2/autoace";

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

  &run_command("/bin/gzip $out");

  $runtime = &runtime;
  $log->write_to("$runtime: Finished making files\n\n");  
  
}

################################################################################
# make list of WBGene IDs to CGC name and Sequence name
################################################################################

sub make_geneID_list {

  $runtime = &runtime;
  $log->write_to("$runtime: making Gene ID list\n");
  # For each 'live' Gene object, extract 'CGC_name' and 'Sequence_name' fields (if present)

  my $tace    = &tace;
  my $command = "Table-maker -p /wormsrv2/autoace/wquery/gene2cgc_name_and_sequence_name.def\nquit\n";
  my $dir     = "/wormsrv2/autoace";
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

  &run_command("/bin/gzip $out");
  
  $runtime = &runtime;
  $log->write_to("$runtime: Finished making list\n\n");
}


################################################################################
# make list of PCR_product connections to CDS and Gene ID plus CGC name
################################################################################

sub make_pcr_list {

  $runtime = &runtime;
  $log->write_to("$runtime: making PCR product 2 gene list list\n");

  my $tace    = &tace;
  my $command = "Table-maker -p /wormsrv2/autoace/wquery/pcr_product2gene.def\nquit\n";
  my $dir     = "/wormsrv2/autoace";
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

  &run_command("/bin/gzip $out");

  $runtime = &runtime;
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


  my $blat_dir  = "/wormsrv2/autoace/BLAT";
  my $blast_dir = "/wormsrv2/wormbase/ensembl_dumps";

  # does this to tidy up???? Not sure why these lines are here, krb
  &run_command("/bin/gzip -f $blast_dir/worm_pep_blastp.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_blastp.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_dna_blastx.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_pep_motif_info.ace");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_motif_info.ace");

  # compress best blast hits files and then copy to FTP site
  &run_command("/bin/gzip -f $blast_dir/worm_pep_best_blastp_hits");
  &run_command("/bin/gzip -f $blast_dir/worm_brigpep_best_blastp_hits");
  &run_command("scp $blast_dir/worm_pep_best_blastp_hits.gz      $targetdir/$WS_name/best_blastp_hits.$WS_name.gz");
  &run_command("scp $blast_dir/worm_brigpep_best_blastp_hits.gz  $targetdir/$WS_name/best_blastp_hits_brigpep.$WS_name.gz");


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
  $runtime = &runtime;
  $log->write_to("$runtime: running $command\n");
  my $status = system($command);
  if($status != 0){
    $log->error;
    $log->write_to("ERROR: $command failed\n");
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



#!/usr/local/bin/perl5.8.0 -w
#
# cgc_names_for_worm_genes.pl
#
# a script to be able to feedback CGC names to CDS, Transcript, and Pseudogenes
# in camace and stlace (so then they can go into EMBL/GenBank)
#
# written by Keith Bradnam (based on locus2seq.pl)
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-04-26 10:36:39 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;
 
##############################
# command-line options       #
##############################

my ($help, $debug,$update_camace);
my $maintainers = "All";

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "update_camace" => \$update_camace);

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


########################################
# Set up and initialise variables      #
########################################

our $log; # for error email to go to Sanger crew
our $stlouis_log; # for email to go to St. Louis
our $cam_out = "/wormsrv2/autoace/acefiles/CAM_cgc_names_for_genes.ace";
our $stl_out = "/wormsrv2/autoace/acefiles/STL_cgc_names_for_genes.ace";
our $all_out = "/wormsrv2/autoace/acefiles/ALL_cgc_names_for_genes.ace";
our $camace_dir = "/wormsrv1/camace";

# use /wormsrv2/geneace as this is frozen, safer to do this than use /wormsrv1/geneace
my $geneace_dir = "/wormsrv2/geneace";


################################
# 
#       Main subroutines
#
################################


# set up log files
&create_log_files;


# grab CGC_name->CDS,Transcript,Pseudogene using common data hash (info comes from /wormsrv2/geneace)
my %worm_gene2cgc = &FetchData('worm_gene2cgc_name');

# look up lab info for each gene to know how to split up the data for the output files
&get_lab_details;


# remove existing camace connections and replace with new ones
# First check for write access to $camace_dir
my $write_access = check_write_access($camace_dir);

&update_camace if (($update_camace) && ($write_access eq "yes")); 




################################
# 
#       Tidy up
#
################################

#copy the ace files to the FTP site

my $runtime = &runtime;
print LOG "\n$runtime: Gzipping files and moving across to ftp site\n";
system("gzip -f /$stl_out") && print LOG "ERROR: gzip failed on STL";
system("gzip -f /$cam_out") && print LOG "ERROR: gzip failed on CAM";
system("gzip -f /$all_out") && print LOG "ERROR: gzip failed on ALL";

system("mv $cam_out.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/") && print LOG "Couldn't move camace file to ftp site\n"; 
system("mv $stl_out.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/") && print LOG "Couldn't move stlouis ace file to ftp site\n"; 
system("mv $all_out.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/") && print LOG "Couldn't move all file to ftp site\n"; 


#inform any interested parties
$runtime = &runtime;
print LOG "\n$runtime: Sending ERROR email to Sanger and emailing St. Louis about updates\n";
close (LOG);


# send main error log email to Sanger
&mail_maintainer($0,"$maintainers",$log);

# Send reminder to St. Louis (unless in debug mode)
my $notify = "jspieth\@watson.wustl.edu,dblasiar\@watson.wustl.edu,krb\@sanger.ac.uk";
if($debug){
  &mail_maintainer($0,"$maintainers",$stlouis_log);
}
else{
  &mail_maintainer($0,"$notify",$stlouis_log);
}

# hasta luego

exit(0);






##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log          = "/wormsrv2/logs/$script_name.$rundate.$$";
  $stlouis_log  = "/wormsrv2/logs/$script_name.st_louis_data.txt";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";


  # Admittedly it's possibly overkill to write a log file each time for such invariant
  # information!

  open(STLOUIS, ">$stlouis_log") || die "Couldn't open $stlouis_log\n";
  print STLOUIS "Updated info linking CDSs, Transcripts, and Pseudogene to their corresponding Gene object\n"; 
  print STLOUIS "(containing CGC approved names) is available from\n";
  print STLOUIS "ftp-wormbase\/pub\/data\/updated_locus2seq\/ in the three files:\n";
  print STLOUIS "CAM_cgc_names_for_genes.ace: loci in Hinxton CDSs.\n";
  print STLOUIS "STL_cgc_names_for_genes.ace: loci in St Louis CDSs.\n";
  print STLOUIS "ALL_cgc_names_for_genes.ace: all loci.\n\n";
  print STLOUIS "The files contain ?Gene objects with approved cgc names and that connect to a\n";
  print STLOUIS "CDS, Transcript, or Pseudogene object.\n";
  close(STLOUIS);

}

########################################################

sub get_lab_details{  
  print LOG "Getting associated lab info from data hash:\n\n";  

  my %worm_gene2lab = &FetchData('worm_genes2lab');   
   
  open (CAMOUT,">$cam_out") || die "cant open CAMOUT";
  open (STLOUT,">$stl_out") || die "cant open STLOUT";
  open (ALLOUT,">$all_out") || die "cant open ALLOUT";


  my $CAMcount = 0;  
  my $STLcount = 0;  
  my $PROBcount = 0;  

  foreach my $worm_gene(keys %worm_gene2cgc){  

    my $data = $worm_gene2cgc{$worm_gene};
    my @details = split(/\s/,$data);  
    my $class    = $details[0];  
    my $cgc_name = $details[1];  
    my $gene_ID  = $details[2];  
    
    my $lab = $worm_gene2lab{$worm_gene};

    if(!defined($lab)){
      print LOG "ERROR: $class: $worm_gene ($cgc_name) has missing From_laboratory tag\n";  
      $PROBcount++;  
    }      
    elsif($lab eq "HX"){  
      print CAMOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
      $CAMcount++;  
      print ALLOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
    }  
    elsif($lab eq "RW"){  
      print STLOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
      $STLcount++;  
      print ALLOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
    }  
    elsif($lab eq "DRW"){  
      print CAMOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
      $CAMcount++;  
      print ALLOUT "Gene : \"$gene_ID\"\n$class\t\"$worm_gene\"\nCGC_name \"$cgc_name\"\n\n";  
    }  
    else{  
      print LOG "ERROR: $class: $worm_gene ($cgc_name) has unknown From_laboratory tag: $lab\n";  
      $PROBcount++;  
    }
 
  }  

  close (CAMOUT);
  close (STLOUT);
  close (ALLOUT);
  
  print LOG "\nFound $CAMcount CGC names on Hinxton CDSs/Transcripts/Pseudogenes.\n";  
  print LOG "Found $STLcount CGC names on StLouis CDSs/Transcripts/Pseudogenes.\n";  
  print LOG "Found $PROBcount CGC names that have problems (see above)\n\n";  
  
  print LOG "Wrote output ACE files to /wormsrv2/autoace/acefiles\n";  
  
}  
   
##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##############################################################
# update camace with new locus->CDS connectons
###############################################################

sub update_camace{


  # Delete existing locus connections in camace and then load new ones in
  my $runtime = &runtime;
  print LOG "\n$runtime: Starting to remove existing worm_genes->Gene connections in camace and replace with new ones\n";
  my $command;
  $command  = "query find worm_genes\n";
  $command .= "eedit -D Gene\n";
  $command .= "save\n";
  $command .= "pparse /wormsrv2/autoace/acefiles/CAM_cgc_names_for_genes.ace\n";
  $command .= "save\nquit\n";

  my $tace = &tace ;

  open (WRITEDB, "| $tace -tsuser cgc_names_for_worm_genes $camace_dir |") || die "Couldn't open pipe to $camace_dir";
  print WRITEDB $command;
  close WRITEDB;

  $runtime = &runtime;
  print LOG "$runtime: Finished updating camace\n";
}


__END__

=pod

=head2 NAME - cgc_names_for_worm_genes.pl

=head1 USAGE

=over 4

=item cgc_names_for_worm_genes.pl  [-options]

=back

This script gets a list of current Gene->CDS, Transcript, Pseudogene connections from a stored data
hash in /wormsrv2/autoace/COMMON_DATA.  It then makes a dump of these to the FTP site, making separate files 
for just St. Louis and Sanger CDSs and also a combined file.  Lab information comes from another data hash
file, such that this script doesn't need to interact with any database to get the data.  A -update_camace
option allows the information to be loaded back into camace, first removing existing ?worm_genes -> Gene
connections. This information will help keep stlace and camace synchronised with changes in geneace.  The
?Gene objects that are loaded to camace contain just the CGC_name field.

cgc_names_for_worm_genes.pl MANDATORY arguments:

=over 4

=item none

=back

cgc_names_for_worm_genes.pl  OPTIONAL arguments:

=over 4

=item -update_camace

This option will specify that existing Gene->CDS/Transcript/Pseudogene connections in camace should 
be removed and replaced with the new ones.

=item -help, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk) but based heavily on a script by Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

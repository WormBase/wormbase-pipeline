#!/usr/local/bin/perl5.8.0 -w
#
# locus2seq.pl
#
# written by Anthony Rogers (ar2@sanger.ac.uk)
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-10-31 15:32:20 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;
use Getopt::Long;
 
##############################
# command-line options       #
##############################

my ($help, $debug,$update,$camace,$geneace);
my $maintainers = "All";

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "geneace=s"  => \$geneace,
	    "camace"     => \$camace,
	    "update"     => \$update
    );


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

my (%seq_locus);
our $log; # for error email to go to Sanger crew
our $stlouis_log; # for email to go to St. Louis
our $cam_out = "/wormsrv2/autoace/acefiles/CAM_locus_seq.ace";
our $stl_out = "/wormsrv2/autoace/acefiles/STL_locus_seq.ace";
our $all_out = "/wormsrv2/autoace/acefiles/ALL_locus_seq.ace";
our $currentDB  = "/nfs/disk100/wormpub/DATABASES/current_DB";
our $camace_dir = "/wormsrv1/camace";
my $geneace_dir;

if ($geneace){
  $geneace_dir = $geneace;
}
else{
  $geneace_dir = "/wormsrv2/geneace";
}
print "\nUsing $geneace_dir as target geneace database\n";


################################
# 
#       Main subroutines
#
################################


# set up log files
&create_log_files;


# grab Locus->Sequence, Locus->Transcript, and Locus->Pseudogene
# connections from geneace
&get_sequence_connections;

# Compare this to current_DB to work out which genes are Sanger / St. Louis
&compare_with_current_DB;



# remove existing camace connections and replace with new ones
# First check for write access to $camace_dir
my $write_access = check_write_access($camace_dir);

&update_camace if ( ($camace) && ($update) && ($write_access eq "yes")); 




################################
# 
#       Tidy up
#
################################


close (CAMOUT);
close (STLOUT);
close (ALLOUT);


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
  print STLOUIS "Updated info linking Loci to Sequences, Transcripts, and Pseudogene is available from\n";
  print STLOUIS "ftp-wormbase\/pub\/data\/updated_locus2seq\/ in the three files:\n";
  print STLOUIS "CAM_locus_seq.ace: loci in Hinxton sequence.\n";
  print STLOUIS "STL_locus_seq.ace: loci in St Louis sequence.\n";
  print STLOUIS "ALL_locus_seq.ace: all loci.\n\n";
  print STLOUIS "These are loci with approved cgc names and that connect to a valid\n";
  print STLOUIS "Sequence (Predicted_gene), Transcript, or Pseudogene object.\n";
  close(STLOUIS);


  open (CAMOUT,">$cam_out") || die "cant open CAMOUT";
  open (STLOUT,">$stl_out") || die "cant open STLOUT";
  open (ALLOUT,">$all_out") || die "cant open ALLOUT";



}


#######################################################################
# Grab locus->sequence and Locus-> Transcript connections from geneace
#######################################################################

sub get_sequence_connections{
 
  print LOG "Getting Locus->Sequence/Transcript/Pseudogene connections from $geneace_dir:\n\n";

  # get locus with confirmed CGC names and the corresponding sequences, transcripts
  # and Pseudogenes.  This uses 3 table_maker queries exported from xace

  my $command1 = "Table-maker -p \"${geneace_dir}/wquery/locus_seq.def\"\nquit\n";
  my $command2 = "Table-maker -p \"${geneace_dir}/wquery/locus_transcript.def\"\nquit\n";
  my $command3 = "Table-maker -p \"${geneace_dir}/wquery/locus_pseudogene.def\"\nquit\n";

  my @commands = ($command1, $command2, $command3);
  my $type = "Sequence";

  # Counter for tracking which iteration you are in
  my $counter = 0;

  foreach my $command (@commands){

    open (GENEACE, "echo '$command' | tace $geneace_dir | ") || die "Couldn't open pipe to $geneace_dir\n";
    while (<GENEACE>){ 
      # skip acedb related lines, i.e. not actual data
      next if ((m/acedb>/) || (m/\/\//));

      my @entry = split(/\s+/,$_);

      if($entry[0] && $entry[1]){
	my $locus = $entry[0]; $locus =~ s/\"//g;
	my $sequence = $entry[1]; $sequence =~ s/\"//g;
	
	# Add to hash with appropriate type prefix
	$seq_locus{"Sequence:".$sequence}   .= "$locus " if ($type eq "Sequence");
	$seq_locus{"Transcript:".$sequence} .= "$locus " if ($type eq "Transcript");
	$seq_locus{"Pseudogene:".$sequence} .= "$locus " if ($type eq "Pseudogene");
       
      }
    }
    close(GENEACE);
    
    # Change type to Transcript/Pseudogene in 2nd/3rd run through loop
    $counter++;
    $type = "Transcript" if ($counter == 1);
    $type = "Pseudogene" if ($counter == 2);
  }
}

########################################################

sub compare_with_current_DB{  
  print LOG "Getting associated lab info from $currentDB:\n\n";  
   
  my $autoace = Ace->connect($currentDB) || die "cant open $currentDB\n";  
   
  my $CAMcount = 0;  
  my $STLcount = 0;  
  my $PROBcount = 0;  
   
  my $count;  
  foreach my $entry(keys %seq_locus){  
    
    # work out what to get (i.e. ?Sequence, ?Transcript, or ?Pseudogene)  
    my @details = split(/\:/,$entry);  
    my $class  = $details[0];  
    my $object = $details[1];  
    
    my $seq = $autoace->fetch($class => "$object");  
    
    if (defined($seq)){  
      my @lab = $seq->at('Origin.From_laboratory');  
      #extract any cases where a sequence contains two loci  
      my @loci = split(/\s+/,"$seq_locus{$entry}");  
      
      foreach my $locus (@loci){  
	if(defined($lab[0])){  
	  my $tag = $class;  
	  $tag = "Genomic_sequence" if ($class eq "Sequence");  
	  if($lab[0] eq "HX"){  
	    print CAMOUT "Locus : \"$locus\"\n$tag\t\"$object\"\n\n";  
	    $CAMcount++;  
	    print ALLOUT "Locus : \"$locus\"\n$tag\t\"$object\"\n\n";  
	  }  
	  elsif($lab[0] eq "RW"){  
	    print STLOUT "Locus : \"$locus\"\n$tag\t\"$object\"\n\n";  
	    $STLcount++;  
	    print ALLOUT "Locus : \"$locus\"\n$tag\t\"$object\"\n\n";  
	  }  
	  else{  
	    print LOG "ERROR: $class: $object ($locus) has unknown From_laboratory tag: $lab[0]\n";  
	    $PROBcount++;  
	  }  
	}  
	else{  
	  print LOG "ERROR: $class: $object ($locus) has no From_laboratory tag in $currentDB\n";  
	  $PROBcount++;  
	}  
      }  
    }  
    else{  
      print LOG "ERROR: $class: $object ($seq_locus{$entry}) not found in $currentDB\n";  
      $PROBcount++;  
    }  
  }  
  $autoace->close;  
  
  print LOG "\nFound $CAMcount loci on Hinxton sequences/transcripts/pseudogene.\n";  
  print LOG "Found $STLcount loci on StLouis sequences/transcripts/pseudogene.\n";  
  print LOG "Found $PROBcount loci that have problems (see above)\n\n";  
  
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
# update camace with new locus->Sequence connectons
###############################################################

sub update_camace{


  # Delete existing locus connections in camace and then load new ones in
  my $runtime = &runtime;
  print LOG "\n$runtime: Starting to remove existing locus->sequence connections in camace and replace with new ones\n";
  my $command;
  $command  = "query find Predicted_gene\n";
  $command .= "eedit -D Locus_genomic_seq\n";
  $command .= "query find Transcript\n";
  $command .= "eedit -D Locus\n";
  $command .= "query find Pseudogene\n";
  $command .= "eedit -D Locus\n";
  $command .= "save\n";
  $command .= "pparse /wormsrv2/autoace/acefiles/CAM_locus_seq.ace\n";
  $command .= "save\nquit\n";

  my $tace = &tace ;

  open (WRITEDB, "| $tace -tsuser locus2seq $camace_dir |") || die "Couldn't open pipe to $camace_dir";
  print WRITEDB $command;
  close WRITEDB;

  $runtime = &runtime;
  print LOG "$runtime: Finished updating camace\n";
}


__END__

=pod

=head2 NAME - locus2seq.pl

=head1 USAGE

=over 4

=item locus2seq.pl  [-options]

=back

This script makes a list of current Locus->Sequence (and Locus->Transcript + Locus->Pseudogene)
connections which are valid and makes a dump of these to the FTP site, making separate files 
for just St. Louis and Sanger sequences and also a combined file.  This information will help 
keep stlace and camace synchronised with changes in geneace.

locus2seq.pl MANDATORY arguments:

=over 4

=item none

=back

locus2seq.pl  OPTIONAL arguments:

=over 4

=item -geneace <path>, location of geneace database

By default this script will compare ~wormpub/DATABASES/current_DB to /wormsrv2/geneace.  The
-d flag allows you to compare against another copy of geneace (i.e. /wormsrv1/geneace). However
the Table maker definition files used by this script are only kept in /wormsrv2/geneace/wquery
at the moment.

=item -camace, update camace

This option will specify that existing Locus->Sequence/Transcript/Pseudogene connections should 
be removed and replaced with the new ones.

=item -help, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

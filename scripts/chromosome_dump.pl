#!/usr/local/bin/perl5.8.0 -w
#
# chromosome_dump.pl 
#
# by Keith Bradnam
#
# A script for dumping dna and/or gff files for chromosome objects in autoace
# see pod for more details
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-09-23 08:15:54 $      


use strict;
use lib '/wormsrv2/scripts/';
use Wormbase;
use Getopt::Long;
use IO::Handle;

$|=1;

######################################################
# Script variables and command-line options          #
######################################################

our $tace   = &tace;
our $giface = &giface;
our $maintainers = "All";
our ($log, $help, $debug, $dna, $gff, $zipdna, $zipgff, $composition, $database, $dump_dir);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "dna"         => \$dna,
	    "gff"         => \$gff,
	    "zipdna"      => \$zipdna,
	    "zipgff"      => \$zipgff,
	    "composition" => \$composition,
	    "database=s"  => \$database,
	    "dump_dir=s"  => \$dump_dir
	   );


##########################################################

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


# Sanity checks
if($database && !$dump_dir){
  die "You have specified a database (-database flag) but not a destination directory\nto dump to (-dump_dir flag).\n";
}
if(!$database && $dump_dir){
  die "You have specified a destination directory to dump to (-dump_dir flag) but not\na source database (-database flag).\n";
}
if(!$gff && !$dna && !$composition && !$zipgff && !$zipdna){
  die "No major option (-dna, -gff, -composition, -zipdna, or -zipgff) has been specified, try again\n";
}

# revert to defaults if no arguments specified
$database = "/wormsrv2/autoace"             if (!defined($database));
$dump_dir = "/wormsrv2/autoace/CHROMOSOMES" if (!defined($dump_dir));


#####################################################
# Main subroutines
#####################################################

&create_log_files;

&dump_dna    if ($dna);
&dump_gff    if ($gff);
&composition if ($composition);
&zip_files   if ($zipdna || $zipgff);

# say goodnight Barry
close(LOG);

exit(0);




#############################################################################
# Subroutines
#############################################################################


#########################
# dump dna files
#########################

sub dump_dna{
  my $command=<<END;
find sequence CHROMOSOME_I
dna -f $dump_dir/CHROMOSOME_I.dna
find sequence CHROMOSOME_II
dna -f $dump_dir/CHROMOSOME_II.dna
find sequence CHROMOSOME_III
dna -f $dump_dir/CHROMOSOME_III.dna
find sequence CHROMOSOME_IV
dna -f $dump_dir/CHROMOSOME_IV.dna
find sequence CHROMOSOME_V
dna -f $dump_dir/CHROMOSOME_V.dna
find sequence CHROMOSOME_X
dna -f $dump_dir/CHROMOSOME_X.dna
find sequence CHROMOSOME_MtDNA
dna -f $dump_dir/CHROMOSOME_MtDNA.dna
quit
END

  &execute_ace_command($command,$tace,$database);
  print LOG "Finished dumping DNA\n\n";
}



#########################
# dump gff files
#########################

sub dump_gff{
  my $command=<<END;
gif seqget CHROMOSOME_I ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_I.gff
gif seqget CHROMOSOME_II ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_II.gff
gif seqget CHROMOSOME_III ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_III.gff
gif seqget CHROMOSOME_IV ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_IV.gff
gif seqget CHROMOSOME_V ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_V.gff
gif seqget CHROMOSOME_X ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_X.gff
gif seqget CHROMOSOME_MtDNA ; seqfeatures -version 2 -file $dump_dir/CHROMOSOME_MtDNA.gff
quit
END

  &execute_ace_command($command,$giface,$database);
  print LOG "Finished dumping GFF files\n\n";
}


###################################
# produce dna composition files
###################################

sub composition{

  print LOG "Generating composition.all\n";	

  chdir $dump_dir;
  system("/bin/cat CHROMOSOME_I.dna CHROMOSOME_II.dna CHROMOSOME_III.dna CHROMOSOME_IV.dna CHROMOSOME_V.dna CHROMOSOME_X.dna | /nfs/disk100/wormpub/bin.ALPHA/composition > composition.all") && die "Couldn't create composition file\n";
  print LOG "Generating totals file\n";
  my $total = 0;
  my $final_total = 0;
  my $minus = 0;
  open(IN,"$dump_dir/composition.all") || die "Couldn't open composition.all\n";
  while(<IN>){
    if(/.*, (\d*) total/){
      $total = $1;
      next;
    }			
    if(/.*\- (\d*)/){
      $minus = $1; 
      last;
    }
  }
  close(IN);
  $final_total = $total - $minus;
  system("echo $total $final_total > totals") && die "Couldn't create totals file\n";
  
  &release_composition;
}

##########################
# zip up files
###########################

sub zip_files{
  foreach my $chr ("I", "II", "III", "IV", "V", "X", "MtDNA"){
    my $dna_file = "$dump_dir"."/CHROMOSOME_".$chr.".dna";
    my $gff_file = "$dump_dir"."/CHROMOSOME_".$chr.".gff";
    if ($zipdna){
      print LOG "Compressing $dna_file\n";
      system ("/bin/gzip -f $dna_file");
    }
    if ($zipgff){
      print LOG "Compressing $gff_file\n";
      system ("/bin/gzip -f $gff_file");
    }
  }	
}


#####################################################
# execute ace command: tace or giface
#####################################################

sub execute_ace_command {
  my ($command,$exec,$dir)=@_;
  open (WRITEDB,"| $exec $dir >> $log") or do {
    print LOG "execute_ace_command failed\n";
    close(LOG);
    die();
  };
  print WRITEDB $command;
  close (WRITEDB);
}

######################################################


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

#####################################################

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

=head1 NAME - chromosome_dump.pl

=head2 USAGE

chromosome_dump.pl is a replacement script for the two existing
shell scripts: chrom_dump_3.0 and gff_dump.  This script can dump
chromosome-length DNA sequences for entire chromsomes in the autoace
database.  It can additionally generate chromosome GFF files, and
finally it can compress these files using gzip.

A log file is written to /wormsrv2/logs/

All dumped files are written to /wormsrv2/autoace/CHROMOSOMES/

chromosome_dump.pl arguments:

=over 4

=item -dna

Dump dna files from specified database (see -p flag), dumps one file for each of 
the six nuclear chromosomes, plus one file for the mitochondrial chromosome.

=back


=over 4

=item -gff

Dump gff files, dumps one file for each chromosome in the database.

=back


=over 4

=item -composition (optional)

Calculates composition statistics for any dna files that are generated.
Use in combination with -d option (see above).  Excludes mitochondrial 
chromosome.

=back

=item -database <database>

Specify database that you wish to dump dna/gff files from.  If -database is not specified
the script will dump from /wormsrv2/autoace by default

=back


=over 4

=item -dump_dir <destination directory for dump files>

Specify destination of the dna and/or gff dump files generated from the -dna or -gff options.
If -dump_dir is not specified, dump files will be written to /wormsrv2/autoace/CHROMOSOMES by default

=back


=over 4

=item -zipdna (optional)

Compresses any dna files using gzip (will remove any existing files first).  The -composition
option will be run before this stage if -composition is specified.

=back

=over 4

=item -zipgff (optional) 

Compresses any gff files using gzip (will remove any existing files first).      

=back

=over 4

=item -help

Show these help files.

=back


=over 4

=item -debug <user>

Only email log file to specified user

=back


=over 4

=over 4

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

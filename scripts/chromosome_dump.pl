#!/usr/local/bin/perl5.6.0 -w

# chromosome_dump.pl 
# by Keith Bradnam aged 12 and a half,  10/08/01
#
# A script for dumping dna and/or gff files for chromosome objects in autoace
# see pod for more details

use strict;
use Getopt::Std;
use lib '/wormsrv2/scripts/';
use Wormbase;

$|=1;


##############################
# Script variables           #
##############################

my $cvs_version = &get_cvs_version("$0");

our $tace   = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";
our $giface = "/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/giface";
our $database = "/wormsrv2/autoace";
our $dump_dir = "/wormsrv2/autoace/CHROMOSOMES";
our ($opt_d,$opt_g,$opt_e,$opt_h);
getopts("dgeh");


#############################
# display help if required  #
#############################

&show_help if (!$opt_d && !$opt_e && !$opt_g && !$opt_h);

##################################################
# Open logfile                                   #
##################################################

my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

our $logfile = "/wormsrv2/logs/chromosome_dump.${rundate}.$$";

open (LOGFILE,">$logfile") || die "Couldn't create $logfile\n";
print LOGFILE "# chromosome_dump.pl\n\n";     
print LOGFILE "# version        : $cvs_version\n";
print LOGFILE "# run details    : $rundate $runtime\n";
print LOGFILE "\n\n";


#####################################################
# Main three subroutines
#####################################################

&dump_dna  if ($opt_d);
&dump_gff  if ($opt_g);
&zip_files if ($opt_e || $opt_h);


# say goodnight Barry

close(LOGFILE);

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
quit
END

  &execute_ace_command($command,$tace,$database);
  print LOGFILE "Finished dumping DNA\n\n";
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
quit
END

  &execute_ace_command($command,$giface,$database);
  print LOGFILE "Finished dumping GFF files\n\n";
}

##########################
# zip up files
###########################

sub zip_files{
	foreach my $chr ("I", "II", "III", "IV", "V", "X"){
		my $dna_file = "$dump_dir"."/CHROMOSOME_".$chr.".dna";
		my $gff_file = "$dump_dir"."/CHROMOSOME_".$chr.".gff";
		if ($opt_e){
			if (-e $dna_file."gz"){
				print LOGFILE "Removing existing *.dna.gz file\n";
				system ("rm -f $dna_file") && die "Couldn't remove files\n";
			}
			print LOGFILE "Compressing $dna_file\n";
			system ("/bin/gzip $dna_file") if ($opt_e);
		}
		if ($opt_h){
 			if (-e $gff_file."gz"){
                        	print LOGFILE "Removing existing *.gff.gz file\n";
                        	system ("rm -f $gff_file") && die "Couldn't remove files\n";
	                }
			print LOGFILE "Compressing $gff_file\n";
                        system ("/bin/gzip $gff_file") if ($opt_h);
		}
	}	
}


#####################################################
# execute ace command: tace or giface
#####################################################

sub execute_ace_command {
  my ($command,$exec,$dir)=@_;
  open (WRITEDB,"| $exec $dir >> $logfile") or do {
    print LOGFILE "execute_ace_command failed\n";
    close(LOGFILE); 
    die();
  };
  print WRITEDB $command;
  close (WRITEDB);
}

######################################################

sub show_help {
  system ('perldoc',$0) && die "Couldn't execute perldoc\n";
  exit(0);
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

=item *

<-d> dumps one dna file for each of the six chromosomes.

=back

=over 4

=item *

<-g> dumps one gff file for each of the six chromosomes

=back

=over 4

=item *

<-e> (optional) compresses any dna files using gzip (will remove 
any existing files first).

=back

=item *
                        
<-h> (optional) compresses any gff files using gzip (will remove        
any existing files first).      

=back

=over 4

=item *

<-h> help page (what you are reading now).

=back


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

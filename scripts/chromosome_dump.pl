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
our $database = "/wormsrv2/autoace";
our $dump_dir = "/wormsrv2/autoace/CHROMOSOMES/";
our ($opt_d,$opt_g,$opt_z,$opt_h);
getopts("dgzh");


#############################
# display help if required  #
#############################

&show_help if ($opt_h);
&show_help if (@ARGV = "");



##################################################
# Open logfile                                   #
##################################################

my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

our $logfile = "/wormsrv2/logs/chromosome_dump.${rundate}.$$";
system ("/bin/touch $logfile") || die "Couldn't create $logfile\n";
open (LOGFILE,">>$logfile") || die "Couldn't append to $logfile\n";

print LOGFILE "# chromosome_dump.pl\n\n";     
print LOGFILE "# version        : $cvs_version\n";
print LOGFILE "# run details    : $rundate $runtime\n";
print LOGFILE "\n\n";


#####################################################
# Main three subroutines
#####################################################

&dump_dna  if ($opt_d);
&dump_gff  if ($opt_g);
&zip_files if ($opt_z);

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

  &execute_tace_command($command,$tace,$database);
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

  &execute_tace_command($command,$tace,$database);
}


#####################################################
# execute tace command
#####################################################

sub execute_tace_command {
  my ($command,$exec,$dir)=@_;
  open (WRITEDB,"| $exec $dir >> $logfile") or do {
    print LOGFILE "execute_tace_command failed\n";
    close(LOGFILE); 
    die();
  };
  print WRITEDB $command;
  close (WRITEDB);
}


sub show_help {
  system ('perldoc',$0);
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

-d dumps one dna file for each of the six chromosomes.

=back

=over 4

=item *

-g dumps one gff file for each of the six chromosomes

=back

=over 4

=item *

-z (optional) compresses all dna and gff files using gzip (will remove 
any existing files first).

=back

=over 4

=item *

-h help page (what you are reading now).

=back


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

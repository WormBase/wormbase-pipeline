#!/usr/local/bin/perl5.6.0 -w

# chromosome_dump.pl 
# by Keith Bradnam aged 12 and a half,  10/08/01
#
# A script for dumping dna and/or gff files for chromosome objects in autoace
# see pod for more details
#
# v1.13 :  dl : Added a chdir command to the -c option. This ensures that you move to the dump directory
#               prior to calculating composition/totals


use strict;
use Getopt::Std;
use IO::Handle;
use lib '/wormsrv2/scripts/';
use Wormbase;

$|=1;


##############################
# Script variables           #
##############################

my $cvs_version = &get_cvs_version("$0");

our $tace   = "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/tace";
#our $giface = "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/giface";
#our $giface = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/giface";
our $giface = "/nfs/disk100/acedb/RELEASE.2002_07_28/bin.ALPHA_4/giface";


our ($opt_d,$opt_g,$opt_e,$opt_h,$opt_c, $opt_p, $opt_q, $opt_t);
getopts("dgehcp:q:t");
our $database;
our $dump_dir;




if($opt_p && !$opt_q){
  die "You have specified a database (-p flag) but not a destination directory\nto dump to (-q flag).\n";
}
if(!$opt_p && $opt_q){
  die "You have specified a destination directory to dump to (-q flag) but not\na source database (-p flag).\n";
}


if (defined($opt_p)){
  $database = $opt_p;
}
else{
  $database = "/wormsrv2/autoace";
}
    
if (defined($opt_q)){
  $dump_dir = $opt_q;
}
else{
  $dump_dir = "/wormsrv2/autoace/CHROMOSOMES";
}


#############################
# display help if required  #
#############################

&show_help if ((!$opt_d && !$opt_e && !$opt_g && !$opt_c && !$opt_t) || $opt_h);




##################################################
# Open logfile                                   #
##################################################

my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

our $logfile = "/wormsrv2/logs/chromosome_dump.${rundate}.$$";

open (LOGFILE,">$logfile") || die "Couldn't create $logfile\n";
LOGFILE->autoflush();
print LOGFILE "# chromosome_dump.pl\n\n";     
print LOGFILE "# version        : $cvs_version\n";
print LOGFILE "# run details    : $rundate $runtime\n";
print LOGFILE "\n\n";


#####################################################
# Main three subroutines
#####################################################

&dump_dna    if ($opt_d);
&dump_gff    if ($opt_g);
&composition if ($opt_c);
&zip_files   if ($opt_e || $opt_t);

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


###################################
# produce dna composition files
###################################

sub composition{

  print LOGFILE "Generating composition.all\n";	

  chdir $dump_dir;
  system("/bin/cat *.dna | /nfs/disk100/wormpub/bin.ALPHA/composition > composition.all") && die "Couldn't create composition file\n";
  print LOGFILE "Generating totals file\n";
  my ($total, $minus, $final_total);
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
  foreach my $chr ("I", "II", "III", "IV", "V", "X"){
    my $dna_file = "$dump_dir"."/CHROMOSOME_".$chr.".dna";
    my $gff_file = "$dump_dir"."/CHROMOSOME_".$chr.".gff";
    if ($opt_e){
      print LOGFILE "Compressing $dna_file\n";
      system ("/bin/gzip -f $dna_file") if ($opt_e);
    }
    if ($opt_t){
      print LOGFILE "Compressing $gff_file\n";
      system ("/bin/gzip -f $gff_file") if ($opt_t);
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

=item -d

Dump dna files from specified database (see -p flag), dumps one file for each of the six chromosomes.

=back


=over 4

=item -c (optional)

Calculates composition statistics for any dna files that are generated.
Use in combination with -d option (see above)

=back


=over 4

=item -h

Show these help files.

=back

=over 4

=item -g

Dump gff files, dumps one file for each chromosome in the database.

=back


=over 4

=item -p <database>

Specify database that you wish to dump dna/gff files from.  If -p is not specified
the script will dump from /wormsrv2/autoace by default

=back


=over 4

=item -q <destination directory for dump files>

Specify destination of the dna and/or gff dump files generated from the -d or -g options.
If -q is not specified, dump files will be written to /wormsrv2/autoace/CHROMOSOMES by default

=back



=over 4

=item -e (optional)

Compresses any dna files using gzip (will remove any existing files first).  The -c
option will be run before this stage if -c is specified.

=back


=over 4

=item -t (optional) 

Compresses any gff files using gzip (will remove any existing files first).      

=back


=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut

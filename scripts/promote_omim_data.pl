#!/software/bin/perl -w
#
# promote_omim_data.pl                           
# 
# by Paul Davis                         
#
# This script promoted the OMIM disease data to the level of the gene.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2012-09-21 08:13:07 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $noload);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noload"     => \$noload,
	   );

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

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################
my ($ace_object,%objects,%acedata);
my ($gene,$database,$id,$type);
my $outfile = $ace_dir."/acefiles/omim_db_data.ace";
my $count = "0";
my $countgenes = "0";
open (OUT, ">$outfile") or $log->log_and_die("Can't write ace file $outfile");

if (-e "$basedir/autoace/wquery/SCRIPT:omim.def") {
  my $def = "$basedir/autoace/wquery/SCRIPT:omim.def";
  my $command = "Table-maker -p $def\nquit\n";
  $log->write_to("\nRetrieving OMIM, using Table-maker...\n");
  
  open (TACE, "echo '$command' | $tace $ace_dir | ") || die "Cannot query acedb. $command  $tace\n";
  while (<TACE>) {
    $count++;
    chomp;
    s/\"//g;
    unless (/(WBGene\d+)\s+(OMIM)\s+(\S+)\s+(\d+)/) {next;}
    if (/(WBGene\d+)\s+(OMIM)\s+(\S+)\s+(\d+)/){
      $gene = $1;
      $database = $2;
      $type = $3;
      $id = $4;
      $ace_object = $1;
      if ($type eq "gene") {
	$objects{$ace_object}++;
	push (@{$acedata{$ace_object}{gene}},$id);
      }
      elsif ($type eq "disease") {
	$objects{$ace_object}++;
	push (@{$acedata{$ace_object}{disease}},$id);
      }
    }
  }
  foreach my $obj (keys %objects) {
    print OUT "\n// $obj\n\n";
    my $line;
    print OUT "Gene : \"$obj\"\n";
    foreach $line (@{$acedata{$obj}{gene}}) {
      print OUT "Database OMIM gene $line\n";
    }
    foreach $line (@{$acedata{$obj}{disease}}) {
      print OUT "Database OMIM disease $line\n";
    }
  }
  $countgenes = (keys %objects);
}
else {die "Cannot query acedb. as $basedir/autoace/wquery/SCRIPT:omim.def does not exist\n";}
close OUT;

if ($noload){
  $log->write_to("Output NOT loaded into ".$wormbase->autoace."\n");
}

unless ($noload) {
  $log->write_to("loading $outfile to ".$wormbase->autoace."\n");
  $wormbase->load_to_database($wormbase->autoace,$outfile,'promote_omim_data.pl', $log);
}

# Close log files and exit
$log->write_to("\nOutput can be found here $outfile\n");
$log->write_to("----------\n");
$log->write_to("Processed $count lines from Table-maker touching $countgenes genes\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

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

=head2 NAME - promote_omim_data.pl

=head1 USAGE

=over 4

=item promote_omim_data.pl  [-options]

=back

This script queries for C. elegans genes that have a Ortholog_other that is an ensembl protein with attached OMIM disease IDs and promotes them to the Database line of the elegans gene.

promote_omim_data.pl MANDATORY arguments:

=over 4

=item None at present.

=back

promote_omim_data.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script in the test build env.

=back

=over 4
    
=item -verbose, prints to screen additional information.

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@ebi.ac.uk)

=back

=cut

#!/usr/bin/env perl 
#
# get_GI.pl
# 
# A script to generate C.elegans CDS/Protein to NCBI_GI number XREFs
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-10-14 09:54:27 $
#
#==================

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Net::FTP;

##############################
# command-line options       #
##############################

my ($test,$debug,$wormbase,$store,$noupdate,$noload,$database, $acefile);

GetOptions (
  "debug:s"     => \$debug,
  "store=s"     => \$store,
  "test"        => \$test, 
  "database:s"  => \$database,
  "noupdate"    => \$noupdate,
  "noload"      => \$noload,
  "acefile=s"   => \$acefile,
	   );

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                           );
}

my $log = Log_files->make_build_log($wormbase);

my $NCBI_GI_FTP_HOST = "ftp.ncbi.nih.gov";
my $NCBI_GI_FTP_DIR = "genbank/livelists";
my $NCBI_GI_FTP_FILE_PREFIX = "GbAccList";
my $LOCAL_DIR =  $wormbase->misc_static . "/XREFS/NCBI_GI";

$database = $wormbase->autoace if not defined $database;
$acefile = $wormbase->acefiles . "/CDS_GI_numbers.ace" if not defined $acefile;

my ($current_file) = glob("$LOCAL_DIR/GbAccList.*");

unless ($noupdate) {
  my $ftp = Net::FTP->new($NCBI_GI_FTP_HOST, Debug => 0, Timeout => 600) 
      or $log->log_and_die("Cannot connect to $NCBI_GI_FTP_HOST: $@");
  $ftp->login("anonymous",'-anonymous@') 
      or $log->log_and_die("Cannot login ", $ftp->message);
  $ftp->cwd($NCBI_GI_FTP_DIR);

  my ($latest_file_on_ftp, $latest_file_on_ftp_date);

  foreach my $file ($ftp->ls()) {
    if ($file =~ /$NCBI_GI_FTP_FILE_PREFIX.(\d{2})(\d{2})\.(\d{4})/) {
      my $file_date = "$3$1$2";
      if (not defined $latest_file_on_ftp_date or $file_date > $latest_file_on_ftp_date) {
        $latest_file_on_ftp = $file;
        $latest_file_on_ftp_date = $file_date;
      }
    }
  }
  if (defined $latest_file_on_ftp) {
    my $download = 0;

    if (not defined $current_file) {
     $log->write_to("No current $NCBI_GI_FTP_FILE_PREFIX file - will download\n");
      $download = 1;
    } elsif ($current_file =~ /$NCBI_GI_FTP_FILE_PREFIX\.(\d{2})(\d{2})\.(\d{4})/) {
      my $current_file_date = "$3$1$2";
      if ($current_file_date >= $latest_file_on_ftp_date) {
        $log->write_to("Current GbAcclist file is up-to-date, will not download\n");
      } else {
        $log->write_to("Current GbAcclist file is out-of-date, will download\n");
        $download = 1;
      }
    } else {
      $log->write_to("Could not get data info from current file - will download\n");
    }

    if ($download) {
      $ftp->binary();
      $ftp->get($latest_file_on_ftp, "$LOCAL_DIR/$latest_file_on_ftp") 
          or $log->log_and_die("Could not successfully download $latest_file_on_ftp\n");
      unlink $current_file if defined $current_file and -e $current_file;
      $current_file = "$LOCAL_DIR/$latest_file_on_ftp";
    }

    $ftp->quit;
  } else {
    die "Cannot find anyfiles on FTP site. Stopping\n";
  }
}

$log->write_to("Fetching protein_id info from CDSs...\n");
my (%cds_info);

my $def = &get_tm_def();
my $tm_query = $wormbase->table_maker_query($database,$def);
while(<$tm_query>) {
  chomp;
  next if (/>/ or /\/\//);
  s/\"//g;

  my($cds,$pid,$pidver) = split(/\t/, $_);
  next unless $pidver;
  my $hkey = "${pid}.${pidver}";
  $cds_info{$hkey} = $cds;
}
unlink $def;

$log->write_to("Reading GbAcclist...\n");
my $entries = 0;

open(my $out_fh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
open(my $fh, "gunzip -c $current_file |") or $log->log_and_die("Could not open gunzip stream to $current_file\n");
while(<$fh>) {
  chomp;
  my ($pid, $pidver, $gi) = split(/,/, $_);

  my $hkey = "${pid}.${pidver}";
  if (exists $cds_info{$hkey}) {
    my $cds = $cds_info{$hkey};
    print $out_fh "\nCDS : \"$cds\"\n";
    print $out_fh "Database NCBI_protein GI_number \"GI:$gi\"\n";
  }

  $entries++;
  if ($entries % 50000000 == 0) {
    $log->write_to("Read $entries entries...\n");
  }
}

if (not $noload) {
  $wormbase->load_to_database($database,$acefile,'get_GI', $log);
}

$log->mail();
exit(0);

##########################################
sub get_tm_def {

  my $tm_def = "/tmp/get_GI.$$.tm_def";
  open(my $outfh, ">$tm_def") 
      or $log->log_and_die("Could not open $tm_def for writing\n");

  my $query =  <<"EOF";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class CDS 
From 1 
 
Colonne 2 
Width 12 
Mandatory
Hidden 
Class 
Class Sequence 
From 1 
Tag Protein_id             
 
Colonne 3 
Width 12 
Optional 
Visible 
Text 
Right_of 2 
Tag HERE              
 
Colonne 4 
Width 12 
Optional 
Visible 
Integer 
Right_of 3 
Tag HERE              

EOF

  print $outfh $query;
  close($outfh);

  return $tm_def;
}


__END__

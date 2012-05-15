#!/software/bin/perl -w
#
# get_GI.pl
# 
# A script to generate C.elegans CDS/Protein to NCBI_GI number XREFs
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2012-05-15 12:12:01 $
#
#==================

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Net::FTP;

##############################
# command-line options       #
##############################

my ($test,$debug,$wormbase,$store,$noftp,$inputfile,$noclean,$database,$usebins,$load,$noupdate);

GetOptions (
	    "debug:s"     => \$debug,     # debug emails only the specified user
	    "store"       => \$store,     # provide a store object
	    "test"        => \$test,      # use the WS666 test env.
	    "noftp"       => \$noftp,     # noftp testing code that allows one to just retrieve a local backup of ncbi file.
	    "file:s"      => \$inputfile, # 
	    "noclean"     => \$noclean,   # don't remove the .dat files, useful for debugging.
	    "usebins:s"   => \$usebins,   # debugging option to use bins that are already in place, specify the number of bins.
	    "database:s"  => \$database,  # a way of specifying the database to use whenm testing or not using the build.
	    "load"        => \$load,      # Load the data into autoace.
	    "noupdate"    => \$noupdate,  # don't update the file links at the end.
	   );

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                           );
}

my $WS_version = $wormbase->get_wormbase_version;
my $wget = "/usr/bin/wget --force-html --no-clobber";

my $log = Log_files->make_build_log($wormbase);

if (defined $database){
  $log->write_to ("Getting new GI numbers from NCBI based on $database\n");
}

else {
  $log->write_to ("Getting new GI numbers for $WS_version\n");
}

# basic paths etc.
#ebi my $basedir = "/net/isilon3/production/panda/ensemblgenomes/wormbase/analysis/GI_numbers/";
my $basedir = $wormbase->wormpub."/analysis/GI_numbers/";
my $ncbiftp = "ftp.ncbi.nih.gov";
my $genbankftp = "ftp://ftp.ncbi.nih.gov/genbank/livelists/";
my $path = "/genbank/livelists/";
my $regex = "GbAccList*";

my $counter;
my $newfile;
my $GI_local;
my $GI_zip;
my $lock_file = $basedir.'PROCESS_LOCK';

#Main Body of Script
# Setup database path if one hasn't been supplied.
unless (defined $database) {
  $database = $wormbase->autoace;
}

my $datfile;
unless ($usebins) { # skip all of the ftp and GI .dat file generation for debugging buit you need the .dat files.
  if ($noftp){
    #testing code to maintain a copy of the file on disk.
    $GI_local = $basedir."GbAccList.1225.2011";
    $GI_zip = $GI_local.".gz";
    $wormbase->run_command("cp ${GI_zip}_bk $GI_zip", $log);
  }

  # Establish ftp connection and get newest file.
  unless ($noftp) {
    my $ftp;
    $ftp = Net::FTP->new("$ncbiftp", Debug => 0) or die "Cannot connect to $ncbiftp: $@";
    $ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;
    $ftp->cwd($path) or die "Cannot change working directory ", $ftp->message;
    my @ftpfiles = $ftp->dir($regex);
    $ftp->binary();
    $newfile = (split(/\s+/,$ftpfiles[-1]))[-1];
    $GI_zip = $basedir.$newfile;
    if ($GI_zip =~ /(\S+)\.gz/) { $GI_local = $1;}
    $ftp->quit;

    # Do we actually need to do anythig?
    if (-e $basedir.$newfile) {
      $log->write_to ("The most up to date file already exists on disk\n");
      print "The most up to date file already exists on file\n";
      if (!-e $lock_file) {
	$log->write_to ("exiting as no work is needed\n");
	print "exiting as no work is needed\n";
	$log->mail();
	exit;
      }
      else {
	$log->write_to ("The most up to date file already exists but it appears the last processing failed to finish, continue\n");
	print "The most up to date file already exists but it appears the last processing failed to finish, continue\n";
	
      }
    }
    else {
      # Using wget as NET:FTP was causing the end of file to be missing??
      $wormbase->run_command("wget -q -O $GI_zip ${genbankftp}$newfile", $log) && die "Couldnt get $newfile $0\n";
    }
  } #unless noftp
  
  #touch the lock file as you are starting to process data.
  $wormbase->run_command("touch $lock_file");
  #Unzip the data file
  $wormbase->run_command("gunzip -f $GI_zip", $log) && die "unzip Failed for ${GI_zip}\n";

  $counter  = 1;
  my $count    = 0;
  my $bin_size = 2500000;
  my $proteinID;

  # Generate 1st bin and open it.
  $datfile = $basedir."GI_number_${counter}.dat";


  open (OUTPUT, ">$datfile");
  print OUTPUT "%GI_number = (\n";
  
  #open the NCBI GI file and then create the .dat file for later processing
  $log->write_to("Opening file ".$GI_local."\n");
  open (GI,"<$GI_local") or die "cant open $GI_local\n";
  $log->write_to("Creating $datfile...\n");

  while (<GI>) {
    chomp;
    #AACY024124488,1,129566182
    my @f = split /\,/;
    $proteinID = $f[0] . "." . $f[1];
    $count ++;
    
    # if the bin is full close the lid and open a new one.
    if ($count == $bin_size) {
      print OUTPUT "\t'$proteinID' => '$f[2]',\n";
      print OUTPUT "             );\n";
      close OUTPUT;
      $counter++;
      $count = 0;
      $datfile = $basedir."GI_number_${counter}.dat";
      $log->write_to("Creating $datfile...\n");
      open (OUTPUT, ">$datfile");
      print OUTPUT "%GI_number = (\n";
    }
    # if the bin is part full just write out the line.
    else {
      print OUTPUT "\t'$proteinID' => '$f[2]',\n";
    }
  }
  # write the final line to the final file.
  print OUTPUT "             );\n";
  # Tidy up this phase.
  close OUTPUT;
  close GI;
  $log->write_to("\nCreated $counter bins\n\n");
} #unless ($usebins)

elsif ($usebins){
  $counter = $usebins;
}

#generate the .ace output file by parsing protein/CDS data and pulling IDs from the .dat files.
$log->write_to("Generate the .ace output file\n\n");
my %GI_number = ();
my $ace = ${basedir}."GI_numbers_WS".${WS_version}.".ace"; #debugging code

my %cds_info;
my $def = $wormbase->basedir."/wquery/SCRIPT:make_wormpep.def";
print STDERR "parsing $def\n";
my $tm_query = $wormbase->table_maker_query($database,$def);
while(<$tm_query>) {
  next if (/>/ or /\/\//);
  s/\"//g;#"
  s/\n//g;
  my($cds,$pid,$pidver,$unip,$unip_ac,$gene,$cgc,$status,$brief) = split(/\t/,$_);
  next unless ($cds and $gene);
  $cds_info{$cds}->{'gene'} = $gene;
  $cds_info{$cds}->{'pid'} = "$pid.$pidver" if ($pid and $pidver);
  $cds_info{$cds}->{'unip'} = $unip_ac if $unip_ac;
  $cds_info{$cds}->{'cgc'} = $cgc if $cgc;
}


open (ACE, ">$ace");
my $i;
my $processed;
my $matchID;
my $xrefs;
my $noGI;
my $datcheck;
my $proteinID;
my $noID;

for ($i = 1; $i <$counter+1; $i++) {
  $processed = "";
  # dat file name generation
  $datfile = $basedir."GI_number_${i}.dat";
  print "Processing $datfile\n";
  $datcheck++;
  last unless (-e $datfile);
  &recreate_hashes;
  
  foreach my $cds (sort keys %cds_info) {
    $processed ++;
    $proteinID = $cds_info{$cds}->{'pid'};
    if (!defined $proteinID) {
      $noID ++;
      next;
    }
    else { 
      $matchID ++;
      #CAA16341,2,26985906
      if (defined $GI_number{$proteinID}) {
	#'CAA16341.2' => 'K01G5.9',
	print ACE "CDS : \"$cds\"\n";
	print ACE "Database NCBI_protein  GI_number GI:$GI_number{$proteinID}\n\n";
	$xrefs ++;
      }
      else {
	#capture an idea of the number of lines where no GI number could be found.
	$noGI ++;
      }
    }
  }
  print "Processed $datfile\n";
}
$log->write_to("Processed $datcheck files\n");
#$log->write_to("\nResults:
#------------------------------------------
#Processed\t\t${processed} CDS objects
#Protein ID \t\t$matchID ProteinID::CDS matches
#Protein ID\t\t$noID CDSs with no protein ID
#GI XREFs\t\t$xrefs CDS/Proteins with a GI mapping
#GI XREFs\t\t$noGI CDS/Proteins with NO GI mapping
#------------------------------------------\n");

# Remove the symbolic link and update to point to new file.
unless ($noupdate) {
  $wormbase->run_command("rm -f GI_numbers.ace", $log);
  $wormbase->run_command("ln -s $ace GI_numbers.ace", $log);
}

# load the data into autoace
if ($load) {
$wormbase->run_script('autoace_builder.pl -load $ace -tsuser "gi_number"', $log);
}

# debugging code to leave .dat files on disk
if ($noclean) {
  $log->write_to("Leaving all .dat files for future debugging\n");
}
else {
  #clean up .dat files and the large NCBI file
  
  $log->write_to("Removing all .dat files from $basedir\n");
  $wormbase->run_command("gzip -9 -f $GI_local", $log) && die "zip Failed for ${GI_local}\n";
  $wormbase->run_command("rm $basedir/GI*.dat", $log);
}

$wormbase->run_command("rm $lock_file", $log);

# subroutines

sub recreate_hashes {
  my $data;
  open (FH, "<$datfile") or die "$datfile : $!";
  undef $/;
  $data = <FH>;
  eval $data;
  die if $@;
  $/ = "\n";
  close FH;
}



$log->mail();
exit(0);






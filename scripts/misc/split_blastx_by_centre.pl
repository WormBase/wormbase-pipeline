#!/usr/bin/perl5.8.0 -w

use strict;
my $scriptdir =  $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use Wormbase;
use File::Basename;
use Log_files;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($debug, $version,$wormbase,$test,$update,$noprpcess,$store,$files,$noprocess);

GetOptions(
	   'test'          => \$test,       #test build
           'debug:s'       => \$debug,
	   'version=s'     => \$version,
	   'wormbase'      => \$wormbase,
	   'files:s'       => \$files,      # stores a list of comma seperated filesnames for processing.
	   'update'        => \$update,     # Load the camace file into the canonical database.
	   'noprocess'     => \$noprocess,  # Don't process the files.
	   "store:s"       => \$store,
	  );



if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                             );
}

if (!$version) {
  $version = "666";
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

#establish global Scalars/arrays and hashes
my $next_build = ($version + 1);

#my $blast_dir = $wormbase->acefiles;
my $blast_dir = "/acari/work2a/wormpipe/dumps/blastx";
my $temp_dir = "/nfs/disk100/wormpub/camace_orig/WS$version-WS$next_build/tmp";
$wormbase->run_command("mkdir $temp_dir", $log);
my $outdir = "/nfs/disk100/wormpub/camace_orig/WS$version-WS$next_build";
print "/nfs/disk100/wormpub/camace_orig/WS$version-WS$next_build\n\n\n" if ($test);

#copy old COMMON_DATA file over to COMMON_DATA dir
$wormbase->run_command("cp " .$wormbase->misc_static."/clone2centre.dat " .$wormbase->common_data."/", $log);
print "cp " .$wormbase->misc_static."/clone2centre.dat " .$wormbase->common_data."/\n\n" if ($test);

my %clone2centre = $wormbase->FetchData("clone2centre") or $log->log_and_die ("ERROR: Cannot find COMMON_DATA/clone2centre.dat\n");
my $out;
# tace executable path
my $tace = $wormbase->tace;

my $stl_out;
my $cam_out;

if (!$noprocess){
  my @files2split;
  if ($files) {
    @files2split = split (/\,/,$files);
  } else {

#TEC_RED
    push (@files2split,"briggsae_blastx.ace");
    push (@files2split,"fly_blastx.ace");
    push (@files2split,"human_blastx.ace");
    push (@files2split,"remanei_blastx.ace");
    push (@files2split,"slimswissprot_blastx.ace");
    push (@files2split,"slimtrembl_blastx.ace");
    push (@files2split,"worm_blastx.ace");
    push (@files2split,"yeast_blastx.ace");
#RepeatMasker
#waba
#TRF.ace missing from WS162
  }
  open ($stl_out,">$outdir/STL_blastx.ace") || die "ERROR Can\'t open STL output $outdir/STL_blastx.ace\n";
  open ($cam_out,">$outdir/CAM_blastx.ace") || die "ERROR Can\'t open CAM output $outdir/CAM_blastx.ace\n";
  $log->write_to("\t\tPROCESSING DATA\n\t\t================================================================\n");
  foreach my $file ( @files2split ){
    $wormbase->run_command("scp ecs4:$blast_dir/$file $temp_dir/", $log);
    my $clone = "";
    my $type = "";
    $log->write_to("\nProcessing $blast_dir/$file\n");
    if ($file =~ (/(\S+)_blastx.ace/)){
      $type = $1;
#print "$type\n\n";
    }
    my $blast_file;
    open ($blast_file, "<$temp_dir/$file" ) or $log->log_and_die ("can\'t open input file $file\t$!\n");
    #Sequence : "cTel33B"
    $out = $cam_out;
    while (<$blast_file>) {
      if( /Sequence \: \"(.*)\"/) {
	$clone = $1; 
	if ($clone) {
	  if($clone2centre{$clone} eq "HX") {
	    $out = $cam_out;
	  }
	  elsif( $clone2centre{$clone} eq "RW") {
	    $out = $stl_out;
	  }
	  else {
	    $log->write_to("ERROR: NO CENTRE FOR $clone\n");
	    $out = *STDOUT;
	  }
	}
      }

      #test clone C30G7
      if ($test) {
	if ($clone =~ ("F58E6") ) {
	  if ($_ =~ /^Sequence : \"/) {
	    print $out "-D Homol_data $clone:wublastx_$type\n\n$_";
	  }
	  else {
	    print $out $_;
	  }
	}
      }
      else {
	if ($_ =~ /^Sequence : \"/) {
	  print $out "-D Homol_data $clone:wublastx_$type\n\n$_";
	}
	else {
	  print $out $_;
	}
      }
    }
    $wormbase->run_command("rm $temp_dir/$file\n\n", $log);
  }

  $log->write_to("\nSplit Output file can be found under $outdir/CAM_blastx.ace - STL_blastx.ace\n");
}

&load_data if ($update);

if ($update){$log->write_to("\n Update was selected so the CAM file has been loaded into the canonical camace\n");}

# Close all filehandles.
close $out;
close $cam_out;
close $stl_out;

$log->mail("pad\@sanger.ac.uk");

exit(0);

  sub load_data {
    print "\nWARNING: You are loading old data as you have chosen not to re-process the blastx data.\n\n" if ($noprocess);
    print "\nLoading data\n\n";
    my $command = "pparse $outdir/CAM_blastx.ace\nsave\nquit\n";
    print "This is the tace command being used:$tace ".$wormbase->database('camace') ." -tsuser merge_split_$next_build\n\n";
    open (DB, "| $tace ".$wormbase->database('camace') ." -tsuser merge_split_$next_build") || die "ERROR: Couldn't open database for $wormbase->camace\n";
    print DB $command;
    close DB;
  }
__END__

=pod

=item This will split the blastx data in to CAM and STL based on clone name.  Clone to centre assignments come from the COMMON_DATA file.

=cut

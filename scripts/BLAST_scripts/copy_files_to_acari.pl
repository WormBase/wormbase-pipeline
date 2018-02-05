#!/software/bin/perl
#
# copy_files_to_acari.pl
#
# dl
#
# copy the latest dna and wormpep files from the wormpub build structure to a machine available to the farm
# This machine used to be called acari - hence the script name :)
#

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;

my ($wormpep, $brigpep, $remapep, $brepep, $chroms, $test, $debug, $store, $ppapep, $jap, $brugpep, $ovolpep, $srapep,$tmupep);
GetOptions (
            'wormpep'   => \$wormpep,
            'brigpep'   => \$brigpep,
            'remapep'   => \$remapep,
            'ppapep'    => \$ppapep,
            'brepep'    => \$brepep,
            'chrom'     => \$chroms,
            "debug=s"   => \$debug,
            "test"      => \$test,
            "store:s"   => \$store,
            'jappep'    => \$jap,
            'brugpep'   => \$brugpep,
            'ovolpep'   => \$ovolpep,
            'srapep'    => \$srapep,
            'tmupep'    => \$tmupep,
	   );

my $wormbase;
if ( $store ) {
  $wormbase = Storable::retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $dir = $ENV{'PIPELINE'} . '/BlastDB';
my $log = Log_files->make_build_log($wormbase);

if ( $chroms ) {
    my @CHROMOSOME = ('I','II','III','IV','V','X');
    print "-------------------------------------------------------\nCopying new Chromosomal DNA\n\n";
    $log->write_to("Removing old DNA files . . . ");
    $wormbase->run_command( "rm -f $dir/*dna",$log) ;
    
    foreach my $chrom (@CHROMOSOME) {
      print "about to get Chromosome $chrom\n";
      $wormbase->run_command( "cp ".$wormbase->chromosomes."/CHROMOSOME_$chrom.dna $dir/" ,$log);
      $log->write_to("chromosome $chrom complete\n");
    }

    $log->write_to("Finished copying DNA to acari\n\n-------------------------------------------------------\n");
  }

&copy_worm_proteins($wormbase) if $wormpep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'briggsae')) if $brigpep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'remanei')) if $remapep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'pristionchus')) if $ppapep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'japonica')) if $jap ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'brenneri')) if $brepep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'brugia')) if $brugpep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'ovolvulus')) if $ovolpep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'sratti')) if $srapep ;
&copy_worm_proteins(Wormbase->new(-debug => $debug,-test => $test, -organism => 'tmuris')) if $tmupep ;

$log->mail;
exit(0);

#this is to copy a worm protein set that we update through curation eg elegans
sub copy_worm_proteins {
	my $species = shift;  # Species object
	my $WS_version = $species->get_wormbase_version;
  	my $wp_file = $species->wormpep .'/'. $species->pepdir_prefix .'pep'. $WS_version;
    
	my $new_file= $species->pepdir_prefix .'pep'. $WS_version .'.pep';
	my $old_file= $species->pepdir_prefix .'pep'. ($WS_version-1) .'.pep';

	$log->write_to("Copying new version of $new_file . . .\n");
        $wormbase->run_command("cp $wp_file $dir/$new_file",$log);
	chmod 0777, "$dir/$new_file";

	if (-e "$dir/$old_file"){
	        $log->write_to("Removing old version of $old_file . . .\n");
        	$wormbase->run_command("rm -f $dir/${old_file}*",$log) ;
	}

        $log->write_to("$new_file copied and ready for action\n\n-------------------------------------------------------\n") if -e "$dir/$new_file";
}


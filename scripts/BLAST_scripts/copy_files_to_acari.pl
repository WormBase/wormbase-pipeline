#!/software/bin/perl
#
# copy_files_to_acari.pl
#
# dl
#
# copy the latest dna and wormpep files from the wormpub build structure to a machine available to the farm
# This machine used to be called acari - hence the script name :)

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;

my ($wormpep, $brigpep, $chroms, $test, $debug, $store);
GetOptions (
	    'wormpep:i' => \$wormpep,
	    'brigpep'   => \$brigpep,
	    'chrom'     => \$chroms,
	    "debug=s"   => \$debug,
	    "test"      => \$test,
	    "store:s"   => \$store
	   );

my $wormbase;
if ( $store ) {
  $wormbase = Storable::retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( 'debug'   => $debug,
                             'test'    => $test,
			     );
}

my $dir = glob("~wormpipe/BlastDB");
my $log = Log_files->make_build_log($wormbase);

my @CHROMOSOME = ('I','II','III','IV','V','X');
if ( $chroms )
  {
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

&copy_worm_proteins('wormpep') if $wormpep ;
&copy_worm_proteins('brigpep') if $brigpep ;


$log->mail;
exit(0);

#this is to copy a worm protein set that we update through curation eg elegans
sub copy_worm_proteins {
	my $pep_set = shift;  # wormpep, brigpep etc
	my $WS_version = $wormbase->get_wormbase_version;
  	my $wp_file = "$pep_set" . $WS_version . ".pep";
   my $wp_old  = "$pep_set" . ($WS_version - 1) . ".pep";
    
   $wormbase->run_command("cp ".$wormbase->$pep_set."/$pep_set$WS_version $dir/$pep_set${WS_version}.pep",$log);

	if($pep_set eq "brigpep") {
		open (BRIG,"<$dir/$pep_set${WS_version}.pep") or $log->log_and_die("cant open $dir/$pep_set${WS_version}.pep");
		open (NEW,">$dir/$pep_set${WS_version}.pep_refmt") or $log->log_and_die("cant open $dir/$pep_set${WS_version}.pep_refmt :\n$!\n\n");
		while (<BRIG>) {
			if(/>CBG\d+\s+(CBP\d+)/) {
				print NEW ">$1\n";
			}
			else { print NEW;}
		}

		$wormbase->run_command("mv -f $dir/$pep_set${WS_version}.pep_refmt $dir/$pep_set${WS_version}.pep",$log)
	}
	
   $log->write_to("Removed old version of $pep_set . . ");
   $wormbase->run_command("rm -f $dir/${wp_old}*",$log) ;

   $log->write_to("$pep_set copied and ready for action\n\n-------------------------------------------------------\n");
   print "finished\n";
}

 ;

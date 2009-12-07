#!/usr/bin/perl -w
#===============================================================================
#         FILE:  build_webdb.pl
#
#        USAGE:  ./build_webdb.pl 
#
#  DESCRIPTION:  build Todd's GFF databases
#
#      $AUTHOR:$
#      COMPANY:  WormBase
#      CREATED:  11/27/09 10:10:08 GMT
#      CHANGED: $Date: 2009-12-07 10:29:27 $
#    $Revision: 1.1 $
#===============================================================================

# need to set the PERl5LIb  to pick up bioperl for the load_gff thing ... maybe have to hardcode it

use threads;
use threads::shared;

use Getopt::Long;
use lib $ENV{CVS_DIR};
use Wormbase;
use strict;

my @species;
my %tmpfiles;
my ($debug,$compress,$currentDB);

GetOptions(
      'species:s' => \@species,
      'debug:s'   => \$debug,
      'compress'  => \$compress,
      'currentDB' => \$currentDB, # use current_DB instead of the build databases
)||die(@!);

# postprocess the collected GFF files per species

@species = qw(elegans) if $currentDB;

## for each species
foreach my $worm(@species){
	my $wormbase = Wormbase->new(
		-species => $worm,
		-debug   => $debug,
	);

	if ($currentDB){
	   $wormbase = Wormbase->new(
		-species => $worm,
		-debug   => $debug,
		-autoace => glob('~wormpub/DATABASES/current_DB'),
	   );
	}

	my $log = Log_files->make_build_log($wormbase);
	$$wormbase{log}=$log; # crude, but works

	$log->write_to("... creating MySQL server\n");
	
	&tmpfile_hook("/tmp/${\$wormbase->species}.cnf");
        my ($mysql)=threads->create('start_mysql',$wormbase->species); 

        # munge the GFF collection
	$log->write_to("... creating munged GFFs\n");
	&process_worm($wormbase);

	# munge fasta files (remove the CHROMOSOME_)
	$log->write_to("... creating munged FASTA\n");
	&munge_fasta($wormbase);

        # create database for sucky gff loader
	$log->write_to("... creating the database\n");
	&create_database($wormbase);

        # load the files
	$log->write_to("... loading data into mysql\n");	
	&load_db($wormbase);	

        # take down the server
	$log->write_to("... stopping MySQL\n");
        &stop_mysql($wormbase);  

	# harvest the mysql thread
	$mysql->join();
	$mysql->detach();

        # compress the tables
	$log->write_to("... compressing tables\n");
	&compress_tables($wormbase) if $compress;
	
	# tar it up and copy it over
	$log->write_to("... tar and copy it\n");
	&tar_and_feather($wormbase);

	# deletes the tmp directory
	&cleanup_mysql($wormbase);
        &clean_tmpfiles;

        $log->mail();
}


print "Hasta Luego\n";

# tars up the database directory and copies it to the ftp server
sub tar_and_feather {
	my ($wb)=@_;
	my $sp = $wb->species;
	my $gspecies = $wb->full_name(-g_species => 1);
	my $ftpDir = $wb->ftp_site.'/web_data/';
	system("cd /tmp/${sp}_datadir/ && tar cvf $gspecies.tar") && die(@!);
	system("pbzip2 -9 /tmp/${sp}_data/$gspecies.tar") && die(@!);
	system("mv /tmp/${sp}_datadir/$gspecies.tar.bz2 $ftpDir") && die(@!);
	system("md5sum $ftpDir/$gspecies.tar.bz2 > /$gspecies.tar.md5") && die(@!);
}

# compresses the myisam tables
sub compress_tables {
	my ($wb)=@_;

	my $species=$wb->species;
        my $dbdir= "/tmp/${species}_datadir/$species/";
	
	foreach my $i (glob("$dbdir/*.MYI")){
	   $$wb{log}->write_to("...... compressing $i\n");
	   system("/software/worm/mysql/bin/myisampack $i"); # the small tables throw errors :-(
	   system("/software/worm/mysql/bin/myisamchk -rq --sort-index --analyze $i") && die(@!);
        }
}

# load gff files + dna into the gff database
sub load_db {
	my ($wb)=@_;
	my $species=$wb->species;
	system('/software/worm/perl_510/bin/perl '.
	       '/software/worm/ensembl/bioperl-live/scripts/Bio-DB-GFF/bulk_load_gff.PLS '.
	       "--create --user=root --fasta /tmp/$species.dna --local ".
	       "--database dbi:mysql:$species:localhost:mysql_socket=/tmp/$species.sock ".
	       "--defaults /tmp/$species.cnf") && die(@!);
}

# create mysql database directories
sub create_database {
	my ($wb) = @_;
	my $species=$wb->species;
	system("/software/worm/mysql/bin/mysql --defaults-file=/tmp/$species.cnf ".
	       "-uroot -e 'CREATE DATABASE $species'") && die(@!);
}

# concatenate the fasta files and remove the CHROMOSOME_
sub munge_fasta {
	my($wb)=@_;
	&tmpfile_hook("/tmp/${\$wb->species}.dna");
	system("cat ${\$wb->chromosomes}/*.dna|sed s/CHORMOSOME_//>/tmp/${\$wb->species}.dna");
}

# takes GFFs from $wormbase and creates a big one in /tmp/$species.gff
sub process_worm {
	my ($wb)=@_;

	my $buildDataDir = '/nfs/disk100/wormpub/BUILD_DATA/SUPPLEMENTARY_GFF';
	my @raw_gffs = glob($wb->chromosomes.'/*.gff');
	my @gz_gffs  = glob($wb->chromosomes.'/*.gff.gz');
        push @gz_gffs, glob("$buildDataDir/*.gff2.gz");
	push @raw_gffs,glob("$buildDataDir/*.gff");

	$$wb{log}->write_to("processing GFF files from ${\$wb->chromosomes} and $buildDataDir\n");
	$$wb{log}->write_to("... processing @gz_gffs\n");
	$$wb{log}->write_to("... processing @raw_gffs\n");

        my $gffile = concatenate_gff($wb,\@gz_gffs,\@raw_gffs);

	process_gff($wb->species,$gffile);
}

# concatenate the files
sub concatenate_gff {
	my ($wb,$gz,$gff)=@_;
        my $outfile = "/tmp/${\$wb->species}.gff_inf";

	tmpfile_hook($outfile);

	# gzipped ones
	$gz = join(' ',@$gz);
	system("zcat $gz>$outfile") && die(@!);
        

	# normal ones
	$gff = join(' ',@$gff);
	system("cat $gff>>$outfile") && die(@!);

	return $outfile;
}

# just the processing bit
sub process_gff {
	my($species,$file)=@_;

	&tmpfile_hook("/tmp/$species.gff");

        print STDERR "/software/bin/perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl ".
	       ">/tmp/$species.gff <$file\n";
	system("/software/bin/perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl ".
	       ">/tmp/$species.gff <$file") && die(@!);
}

# unlink the file and add it to the final cleanup
sub tmpfile_hook {
	my ($tmpfile)=@_;
	unlink $tmpfile if -e $tmpfile;
	$tmpfiles{$tmpfile}=1;
}

# that one should get threaded
sub start_mysql {
	my ($sp)=@_;
	# create config file
	# * file should be /tmp/$species.mysql.conf
	# * socket should become /tmp/$species.sock
	# * datadir should become /tmp/$species_datadir/
	# * needs only to listen to the socket and not TCP/IP
	# * needs only MyISAM table support
	
	# config file creation
	
	print STDERR "sed s/SPECIES/$sp/g $ENV{CVS_DIR}/web_data/maria_build.cnf > /tmp/$sp.cnf\n";
	system("sed s/SPECIES/$sp/g $ENV{CVS_DIR}/web_data/maria_build.cnf > /tmp/$sp.cnf") && die(@!);

        sleep(1);

	die("OMG the darn  /tmp/$sp.cnf does not exist!!!\n") unless -e "/tmp/$sp.cnf";

	# create database
	system("/software/worm/mysql/bin/mysql_install_db --defaults-file=/tmp/$sp.cnf") && die(@!);

	# start the thing
	system("/software/worm/mysql/bin/mysqld_safe --defaults-file=/tmp/$sp.cnf"); 
	# can't die that one, as the mysqladmin sends it some funky signal which it returns
}

sub stop_mysql {
        my ($wb)=@_;
	system("/software/worm/mysql/bin/mysqladmin --defaults-file=/tmp/${\$wb->species}.cnf -uroot shutdown");
	sleep(1);
}

sub cleanup_mysql {
	my ($wb)=@_;
	system("rm -rf /tmp/${\$wb->species}_datadir") && die(@!);
}

# pack it up for Todd

sub clean_tmpfiles {
	foreach my $key (keys %tmpfiles){
		unlink($key) if (-e $key);
	}
}

# just in case, you never know
END {
	&clean_tmpfiles();
}

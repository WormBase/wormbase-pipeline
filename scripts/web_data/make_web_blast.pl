#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  make_web_blast.pl
#
#        USAGE:  ./make_web_blast.pl 
#
#  DESCRIPTION:  create web blast databases for OICR
#
#      OPTIONS:  -species SPECIES
# REQUIREMENTS:  needs chromosome files
#         BUGS:  ---
#        NOTES:  ---
#       Author: $Author: mh6 $
#      COMPANY: WormBase
# LAST CHANGED: $Date: 2010-08-25 14:20:35 $
#     REVISION: $Revision: 1.3 $
#===============================================================================

# find the dna files
# blastify them
# stage them on ftp

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;


my ($species,$test,$debug,$store,$wormbase);

GetOptions(
     'species:s' => \$species,
     'test'      => \$test,
     'debug:s'   => \$debug,
     'store:s'   => \$store,
)||die(@!);

if ($store){
    $wormbase = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store");}
else{
 $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
    -organism=> $species
 );
}

$species||=$wormbase->species;

my $log = Log_files->make_build_log($wormbase);

$log->write_to("creating blast databases for $species\n");

# CHROMOSOME directory
my @infiles;
if (-e $wormbase->chromosomes . '/supercontigs.fa'){
	push @infiles,$wormbase->chromosomes . '/supercontigs.fa'
}
else {
	@infiles = glob($wormbase->chromosomes . '/*.dna') unless $infiles[0]
}

my $infile=join(' ',@infiles);

# that one is for Todd, as he does not like the CHROMOSOME_
if ($wormbase->species eq 'elegans'){
	system("cat $infile |sed s/CHROMOSOME_//>/tmp/elegans.dna") && die(@!);
	$infile = '/tmp/elegans.dna';
}

# to get a x_abc species name
my $outfile=$wormbase->ftp_site . "/web_data/blastdb/${\$wormbase->full_name(-g_species => 1)}";



# WORMPEP file
my $pepfile = $wormbase->wormpep.'/'.$wormbase->pepdir_prefix.'pep'.$wormbase->version;
system("cat $pepfile ".'| perl -pne "s/\t/ /g" > /tmp/'."$species.pep")&& die(@!);

$log->write_to("writing to $outfile\n");
system(qq(/software/worm/bin/ncbi_blast/formatdb -t "$species DNA" -p F -o T -i "$infile" -n $outfile -l /dev/null)) && die(@!);
system(qq(/software/worm/bin/ncbi_blast/formatdb -t "$species Proteins" -p T -o T -i /tmp/$species.pep -n $outfile -l /dev/null)) && die(@!);

$log->write_to("packing up $outfile files\n");

# that is to get the truncated file names for tar
my @blast_files = map {$_=~s/\/.*\///;$_} glob("$outfile.n* $outfile.p*");
my $blast_file  = join(' ',@blast_files);

system("tar -c -C ${\$wormbase->ftp_site}/web_data/blastdb/ -f $outfile.tar $blast_file")
  && die(@!);
system("rm $outfile.n* $outfile.p*")&& die(@!);

# cleanup
unlink('/tmp/elegans.dna') if (-e '/tmp/elegans.dna');
unlink("/tmp/$species.pep") if (-e "/tmp/$species.pep");
unlink("$outfile.tar.bz2") if (-e "$outfile.tar.bz2");
unlink("$outfile.tar.md5") if (-e "$outfile.tar.md5");

system("bzip2 -9 $outfile.tar")&& die(@!);
system("md5sum $outfile.tar.bz2 > $outfile.tar.md5")&& die(@!);

$log->mail();

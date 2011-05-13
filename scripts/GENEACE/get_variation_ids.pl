#!/software/bin/perl -w
# a script to batch request variation ids based on lists of public_names
# Last change by $Author: mh6 $ on $Date: 2011-05-13 10:38:15 $
# usage: perl get_variation_ids.pl -species elegans -user me -pass me -in one_public_id_per_line -out varId_pubId_per_line


use lib $ENV{'CVS_DIR'};
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib '/software/worm/lib/perl';
use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;
use IO::File;
use strict;


my ($PASS,$USER); # mysql ones
my $DOMAIN  = 'Variation'; # hardcoded to variation
my ($debug, $test, $store, $species, $maria,$infile,$outfile,$weak,$nocheck);

GetOptions (
	    "debug=s"    => \$debug,   # send log emails only to one user
	    "test"       => \$test,    # run against the test database on mcs4a
	    "store:s"    => \$store,   # if you want to pass a Storable instead of recreating it
	    "species:s"  => \$species, # elegans/briggsae/whateva .. needed for logging
	    "user:s"	 => \$USER,    # mysql username
	    "pass:s"	 => \$PASS,    # mysql password
	    'maria'      => \$maria,   # try the more experimantal MariaDB for testing
	    'infile:s'   => \$infile,  # input file name
	    'outfile:s'  => \$outfile, # output file name
	    'weak'       => \$weak,    # weak check public_names
	    'nocheck'    => \$nocheck, # don'tcheck public_names
	   )||die(@!);


my $wormbase;

die "Species option is mandatory\n" unless $species; # need that for NameDB to fetch the correct regexps

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $DB = $wormbase->test ? 'test_wbgene_id;mcs4a:3307' : 'wbgene_id;shap:3303';
$DB='namedb;wormbase2:3307' if $maria; # there exist only one namedb database on the mirror

my $db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/data');
$db->setDomain($DOMAIN);

# there should be *only* one public_name per line in the input_file
# like:
# abc123
# abc234
my $inHandle = new IO::File $infile,'r' ||die(@!);
my $outHandle = new IO::File $outfile,'w' ||die(@!);
while (<$inHandle>) {
	chomp;
	
        if(&_check_name($_)||$nocheck) {
           my $id = $db->idCreate;
           $db->addName($id,'Public_name'=>$_);
           print $outHandle "$id\t$_\n";
	}else{print STDERR "ERROR: cannot create $_\n"}
}

$inHandle->close;
$outHandle->close;

$log->mail;

# small function to check the variation public_name for sanity
sub _check_name {
    my $name = shift;
    my $var = $db->idGetByTypedName('Public_name'=>$name)->[0];
    if($var) {  
        print STDERR "ERROR: $name already exists as $var";
        return undef;
    }
    elsif ($name =~ /(\w+)[a-z]+$/) {
	unless ($weak){
         my $short_var = $db->idGetByTypedName('Public_name'=>$1)->[0];
         if($short_var) {        
            print STDERR "ERROR: $name looks like a poorly name version of $short_var";
            return undef;
         }
        }
    }
    return 1;
}


#!/usr/bin/env perl
#
# a script to batch request variation ids based on lists of public_names
# Last change by $Author: mh6 $ on $Date: 2015-03-31 10:16:20 $
# usage: perl get_variation_ids.pl -species elegans -user me -pass me < file containing varId_pubId_per_line


use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";

use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;
use IO::File;
use strict;


my ($PASS,$USER, $DB); # mysql ones
my $DOMAIN  = 'Variation'; # hardcoded to variation

my ($wormbase, $debug, $test, $store, $species,$infile,$outfile,$nocheck, 
    $ace_file_template, $input_is_public_names, $input_is_other_names);

GetOptions (
  "debug=s"       => \$debug,   # send log emails only to one user
  "test"          => \$test,    # run against the test database on utlt-db
  "store:s"       => \$store,   # if you want to pass a Storable instead of recreating it
  "species:s"     => \$species, # elegans/briggsae/whateva .. needed for logging
  "user:s"	  => \$USER,    # mysql username
  "pass:s"	  => \$PASS,    # mysql password
  'nocheck'       => \$nocheck, # don't check public_names
  'acetemplate=s' => \$ace_file_template,
    )||die(@!);
    
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

$DB = $wormbase->test ? 'test_wbgene_id;utlt-db:3307' : 'wbgene_id;shap:3303';

my $db = NameDB_handler->new($DB,$USER,$PASS,$wormbase->wormpub . '/DATABASES/NameDB');
$db->setDomain($DOMAIN);

# there should be *only* one public_name per line in the input_file
# like:
# abc123
# abc234

if (defined $ace_file_template and -e $ace_file_template) {
  my $tmp_str = "";

  open(my $atfh, $ace_file_template) 
      or $log->log_and_die("Could not open $ace_file_template for reading\n");
  while(<$atfh>) {
    /\S/ and $tmp_str .= $_;
  }
  
  $ace_file_template = $tmp_str;
}

my @entries;
while (<>) {
  chomp;

  my @names = split(/\s+/, $_);
  if (scalar(@names) == 2 or scalar(@names) == 1) {
    push @entries, \@names,
  } else {
    $log->log_and_die("Bad line:$_\n");
  }
}

foreach my $entry (@entries) {
  my ($other_name, $public_name) = @$entry;

  if (not $nocheck and defined $public_name and not defined &_check_name($public_name)) {
    next;
  }

  my $id = $db->idCreate;
  
  $public_name = $id if not defined $public_name;

  $db->addName($id,'Public_name'=>$public_name);

  if (defined $ace_file_template) {
    print "Variation : \"$id\"\n";
    print "Public_name $public_name\n";
    print "Other_name $other_name\n";
    print "$ace_file_template\n";
  } else {
    print join("\t", $id, $public_name, $other_name), "\n"; 
  } 
}

$log->mail;

# small function to check the variation public_name for sanity
sub _check_name {
  my $name = shift;
  my $var = $db->idGetByTypedName('Public_name'=>$name)->[0];
  if($var) {  
    print STDERR "ERROR: $name already exists as $var\n";
    return undef;
  }
  return 1;
}

__END__

=pod

=head1 NAME - get_variation_ids.pl

=head1 DECRIPTION

This script takes a list of submitter-supplied variation names  and generates new WBVar identifiers in the name server, 
writing the results to the given output file. 

The input file is expected to contain a list of submitter-supplied variation names, one per line. Optionally, each 
line may also contain an additional name, which is added to the NameDB as the public name. If no second name is 
defined, the generated WBVar identifier will be used as the public name

By default, the script will refuse to create new WBVar identifiers for public names that already exist in the
database; this behaviour can be switched off with the -nocheck option. 

-head1 OUTPUT

The default output is a 3-column table of 

WBVarid<tab>Public name<tab>Submitter-supplied name

However, if the -acetemplate  option is used, Ace format can be generated (see -acetemplate below). With Ace format, 
the submitter-supplied name is used to populate the Other_name tag. 


=head1 USAGE EXAMPLES

get_variation_ids.pl -species briggsae input_id_-user me -pass me file.txt > output_id_file.txt

get_variation_ids.pl -species elegans -user me -pass me -acetemplate wild_template.ace input_id_file.txt > output_file.ace

=head1 OPTIONS

-species      WormBase species name, required for logging (defaults to elegans)

-user         MySQL NameDB account user name

-pass         MySQL NameDB account password

-acetemplate  Write Ace format, using the given file as a template; the contents will be tagged to the end of each entry
 
-nocheck      Do not check for public name clashes

=cut


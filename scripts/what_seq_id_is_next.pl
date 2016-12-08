#!/usr/bin/env perl
#
# a script to return the next sequence_id available in the nameserver for a given clone or list of clones.
# Last change by $Author: pad $ on $Date: 2015-07-02 13:48:58 $
# usage: perl what_seq_id_is_next -species elegans -user me -pass me -in one_clone_per_line -out seqId_dash_clone_per_line


use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";

use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;
use IO::File;
use strict;
# PERL MODULES WE WILL BE USING
use DBI;
use DBD::mysql;

my ($wormbase, $debug, $test, $store, $species, $in,$out,$pw,$user,$seq);

GetOptions (
	    'debug=s'     => \$debug,   # send log emails only to one user
	    'test'        => \$test,    # run against the test database on mcs4a
	    'store:s'     => \$store,   # if you want to pass a Storable instead of recreating it
	    'user:s'	  => \$user,    # mysql username
	    'pass:s'	  => \$pw,    # mysql password
	    'in:s'        => \$in,
	    'out:s'       => \$out,
	    'clone:s'     => \$seq,
	   )||die(@!);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
      );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $version = $wormbase->version;

# there should be *only* one clone per line + a gene value in the input_file
# like:
# AH6 unc-20
# Y52D7 lev-10

unless (defined$out){
  if (defined $in){
    $out = "${in}_output";
  }
  elsif(defined $seq){
    $out = "${seq}_output";
  }
}

open (OUT, "> $out") or $log->log_and_die("Could not open $out output file\n");


my (@clones, @genes);
if (defined$in){

  open(IN, $in) or $log->log_and_die("Could not open $in file for reading\n");
  $log->write_to ("Processing $in file\nResults will be stored in $out\n\n");

  while (<IN>) {
    chomp;
    my @names = split(/\s+/, $_);
    if (scalar(@names) == 2) {
      push @clones, $names[0];
      push @genes, $names[1];
    } else {
      $log->log_and_die("Bad line:$_\n");
    }
  }
}
elsif (defined$seq){
  push @clones, $seq;
}

else {
  $log->log_and_die("No input file specified or clone given on command line, please use -in <File> or -clone BLAH\n");
}

# MYSQL CONFIG VARIABLES
my ($port,$host,$database,);
if ($test) {
 $host = 'utlt-db';
 $database = 'test_wbgene_id';
 $port = 3307;
}else {
 $host = 'web-wwwdb-core-02';
 $database = 'nameserver_live';
 $port = 3449;
}


# DATA SOURCE NAME
my $dsn = "dbi:mysql:$database:$host:$port" or die "Unable to connect: $DBI::errstr\n";

# PERL DBI CONNECT
my $connect = DBI->connect($dsn, $user, $pw);


my (%storedIDs,$myquery, $execute,);
my $counter = "-1";
foreach my $clone (@clones) {
  $counter ++;
  $log->write_to ("$clone $genes[$counter] - ID assigned = ");
  print OUT '-R Transcript ',$genes[$counter]," ";
  my ($query_handle,$name);
  $myquery = "select distinct(object_name) from secondary_identifier where object_name LIKE \"$clone.%\" order by object_id";
  $query_handle = $connect->prepare($myquery);
  $query_handle->execute(); 
  $query_handle->bind_columns(\$name);
  my $ID;
  my $root;
  my $test = "0";
  while($query_handle->fetch()) {
    if ($name =~ /(\S+)\.(\d+)/) {
      $root = $1;
      if ($2 > $test) {
	$test = $2
      }
    }
  }
  if (!defined$root) {
    print "Warning: This appears to be the 1st gene model on the clone\nAssigning ${clone}.1\n\n";
    $root = ${clone};
    $test = "0";
  }
  $ID = "$root\."."$test";
  my ($digit,$new_id,);
  if ($ID =~ /(\S+\.)(\d+)/) {
    my $digit = $2;
    $digit ++;
    $new_id = "${1}$digit";
    $log->write_to ("1st ID assigned = $new_id\n");
    # check to see if the new name exists in the name server.
    my ($check_test,$capture);
    my $mycheck = "select distinct(object_name) from secondary_identifier where object_name = \'$new_id\'";
    $check_test = $connect->prepare($mycheck);
    $check_test->execute();
    my (@row,$capture);
    while (@row = $check_test->fetchrow_array()) {
      $capture = $row[0];
      $log->write_to ("ERROR: $capture already exists in the database\n");
    }

    #check to see if it has been assigned in this round already
    my $itterator = "0";
  CHECK_ID: while ($storedIDs{$new_id}){
	$log->write_to ("$new_id - Already assigned this ID itteration ${itterator}.....\n");
	$digit ++;
	$new_id = "${root}.$digit";
	$log->write_to ("Next ${itterator} ID assigned $new_id\n");
	$itterator ++;
	last CHECK_ID if (!$storedIDs{$new_id});
	$storedIDs{$new_id} = 1;
      }
    $log->write_to("$new_id passes internal use test\n");

    # check to see if the ID has ever been used for a CDS
    my $new      = `fgrep -c $root"\.$digit" /nfs/wormpub/BUILD/WORMPEP/wormpep_${version}./wormpep.history*`;
    if ($new > 0){
      $log->write_to("ROUND1 $new_id fails wormpep test trying again.......\n");
      $storedIDs{$new_id} = 1;
      $digit ++;
      $new_id = "${root}.$digit";
      my $new2      = `fgrep -c $root"\.$digit" /nfs/wormpub/BUILD/WORMPEP/wormpep_current/wormpep.history*`;
      if ($new2 > 0){
	$log->write_to("ROUND2 $new_id fails wormpep test trying again.......\n");
	$storedIDs{$new_id} = 1;
	$digit ++;
	$new_id = "${root}.$digit";
	my $new3      = `fgrep -c $root"\.$digit" /nfs/wormpub/BUILD/WORMPEP/wormpep_current/wormpep.history*`;
	if ($new3 > 0){	
	  $log->write_to ("ROUND3 $new_id fails wormpep tests <Manually assign>\n");
	}
	else {
	  $log->write_to ("$new_id //assigned in Round 3\n");
	  $storedIDs{$new_id} = 1;
	}
      }
      else {
	$log->write_to ("$new_id //assigned in Round 2\n");
	$storedIDs{$new_id} = 1;
      }
    }
    else {
      $log->write_to ("$new_id //assigned\n");
      $storedIDs{$new_id} = 1;
    }
    print "Transcript $new_id\n";
    print "CDS $new_id\n";
    print "Gene_name $new_id\n";
    print "Pseudogene $new_id\n";
    print "Sequence_name $new_id\n";
  }
  print OUT "$new_id\n";
}
$counter ++;
$log->write_to ("Processed $counter clones from $in\n\n");
$log->mail;


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

get_variation_ids.pl -species briggsae input_id_file.txt > output_id_file.txt

get_variation_ids.pl -species elegans -acetemplate wild_template.ace input_id_file.txt > output_file.ace

=head1 OPTIONS

-species      WormBase species name, required for logging (defaults to elegans)

-user         MySQL NameDB account user name

-pass         MySQL NameDB account password

-acetemplate  Write Ace format, using the given file as a template; the contents will be tagged to the end of each clone
 
-nocheck      Do not check for public name clashes

=cut


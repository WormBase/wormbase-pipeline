#!/usr/bin/env perl
#
# a script to batch request variation ids based on lists of public_names
# Last change by $Author: mh6 $ on $Date: 2015-03-31 10:16:20 $
# usage: perl get_NS_ids.pl -species elegans -user me -pass me < file containing varId_pubId_per_line


use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";

use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;
use IO::File;
use strict;


my ($PASS,$USER, $DB); # mysql ones
my $DOMAIN;

my ($wormbase, $debug, $test, $store, $species,$infile,$outfile,$nocheck, 
    $ace_file_template, $input_is_public_names, $input_is_other_names, $input, $output,$public,$varname);

GetOptions (
    "debug=s"   => \$debug,   # send log emails only to one user
    "test"      => \$test,    # run against the test database on utlt-db
    "store:s"   => \$store,   # if you want to pass a Storable instead of recreating it
    "species:s" => \$species, # elegans/briggsae/whateva .. needed for logging
    'nocheck'   => \$nocheck, # don't check public_names
    'input:s'   => \$input,   # File containing new public_names for new Variations
    "output:s"  => \$output,  # File to capture the new WBVar and name associations
    'acetemplate=s' => \$ace_file_template, #unknown template file that is also printed along with the IDs retrieved by the script
    'varname'   => \$varname,
    'public'    => \$public,   # The file only contains public_names so these should be used (This should be used for all types unless Variation with other names).
    'domain:s'  => \$DOMAIN, # specify Gene, Variation, Feature or Strain to retrieve this type of ID.
    )||die(@!);
    
die "Species option is mandatory\n" unless $species; # need that for NameDB to fetch the correct regexps
die "No domain specified which is mandatory\n" unless $DOMAIN; # need that for NameDB to assign correct ID type
#die "Must specify -public if domain isn't Variation\n" unless (($public) && ($DOMAIN ne "Variation"));
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

my $db = NameDB_handler->new($wormbase);

# there should be *only* one public_name per line in the input_file
# like:
# abc123
# abc234

if (defined $ace_file_template) {
  my $tmp_str = "";

  open(my $atfh, $ace_file_template) 
      or $log->log_and_die("Could not open $ace_file_template for reading\n");
  while(<$atfh>) {
    /\S/ and $tmp_str .= $_;
  }
  
  $ace_file_template = $tmp_str;
}

my @entries;
open (OUT, ">$output") or $log->log_and_die("Could not open $output for reading\n");
open (IN, "<$input") or $log->log_and_die("Could not open $input for reading\n");
while (<IN>) {
    chomp;
    s/\"//g;
    my @names;
    @names = $_;

    my $id;
    foreach my $entry (@names) {
        if ($DOMAIN =~ /Strain/){
            my $strain_ref = $db->info_strain($entry);
            unless (defined $strain_ref) {
                print OUT "\/\/Cannot find $entry\n";
                next;
            }
            if (defined $strain_ref) {
                my $WBstrain_id = $strain_ref->{'id'};
                print OUT "-R $DOMAIN $entry $WBstrain_id\n";
            }
        }
    
        
        if ($DOMAIN =~ /Variation/){
            my $var_ref = $db->info_variation($entry);
            my $WBvar_id;
            unless (defined $var_ref) {
                print OUT "\/\/Cannot find $entry\n";
                next;
            }
            if (defined $var_ref) {
                if ($entry =~ /WBVar/) {
                    $WBvar_id = $var_ref->{'name'};
                }
                else {
                    $WBvar_id = $var_ref->{'id'};
                }
                print OUT "-R $DOMAIN $entry $WBvar_id\n";
            }
        }
        
        
        if ($DOMAIN =~ /Gene/){
            my $gene_ref = $db->find_genes($entry);
            unless (defined $gene_ref->[0]{'id'}) {
                print OUT "\/\/Cannot find $entry\n";
                next;
            }
            if (defined $gene_ref) {
                my $WBgene_id = $gene_ref->[0]{'id'};
                print OUT "-R $DOMAIN $entry $WBgene_id\n";
            }
        }
    }
}

$log->mail;
# small function to check the strain public_name for sanity



__END__

=pod

=head1 NAME - show_NS_ids.pl

=head1 DECRIPTION

This script takes a list of submitter-supplied <single_domian> public_names 

-head1 OUTPUT

The default output is a 3-column table of 

WBVarid<tab>Public name<tab>Submitter-supplied name

However, if the -acetemplate  option is used, Ace format can be generated (see -acetemplate below). With Ace format, 
the submitter-supplied name is used to populate the Other_name tag. 


=head1 USAGE EXAMPLES

get_variation_ids.pl -species briggsae input_id_-user me -pass me input file.txt -output output_id_file.txt

get_variation_ids.pl -species elegans -user me -pass me -acetemplate wild_template.ace -input input_id_file.txt -output output_file.ace

=head1 OPTIONS

-species      WormBase species name, required for logging (defaults to elegans)

-user         MySQL NameDB account user name

-pass         MySQL NameDB account password

-acetemplate  Takes the file name of a file containing Ace template format which will be printed to the end of the results to create a fuller entry.
 
-nocheck      Do not check for public name clashes

-input        Takes the file name that you have generated as input

-output       Takes the name of a file that the results from the nameserver are written to. 

=cut


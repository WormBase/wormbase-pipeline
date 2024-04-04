#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;


use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

use lib "$ENV{CVS_DIR}/ONTOLOGY";
use GAF;

my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $include_iea, $ws_version, $outfh);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "electronic"  => \$include_iea,
  "wsversion=s" => \$ws_version,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = &get_GAF_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

if ($outfile) {
  open($outfh, ">$outfile") or die("cannot open $outfile : $!\n");  
} else {
  $outfh = \*STDOUT;
}


my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my ( $count, $it);

&print_DAF_header($outfh);
&print_DAF_column_headers($outfh);

$it = $db->fetch_many(-query=>'find Gene Disease_info');

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
  next unless $obj->Species->name eq $full_name;

  my $g = $obj->name;

  if ($include_iea) {
    foreach my $doterm ($obj->Potential_model) {
      
      my $meth = $doterm->right->right;
      if ($meth->name eq 'Inferred_automatically') {
        my $text = $meth->right->name;
        my ($with_from_list) = $text =~ /\((\S+)\)/;
        
        my @ens = map { "ENSEMBL:$_" } grep { $_ =~ /ENSG\d+/ } split(/,/, $with_from_list);
        my @omim = grep { $_ =~ /O?MIM:/ } split(/,/, $with_from_list);
        
        my $obj =  {
          object_type => "gene",
          object_id => $g,
          inferred_ga => "WB:$g",
          object_symbol =>  $obj->Public_name->name,
          do_id => $doterm,
          reference => "PMID:19029536",  # this is the reference for Ensembl Compara
          evidence => "IEA", 
          with => join("|", @ens,@omim),
          date => "",
        };
        
        &print_DAF_line($outfh, $obj);
      }
    }
  }

  foreach my $doterm ($obj->Experimental_model) {
    my (@papers, $evi_date);
    foreach my $evi ($doterm->right->col) {
      if ($evi->name eq 'Paper_evidence') {
        foreach my $paper ($evi->col) {
          my $pmid;
          foreach my $db ($paper->Database) {
            if ($db->name eq 'MEDLINE') {
              $pmid = $db->right->right->name;
              last;
            }
          }
          push @papers, $pmid if $pmid;
        }
      } elsif ($evi->name eq 'Date_last_updated') {
        $evi_date = $evi->right;
        $evi_date->date_style('ace');
        my ($y, $m, $d) = split(/\-/, $evi_date);
        $evi_date = sprintf("%4d%02d%02d", $y, $m, $d);
      }
    }
    foreach my $paper (@papers) {
      my $obj = {
        object_type => "gene",
        object_id => $g,
        inferred_ga => "WB:$g",
        object_symbol =>  $obj->Public_name,
        do_id => $doterm,
        evidence => "IMP", 
        reference => "PMID:$paper", 
        assoc_type => 'causes_condition',
        sex => "hermaphrodite",
        date => $evi_date,
      };
      &print_DAF_line($outfh, $obj);  
    }
  }
}

$db->close;

exit(0);

	
###########################################################3
sub print_DAF_header {
  my ($fh) = @_;

  print $fh "!daf-version 0.1\n";
  print $fh "!Date: $date\n";
  print $fh "!Project_name: WormBase (WB) Version $ws_version\n";
  print $fh "!URL: http://www.wormbase.org/\n";
  print $fh "!Contact Email: wormbase-help\@wormbase.org\n";
  print $fh "!Funding: NHGRI at US NIH, grant number U41 HG002223\n";

}


###########################################################3
sub print_DAF_column_headers {
  my ($fh) = @_;

  my @col_headers  = (
    'Taxon',
    'DB Object Type',
    'DB',
    'DB Object ID',
    'DB Object Symbol',
    'Inferred gene association',
    'Gene Product Form ID',
    'Experimental conditions',
    'Association type',
    'Qualifier',
    'DO ID',
    'With',
    'Modifier - assocation type',
    'Modifier - Qualifier',
    'Modifier - genetic',
    'Modifier - experimental conditions',
    'Evidence Code',
    'genetic sex',
    'DB:Reference',
    'Date',
    'Assigned By');

  print $fh join("\t", @col_headers), "\n";
}


###########################################################3
sub print_DAF_line {
  my ($fh, $obj) = @_;

  printf($fh "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         "taxon:$taxid",
         $obj->{object_type}, 
         "WB", 
         $obj->{object_id},
         $obj->{object_symbol},
         (exists $obj->{inferred_ga}) ? $obj->{inferred_ga} : "",
         (exists $obj->{product_id}) ? $obj->{product_id} : "",
         (exists $obj->{exp_cond}) ? $obj->{exp_cond} : "",
         $obj->{assoc_type},
         (exists $obj->{qualifier}) ? $obj->{qualifier} : "",
         $obj->{do_id},
         (exists $obj->{with}) ? $obj->{with} : "",
         (exists $obj->{mod_assoc_type}) ? $obj->{mod_assoc_type} : "",
         (exists $obj->{mod_qualifier}) ? $obj->{mod_qualifier} : "",
         (exists $obj->{mod_genetic}) ? $obj->{mod_genetic} : "",
         (exists $obj->{mod_exp_cond}) ? $obj->{mod_exp_cond} : "",
         $obj->{evidence},
         $obj->{sex},
         $obj->{reference},
         ($obj->{date}) ? $obj->{date} : $date,
         "WB");

}

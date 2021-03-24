#!/usr/bin/env perl

use strict;
use Ace;
use Getopt::Long;
use Storable;
          
use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

use lib "$ENV{CVS_DIR}/ONTOLOGY";
use GAF;

my ($help, $debug, $test, $store, $wormbase,$tace);
my ($output, $acedbpath, $rnai, $variation, $skiplist);

my %opts=();
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "rnai"   	 => \$rnai,
	    "variation"  => \$variation,
	    "output:s"   => \$output,
	    "skiplist:s" => \$skiplist,
	    )||die(@!);

if ($help) {
    print "usage: $0 [options] -output output -database database\n";
    print "       -help            help - print this message\n";
    print "       -output <output>   output file\n";
    print "       -database <database> path to database\n";
    print "       -rnai            generate RNAi2GO mapping; default no\n";
    print "       -variation       generate variation2GO mapping; default no\n";
    print "       -skiplist <papers_to_skip> list of papers to skip when generating RNAi2GO mappings\n";
    exit;
}

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my %aspect=(molecular_function=>'F',
	    biological_process=>'P',
	    cellular_component=>'C'
	    );
	    
$acedbpath = $wormbase->autoace if not defined $acedbpath;
$tace      = $wormbase->tace;

my $date = &get_GAF_date();
my $full_name = $wormbase->full_name;
my $taxid = $wormbase->ncbi_tax_id;
my $gene_info = &get_gene_info( $acedbpath, $wormbase, $full_name );


my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);

my (%paper_fields, %paper_accno, %papers_to_skip, %go_term_types, %phen2go);

my @aql_results = $db->aql('select a, a->Database[2], a->Database[3]  from a in class paper');
foreach (@aql_results) {
  $paper_fields{$_->[0]} = $_->[1]; # database field e.g. 'PMID'
  $paper_accno{$_->[0]}  = $_->[2];  # database accession number e.g. '12393910'
}

@aql_results = $db->aql('select a, a->Type[1] from a in class GO_term');
foreach my $res (@aql_results) {
  $go_term_types{$res->[0]} = $res->[1];
}

if ($skiplist) {
  open (my $skip_fh, "<$skiplist") or $log->log_and_die("cannot open $skiplist : $!\n");
  while (<$skip_fh>) {
    chomp;
    s/\s+//g;
    next unless /WBPaper/;
    $papers_to_skip{$_}=1;
  }
  $log->write_to( scalar(keys %papers_to_skip) . " papers will be skipped\n");
}

$output = $wormbase->ontology."/phenotype2go.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open(my $out, ">$output") or $log->log_and_die("cannot open $output : $!\n");

&print_wormbase_GAF_header($out, $wormbase->get_wormbase_version_name);

if ($rnai or $variation) {
  my @aql_results = $db->aql("select p, p->GO_term from p in class Phenotype where exists p->GO_term");
  foreach(@aql_results) { 
    if (defined $_->[1]) { 
      $phen2go{$_->[0]}{$_->[1]}=1;
    }
  }
  $log->write_to( scalar(keys %phen2go) . " phenotypes read\n");
}

if ($rnai) {
  
  my (@phen_array);
  
  my $it=$db->fetch_many(-query=>'find phenotype go_term; follow rnai');
    
  while (my $obj=$it->next) {
    next unless $obj->isObject();

    my $ref = $obj->Reference||'';
    next if $papers_to_skip{$ref};
    
    my @genes=();
    foreach my $g ($obj->Gene) {
      if ($g->right(2) eq 'RNAi_primary') {
        push @genes, $g;
      }
    }

    @genes = grep { exists $gene_info->{$_} and $gene_info->{$_}->{status} ne 'Dead' } @genes;

    foreach my $gene (@genes) {
      my $pub_name = $gene_info->{$gene}->{public_name};
      my $seq_name = $gene_info->{$gene}->{sequence_name};
      
      foreach my $phen ($obj->Phenotype) {
        next if not $phen2go{$phen};
        
        my $ref_field = '';
        if ($paper_fields{$ref}) {
          $ref_field="WB_REF:$ref|$paper_fields{$ref}:$paper_accno{$ref}";
        }
        elsif ($ref) {
          $ref_field="WB_REF:$ref";
        } else {
          # if the RNAi object does not have a Reference, give a warning and skip this RNAi because the $ref_field is a mandatory term
          $log->write_to("No reference for $obj so skipping\n");
          next;
        }
        my $with="WB:$obj|$phen";
        
        foreach my $term (keys %{$phen2go{$phen}}) {
          my $go_type = $go_term_types{$term};
          my $a = $aspect{lc $go_type};
          
          &print_wormbase_GAF_line($out,
                                   $gene,
                                   $pub_name,
                                   "",
                                   $term,
                                   $ref_field,
                                   "IEA",
                                   $with,
                                   $a,
                                   $seq_name,
                                   $taxid,
                                   $date);
          
        }
      }
    }
  }
}


if ($variation) {
  
  my $it=$db->fetch_many(-query=>'find phenotype go_term; follow variation');

  while (my $obj=$it->next) {
    next unless $obj->isObject();
  
    my @genes = $obj->Gene;
    @genes = grep { exists $gene_info->{$_} and $gene_info->{$_}->{status} ne 'Dead' } @genes;

    next if not @genes;

    my %phen_hash=();
    foreach my $phen ($obj->Phenotype) {
      next if not $phen2go{$phen};

      if ($phen->at("Variation_effect")) {
        my $gain_of_function = 0;
        foreach my $effect ($phen->at("Variation_effect")) {
          if ($effect =~ /gain_of_function/i) {
            $gain_of_function = 1;
          }
        }
        if ($gain_of_function) {
          $log->write_to("Skipping $obj<->$phen because annotated as Gain-of-function\n");
          next;
        }
      }

      $phen_hash{$phen}{count}++;
      if ($phen->at("Paper_evidence")) {
        $phen_hash{$phen}{Paper_evidence} = $phen->at("Paper_evidence[1]");
      }
      if ($phen->at("Person_evidence")) {
        $phen_hash{$phen}{Person_evidence} = $phen->at("Person_evidence[1]")=~/(WBPerson\d+)/ ? $1 : '';
      }
      if ($phen->at("Curator_confirmed")) {
        $phen_hash{$phen}{Curator_confirmed} = $phen->at("Curator_confirmed[1]")=~/(WBPerson\d+)/ ? $1 : '';
      }
    }
  
    foreach my $gene (@genes) {
      my $pub_name = $gene_info->{$gene}->{public_name};
      my $seq_name = $gene_info->{$gene}->{sequence_name};

      foreach my $phen (sort {$a cmp $b} keys %phen_hash) {

        my $ref_field='';
        if ($phen_hash{$phen}{Paper_evidence}) {
          if ($paper_fields{$phen_hash{$phen}{Paper_evidence}}) {
            $ref_field = join("|", 
                              "WB_REF:$phen_hash{$phen}{Paper_evidence}", 
                              "$paper_fields{$phen_hash{$phen}{Paper_evidence}}:$paper_accno{$phen_hash{$phen}{Paper_evidence}}"); 
          } else {
            $ref_field = "WB_REF:$phen_hash{$phen}{Paper_evidence}";
          }
        } elsif ($phen_hash{$phen}{Person_evidence}) {
          $ref_field = "WB:$phen_hash{$phen}{Person_evidence}";
        } elsif ($phen_hash{$phen}{Curator_confirmed}) {
          $ref_field = "WB:$phen_hash{$phen}{Curator_confirmed}";
        }
        my $with = "WB:$obj|$phen";

        foreach my $term (keys %{$phen2go{$phen}}) {
          my $go_type = $go_term_types{$term};
          my $a = $aspect{lc $go_type};
          
          &print_wormbase_GAF_line($out,
                                   $gene,
                                   $pub_name,
                                   "",
                                   $term,
                                   $ref_field,
                                   "IEA",
                                   $with,
                                   $a,
                                   $seq_name,
                                   $taxid,
                                   $date);
          
        }
      }
    }
  }
}

close($out);
$db->close;
$log->mail;

exit(0);

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
my ($output, $acedbpath, $skiplist);

my %opts=();
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "output:s"   => \$output,
	    "skiplist:s" => \$skiplist,
	    )||die(@!);

if ($help) {
    print "usage: $0 [options] -output output -database database\n";
    print "       -help            help - print this message\n";
    print "       -output <output>   output file\n";
    print "       -database <database> path to database\n";
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
my $gene_info = &get_gene_info( $acedbpath, $wormbase );

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);
$db->date_style('ace');

my (%paper_fields, %paper_accno, %papers_to_skip, %go_term_types);

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

$output = $wormbase->ontology."/gene_association.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open(my $out, ">$output") or $log->log_and_die("cannot open $output : $!\n");

&print_wormbase_GAF_header($out);

my $it=$db->fetch_many(-query=>'FIND Gene GO_term AND NOT Dead');

while (my $obj = $it->next) {
  next unless $obj->isObject();
  
  my $wbgene = $obj->name;
  next unless exists $gene_info->{$wbgene};
  
  my $public_name = $gene_info->{$wbgene}->{public_name};
  my $seq_name    = $gene_info->{$wbgene}->{sequence_name};
  my $taxid       = $gene_info->{$wbgene}->{taxid};
  
  foreach my $term ($obj->GO_term){
    my $go_type = $go_term_types{$term};
    my $evi_code = $term->right;
    my $ec = $evi_code->name;

    my ($ref,$with) = ("", "");
    my %attr;
    
    foreach my $subtree ($evi_code->col()) {
      my @l = $subtree->col();
      if ($subtree->name ne 'Inferred_automatically') {
        $attr{$subtree->name} = $l[0]->name;
      } else {
        $attr{$subtree->name} = [map { $_->name } @l];
      }
    }
    
    if (exists $attr{Paper_evidence}) {
      my $pap =  $attr{Paper_evidence};
      $ref = sprintf("WB_REF:%s", $pap);
      if ($paper_fields{$pap}) {
        # type of database field and the accession_number e.g. 'PMID:12393910'
        $ref.="|$paper_fields{$pap}:$paper_accno{$pap}"; 
      }
    } elsif (exists $attr{Inferred_automatically}) {
      $ref = "GO_REF:0000002";
      $with = join("|", @{$attr{Inferred_automatically}});
      $with =~ s/INTERPRO/InterPro/g; 
    } else {
      next;
    }
    
    my $this_date = $date;
    if (exists $attr{Date_last_updated}) {
      $this_date = $attr{Date_last_updated};
      $this_date =~ s/\-//g;
    }

    my $a = $aspect{lc $go_type};
    
    &print_wormbase_GAF_line($out,
                             $obj,
                             $public_name,
                             "",
                             $term,
                             $ref,
                             $ec,
                             $with,
                             $a,
                             $seq_name,
                             $taxid,
                             $this_date);
  }
}

close($out);
$db->close;

#
#separate species for gene associations
#
my %coreSpecies = $wormbase->species_accessors;
foreach my $species ($wormbase, values %coreSpecies){
  my $out_file = $output.".".$species->full_name(-g_species => 1);  
  my $taxid = $species->ncbi_tax_id;

  open(my $full_out, $output) or $log->log_and_die("Could not open $output when trying to make per-species files\n");
  open(my $spec_out, ">$out_file") or $log->log_and_die("Could not open $out_file for writing\n");

  &print_wormbase_GAF_header($spec_out);
  while(<$full_out>) {
    next if /^\!/;

    print $spec_out $_ if /\s+taxon:$taxid\s+/;
  }

  close($spec_out) or $log->log_and_die("Could not close $out_file after writing\n");
}

$log->mail;
exit(0);

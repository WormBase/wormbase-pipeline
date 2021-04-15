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
my ($outputdir, $acedbpath, @all_gaf_lines, @nr_gaf_lines, @per_species_files, %dead_genes, $count);

my %opts=();
GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "database:s" => \$acedbpath,
  "outdir:s"   => \$outputdir,
	    )||die(@!);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}


my $log = Log_files->make_build_log($wormbase);
	    
$acedbpath = $wormbase->autoace if not defined $acedbpath;
$tace      = $wormbase->tace;
$outputdir = $wormbase->ontology if not defined $outputdir;

my $date      = &get_GAF_date();

$log->write_to("Fetching gene data...\n") if $debug;
my $gene_info = &get_gene_info( $acedbpath, $wormbase );

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) 
    or $log->log_and_die("Connection failure: ". Ace->error);
$db->date_style('ace');

my (%paper_links, %go_aspect);

$log->write_to("Fetching Paper data...\n") if $debug;
my @aql_results = $db->aql('select a, a->Name, a->Database[2], a->Database[3]  from a in class paper');
foreach my $res (@aql_results) {
  my ($wbpap, $nm, $db, $dbf) = @$res;
  if ($nm and $nm =~ /doi(\S+)/) {
    $paper_links{$wbpap}->{DOI} = $1;
  }
  if ($db and ($db eq 'PMID' or $db eq 'DOI')) {
    $paper_links{$wbpap}->{$db} = $dbf;
  }
}

$log->write_to("Fetching GO type data...\n") if $debug;
@aql_results = $db->aql('select a, a->Type[1] from a in class GO_term');
foreach my $res (@aql_results) {
  next if not $res->[1];
  if (lc($res->[1]) eq 'molecular_function') {
    $go_aspect{$res->[0]} = 'F';
  } elsif (lc($res->[1]) eq 'biological_process') {
    $go_aspect{$res->[0]} = 'P';
  } elsif (lc($res->[1]) eq 'cellular_component') {
    $go_aspect{$res->[0]} = 'C';
  } else {
    $log->log_and_die("Unknown type for $res->[0] : $res->[1]\n");
  }
}

$db->close();

$log->write_to("Processing annotations...\n") if $debug;
# AcePerl/tace combo leaks memory over time, so best to do the work
# in batches, reconnecting for each batch
foreach my $suf (0..9) {
  $log->write_to("Fetching annotations with suffix '${suf}'...\n") if $debug;

  $db = Ace->connect(-path => $acedbpath,  -program => $tace) 
      or $log->log_and_die("(Re)Connection failure: ". Ace->error);
  $db->date_style('ace');
  
  my $it=$db->fetch_many(-query=> "FIND GO_annotation *${suf}");
  
  while (my $obj = $it->next) {
    next unless $obj->isObject();
    print "$obj\n" if ($debug);

    unless ($obj->Gene) {
        $log->write_to("Skipping annotation $obj because there is no Gene connection\n");
        next;
    }

    my $wbgene = $obj->Gene->name;
    if (not exists $gene_info->{$wbgene} or $gene_info->{$wbgene}->{status} eq 'Dead') {
      if (not exists $dead_genes{$wbgene}) {
        $log->write_to("Skipping annotation for $wbgene because could not find a live gene with that id in db\n");
        $dead_genes{$wbgene} = 1;
      }
      next;
    }
    
    my $gaf_line = {
      annot_num   => $obj->name,
      annot_nums  => [$obj->name],
      db          => "WB", 
      gene        => $wbgene,
      go_term     => $obj->GO_term->name,
      go_code     => $obj->GO_code->name,
      aspect      => $go_aspect{$obj->GO_term->name},
      date        => ($obj->Date_last_updated) ? &get_GAF_date($obj->Date_last_updated->name) : $date,
      object_type => "gene",
      symbol      => $gene_info->{$wbgene}->{public_name},
      species     => $gene_info->{$wbgene}->{species},
      isoform     => "",
    };
    
    # DB synonym
    my @synonyms;
    foreach my $name ($gene_info->{$wbgene}->{sequence_name}, @{$gene_info->{$wbgene}->{other_name}}) {
      if ($name ne $gaf_line->{symbol} and $name !~ /^CELE/) {
        push @synonyms, $name;
      }
    }
    $gaf_line->{synonym} = join("|", @synonyms);
    
    # Taxon id
    my @taxids;
    if ($obj->Interacting_species) {
      my ($species_level, $strain_level);
      
      $species_level = $obj->Interacting_species->NCBITaxonomyID;
      if ($obj->Interacting_species(1)) {
        $strain_level = $obj->Interacting_species(1)->NCBITaxonomyID;
      }
      
      if (defined $strain_level) {
        push @taxids, $strain_level->name;
      } elsif (defined $species_level) {
        push @taxids, $species_level->name;
      }
    }  
    $gaf_line->{taxon} = \@taxids;
    
    # qualifier (Annotation_relation)
    my @annot_rel;
    # interpro defaults based on Kimberly
    if ($obj->Motif){
	    my $type;
	    if ("${\$obj->GO_term->Type}" eq 'Molecular_function'){
		    $type = 'enables';
	    }elsif("${\$obj->GO_term->Type}" eq 'Biological_process'){
		    $type = 'involved_in';
	    }elsif("${\$obj->GO_term->Type}" eq 'Cellular_component'){
                    if (grep {"$_" eq 'GO:0032991'} $obj->GO_term->Ancestor){
			    $type = 'part_of';
		    }else{
			    $type = 'located_in';
		    }
	    }
	    $gaf_line->{qualifier} = $type if $type;
    }
    
    if ($obj->Annotation_relation_not) {
      push @annot_rel, "NOT";
      my $al = $obj->Annotation_relation_not->Name;
      $al =~ s/\s+/_/g;
      $al =~ s/,//g; # hopefully only needed for 278 and 279
      if ($al eq 'colocalizes_with' or
          $al eq 'contributes_to') {
        push @annot_rel, "$al";
      }
      $gaf_line->{qualifier} = 'NOT|'.$al;
    }
    if ($obj->Annotation_relation) {
      my $al = $obj->Annotation_relation->Name;
      $al =~ s/\s+/_/g;
      $al =~ s/,//g; # hopefully only needed for 278 and 279
      if ($al eq 'colocalizes_with' or
          $al eq 'contributes_to') {
        push @annot_rel, $al;
      }
      $gaf_line->{qualifier} = "$al";
    }

    # special for WS278, data should be fixed for release WS279+
    unless($gaf_line->{qualifier}){
	    if ($obj->GO_code->name =~/^(IMP|IGI)$/){
		    $gaf_line->{qualifier} = 'acts_upstream_of_or_within';
	    }elsif($obj->GO_code->name =~/^(IDA|ND)$/){
		    $gaf_line->{qualifier} = 'located_in';
	    }
    }

    # that is another bit, that exists only temporary until upstream sorts out the annotation
    if ($obj->GO_term){
	    my $type;
	    if("${\$obj->GO_term->Type}" eq 'Cellular_component'){
                    if (grep {"$_" eq 'GO:0032991'} $obj->GO_term->Ancestor){
			    $type = 'part_of';
		    }else{
			    $type = 'located_in';
		    }
	    }
	    $gaf_line->{qualifier} = $type if $type;
    }
 
    # Reference
    my @reference;
    if ($obj->Reference) {
      my $wbpaper = $obj->Reference->name;
      @reference = ("WB_REF:$wbpaper");
      if (exists $paper_links{$wbpaper}) {
        if (exists $paper_links{$wbpaper}->{PMID}) {
          push @reference, "PMID:" . $paper_links{$wbpaper}->{PMID};
        } elsif (exists $paper_links{$wbpaper}->{DOI}) {
          push @reference, "DOI:" . $paper_links{$wbpaper}->{DOI};
        }
      }
    } elsif ($obj->GO_reference) {
      my ($prefix) = $obj->GO_reference(2);
      my ($suffix) = $obj->GO_reference(3);
      
      push @reference, "${prefix}:${suffix}";
    }
    $gaf_line->{reference} = join("|", @reference);
    
    # 
    # With/from
    #
    my %with_from;
    foreach my $line ($obj->Annotation_made_with) {
      # WormBase entities
      if ($line->name eq 'Interacting_gene' or
          $line->name eq 'RNAi_result' or 
          $line->name eq 'Variation') {
        foreach my $ref ($line->col) {
          $with_from{join(":", "WB", $ref->name)} = 1;
        }
      } elsif ($line->name eq 'Phenotype' or 
               $line->name eq 'Inferred_from_GO_term' or 
               $line->name eq 'Motif') {
        foreach my $ref ($line->col) {
          my $txt = $ref->name;
          $txt =~ s/INTERPRO/InterPro/;
          $with_from{$txt} = 1;
        }
      } elsif ($line->name eq 'Database') {
        foreach my $dblink ($line->col) {
          my $dbname = $dblink->name;
          my $dbfield = $dblink->right->name;
          foreach my $dbv ($dblink->right->col) {
            my $dbvalue = $dbv->name;

            if ($dbname eq 'KEGG') {
              $dbname = "EC";
            } elsif ($dbname eq 'Panther') {
              $dbname = 'PANTHER';
            } elsif ($dbfield eq 'UniProtKB-KW' or 
                     $dbfield eq 'UniProtKB-SubCell' or 
                     $dbfield eq 'UniRule') {
              $dbname = $dbfield;
            } elsif ($dbname =~ /UniProt/) {
              $dbname = "UniProtKB";
            } else {
              # do nothing; dbname is already correct
            }
            $with_from{join(":", $dbname, $dbvalue)} = 1;
          }
        }
      }
    }
    $gaf_line->{with_from} = join("|", sort keys %with_from);
    
    #
    # Contributed_by
    #
    my $cb = 'WB';
    if ($obj->Contributed_by) {
      $cb = $obj->Contributed_by->name;
      if ($cb eq 'WormBase'){$cb = 'WB'}
      elsif($cb=~/UniProt/){$cb = 'UniProt'}
    }else{
	    $cb = 'InterPro' if $obj->Motif && $obj->Motif->name=~/INTERPRO/;
    }
    $gaf_line->{contributor} = $cb;
    
    #
    # Annotation_extention
    #
    my @annot_ext;
    foreach my $line ($obj->Annotation_extension) {
      my $relation = $line->right->name;
      my $subject = $line->right->right;
      my $prefix = "";
      if ($line->name eq 'Gene_relation') {
        $prefix = "WB:";
      } elsif ($line->name eq 'Molecule_relation') {
        my $mol = $subject->fetch();
        $prefix = "WB:";
        foreach my $dblink ($mol->Database) {
          if ($dblink->name eq 'ChEBI') {
            $prefix = "ChEBI:";
            $subject = $dblink->right->right->name;
          }
        }
        $mol->DESTROY();
      }
      push @annot_ext, sprintf("%s(%s%s)", $relation, $prefix, $subject);
    }
    $gaf_line->{annotation_extension} = (@annot_ext) ? [join(",", @annot_ext)] : [];
    
    
    push @all_gaf_lines, $gaf_line;
    
    $log->write_to(sprintf("  processed %d annotations...\n", scalar(@all_gaf_lines)))
        if scalar(@all_gaf_lines) % 10000 == 0 and $debug;
    
    $obj->DESTROY();
  }
  $db->close;
}

$log->write_to("Merging annotations...\n") if $debug;
#
# Now need to sort and merge redundant lines
#
@all_gaf_lines = sort {
  $a->{gene} cmp $b->{gene} or
      $a->{go_term} cmp $b->{go_term} or
      $a->{contributor} cmp $b->{contributor} or
      $a->{reference} cmp $b->{reference} or
      $a->{go_code} cmp $b->{go_code} or
      $a->{with_from} cmp $b->{with_from} or
      $a->{qualifier} cmp $b->{qualifier} or
      $a->{date} cmp $b->{date}
} @all_gaf_lines;

foreach my $gaf (@all_gaf_lines) {
  if (@nr_gaf_lines and
      $nr_gaf_lines[-1]->{gene} eq $gaf->{gene} and
      $nr_gaf_lines[-1]->{go_term} eq $gaf->{go_term} and
      $nr_gaf_lines[-1]->{reference} eq $gaf->{reference} and
      $nr_gaf_lines[-1]->{contributor} eq $gaf->{contributor} and
      $nr_gaf_lines[-1]->{go_code} eq $gaf->{go_code} and
      $nr_gaf_lines[-1]->{with_from} eq $gaf->{with_from} and
      $nr_gaf_lines[-1]->{qualifier} eq $gaf->{qualifier} and
      $nr_gaf_lines[-1]->{date} eq $gaf->{date}) {
    
    if ($debug) {
      $log->write_to(sprintf("Merging redundant annotations: %s into %s:\n", 
                             $gaf->{annot_num}, 
                             $nr_gaf_lines[-1]->{annot_num}));
      $log->write_to(" " . &get_gaf_line($nr_gaf_lines[-1]));
      $log->write_to(" " . &get_gaf_line($gaf));
    }

    push @{$nr_gaf_lines[-1]->{annotation_extension}}, @{$gaf->{annotation_extension}};
    push @{$nr_gaf_lines[-1]->{annot_nums}}, $gaf->{annot_num};
  } else {
    push @nr_gaf_lines, $gaf;
  }
}

$log->write_to("Writing annotations...\n") if $debug;

#
# separate file for each species, and collated file with all species
#
my $collated_file = "$outputdir/gene_association." . $wormbase->get_wormbase_version_name. ".wb" ;
open (my $colfh, ">$collated_file") or $log->log_and_die("Could not open $collated_file for writing\n");
&print_wormbase_GAF_header($colfh, $wormbase->get_wormbase_version_name, 'GO', '2.2');
for my $gaf (@nr_gaf_lines) {
    my $line = &get_gaf_line($gaf);
    print $colfh $line;
}
close($colfh) or $log->log_and_die("Could not close $outputdir/$collated_file after writing\n");

&make_species_files($wormbase, $collated_file);

$log->mail;
exit(0);


#####
sub get_gaf_line {
  my ($gaf) = @_;

  my $line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                     $gaf->{db},
                     $gaf->{gene},
                     $gaf->{symbol},
                     $gaf->{qualifier},
                     $gaf->{go_term},
                     $gaf->{reference},
                     $gaf->{go_code},
                     $gaf->{with_from},
                     $gaf->{aspect},
                     "",
                     $gaf->{synonym},
                     $gaf->{object_type},
                     join("|", map { "taxon:$_" } @{$gaf->{taxon}}),
                     $gaf->{date},
                     $gaf->{contributor},   # was renamed to "Assigned By"
                     join("|", @{$gaf->{annotation_extension}}),
                     $gaf->{isoform});
  
  return $line;
}

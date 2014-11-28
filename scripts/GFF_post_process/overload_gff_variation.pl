#!/usr/bin/env perl
#
# overload_gff_variation.pl
#
# Overloads Variation lines with extra info (consequence etc)
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-11-28 14:30:09 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ( $debug, $test, $store, $wormbase,$database);
my ($datavase, $species, $gff3, $infile, $outfile, %var, $changed_lines, $added_lines);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "database:s" => \$database,
  "gff3"       => \$gff3,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$database = $wormbase->autoace if not defined $database;
$species = $wormbase->species;
my $sp_full_name = $wormbase->full_name;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

&get_molecular_consequences();

my $db = Ace->connect(-path => $database);

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while (<$gff_in_fh>) {
  if (/Variation \"(WBVar\d+)\"/ or 
      /Variation:(WBVar\d+)/ or 
      /variation=(WBVar\d+)/) {
    my $allele = $1; 
    
    my $is_putative_change_of_function_allele = 0;

    my @new_els;
    my @current_els = split(/\t/, $_);
    pop @current_els;
    push @new_els, ['Variation', $allele];

    my $variation = $db->fetch(Variation => $allele);
    
    my @var_types = $variation->at('Variation_type');
    my @other_names = $variation->at('Name.Other_name');
    my ($public_name) = $variation->at('Name.Public_name');
    my $natural_variant = 0;

    if ($public_name) {
      push @new_els, ['Public_name', $public_name];
    }
    if (@other_names) {
      push @new_els, ['Other_name', \@other_names];
    }
    
    if ($variation->Strain) {
      my @strains = $variation->Strain;
      push @new_els, ['Strain', \@strains];
    }
    
    foreach my $tp (@var_types) {
      push @new_els, ['Status', 'Confirmed'] if $tp =~ /confirmed/i;

      push @new_els, ['RFLP'] if $tp =~ /RFLP/;
      if ($tp eq 'Natural_variant') {
        $natural_variant = 1;
        if ($current_els[2] ne 'transposable_element_insertion_site' and 
            $current_els[1] ne 'Polymorphism') {
          $current_els[1] .= "_Polymorphism";
        }
        push @new_els, ['Polymorphism'];
      }
    }
    
    my @types = $variation->at('Sequence_details.Type_of_mutation');
    foreach my $type (@types) {
      if ($type eq 'Substitution' and defined $variation->Substitution) {
        my @substitution = $variation->Substitution->row;
        #print NEW " ; Substitution \"$substitution[0]/$substitution[1]\"";
        push @new_els, ['Substitution', join("/", uc($substitution[0]), uc($substitution[1]))];
      } elsif ($type eq 'Insertion' and defined $variation->Insertion) {
        my ($insert) = $variation->Insertion->row;
        push @new_els, ['Insertion', uc($insert)];
      }
    }

    if ($current_els[2] eq 'sequence_alteration') {
      # general term used when set contains a mixture of substitutions and indels.
      # Try to make it more specific here
      my $new_term = "sequence_alteration";
      if (scalar(@types) == 1) {
        my ($tp) = @types;
        
        if ($tp eq 'Insertion' or $tp eq 'Deletion' or $tp eq 'Inversion') {
          $new_term = lc($tp);
          $new_term .= "_site" if $new_term eq 'insertion';
        } elsif ($tp eq 'Substitution') {
          if ($current_els[3] == $current_els[4]) {
            # single nucleotide
            if ($natural_variant) {
              $new_term = 'SNP';
            } else {
              $new_term = 'point_mutation';
            }
          } else {
            $new_term = 'substitution';
          }
        } elsif ($tp eq 'Tandem_duplication') {
          $new_term = 'tandem_duplication';
        }
      } elsif (scalar(@types) > 1) {
        if (grep { $_ eq 'Tandem_duplication' } @types) {
          # some tandem duplications are assoiated with micro-insertions/deletions at the flanks, 
          # so will have multiple Type_of_mutation. We still want to classify them as tandem_duplications
          $new_term = "tandem_duplication";
        } else {
          # more than one type, so just put complex_substitution
          $new_term = "complex_substitution";
        }
      }
      
      $current_els[2] = $new_term;
    }
    
    if (exists $var{$allele}->{mol_change}) {
      my $cons = $var{$allele}->{mol_change};

      if ($cons eq 'Frameshift' or
          $cons eq 'Missense' or
          $cons eq 'Nonsense' or
          $cons eq 'Readthrough' or
          $cons eq 'Coding_exon') {
        if ($current_els[2] ne 'transposable_element_insertion_site' and 
            $current_els[2] ne 'tandem_duplication') {
          $is_putative_change_of_function_allele = 1;
        }
      }

      push @new_els, ['Consequence', $cons];
      
      my @consequences = $variation->at('Affects.Predicted_CDS[2]');
      
      if (grep { $_ =~ /Missense|Nonsense|Readthrough/ } @consequences){
        # 2.) the missense
        my @missense = $variation->at('Affects.Predicted_CDS[4]');
        @missense = grep {/\sto\s/} @missense;
        #print NEW " ; AAChange \"$missense[0]\"";
        push @new_els, ['AAChange', $missense[0]];
      }
    }
    
    my @new_el_strings;
    foreach my $el (@new_els) {
      if (scalar(@$el) == 1) {
        if ($gff3) {
          push @new_el_strings, sprintf("%s=1", lc($el->[0]));
        } else {
          push @new_el_strings, $el->[0];
        }
      } else {
        my ($k, $v) = @$el;
        if ($gff3) {
          if (ref($v) eq 'ARRAY') {
            $v = join(",", @$v);
          } 
          push @new_el_strings, join("=", lc($k), $v);
        } else {
          if (ref($v) eq 'ARRAY') {
            foreach my $vv (@$v) {
              push @new_el_strings, sprintf("%s \"%s\"", $k, $vv);
            }
          } else {
            push @new_el_strings, sprintf("%s \"%s\"", $k, $v);
          }
        }
      }
    }
    
    my $join_str = ($gff3) ? ";" : " ; ";
    my $group = join($join_str, @new_el_strings); 
    print $gff_out_fh join("\t", @current_els, $group), "\n";

    $changed_lines++;
        
    # duplicate putative change-of-function alleles, with a different source, so that they
    # can be drawn in a separate track
    if ($is_putative_change_of_function_allele) {
      $current_els[1] = "PCoF_" . $current_els[1];
      print $gff_out_fh join("\t", @current_els, $group), "\n";
      $added_lines++;
    }
    
    if ($changed_lines % 50000 == 0) {
      $db->close();
      $db = Ace->connect(-path => $database);
    }
  } else {
    print $gff_out_fh "$_";
  }
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$db->close();

$log->write_to("Finished processing : $changed_lines lines modified, $added_lines added\n");
$log->mail();
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub get_molecular_consequences {
  my %interested = ('Genomic_neighbourhood' => 1,
                    'Regulatory_feature'    => 2,
                    'Promoter'              => 3,
                    'UTR_5'                 => 4,
                    'UTR_3'                 => 5,
                    'Intron'                => 6,
                    'Coding_exon'           => 7,
                    'Silent'                => 8,
                    'Splice_site'           => 9,
                    'Readthrough'           => 10,
                    'Nonsense'              => 11,
                    'Frameshift'            => 12,
                    'Missense'              => 13,
		      );

  my $table = $wormbase->table_maker_query($database, &write_mol_change_def_file);
  
  while(<$table>) {
    chomp;
    s/\"//g; #"
    next if (! defined $_);
    next if (/acedb/ or /\/\//);
    my ($var_name, @mut_affects) = split(/\s+/,$_);
    next if not defined $var_name;
    
    foreach my $mut_affects (@mut_affects) {
      
      if($mut_affects and $interested{$mut_affects}){
        if( not exists $var{$var_name}->{mol_change} or 
            $interested{$mut_affects} > $interested{ $var{$var_name}->{mol_change} }) {
          $var{$var_name}->{mol_change} = $mut_affects;
        }
      }
    }
  }
  close($table);
}


sub write_mol_change_def_file {
  my $def = '/tmp/overload_SNP_GFF_mol_chng.def';
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Variation 
From 1 

Colonne 2
Width 12
Mandatory
Hidden
Show_Tag
From 1
Tag Sequence

Colonne 3
Width 12 
Mandatory 
Hidden 
Class 
Class CDS 
From 1 
Tag Predicted_CDS 
 
Colonne 4 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Missense 
 
Colonne 5 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Nonsense 
 
Colonne 6 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Splice_site 
 
Colonne  7
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Frameshift 
 
Colonne 8 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Intron 
 
Colonne 9 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Coding_exon 
 
Colonne 10 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Promoter 
 
Colonne 11 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # UTR_3 
 
Colonne 12 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # UTR_5 
 
Colonne 13 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Regulatory_feature 
 
Colonne 14 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3  
Tag  HERE  # Genomic_neighbourhood 
 
Colonne 15 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3  
Tag  HERE  # Silent 

Colonne 16
Width 12
Optional
Visible
Show_Tag
Right_of 3
Tag  HERE  # Readthrough

END

  print TMP $txt;
  close TMP;
  return $def;
}

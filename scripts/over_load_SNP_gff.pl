#!/software/bin/perl -w                 
#
# This scripts decorates the variations in the GFF with all sorts of extra useful info
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-01-14 14:11:14 $      

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

my ( $debug, $test, $verbose, $store, $wormbase,$database);
my ($species, $gff_file, $gff_v3, %var);

GetOptions (
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database,
	    "file:s"     => \$gff_file,
            "gff3"       => \$gff_v3,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$wormbase->{autoace} = $database if $database;
$species = $wormbase->species;
my $sp_full_name = $wormbase->full_name;

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################


#
# Get molecular consequences for all vars
#
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

my $table = $wormbase->table_maker_query($wormbase->autoace, &write_mol_change_def_file);

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


my $stat = 0;

my @gff_files;
if ($gff_file) {
  if (not -e $gff_file or -z $gff_file) {
    $log->log_and_die("Non-existent or zero length GFF file");
  }
  @gff_files = ($gff_file);
} else {
  if ($wormbase->assembly_type eq 'contig'){
    @gff_files = ($wormbase->species);
  } else {
    @gff_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
  }
  for(my $i=0; $i < @gff_files; $i++) {
    $gff_files[$i] = sprintf("%s/%s.gff", $wormbase->gff, $gff_files[$i]);
    if (not -e $gff_files[$i] or -z $gff_files[$i]) {
      $log->log_and_die("Non-existent or zero-length GFF file $gff_files[$i]");
    }
  }
}

foreach my $file (@gff_files) {
  # open and close db connection for each chrom, since AcePerl/giface seems to 
  # hang on to memory...
  my $db = Ace->connect(-path => $wormbase->autoace);

  open(GFF,"<$file") or $log->log_and_die("cant open $file");
  open(NEW,">$file.tmp") or $log->log_and_die("cant open $file tmp file\n");
  while( <GFF> ) {
    if (/Variation \"(\S+)\"/ or /Variation:(\S+)/) {
      my $allele = $1;
      
      my @new_els;
      my @current_els = split(/\t/, $_);
      pop @current_els;
      push @new_els, ['Variation', $allele];

      my $variation = $db->fetch(Variation => $allele);

      my @var_types = $variation->at('Variation_type');
      my @other_names = $variation->at('Name.Other_name');
      my @public_names = $variation->at('Name.Public_name');
      my $natural_variant = 0;

      foreach my $tp (@var_types) {
        push @new_els, ['Status', 'Confirmed'] if $tp =~ /confirmed/i;
        push @new_els, ['Status', 'Predicted'] if $tp =~ /predicted/i;
        push @new_els, ['RFLP'] if $tp =~ /RFLP/;
        if ($tp eq 'Natural_variant') {
          $natural_variant = 1;
        }
      }

      foreach my $pb (@public_names) {
        #print NEW " ; Public_name \"$pb\"";
        push @new_els, ['Public_name', $pb];
      }
      foreach my $on (@other_names) {
        #print NEW " ; Other_name \"$on\"";
        push @new_els, ['Other_name', $on];
      }

      my @types = $variation->at('Sequence_details.Type_of_mutation');
      foreach my $type (@types) {
        if ($type eq 'Substitution' and defined $variation->Substitution) {
          my @substitution = $variation->Substitution->row;
          #print NEW " ; Substitution \"$substitution[0]/$substitution[1]\"";
          push @new_els, ['Substitution', join("/", uc($substitution[0]), uc($substitution[1]))];
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
          }
        } elsif (scalar(@types) > 1) {
          # more than one type, so just put complex_substitution
          $new_term = "complex_substitution";
        }

        $current_els[2] = $new_term;
      }
        
      if (exists $var{$allele}->{mol_change}) {
        push @new_els, ['Consequence', $var{$allele}->{mol_change}];
        
        my @consequences = $variation->at('Affects.Predicted_CDS[2]');
        
        if (grep { $_ =~ /Missense|Nonsense|Readthrough/ } @consequences){
          # 2.) the missense
          my @missense = $variation->at('Affects.Predicted_CDS[4]');
          @missense = grep {/to/} @missense;
          #print NEW " ; AAChange \"$missense[0]\"";
          push @new_els, ['AAChange', $missense[0]];
        }
      }
    
      # 3.) the strain
      #print NEW " ; Strain \"${\$variation->Strain}\"" if $variation->Strain;
      push @new_els, ['Strain', $variation->Strain] if $variation->Strain;

      my @new_el_strings;
      foreach my $el (@new_els) {
        if (scalar(@$el) == 1) {
          push @new_el_strings, $el->[0];
        } else {
          if ($gff_v3) {
            push @new_el_strings, join(":", @$el);
          } else {
            push @new_el_strings, sprintf("%s \"%s\"", @$el);
          }
        }
      }

      my $join_str = ($gff_v3) ? ";" : " ; ";
      my $group = join($join_str, @new_el_strings); 
      print NEW join("\t", @current_els, $group), "\n";

      $stat++;

      if ($stat % 50000 == 0) {
        $db->close();
        $db = Ace->connect(-path => $wormbase->autoace);
      }
    } else {
      print NEW "$_";
    }
  }
  $wormbase->run_command("mv -f $file.tmp $file", $log);

  $db->close();
}

##################
# Check the files
##################

if($gff_file) {
  $log->write_to("Not checking ad hoc file\n");
} else { 
  foreach my $file (@gff_files) {
    my $minsize = ($file=~/random|un/)?170000:1500000;
    $wormbase->check_file($file, $log,
                          minsize => $minsize,
                          lines => ['^##',
                                    "^\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
        );
  }
}
  

# Close log files and exit
$log->write_to("\n\nChanged $stat lines\n");
$log->write_to("----------\n\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

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

__END__
# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


=pod

=head2 NAME - over_load_SNP_gff.pl

=head1 USAGE

=over 4

=item over_load_SNP_gff.pl  [-options]

=back

This script adds the status (confirmed / predicited) of SNPs and whether are RFLPs to GFF lines.  This is so that the web
page can easily distinguish between them

=over 4

=item None at present.

=back

$0.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

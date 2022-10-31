#!/usr/bin/env perl

=pod

=head1 SYNOPSIS

  prepare_conf_and_fasta.pl [yaml_file]

=head1 DECRIPTION

Examines GFF3 and FASTA in an assembly directory created by 
C<prepare_data_folder.pl>. Checks GFF3 sources are correct (WIP: could add more 
validation/verification). Splits FASTA and creates AGP file if required.

Also writes a YAML configutation file suitable for data import scripts, 
including C<worm_lite.pl>.  This configutation is similar to the input, but with 
some additional parameters that describe the GFF3, FASTA and (if applicable) AGP 
data.

Input configuration is provided as YAML.  This is in the format written to 
standard output by C<prepare_data_folder.pl>.  This must have a single top-level 
element that has the same name as the assembly data directory.

Checks the newly created YAML configuration file is complete, and if so
concatenates it to standard output.

Conventions:

=over

=item *

The root of the assembly data directories is defined by PARASITE_DATA in  the user
environment

=item *

The directory for this assembly is named according to the ParaSite naming convention:
C<genus_species_bioproject>

=item *

The FASTA file is also named according to the ParaSite naming convention:
C<genus_species_bioproject.fa>

=back

These will all be fine if you used C<prepare_data_folder.pl> to create the 
assembly directory.

=head1 OPTIONS AND ARGUMENTS

Reads assembly metadata as YAML from standard input or named file.

=over

=item force

force writing of new configuration file, even if it already exists, and (if 
splitting required) rewrting of FASTA file & weriting of AGP file. Same as 
setting REDO_FASTA in the environment

=item split_fasta

Splits FASTA and creates AGP I<when required>. Default B<true>, negatable with 
C<-nosplit_fasta>; can also be negated by setting SKIP_SPLIT_FASTA in the 
environment

=item set

Set a value in the configuration.  Can be repeated to set multiple
values.  For example:

   -set "mitochondrial=mt_seq_id" -set "taxon_id=1234"

(Note that if a configuration file already exists, you need to use
C<-force> to have your new value(s) written to the file).

To set metadata values, see C<-meta>.

=item meta

Set a metadata value in the configuration.  Can be repeated to set multiple
values.  For example:

   -meta "species.nematode_clade=IV" -meta "provider.url=https://foo.bar"

(Note that if a configuration file already exists, you need to use
C<-force> to have your new value(s) written to the file).

=item delete_meta

Deletes a metadata item from the configuration.  Can be repeated to set multiple 
items. It is not an error if the item doesn't exist.
For example:

   -delete_meta "species.strain"

(Note that if a configuration file already exists, you need to use
C<-force> to have your new value(s) written to the file).

=item help

print this message and exit

=back

=head1 TO DO

=over

=item *

Stop using deprecated YAML package (or at least use it as YAML::Old)

=item *

Validation of input

=item *

More thorough verification of data files

=item *

Rationalise -set, -meta and -delete_meta options

=back

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::IO::GFFParser;
use Carp;
use CoreCreation::Config::Utils;
use CoreCreation::Fasta;
use CoreCreation::GFF3::Utils;
use File::Basename;
use File::Spec;
use File::Slurp;
use Getopt::Long;
use Pod::Usage;
use Storable;
use Try::Tiny;
use YAML;

my($force, $split_fasta, @user_supplied, @user_supplied_metadata, @metadata_to_delete, $help);
$split_fasta = 1; # negatable option
GetOptions( 'force'           => \$force,
            'split_fasta!'    => \$split_fasta,
            'set=s'           => \@user_supplied,
            'meta=s'          => \@user_supplied_metadata,
            'delete_meta=s'   => \@metadata_to_delete,
            'help'            => \$help
            )
            || pod2usage({-exitval=>1});
$help && pod2usage({-verbose=>2, -exitval=>0});

my $conf;
# read STDIN or from a named file, in the standard Perl fashion,
# and catch YAML::Load errors
# my $conf = Load(do {local $/; <>} || "--- {}");
# die "usage: PARASITE_DATA=... $0 < conf" unless $conf;
try{
   $conf = Load join('',<>);
} catch {
   croak "Error parsing input YAML at line $.\n".$_;
};
pod2usage(255) unless $conf;

my ($data_dir_name, @others) = keys %{$conf};

croak "Input must contain exactly one data directory name, but all these were found:\n$data_dir_name @others\n"
   unless $data_dir_name and not @others and CoreCreation::Config::Utils::parasite_data_id_is_valid($data_dir_name);
my $data_dir_path = join ("/", $ENV{PARASITE_DATA}, $data_dir_name);
my $conf_path = File::Spec->catfile($data_dir_path, "$data_dir_name.conf");

foreach my $u (@user_supplied) {
   my($k,$v) = split(/=/, $u, 2);
   croak "Badly formed -set argument: \"$u\"" unless $k && defined $v;
   $conf->{$data_dir_name}->{$k} = $v;
}

foreach my $d (@metadata_to_delete) {
   exists $conf->{$data_dir_name}->{meta}->{$d} && delete $conf->{$data_dir_name}->{meta}->{$d};
}

foreach my $m (@user_supplied_metadata) {
   my($k,$v) = split(/=/, $m, 2);
   croak "Badly formed -meta argument: \"$m\"" unless $k && defined $v;
   $conf->{$data_dir_name}->{meta}->{$k} = $v;
}

if( $force or not -s $conf_path or $ENV{REDO_FASTA} ) {
   # less confusing (?) to give this reference a new name rather than reassigning $conf
   # also create copy of hash so we're not rewriting the input params
   # $conf = $conf->{$data_dir_name};
   my $this_assembly = Storable::dclone $conf->{$data_dir_name};

   # check existence of GFF3 file in specified location
   $this_assembly->{gff3} //=  File::Spec->catfile($data_dir_path, "$data_dir_name.gff3");
   croak "Didn't find expected GFF3 file at ".$this_assembly->{gff3} unless -f $this_assembly->{gff3} and File::Spec->file_name_is_absolute($this_assembly->{gff3});

   # validate GFF3
   my($valid_gff3, @validation_errors);
   try {
      ($valid_gff3, @validation_errors) = CoreCreation::GFF3::Utils::validate($this_assembly->{gff3});
   } catch {
      # CoreCreation::GFF3::Utils::validate dies if there's an error in running the validation, not an indication of invalid GFF3
      croak "Unable to validate GFF3:\n".$_;
   };
   # exit on GFF3 validation error
   $valid_gff3 || print "$this_assembly->{gff3} is not valid GFF3:\n".termcap_bold(join("\n",@validation_errors))."\n  ";
   
   # my $check_sources_column = "grep -c $this_assembly->{gff_sources} $this_assembly->{gff3}";
   # die "Failed: $check_sources_column" unless 0 < `$check_sources_column`;
   
   # limited verification of GFF3 data to provide early warnings of major problems
   # (most verification during database load)
   {
      # load features
      my @features;
      try {
         open(GFF3,$this_assembly->{gff3}) || die "can't read file: $!";
         my $gff3_parser =  Bio::EnsEMBL::Utils::IO::GFFParser->new(\*GFF3)
                           || die "failed to create Bio::EnsEMBL::Utils::IO::GFFParser";
         $gff3_parser->parse_header(); # discard headers
         while( my $f = $gff3_parser->parse_next_feature() ) {
            push(@features, $f);
         }
         close(GFF3) || die "error whilst reading file: $!";
      } catch {
         croak "Error whilst parsing GFF3 file $this_assembly->{gff3}:\n".termcap_bold($_);
      };

      # sources: verify they match the config
      my @sources = split(/,/, $this_assembly->{gff_sources});
      FEAT: foreach my $feature (@features){
	foreach my $source (@sources){
	   next FEAT if $feature->{source} eq $source;
	}
	die qq~incorrect GFF source "$feature->{source}" (expected values are : "@sources") at $this_assembly->{gff3} line $.\n~;
      }

      my %gff = (gene=>[], pseudogene=>[]);
      
      # genes and pseudogenes: verify uniqueness
      my %genes = ();
# using type 'gene' instead of pseudogene
#       my %pseudogenes = ();
      foreach my $feature (@features){
         if( 'gene' eq $feature->{type} ) {
            my $hash = { ID=>$feature->{attribute}->{ID}, RNA=>[] };
            die qq~gene ID $feature->{attribute}->{ID} is not unique at $this_assembly->{gff3} line $.\n~ if exists $genes{ $feature->{attribute}->{ID} };
            push( @{$gff{gene}}, $hash);
            $genes{ $feature->{attribute}->{ID} } = $hash;
         } 
# using type 'gene' instead of pseudogene
#          elsif( 'pseudogene' eq $feature->{type} ) {
#             my $hash = { ID=>$feature->{attribute}->{ID}, pseudogenic_transcript=>[] };
#             die qq~pseudogene ID is not unique at $this_assembly->{gff3} line $.\n~ if exists $pseudogenes{ $feature->{attribute}->{ID} };
#             push( @{$gff{pseudogene}}, $hash);
#             $pseudogenes{ $feature->{attribute}->{ID} } = $hash;
#          }
      }
      # RNAs: verify uniqueness and gene parentage, link to parent
      # pseudogenic_transcripts: verify uniqueness and gene parentage, link to parent
      my %RNAs = ();
      my %pseudogenic_transcripts = ();
      foreach my $feature (@features){
         my $hash = {ID=>$feature->{attribute}->{ID}, type=>$feature->{type}, Parent=>$feature->{attribute}->{Parent}, CDS=>[], exon=>[] };
         if( 'mRNA' eq $feature->{type} || 'transcript' eq $feature->{type} || $feature->{type} =~ m/^..?RNA$/ || 'nontranslating_transcript' eq $feature->{type}) {
            die "$feature->{type} ID $feature->{attribute}->{ID} is not unique at $this_assembly->{gff3} line $.\n"
               if exists $RNAs{ $feature->{attribute}->{ID} };
            # parent must exists *and* must be a gene
            die "RNA $feature->{attribute}->{ID} parent doesn't exist or is not a gene in $this_assembly->{gff3}"
               unless exists $genes{ $feature->{attribute}->{Parent} };
            push( @{$genes{ $feature->{attribute}->{Parent} }->{RNA}}, $hash);
            $RNAs{ $feature->{attribute}->{ID} } = $hash;
         }
         if( 'pseudogenic_transcript' eq $feature->{type} ) {
            die "pseudogenic_transcript ID $feature->{attribute}->{ID} is not unique at $this_assembly->{gff3} line $.\n"
               if exists $pseudogenic_transcripts{ $feature->{attribute}->{ID} };
            # parent must exists *and* must be a pseudogene
# using type 'gene' instead of pseudogene
#             die "pseudogenic_transcript $feature->{attribute}->{ID} parent doesn't exist or is not a pseudogene in $this_assembly->{gff3}"
#                unless exists $pseudogenes{ $feature->{attribute}->{Parent} };
#             push( @{$pseudogenes{ $feature->{attribute}->{Parent} }->{pseudogenic_transcript}}, $hash);
            die "pseudogenic_transcript $feature->{attribute}->{ID} parent doesn't exist or is not a gene in $this_assembly->{gff3}"
               unless exists $genes{ $feature->{attribute}->{Parent} };
            push( @{$genes{ $feature->{attribute}->{Parent} }->{pseudogenic_transcript}}, $hash);
            $pseudogenic_transcripts{ $feature->{attribute}->{ID} } = $hash;
         }
      }
      # CDS: verify RNA parentage, link to parent
      # exon: verify RNA *or* pseudogenic transcript parentage, link to whichever is the parent
      foreach my $feature (@features){
         if( 'CDS' eq $feature->{type} ) {
            die "CDS $feature->{attribute}->{ID} references a non-existent parent in $this_assembly->{gff3}"
               unless exists $RNAs{ $feature->{attribute}->{Parent} };
            push( @{$RNAs{ $feature->{attribute}->{Parent} }->{ $feature->{type} }},
                  $feature->{attribute}->{ID});
         } elsif( 'exon' eq $feature->{type} ) {
            if( exists $RNAs{ $feature->{attribute}->{Parent} } ) {
               push( @{$RNAs{ $feature->{attribute}->{Parent} }->{ $feature->{type} }},
                     'exon');
            } elsif( exists $pseudogenic_transcripts{ $feature->{attribute}->{Parent} } ) {
               push( @{$pseudogenic_transcripts{ $feature->{attribute}->{Parent} }->{ $feature->{type} }},
                     'exon');
            } else {
               die "exon $feature->{attribute}->{ID} references a non-existent parent in $this_assembly->{gff3}"
            }
         }
      }
      # at this point the gene -> RNA -> CDS/exon structure should be complete
      # *and* it's each child's parent has been verified
      # => next check that each gene has child(ren)
      foreach my $gene_ID (keys %genes) {
         die "gene $gene_ID has no RNA/pseudogenic transcript in $this_assembly->{gff3}" unless $genes{$gene_ID}->{RNA}->[0] || $genes{$gene_ID}->{pseudogenic_transcript}->[0];
         foreach my $this_RNA (@{$genes{$gene_ID}->{RNA}}) {
            my $RNA_ID = $this_RNA->{ID};
            die "RNA $RNA_ID has no CDS or exon in $this_assembly->{gff3}" unless $RNAs{$RNA_ID}->{CDS}->[0] || $RNAs{$RNA_ID}->{exon}->[0];
         }
      }
# using type 'gene' instead of pseudogene
#       foreach my $pseudogene_ID (keys %pseudogenes) {
#          die "pseudogene $pseudogene_ID has no pseudogenic_transcript in $this_assembly->{gff3}" unless $pseudogenes{$pseudogene_ID}->{pseudogenic_transcript}->[0];
#          foreach my $this_pseudogenic_transcript (@{$pseudogenes{$pseudogene_ID}->{pseudogenic_transcript}}) {
#             my $pseudogenic_transcript_ID = $this_pseudogenic_transcript->{ID};
#             die "pseudogenic transcript $pseudogenic_transcript_ID has a CDS $this_assembly->{gff3}" if $pseudogenic_transcripts{$pseudogenic_transcript_ID}->{CDS}->[0];
#             die "pseudogenic transcript $pseudogenic_transcript_ID has no exon $this_assembly->{gff3}" unless $pseudogenic_transcripts{$pseudogenic_transcript_ID}->{exon}->[0];
#          }
#       }
   }

   # check existence of FASTA file in specified location
   my $fasta_path = File::Spec->catfile($data_dir_path,"$data_dir_name.fa");
   croak "Didn't find expected FASTA file at $fasta_path" unless -f $fasta_path and File::Spec->file_name_is_absolute($fasta_path);
   $this_assembly->{toplevel} = "scaffold";
   my $fasta = CoreCreation::Fasta->new($fasta_path) || croak "Failed to create CoreCreation::Fasta";
   if($split_fasta and $fasta->needs_contig_structure and not $ENV{SKIP_SPLIT_FASTA}){
      (my $split_fasta_path   = $fasta_path) =~ s/.fa/.split.fa/;
      (my $agp_path           = $fasta_path) =~ s/.fa/.toplevel.agp/;
      $fasta->split(fasta => $split_fasta_path, agp => $agp_path);
      $this_assembly->{fasta}    = $split_fasta_path;
      $this_assembly->{agp}      = $agp_path;
      $this_assembly->{seqlevel} = "contig";
   } else {
      $this_assembly->{fasta}    = $fasta_path;
      $this_assembly->{seqlevel} = "scaffold";
   }
   my $mito = $fasta->mito;
   $this_assembly->{mitochondrial} = $mito if $mito;
  
  open(FH, '>', $conf_path) or croak "Can't write to $conf_path: $!";
  print FH Dump({$data_dir_name => $this_assembly});
  close(FH) or croak "Error whilst writing $conf_path: $!";
}

# this check for missing values is a wee bit fragile (any occurence of a
# CoreCreation::Config::Utils::MISSING_METADATA_PATTERN (used to be hardcoded
# as '?')anywhere in the text)
# is isn't extensible for additional validation/verification
# the slurping as text and a separated YAML::Load->YAML::Dump is a bit odd too; is this
# meant to check it is valid YAML?
# my $text = File::Slurp::read_file($conf_path);
# if ($text =~ /\?/) {
#   die "$conf_path: complete the config file and remove all the ?s!";
# } 
# print YAML::Dump (YAML::LoadFile($conf_path)); 

my $new_conf = YAML::LoadFile($conf_path) or croak "YAML parser barfed on $conf_path";

# check configuration for missing values
my $missing = 0;
my $flat = CoreCreation::Config::Utils::flatten_hash(Storable::dclone $new_conf->{$data_dir_name});
while (my ($conf_key, $conf_value) = each %{$flat}) {
   if( CoreCreation::Config::Utils::MISSING_METADATA_PATTERN eq $conf_value ) {
      ++$missing;
      print "ERROR: configuration has a missing value for ".termcap_bold($conf_key)."\n";
   }
}
die "To proceed further, provide the missing value".($missing>1?'s':'')." run again. Tip: you can rerun ".basename($0)." using \$PARASITE_DATA/$data_dir_name/".basename($conf_path)." as input (use --force is providing metadata values with --meta).\n"
   if $missing;

# configuration checked: print
print YAML::Dump $new_conf;

sub termcap_bold {
   my @input = @_;
   my @bold = ( `tput bold`,@input,`tput sgr0` );
   return wantarray ? @bold : "@bold";
}

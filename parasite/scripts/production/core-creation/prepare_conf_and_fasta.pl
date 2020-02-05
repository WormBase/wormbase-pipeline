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

=item help

print this message and exit

=back

=head1 TO DO

=over

=item *

Validation of input

=item *

More thorough verification of data files

=item *

Create package to eliminate code (e.g. species name filters) duplicated
between ParaSite scripts used for NCBI checking & import etc.

=back

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::IO::GFFParser;
use Carp;
use CoreCreation::Fasta;
use File::Basename;
use File::Spec;
use File::Slurp;
use Getopt::Long;
use Pod::Usage;
use Storable;
use Try::Tiny;
use YAML;

use constant GENOME_NAME_FILTER  => qr/^[a-z\d]+_[a-z\d]+_[a-z]+\d+$/;

my($force, $split_fasta, $help);
$split_fasta = 1; # negatable option
GetOptions( 'force'        => \$force,
            'split_fasta!' => \$split_fasta,
            'help'         => \$help
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

croak "Input must contain exactly one data directory name, but all these were found:\n$data_dir_name @others\n" unless $data_dir_name and not @others and $data_dir_name =~ GENOME_NAME_FILTER;
my $data_dir_path = join ("/", $ENV{PARASITE_DATA}, $data_dir_name);
my $conf_path = File::Spec->catfile($data_dir_path, "$data_dir_name.conf");

if( $force or not -s $conf_path or $ENV{REDO_FASTA} ) {
   # less confusing (?) to give this reference a new name rather than reassigning $conf
   # also create copy of hash so we're not rewriting the input params
   # $conf = $conf->{$data_dir_name};
   my $this_assembly = Storable::dclone $conf->{$data_dir_name};

   # check existence of GFF3 file in specified location
   $this_assembly->{gff3} //=  File::Spec->catfile($data_dir_path, "$data_dir_name.gff3");
   croak "Didn't find expected GFF3 file at ".$this_assembly->{gff3} unless -f $this_assembly->{gff3} and File::Spec->file_name_is_absolute($this_assembly->{gff3});

   # WIP: verify GFF3 data 
   # note Bio::EnsEMBL::Utils::IO::GFFParser doesn't validate, but if the file is
   # badly mangled this will throw some sort of error
   # my $check_sources_column = "grep -c $this_assembly->{gff_sources} $this_assembly->{gff3}";
   # die "Failed: $check_sources_column" unless 0 < `$check_sources_column`;
   try {
      open(GFF3,$this_assembly->{gff3}) || die "can't read file: $!";
      my $gff3_parser =  Bio::EnsEMBL::Utils::IO::GFFParser->new(\*GFF3)
                        || die "failed to create Bio::EnsEMBL::Utils::IO::GFFParser";
      $gff3_parser->parse_header(); # discard headers
      while( my $feature = $gff3_parser->parse_next_feature() ) {
         die qq~$feature->{seqid} has incorrect GFF source "$feature->{source}" (expected "$this_assembly->{gff_sources}")\n ~
            if $feature->{source} ne $this_assembly->{gff_sources};
      }
      close(GFF3) || die "error whilst reading file: $!";
   } catch {
      croak "Error whilst parsing GFF3 file $this_assembly->{gff3}:\n".$_;
   };

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

my $text = File::Slurp::read_file($conf_path);
if ($text =~ /\?/) {
  die "$conf_path: complete the config file and remove all the ?s!";
} 
# not sure why this rewrites the file, and if that's required then why
# $text isn't just written to the file rather than using YAML package..?
print YAML::Dump (YAML::LoadFile($conf_path)); 

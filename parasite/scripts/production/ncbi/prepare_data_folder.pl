#!/usr/bin/env perl

=pod

=head1 SYNOPSIS

  prepare_data_folder [options] [yaml_file]

=head1 DESCRIPTION

Creates local directory into which an assembly will be downloaded from NCBI. 
Local files for the ParaSite import are created in the same directory.  The 
directory name is based on species and bioproject names according to the 
ParaSite convention; it is created under the data root directory defined by 
PARASITE_DATA in the user environment.

Input data are provided as YAML.  This is in the format produced by 
doc_for_assembly.pl, which is a scratch XML-to-YAML conversion based on NCBI 
reports for an assembly.

Writes YAML to standard output that provides input data for
prepare_conf_and_fasta.pl.

=head1 OPTIONS AND ARGUMENTS

Reads assembly metadata as YAML from standard input or named file.

=over

=item force

force downloading of assembly from NCBI and creation of local files, even
if they already exist

=item gff3_download

attempt to download GFF3 (as well as FASTA) from NCBI

=item help

print this message and exit

=back

=head1 TO DO

=over

=item * Stop using deprecated YAML package (or at least use it as YAML::Old)

=item * Validation of input

=item * Create package to eliminate code (e.g. species name filters) duplicated
        between ParaSite scripts used for NCBI checking & import etc.
        
=item * More SubmitterOrganization -> URL mappings?
  
=back

=cut


use strict;
use warnings;

use Carp;
use CoreCreation::Config::Utils;
use File::Basename;
use File::Path;
use Getopt::Long;
use IO::Uncompress::Gunzip;
use Net::FTP;
use Pod::Usage;
use Try::Tiny;
use YAML;

use constant NCBI_FTP_SERVER  => 'ftp.ncbi.nlm.nih.gov';

my($force, $gff3_download, $help);
GetOptions( 'force'           => \$force,
            'gff3_download'   => \$gff3_download,
            'help'            => \$help
            )
            || pod2usage({-exitval=>1});
$help && pod2usage({-verbose=>2, -exitval=>0});

my $conf;
# read STDIN or from a named file, in the standard Perl fashion,
# and catch YAML::Load errors
# $conf = Load(do {local $/; <STDIN>} || "--- {}") unless -t STDIN;
try{
   $conf = Load join('',<>);
} catch {
   croak "Error parsing input YAML at line $.\n".$_;
};
pod2usage(255) unless $conf;

my ($species,$spe,$cies) = CoreCreation::Config::Utils::parasite_data_id($conf->{SpeciesName}, $conf->{GB_BioProjects}{Bioproj}{BioprojectAccn});

# this variable is treated as a string, but can be read as a HASH
# (on incidences I have seen are empty hashes, but I suspect this is just an artefact of
# they way the YAML is created running YAML::Load on a structure created by feeding
# NCBI XML to XML::Simple::XMLin...)
# better to treat this as an empty string when this happens
# this may not be the desired effect, but it's hard to see how it can be less wrong than
# doing string matches on hash references :-/
# Tim S. 2020-01-04
my $isolate = $conf->{Biosource}{Isolate};
$isolate = undef if ref($isolate) && ref({}) eq ref($isolate);

# goto considered harmful
# goto PRINT unless $ENV{PARASITE_DATA};
if( $ENV{PARASITE_DATA} ) {
   
   -d $ENV{PARASITE_DATA} || croak "environment variable PARASITE_DATA has been set to a non-existent directory: $ENV{PARASITE_DATA}";
   -w $ENV{PARASITE_DATA} || croak "environment variable PARASITE_DATA set to $ENV{PARASITE_DATA}: you cannot write to that directory";
   my $parasite_data_dir = join("/", $ENV{PARASITE_DATA}, $species);
   File::Path::make_path($parasite_data_dir);

   my ($assembly_path) = glob("$parasite_data_dir/*genomic.fna");
   
   if( !$force and $assembly_path and -s $assembly_path ) {
      warn "$assembly_path already exists: no new download\n";
   } else {
      my $ftp = Net::FTP->new(NCBI_FTP_SERVER)
         or croak "$@";
       
      my $ftp_uri = $conf->{FtpPath_GenBank};
      $ftp_uri =~ m~^ftp://${\NCBI_FTP_SERVER}~ or croak "URI doesn't match the expected host name ${\NCBI_FTP_SERVER}: $ftp_uri";
      try{
         # oh my days
         # (my $uri = $conf->{FtpPath_GenBank}) =~ s{ftp://ftp.ncbi.nlm.nih.gov}{};
         my $ftp_path   = substr($ftp_uri,length('ftp://'.NCBI_FTP_SERVER));

         $ftp->login("anonymous",'-anonymous@') or die "Cannot log in: ".$ftp->message;
         $ftp->cwd($ftp_path) or die "Cannot access directory $ftp_path: ".$ftp->message;
         my @files = $ftp->ls() or die "Cannot get directory listing: ".$ftp->message;
         my ($remote, @others) = grep {$_ !~ /from_genomic\.fna\.gz$/} grep {$_ =~ /genomic\.fna\.gz$/ } @files;
         die "Multiple matches for the expected assembly file name pattern (should be unique): ".  join("\t", $remote, @others) if $others[0];
         die "No match for the expected assembly file name pattern in directory $ftp_path" unless $remote;
         $assembly_path = "$parasite_data_dir/$remote";
         $ftp->binary() or die "Cannot select binary file transfer: ".$ftp->message;
         $ftp->get($remote, $assembly_path) or die "Cannot download file $remote: ".$ftp->message;
         if($gff3_download) {
            # note that GFF3 files at NCBI are named .gff, but are stored locally as .gff3
            my $gff3_remote   = $remote;
            my $gff3_path     = $assembly_path;
            $gff3_remote   =~ s/\.fna\.gz\s*$/\.gff\.gz/;
            $gff3_path     =~ s/\.fna\.gz\s*$/\.gff3\.gz/;
            if( $ftp->get($gff3_remote, $gff3_path) ) {
               my $gunzipped = $gff3_path;
               $gunzipped =~ s/\.gz\s*$//;
               warn "NCBI annotation file $gff3_remote will be stored locally as ".basename($gunzipped)."\n";
               IO::Uncompress::Gunzip::gunzip($gff3_path => $gunzipped) or croak "failed to gunzip $gff3_path: $IO::Uncompress::Gunzip::GunzipError";
               unlink $gff3_path || warn "couldn't remove gzipped source file $gff3_path";
            } else {
               warn "Cannot download file $gff3_remote: ".$ftp->message;
            }
         }
      } catch {
         croak "Error retrieving assembly from $ftp_uri: $_";
      };
      #system("curl --silent ftp://ftp.ncbi.nlm.nih.gov/$ftp_path/$remote > $assembly_path") and die "$@";
      # don't use subshell; IO::Uncompress::Gunzip provide interface
      #system("gunzip $assembly_path") and die "Failed: gunzip $assembly_path";  
      my $gunzipped = $assembly_path;
      $gunzipped =~ s/\.gz$//;
      IO::Uncompress::Gunzip::gunzip($assembly_path => $gunzipped) or croak "failed to gunzip $assembly_path: $IO::Uncompress::Gunzip::GunzipError";
      unlink $assembly_path || warn "couldn't remove gzipped source file $assembly_path";
      $assembly_path = $gunzipped;
   }
   croak "Failed to write local copy of assembly to $assembly_path" unless -s $assembly_path;

   my $assembly_destination_path = join("/", $parasite_data_dir, "$species.fa");
   my $seq_region_synonyms_path  = join("/", $parasite_data_dir, "$species.seq_region_synonyms.tsv");
   if( !$force and -s $assembly_destination_path and -s $seq_region_synonyms_path ) {
      warn "$assembly_destination_path and $seq_region_synonyms_path already exist: will not be recreated\n";
   } else {
   # if( not -s $assembly_destination_path or not -s $seq_region_synonyms_path){
      open(my $assembly_in_fh, "<", $assembly_path) or croak "cannot read $assembly_path: $!";
      open(my $assembly_out_fh, ">", $assembly_destination_path) or croak "cannot write to $assembly_destination_path: $!";
      open(my $seq_region_synonyms_fh, ">", $seq_region_synonyms_path) or croak "cannot write to $seq_region_synonyms_path: $!";
      ASSEMBLY_IN: while(my $line = <$assembly_in_fh>){
         unless($line =~ m/^>/){
            print $assembly_out_fh $line;
            next ASSEMBLY_IN;
         }
         chomp $line;
         # this old code worked, but I don't really like the way $line is incrementally modified
         # as we then lose the original FASTA header as a reference
         # $line =~ s/^>//;
         # my ($scaffold_name_ncbi) = $line =~ /(.*?) /;
         # $line =~ s/^$scaffold_name_ncbi//;
         # $line =~ s/, whole genome shotgun sequence$//;
         # my ($scaffold_name_local) = $line =~ m/.*? (\S+)$/;
         my ($scaffold_name_ncbi) = $line =~ /^>(.*?)\s/;
         my $scaffold_name_local = $line;
         $scaffold_name_local =~ s/^$scaffold_name_ncbi//;
         $scaffold_name_local =~ s/, whole genome shotgun sequence$//;
         $scaffold_name_local =~ s/.*\s//;

         # $scaffold_name_local is just a best guess, and may fail to find a scaffold name/ID
         # genuine names typically have  both letters and digits so that may be a useful sanity check?
         warn "*WARNING* Guessed at local scaffold name \"$scaffold_name_local\" but that doesn't look too good.\n"
            . "(FASTA header was \"$line\")\n"
            . "If it is wrong, you will need to manually change it in $species.fa  $species.seq_region_synonyms.tsv\n"
            unless $scaffold_name_local =~ m/[a-zA-Z]/ && $scaffold_name_local =~ m/[0-9]/;

         # behaviour of this original line wasn't very clear $isolate could be undefined
         # and the die output didn't state what the problem was
         #die "$scaffold_name_local | $_" if grep {$scaffold_name_local =~ /$_/i} ("$spe.$cies", $isolate);
         croak "Found species name $spe $cies in local scaffold name $scaffold_name_local: baled" if $scaffold_name_local =~ /$spe.$cies/i;
         croak "Found isolate name $isolate in local scaffold name $scaffold_name_local: baled" if $isolate && $scaffold_name_local =~ /$isolate/i;
         print $assembly_out_fh ">$scaffold_name_local\n";
         print $seq_region_synonyms_fh join("\t", "toplevel", $scaffold_name_local, $scaffold_name_ncbi, "INSDC")."\n"; 
      }
      # close $_ for ($assembly_in_fh, $assembly_out_fh, $seq_region_synonyms_fh);
      foreach my $fh ($assembly_in_fh, $assembly_out_fh, $seq_region_synonyms_fh) {
         close($fh);
      }
   }
}

# PRINT:

my %urls = (
  "University of Hull"        => "https://www.hull.ac.uk",
  "WTSI"                      => "http://www.sanger.ac.uk",
  "University of Melbourne"   => "http://www.gasserlab.org",
  "Oregon State University"   => "https://oregonstate.edu",
  'University College London' => 'https://www.ucl.ac.uk',
);
my %synonyms = (
  "WTSI"                      => "Wellcome Sanger Institute",
  "ED"                        => "University of Edinburgh",
  'UCL'                       => 'University College London',
  'UNIVERSITY COLLEGE LONDON' => 'University College London',
);
my $this_url = exists $urls{$conf->{SubmitterOrganization}}
             ? $urls{$conf->{SubmitterOrganization}}
             : (  exists $synonyms{$conf->{SubmitterOrganization}} && exists $urls{$synonyms{$conf->{SubmitterOrganization}}}
                  ?  $urls{$synonyms{$conf->{SubmitterOrganization}}}
                  :  CoreCreation::Config::Utils::MISSING_METADATA_PATTERN
                  );

my $output_conf = {  taxon_id                => $conf->{SpeciesTaxid} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                     assembly_version        => $conf->{AssemblyName} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                     gff_sources             => "WormBase_imported",
		     gff_codinggene_logic    => CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
		     gff_noncodinggene_logic => CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
	             PMID                    => CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                     PM_text		     => CoreCreation::Config::Utils::MISSING_METADATA_PATTERN, 
                     meta              => {  "assembly.accession"       => $conf->{AssemblyAccession} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             "provider.name"            => $synonyms{$conf->{SubmitterOrganization}//""} // $conf->{SubmitterOrganization} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             "provider.url"             => $this_url // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             "species.strain"           => $conf->{Biosource}{Isolate} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             "species.biosample"        => $conf->{BioSampleAccn} // CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             "species.nematode_clade"   => $species =~ /meloidogyne|globodera|heterodera/ ? "IV" : $species =~ /pristionchus|caenorhabditis/ ? "V": CoreCreation::Config::Utils::MISSING_METADATA_PATTERN,
                                             }
                     };
# if the source config indicated there mitochondrial DNA in this assembly, flag it by inserting a null value
# (this /could/ get populated with a real value in the following steps, but at least it tells the user that
# there's a mtDNA sequence in here somewhere that needs to be idnetified)
# required end value for this is a comma-separate list of all sequences that are mtDNA
if( grep( $_ eq 'has-mitochondrion', @{$conf->{PropertyList}->{string} || []}) ) {
   $output_conf->{mitochondrial} = CoreCreation::Config::Utils::MISSING_METADATA_PATTERN;
   warn  "*WARNING* there is at least one mitochondrial DNA sequence in this assembly.  These may be\n"
      .  "automatically identified, but if not then the 'mitochondrial' config value must manually be\n"
      .  "set, with a comma-separated list of the mtDNA sequence identifiers\n";
}
                  
print Dump( { $species => $output_conf } );

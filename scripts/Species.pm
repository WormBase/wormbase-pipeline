#
# Species-specific properties 
#
# look at Wormbase::initialize for usage
#
# The TSL sequences are taken from
# PLOS Genetics
# Nov 2006
# Operon Conservation and the Evolution of trans-Splicing in the Phylum Nematoda
# David B. Guiliano, Mark L. Blaxter
# http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.0020198
# Figure 4
#
# All nematodes have the SL1 TSL, but vary in additional sequences.

use strict;
use Carp;

package Species;

our @ISA = qw(Wormbase);

sub flatten_params {
  shift;    # get rid of the class name
  my %param_hash = %{ shift() };
  my @param_array;
  while ( my ( $k, $v ) = each %param_hash ) { push @param_array, $k, $v }
  return @param_array;
}

sub _new {
  my $class = shift;
  my %param = %{ shift(@_) };
  
  my $self = $class->initialize( $class->flatten_params( \%param ) );
  
  bless $self, $class;
}

# use with  (-prefix => 'whatever', -mito => 1)
sub get_chromosome_names {
  my $self= shift;
  my %options = @_;
  my @chromosomes = $self->chromosome_names;
  my $prefix=$options{'-prefix'};
  $prefix = $self->chromosome_prefix if ($prefix && length $prefix < 2);
  push @chromosomes, $self->mt_name if $options{'-mito'} && $self->mt_name;
  return map {$prefix.$_} @chromosomes if $options{'-prefix'};
  return @chromosomes
}

# use with  (-prefix => 'whatever', -mito => 1)
sub get_chromosome_lengths {
  my $self= shift;
  my %options = @_;
  my $chr_lengths = $self->chromosome_lengths;
  my $prefix = $options{'-prefix'};
  
  $prefix = $self->chromosome_prefix if ($prefix && length $prefix < 2);
  if ($prefix) {
    foreach my $k (keys %$chr_lengths) {
      my $pk = $prefix . $k;
      my $val = $chr_lengths->{$k};
      delete $chr_lengths->{$k};
      $chr_lengths->{$pk} = $val;
    }
  }

  return $chr_lengths;
}

sub _seq_names_and_length_cache {
  my ($self) = @_;
  return $self->common_data . "/toplevel_seqs.lst";
}

sub _cache_chromosome_names_and_lengths {
  my ($self) = @_;

  my $cache_file = $self->_seq_names_and_length_cache;

  my $db = Ace->connect(-path => $self->autoace) || die(Ace->error);

  my $species = $db->fetch(Species => $self->full_name);
  die("cannot find species ${\$self->full_name} in ${\$self->autoace}\n") unless $species;
  
  my $assembly;
  
  my @assemblies = $species->Assembly;
  while (my $a = shift @assemblies){
    next unless $a->Status eq 'Live';
    my @hit = grep {$self->ncbi_bioproject eq $_} $a->DB_info->col(3);
    $assembly = $a if $hit[0];
  }
  my @sequences = $assembly->follow(-tag=>'Sequences');
  my %h;
  require Ace::Sequence;
  foreach my $seq (@sequences) {
    my $ace_seq = Ace::Sequence->new($seq);
    $h{$seq->name} = $ace_seq->length;
  }

  if (keys %h) {
    open (my $outf,">$cache_file") || die($!);
    foreach my $k (keys %h) {
      print $outf "$k $h{$k}\n";
    }
    close $outf;
  } else {
    die "Could not find names/lengths for any toplevel sequences\n";
  }
}


sub chromosome_names {
  my ($self) = @_;
  
  my $cache =  $self->_seq_names_and_length_cache; 
  if (not -e $cache) {
    $self->_cache_chromosome_names_and_lengths;
  }

  my @names;
  open(my $genomefh, $cache) or die "Could not open $cache for reading\n";
  while(<$genomefh>){ 
    /^(\S+)/ and push @names, $1;
  }
  close($genomefh);

  return @names;
}

sub chromosome_lengths {
  my ($self) = @_;
  
  my $cache =  $self->_seq_names_and_length_cache; 
  if (not -e $cache) {
    $self->_cache_chromosome_names_and_lengths;
  }

  my %lens;
  open(my $genomefh, $cache) or die "Could not open $cache for reading\n";
  while(<$genomefh>){ 
    /^(\S+)\s+(\d+)/ and $lens{$1} = $2;
  }
  close($genomefh);

  return \%lens;
}

# denotes whether this is the canonical project for this species
sub is_canonical {
  return 1;
}

sub TSL {  
  ('SL1'  => "GGTTTAATTACCCAAGTTTGAG"); # the SL1 sequence is conserved across the majority of nematodes
}

sub full_name {
  my $self = shift;
  my %param = @_ ;
  if($param{'-short'}){
    return $self->short_name;
  } elsif($param{'-g_species'}){
    return $self->gspecies_name;
  }
  else { 
    return $self->long_name;
  }
}

sub chromosome_prefix { return '' }
sub mt_name {return undef}
sub pep_prefix {return undef}
sub wormpep_prefix {return undef}
sub pepdir_prefix{ return undef };
sub ncbi_bioproject{ return undef};
sub repeatmasker_library {return undef}

#########################################
#
# Core species
#
#########################################

package Elegans;

our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ELE/elegans.lib'}
sub chromosome_prefix {'CHROMOSOME_'}
sub chromosome_names {qw(I II III IV V X)}
sub mt_name {'MtDNA'}
sub pep_prefix {'CE'}
sub pepdir_prefix{'worm'};
sub cds_regex{qr/^[A-Z0-9_cel]+\.[1-9]\d{0,3}[A-Za-z]?$/};  #the cel is for telomeric clone CDSs cTel54X.1
sub seq_name_regex{qr/^[A-Z0-9_cel]+\.[1-9]\d{0,3}/};  #to get just the Sequence_name part
sub cds_regex_noend{qr/^[A-Z0-9_cel]+\.[1-9]\d{0,3}[A-Za-z]?/};  # for getting the CDS part of a Transcript name

sub ncbi_tax_id {'6239'};
sub ncbi_bioproject {'PRJNA13758'};
sub bioproject_description {'C.elegans Sequencing Consortium genome project'};
sub short_name {'C. elegans'}
sub gspecies_name {'c_elegans'}
sub long_name {'Caenorhabditis elegans'}
sub wormpep_prefix{'WP'}
sub assembly_type {'chromosome'};
sub seq_db {my $self = shift;return $self->database('camace');}
sub upload_db_name {return ('citace')};
sub TSL {
  (
    'SL1' => 'GGTTTAATTACCCAAGTTTGAG',
    'SL2' => 'GGTTTTAACCCAGTTACTCAAG',
    'SL2a' => 'GGTTTATACCCAGTTAACCAAG',
    'SL2b' => 'GGTTTTAACCCAGTTTAACCAAG',
    'SL2c' => 'GGTTTTAACCCAGTTACCAAG',
    'SL2d' => 'GGTTTTTACCCAGTTAACCAAG',
#    'SL2e' => 'GGTTTAAAACCCAGTTAACAAG', This is erroneous but has been around since the late 90s.
    'SL2f' => 'GGTTTTAACCCAGTTAACCAAG',
    'SL2g' => 'GGTTTTAACCAGTTAACTAAG',
    'SL2h' => 'GGTTTTAACCCATATAACCAAG',
    'SL2i' => 'GGTTTTAACCCAAGTTAACCAAG',
    'SL2j' => 'GGTTTAAAACCCAGTTACCAAG',
    'SL2k' => 'GGTTTTAACCCAGTTAATTGAG',
  );

};

########################################
package Briggsae;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/BRI/briggsae_ws200_filtered.fas'}
sub pep_prefix {'CBP'}
sub pepdir_prefix{'brig'};
sub cds_regex{qr/^CBG\d{5}[a-z]*$/};
sub seq_name_regex{qr/^CBG\d{5}/};
sub cds_regex_noend{qr/^CBG\d{5}[a-z]*/}; # for getting the CDS part of a Transcript name
sub ncbi_tax_id {'6238'};
sub ncbi_bioproject {'PRJNA10731'};
sub bioproject_description {'C. briggsae Sequencing Consortium genome project'};
sub assembly_type {'contig'};
sub seq_db {my $self = shift;return $self->database('briggsae');}
sub short_name {'C. briggsae'}
sub gspecies_name {'c_briggsae'}
sub long_name {'Caenorhabditis briggsae'}
sub wormpep_prefix{'BP'}
sub upload_db_name {'briggsae'}
sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'Cb_SL2' => "GGTTTTAACCCAGTTACTCAAG",
	  'Cb_SL3' => "GGTTTTAACCCAGTTAACCAAG",
	  'Cb_SL4' => "GGTTTTAACCCAGTTTAACCAAG",
	  'Cb_SL10' => "GGTTTTAACCCAAGTTAACCAAG",
	  'Cb_SL13' => "GGATTTATCCCAGATAACCAAG",
	  'Cb_SL14' => "GGTTTTTACCCTGATAACCAAG",
)};


#########################################
package Remanei;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/REM/remanei_ws200_filtered.fas'}
sub short_name {'C. remanei'}
sub gspecies_name{'c_remanei'}
sub long_name{'Caenorhabditis remanei'}
sub wormpep_prefix {'RP'}
sub pep_prefix {'RP'}
sub pepdir_prefix{'rema'};
sub ncbi_tax_id {'31234'};
sub ncbi_bioproject {'PRJNA53967'};
sub bioproject_description {'C.remanei Sequencing Consortium genome project'};
sub cds_regex{qr/CRE\d{5}[a-z]*/};
sub seq_name_regex{qr/^CRE\d{5}/};
sub cds_regex_noend{qr/CRE\d{5}[a-z]*/}; # for getting the CDS part of a Transcript name
sub assembly_type {'contig'};
sub seq_db {my $self = shift;return $self->database('remanei');}
sub upload_db_name {'remanei'};

# The genes CRE33181 and CRE33245 defined the additional sequences
# GGTTTTAACCGATTTAACGAAG and GGTTTATGCCCAGTTAGCCAAG, but these are not
# seen in the short read data, so these genes are likely to be
# pseudogenes.
sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'CRE_SL2' => "GGTTTTAACCCAGTTACTCAAG",
	  'CRE_SL2a' => "GGATTTATCCCAGTTAACCAAG",
	  'CRE_SL2b' => "GGTTTTAACCCAGTTTAACCAAG",
	  'CRE_SL2c' => "GGTTTTAACCCAGTTAACAAAG",
	  'CRE_SL2d' => "GGTTTTTACCCAGTTAACCAAG",
	  'CRE_SL2e' => "GGTTTTAACCCAGTTAACTAAG",
	  'CRE_SL2f' => "GGTTTTAACCCAGTTAACCAAG",
	  'CRE_SL2g' => "GGTTTTAACCCCAGTAACCAAG",
	  'CRE_SL2h' => "GGTTTTTAACCCAGTTACTCAAG",
	  'CRE_SL2i' => "GGTTTTAACCCAAGTTAACCAAG",
	)};


#######################################################
package Brenneri;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/BRE/brenneri_ws200_filtered.fas'}
sub pep_prefix {'CN'}
sub pepdir_prefix{'bre'};
sub ncbi_tax_id {'135651'};
sub ncbi_bioproject {'PRJNA20035'};
sub bioproject_description {'C.brenneri Sequencing Consortium genome project'};
sub short_name {'C. brenneri'}
sub gspecies_name{'c_brenneri'}
sub long_name{'Caenorhabditis brenneri'}
sub assembly_type {'contig'};
sub seq_db {my $self = shift;return $self->database('brenneri');}
sub wormpep_prefix {'CN'}
sub cds_regex{qr/CBN\d{5}[a-z]*/};
sub seq_name_regex{qr/^CBN\d{5}/};
sub cds_regex_noend{qr/CBN\d{5}[a-z]*/}; # for getting the CDS part of a Transcript name
sub upload_db_name {'brenneri'};

# The following genes define some additional SL2 sequences but these
# are seen seen in the short read data at low frequency, so these
# genes are likely to be pseudogenes.
# ggttttaacccagtgtaaccaag	CBN33679
# ggttttaaccctattaccaag	CBN33735
# agttttaacccagtaaccaag	CBN33753
# ggttttaaccctgataaccaag	CBN33765
# ggttttaaccctgataaccaag	CBN33767
# ggtttcaacccattataaccaag	CBN33769
# ggttttaacccttttaccaag	CBN33771
# ggatttatcccagttacccaag	CBN33774
# ggttttaaccctgataaccaag	CBN33798
# ggtttttacccagtataaaccag	CBN33807
# ggttttaaccctttaaccaag	CBN33816
# ggtttttacccagtttaaccaaag	CBN33818
# ggttttaaccctattaccaag		CBN33851
# ggttttaaccctatatccaag		CBN33859
# gatatcattgagatttaaccaag		CBN33918
# ggttttggcccagtaaccaag		CBN33919
# ggttcttgctactattgcgccccag	CBN33921
# agtacaaccccaacgggtgaaacaaaag	CBN33924
# ggccgtgaccctttaatcaaggttag	CBN33946
# gatttcccgtataaactaaattttggaaccg	CBN34178
# ggctttttaccctgtatagccttag	CBN34182

sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'CBN_SL2' => "GGTTTTAACCCAGTTACTCAAG",
	  'CBN_SL2a' => "GGATTTATCCCAGATAACCAAG",
	  'CBN_SL2b' => "GGTTTTAACCCAGTTTAACCAAG",
	  'CBN_SL2c' => "GGTTTTAACCCAGTTACCAAG",
	  'CBN_SL2d' => "GGTTTTTACCCAGTTAACCAAG",
	  'CBN_SL2e' => "GGATTTATCCCAGTTACTCAAG",
	  'CBN_SL2f' => "GGTATTAACCCTGATAACCAAG",
	  'CBN_SL2g' => "GGTTTCAACCCTGATAACCAAG",
	  'CBN_SL2h' => "GGTTTTAACCCAGATAACCAAG",
	  'CBN_SL2i' => "GGTTTTAACCCAAGTTAACCAAG",
	  'CBN_SL2j' => "GGTTTTAACCCAGTAACCAAG",
	  'CBN_SL2k' => "GGTTTTAACCCAGTATAACCAAG",
	  'CBN_SL2l' => "GGTTTTAACCCAGTTACCCAAG",
	  'CBN_SL2m' => "GGTTTTAACCCCAGTTACCAAG",
	  'CBN_SL2n' => "GGTTTTAACCCTATAACAAAG",
	  'CBN_SL2o' => "GGTTTTAACCCTATAACCAAG",
	  'CBN_SL2p' => "GGTTTTAACCCTATAACCTATAACCAAG",
	  'CBN_SL2q' => "GGTTTTTACCCAGTATAACCAAG",
	  'CBN_SL2r' => "GGTTTTTACCCAGTTTAACAAAG",
	  'CBN_SL2s' => "GGTTTTTACCCAGTTTAACCAAG",
	  'CBN_SL2t' => "GGTTTTTACCCCAGTTAACCAAG",
	  'CBN_SL2u' => "GGTTTTTACCCTAATTACCAAG",
	  'CBN_SL2v' => "GGTTTTTACCCTAGTTAACCAAG",
	  'CBN_SL2w' => "GGTTTTTACCCTGTTAACCAAG",
	)};


#######################################################

package Japonica;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/JAP/japonica_ws200_filtered.fas'}
sub pep_prefix {'JA'}
sub pepdir_prefix{'jap'};
sub ncbi_tax_id {'281687'};
sub ncbi_bioproject {'PRJNA12591'};
sub bioproject_description {'C.japonica Sequencing Consortium genome project'};
sub cds_regex{qr/CJA\d{5}[a-z]*/};
sub seq_name_regex{qr/^CJA\d{5}/};
sub cds_regex_noend{qr/CJA\d{5}[a-z]*/}; # for getting the CDS part of a Transcript name
sub assembly_type {'contig'};
sub seq_db {my $self = shift;return $self->database('japonica');}
sub wormpep_prefix {'JA'}
sub short_name {'C. japonica'}
sub gspecies_name{'c_japonica'}
sub long_name{'Caenorhabditis japonica'}
sub upload_db_name {'japonica'};

# The gene CJA48849 defines the additional sequence
# GGTTTCAACCTAGTTAATCAAG, but this is not seen in the short read data,
# so this gene is likely to be a pseudogene.
sub TSL {
  (
    'SL1' => 'GGTTTAATTACCCAAGTTTGAG',
    'CJA_SL2' => 'GGTTTTAACCCAGTTACTCAAG',
    'CJA_SL2a' => 'GGTTTTAACCCAGAAACTAAAG',
    'CJA_SL2b' => 'GGTTTTAACCCAGTAACTAAAG',
    'CJA_SL2c' => 'GGTTTTAACCCAGTTACCCAAG',
    'CJA_SL2d' => 'GGTTTTAACCCCAGTTAACCAAG',
    'CJA_SL2e' => 'GGTTTTAACCCCAGTTAATCAAG',
    'CJA_SL2f' => 'GGTTTTAACCCAGTTAACCAAG',
    'CJA_SL2g' => 'GTTAATACCCCAGTTATCAAG',
    'CJA_SL2h' => 'GTTTTAAACCCAGTTAACCAAG',
    'CJA_SL2i' => 'GGTTTTAACCCAAGTTAACCAAG',
  )};


#######################################################

package Pristionchus;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/PPA/ppa_repeats.fas'}
sub short_name {'P. pacificus'}
sub gspecies_name{'p_pacificus'}
sub long_name{'Pristionchus pacificus'}
sub cds_regex{qr/PPA\d{5}[a-z]*/};
sub seq_name_regex{qr/^PPA\d{5}/};
sub cds_regex_noend{qr/PPA\d{5}[a-z]*/}; # for getting the CDS part of a Transcript name
sub pep_prefix {'PP'}
sub pepdir_prefix{'ppa'};
sub ncbi_tax_id {'54126'};
sub ncbi_bioproject {'PRJNA12644'};
sub bioproject_description { 'Max Planck Institute for Developmental Biology P. pacificus genome project' }
sub assembly_type {'contig'};
sub wormpep_prefix {'PP'}
sub upload_db_name {''}; # we already hold the data in the primary database, it is not downloaded
sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'Pp_SL2a' => "GGTTTTTACCCAGTATCTCAAG",
	  'Pp_SL2b' => "GGTTTTAACCCAGTATCTCAAG",
	  'Pp_SL2c' => "GGTTTATACCCAGTATCTCAAG",
	  'Pp_SL2d' => "GGTTTTTAACCCAGTATCTCAAG",
	  'Pp_SL2e' => "GGTTTTTACTCAGTATCTCAAG",
	  'Pp_SL2f' => "GGTCTTTACCCAGTATCTCAAG",
	  'Pp_SL2g' => "GGTTTTAACCCGGTATCTCAAG",
	  'Pp_SL2h' => "GGTTTTAACCCAGTATCTTAAG",
	  'Pp_SL2i' => "GGTTTTGACCCAGTATCTCAAG",
	  'Pp_SL2j' => "GTTTTATACCCAGTATCTCAAG",
	  'Pp_SL2k' => "GGTTTATACCCAGTATCTCAAG",
	  'Pp_SL2l' => "GGTTTAAACCCAGTATCTCAAG",
)};


####################################################

package Brugia;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/BRU/brugia_ws236_filtered.fas'}
sub short_name {'B. malayi'}
sub gspecies_name{'b_malayi'}
sub long_name{'Brugia malayi'}
sub seq_name_regex{qr/^Bm\d+/}
sub pep_prefix {'BM'}
sub pepdir_prefix{'brug'}
sub cds_regex_noend{qr/Bm\d+[a-z]*/}
sub cds_regex{qr/Bm\d+[a-z]*/}
sub ncbi_tax_id {'6279'}
sub ncbi_bioproject {'PRJNA10729'}
sub bioproject_description { 'University of Pittsburgh B. malayi genome project' }
sub assembly_type {'contig'}
sub seq_db {my $self = shift;return $self->database('brugia');}

sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'Bm_SL1a' => "GGTTTTAATTACCCAAGTTTGAG",
	  'Bm_SL1b' => "GGTTTAATCACCCAAGTTTGAG",
	  'Bm_SL1c' => "GGTTTAACTACCCAAGTTTGAG",
)};

sub wormpep_prefix {'BM'}
sub upload_db_name {'brugia'};

####################################################

package Ovolvulus;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ONV/onchocerca_volvulus_repeatModeller.fa.filtered'}
sub short_name {'O. volvulus'}
sub gspecies_name{'o_volvulus'}
sub long_name{'Onchocerca volvulus'}
sub seq_name_regex{qr/^OVOC\d+/}
sub pep_prefix {'OVP'}
sub pepdir_prefix{'ovol'};
sub cds_regex_noend{qr/OVOC\d+[a-z]*/}
sub cds_regex{qr/OVOC\d+[a-z]*/}
sub ncbi_tax_id {6282};
sub ncbi_bioproject {'PRJEB513'}
sub bioproject_description { 'Wellcome Trust Sanger Institute O. volvulus genome project' }
sub assembly_type {'contig'};
sub seq_db {my $self = shift;return $self->database('ovolvulus');}

sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
)};

sub wormpep_prefix {'OV'}
sub upload_db_name {'ovolvulus'};

######################################################
package Tmuris;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/TMU/TMUE.fa'}
sub short_name {'T. muris'}
sub gspecies_name{'t_muris'}
sub long_name{'Trichuris muris'}
sub ncbi_tax_id {'70415'}
sub ncbi_bioproject {'PRJEB126'}
sub bioproject_description { 'Wellcome Trust Sanger Institute T. muris genome project' }
sub assembly_type {'contig'}

sub cds_regex{qr/^TMUE_(0|1|2|3|M)\d{9}[a-z]*/};
sub seq_name_regex{qr/^TMUE_(0|1|2|3|M)\d{9}/};
sub cds_regex_noend{qr/^TMUE_(0|1|2|3|M)\d{9}[a-z]*/}; # for getting the CDS part of a Transcript name

sub pep_prefix {'TMP'}
sub pepdir_prefix{'tmu'};
sub seq_db {my $self = shift;return $self->database('tmuris');}
sub wormpep_prefix {'TMP'}

# T. muris TSL sequences, taken from:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4125394/
my %TSL = (
           'SL1'  => 'GGTTATTTACCCTGTTAACAAG',
           'SL2A' => 'GGTTAATTACCCAATTTAAAAG',
           'SL2B' => 'GGTAAATTTACCCACAGTGAAG',
           'SL2C' => 'GGTTAAGTTTACCCAATTGAAG',
           'SL2D' => 'GGTTATTTACCCATTGCCACAG',
           'SL2E' => 'GGTTAATTACCCAATTTCAAAG',
           'SL2F' => 'GGTTTTATACCCTCTCACCAAG',
           'SL2G' => 'GGTAAATTTACCCGCAATAAAG',
           'SL2H' => 'GGTACATTTACCCACAGTTAAG',
           'SL2I' => 'GGTTTATTACCCAAATCGAAAG',
           'SL2J' => 'GGTTATTTATCCCCTGACCAAG',
           'SL2K' => 'GGTTAAATTTACCCCTCAAAAG',
           'SL2L' => 'GGTATTTACCCAACGTTGACTG',
);


######################################################
package Sratti;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/SRA/SRAE.repeatLib.fa'}
sub short_name {'S. ratti'}
sub gspecies_name{'s_ratti'}
sub long_name{'Strongyloides ratti'}
sub ncbi_tax_id {'34506'}
sub ncbi_bioproject {'PRJEB125'}
sub bioproject_description { 'Wellcome Trust Sanger Institute S. ratti genome project' }
sub assembly_type {'contig'}

sub seq_name_regex{qr/^SRAE_[\dXM]\d+/};
sub pep_prefix {'SRP'}
sub pepdir_prefix{'sra'};
sub cds_regex_noend{qr/SRAE_[\dXM]\d+[a-z]*/}; # for getting the CDS part of a Transcript name
sub cds_regex{qr/SRAE_[\dXM]\d+[a-z]*/};
sub seq_db {my $self = shift;return $self->database('sratti');}

sub wormpep_prefix {'SRP'}

sub TSL {(
	  'SL1'  => "GGTTTAATTACCCAAGTTTGAG",
	  'Sr_SL1a' => "GGTTTATAAAACCCAGTTTGAG",
	  'Sr_SL1b' => "GGTTTAAAAAACCCAGTTTGAG",
	  'Sr_SL1c' => "GGTTTAAAAACCCAGTTTGAG",
	  'Sr_SL1d' => "GGTTTTAAAACCCAGTTTGAG",
	  'Sr_SL1e' => "GGTTTAAAAACCCAATTTGAG",
	  'Sr_SL1f' => "GGTTTAAATAACCCAGTTTGAG",
	  'Sr_SL1g' => "GGTTTAAATAACCCATATAGAG",
	  'Sr_SL1h' => "GTTTTTTAAATAACCAAGTTTGAG",
	  'Sr_SL1i' => "GGTTTAAGAAAACCCATTCAAG",
	  'Sr_SL1j' => "GGTTTTATAAAACCCAGTTTGAG",
	  'Sr_SL1k' => "GGTTTATAAAACCCAGTTTAAG",
	  'Sr_SL1l' => "GGTTTAAAAACCCGATTTTGAG",
	  'Sr_SL1m' => "GGTTTTAAATAACCCAGTTTGAG",
	  'Sr_SL1n' => "GGTTTATATAACCCAGTTTGAG",
	  'Sr_SL1o' => "GGTTTAAAAACCCAAATTAAA",
	  'Sr_SL1p' => "GGTTTTAAAAACCCAGTTTGAG",
	  'Sr_SL1q' => "GGTTTATACAACCCAGTTTGAG",
	  'Sr_SL1r' => "GGTTTAAGAAACCCTGTTTGAG",
	  'Sr_SL1s' => "GGTTTAAAAAACCCAGTTTAAG",
)};

######################################################
#
# Non-Core species
#
######################################################

package Cangaria;
use Carp;

our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/CAN/caenorhabditis_angaria_RepeatModeller.fa.classified'}
sub short_name {'C. angaria'}
sub gspecies_name{'c_angaria'}
sub long_name{'Caenorhabditis angaria'}
sub ncbi_tax_id {'860376'};
sub ncbi_bioproject {'PRJNA51225'};
sub bioproject_description { 'California Institute of Technology C. angaria genome project' }
sub assembly_type {'contig'};

#######################################################

package Cnigoni;
use Carp;

our @ISA = qw(Species);

sub short_name {'C. nigoni'}
sub gspecies_name{'c_nigoni'}
sub long_name{'Caenorhabditis nigoni'}
sub ncbi_tax_id {1611254};
sub ncbi_bioproject {'PRJNA384657'};
sub bioproject_description { 'Caenorhabditis nigoni strain JU1422, whole genome shotgun sequencing project' }
sub assembly_type {'contig'};
sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/CNI/nigoni_repeats_2016.08.14.fa'}

#######################################################

package Csinica;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/C5/caenorhabditis_species5.fa.classified'}
sub short_name {'C. sinica'}
sub gspecies_name{'c_sinica'}
sub long_name{'Caenorhabditis sinica'}
sub ncbi_tax_id {'1550068'}; # this is the ID of the DRD-2008 strain
sub ncbi_bioproject {'PRJNA194557'};
sub bioproject_description { 'University of Edinburgh Caenorhabditis sinica genome project' }
sub assembly_type {'contig'};

######################################################
package Cafra;
use Carp;
our @ISA = qw(Species);

sub short_name {'C. afra'}
sub gspecies_name{'c_afra'}
sub long_name{'Caenorhabditis afra'}
sub ncbi_tax_id {'1094335'};
sub ncbi_bioproject {'PRJNA51171'};
sub bioproject_description { 'Genome Institute at Washington University Caenorhabditis afra genome project' }
sub assembly_type {'contig'};

######################################################
package Ctropicalis;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/C11/caenorhabditis_species11_repeatModeller.fa.classified'}
sub short_name {'C. tropicalis'}
sub gspecies_name{'c_tropicalis'}
sub long_name{'Caenorhabditis tropicalis'}
sub ncbi_tax_id {'1561998'};
sub ncbi_bioproject {'PRJNA53597'};
sub bioproject_description { 'Genome Institute at Washington University Caenorhabditis sp. 11 genome project' }
sub assembly_type {'contig'};

######################################################
package Panagrellus;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/PRED/panagrellus_redivivus_RepeatModeller.fa.classified_marissa'}
sub short_name {'P. redivivus'}
sub gspecies_name{'p_redivivus'}
sub long_name{'Panagrellus redivivus'}
sub ncbi_tax_id {'6233'};
sub ncbi_bioproject {'PRJNA186477'};
sub bioproject_description { 'California Institute of Technology P. redivivus genome project'};
sub assembly_type {'contig'};

#######################################################

package Remanei_px356;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/REM/remanei_ws200_filtered.fas'}
sub short_name {'C. remanei'}
sub gspecies_name{'c_remanei'}
sub long_name{'Caenorhabditis remanei'}
sub ncbi_tax_id {'31234'};
sub ncbi_bioproject {'PRJNA248909'};
sub bioproject_description { 'University of Oregon University C. remanei PX356 genome project' }
sub assembly_type {'contig'};
sub is_canonical { 0 };

#######################################################
package Remanei_px439;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/REM/remanei_ws200_filtered.fas'}
sub short_name {'C. remanei'}
sub gspecies_name{'c_remanei'}
sub long_name{'Caenorhabditis remanei'}
sub ncbi_tax_id {'31234'};
sub ncbi_bioproject {'PRJNA248911'};
sub bioproject_description { 'University of Oregon University C. remanei PX356 genome project' }
sub assembly_type {'contig'};
sub is_canonical { 0 };

#######################################################

package Elegans_hawaii;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ELE/elegans.lib'}
sub short_name {'C. elegans'}
sub gspecies_name{'c_elegans'}
sub long_name{'Caenorhabditis elegans'}
sub ncbi_tax_id {'6239'};
sub ncbi_bioproject {'PRJNA275000'};
sub bioproject_description { 'University of Washington C. elegans CB4856 genome project' }
sub assembly_type {'contig'};
sub is_canonical { 0 };

#######################################################

package Elegans_vc2010;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ELE/elegans.lib'}
sub short_name {'C. elegans'}
sub gspecies_name{'c_elegans'}
sub long_name{'Caenorhabditis elegans'}
sub ncbi_tax_id {'6239'};
sub ncbi_bioproject {'PRJEB28388'};
sub bioproject_description { 'Caenorhabditis elegans strain VC2010 genome sequencing project' }
sub assembly_type {'chromosomes'};
sub is_canonical { 0 };

#######################################################

package Cinopinata;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ELE/elegans.lib'}
sub short_name {'C. inopinata'}
sub gspecies_name{'c_inopinata'}
sub long_name{'Caenorhabditis inopinata'}
sub ncbi_tax_id {'1978547'};
sub ncbi_bioproject {'PRJDB5687'};
sub bioproject_description { 'University of Miyazaki' }
sub assembly_type {'contig'};
sub is_canonical { 1 };

#######################################################

package Clatens;
use Carp;
our @ISA = qw(Species);

sub repeatmasker_library{my($self)=@_;$self->misc_static.'/REPEATMASKER/ELE/elegans.lib'}
sub short_name {'C. latens'}
sub gspecies_name{'c_latens'}
sub long_name{'Caenorhabditis latens'}
sub ncbi_tax_id {'1503980'};
sub ncbi_bioproject {'PRJNA248912'};
sub bioproject_description { 'University of Oregon' }
sub assembly_type {'contig'};
sub is_canonical { 1 };

#######################################################

package Otipulae;
use Carp;
our @ISA = qw(Species);

sub short_name {'O. tipulae'}
sub gspecies_name{'o_tipulae'}
sub long_name{'Oscheius tipulae'}
sub ncbi_tax_id {141969};
sub ncbi_bioproject {'PRJEB15512'};
sub bioproject_description { 'The genome of Oscheius tipulae' }
sub assembly_type {'contig'};
sub is_canonical { 1 };

1;

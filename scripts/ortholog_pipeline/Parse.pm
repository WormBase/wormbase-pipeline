# Generic parsing subroutines of analysis data
package Parse;

use strict 'vars';
use Carp 'croak','cluck';
#use Rearrange;
use vars qw/@ISA/;
@ISA = '';


# Paths to various analysis files
use constant ELEGANS  => "/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/elegans_positions.out";
use constant BRIGGSAE => "/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/cb25.hybrid.genes2.00.gff";

use constant ORTHOLOGS => "/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/ortholog_genes-current";
use constant ORTHOLOGS_TEMP => "/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/orthologs.all";

use constant ELEGANS_BLAST => '/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/bestgenes_el_brig';
use constant BRIGGSAE_BLAST => '/nfs/team71/worm/ar2/wormbase/scripts/ortholog_pipeline/bestgenes_brig_el';

use constant SW => '/Users/todd/projects/briggsae/data/blastp-current/SW/split/elegans_vs_briggsae.out';

#use constant ELEGANS  => '/Users/todd/projects/briggsae/data/gene_set-current/longest/elegans_positions.out';
#use constant BRIGGSAE => '/Users/todd/projects/briggsae/data/gene_set-current/cb25.hybrid.genes2.00.gff';

#use constant KAKS      => '/Users/todd/projects/briggsae/data/orthologs-2.00/orthologs-kaks-needle-2_00.rnai.dat';
#use constant ORTHOLOGS => '/Users/todd/projects/briggsae/data/orthologs-current';
#use constant ORTHOLOGS_TEMP => '/Users/todd/projects/briggsae/data/orthologs-current/orthologs.all';

#use constant ELEGANS_BLAST => '/Users/todd/projects/briggsae/data/blastp-current/seg-on/elegans_vs_briggs.out';
#use constant BRIGGSAE_BLAST => '/Users/todd/projects/briggsae/data/blastp-current/seg-on/briggs_vs_elegans.out';

#use constant SW => '/Users/todd/projects/briggsae/data/blastp-current/SW/split/elegans_vs_briggsae.out';

########################
##   NEW CONSTRUCTOR  ##
########################
sub new {
  my ($self,@p) = @_;
  my $this = bless {},$self;
  return $this;
}



# Pass chrom to get genes seperated by chromosome,
# or gene to get simple hash keyed by gene.
sub elegans {
  my ($self,$flag) = @_;
  open IN,ELEGANS or die "$! ",ELEGANS,"\n";
  
  my $parsed = {};
  while (<IN>) {
    chomp;
    next if (/^\#/);
    # The columns of the elegans positions file...
    # 0 gene protein clone chrom strand 
    # 5 pep_length spliced_length 
    # 7 start stop chrom_start chrom_stop 
    # 11 genome_start genome_stop 
    # 13 gmap source version method
    
    my ($gene,$protein,$clone,$chrom,$strand,$pep_length,$spliced_length,$start,$stop,$chrom_start,$chrom_stop,$genome_start,$genome_stop,$gmap,$source,$version,$method) 
      = split("\t");
    
    if ($flag eq 'chrom') {
      push (@{$parsed->{$chrom}},{ gene           => $gene,
				   protein        => $protein,
				   chrom          => $chrom,
				   strand         => $strand,
				   pep_length     => $pep_length,
				   spliced_length => $spliced_length,
				   start          => $start,
				   stop           => $stop,
				   chrom_start    => $chrom_start,
				   chrom_stop     => $chrom_stop,
				   genome_start   => $genome_start,
				   genome_stop    => $genome_stop,
				   gmap           => $gmap,
				   source         => $source,
				   version        => $version,
				   method         => $method });
    } else {
      my $key = ($flag eq 'protein') ? $protein : $gene;
      $parsed->{$key} = { gene           => $gene,
			   protein        => $protein,
			   chrom          => $chrom,
			   strand         => $strand,
			   pep_length     => $pep_length,
			   spliced_length => $spliced_length,
			   start          => $start,
			   stop           => $stop,
			   chrom_start    => $chrom_start,
			   chrom_stop     => $chrom_stop,
			   genome_start   => $genome_start,
			   genome_stop    => $genome_stop,
			   gmap           => $gmap,
			   source         => $source,
			   version        => $version,
			   method         => $method };
    }
  }
  return $parsed;
}



# Pass in ng or ml to check for different KaKs values.
sub kaks {
  my ($self,$flag) = @_;
  open IN,KAKS or warn "$!\n";
  my $parsed = {};
  while (<IN>) {
    chomp;
    next if (/^\#/);
    my ($ce,$cb,$method,$rnai,$ng_ka,$ng_ks,$ng_omega,$ml_ka,$ml_ks,$ml_omega)
      = split(",");

    my ($kaks,$ka,$ks);
    if ($flag eq 'ng') {
      $kaks = $ng_omega;
      $ks   = $ng_ks;
      $ka   = $ng_ka;
    } else {
      $kaks = $ml_omega;
      $ks   = $ml_ks;
      $ka   = $ml_ka;
    }

    # Filtering
    #    next if ($ng_ks < 0 || $ng_ks > 4);
    next if ($ng_ks < 0 );

    $parsed->{$ce}->{kaks} = $kaks;
    $parsed->{$ce}->{ks}   = $ks;
    $parsed->{$ce}->{ka}   = $ka;
    $parsed->{$ce}->{cb}    = $cb;
  }
  return $parsed;
}



# Pass 'super' to organize returned hash by genes on the supercontig.
# or 'gene' to simply key the hash by the briggsae gene
sub briggsae {
  my ($self,$flag) = @_;
  open IN,BRIGGSAE or die "$!\n";
  my $parsed = {};
  while (<IN>) {
    chomp;
    my ($super,$source,$method,$start,$stop,$phase,$strand,$junk,$full_gene) = split("\t");
    my $gene = ($full_gene =~ /Sequence\s\"(.*)\"/) ? $1 : $full_gene;
    if ($flag eq 'super') {
      push (@{$parsed->{$super}},{ supercontig => $super,
				   source      => $source,
				   method      => $method,
				   start       => $start,
				   stop        => $stop,
				   phase       => $phase,
				   strand      => $strand,
				   gene        => $gene});
    } else {
      $parsed->{$gene} = { supercontig => $super,
			   source      => $source,
			   method      => $method,
			   start       => $start,
			   stop        => $stop,
			   phase       => $phase,
			   strand      => $strand,
			   gene        => $gene};
    }
  }
  return $parsed;
}





sub orthologs_temp {
  my ($self,$file) = @_;
  if ($file) {
    open IN,$file or warn "$!\n";
  } else {
    open IN,ORTHOLOGS_TEMP or warn "$!\n";
  }
  my $orthos = {};
  while (<IN>) {
    chomp;
    my ($cb,$ce,$e,$conf,$method) = split("\t");
    $orthos->{$cb}->{mate} = $ce;
    $orthos->{$cb}->{eval} = $e;
    $orthos->{$cb}->{meth} = $method;
    $orthos->{$cb}->{conf} = $conf;

    $orthos->{$ce}->{ce}   = $ce;
    $orthos->{$ce}->{cb}   = $cb;
    $orthos->{$ce}->{mate} = $cb;
    $orthos->{$ce}->{eval} = $e;
    $orthos->{$ce}->{meth} = $method;
    $orthos->{$ce}->{conf} = $conf;
  }
  return $orthos;
}


sub orthologs {
  my ($self,$file) = @_;
  # If I'm passing a flag, then I'm trying to parse the
  # extended format pseudo-gff file
  my $orthos = {}; 
  if ($file) {
    open IN,$file or warn "$!\n";
    while (<IN>) {
      chomp;
      my @cols = qw/chrom ce ce_start ce_stop ce_strand
	super cb cb_start cb_stop cb_strand
	  evalue per_id conf meth/;
      my @f = split("\t");
      my $c;
      my $temp = {};
      foreach my $val (@f) {
	$temp->{$cols[$c]} = $val;
	$c++;
      }
      $orthos->{$temp->{ce}} = $temp;
    }
  } else {
    open IN,ORTHOLOGS or warn "$!\n";
    while (<IN>) {
      chomp;
      my ($cb,$ce,$e,$conf,$method) = split("\t");
      $orthos->{$cb}->{mate} = $ce;
      $orthos->{$ce}->{mate} = $cb;
    }
  }
  return $orthos;
}

sub blast {
  my ($self,$flag) = @_;
  if ($flag eq 'elegans') {
    open IN,ELEGANS_BLAST or warn "$!\n";
  } elsif ($flag eq 'briggsae') {
    open IN,BRIGGSAE_BLAST or warn "$!\n";
  } else {
    open IN,$flag or warn "$! $flag\n";
  }
  
  my $parsed = {};
  while (<IN>) {
    chomp;
    s/\w{2}://;
    my @temp = split(/\,/);
    
    my $evalue  = $temp[2];
    if (substr($evalue,0,1) eq 'e') { $evalue = "1".$evalue;}
    my $query   = $temp[0];
    my $subject = $temp[1];
    push (@{$parsed->{$query}},[$subject,$evalue]);
  }
  return $parsed;
}


sub sw {
  # columns are query subject evalue score sw-score bits pid gid alen qlen hlen qstart qend hstart hend
  open IN,SW or die;
  my $hash = {};
  while (<IN>) {
    chomp;
    next if (/^\#/);
    my ($ce,$cb,$evalue,$bits,$per_id,$qs,$qe,$ss,$se) = split("\t");

    push (@{$hash->{$cb}},{cb => $cb,
			   ce => $ce,
			   eval => $evalue,
			   bits => $bits,
			   per_id => $per_id,
			   qs => $qs,
			   qe => $qe,
			   ss => $ss,
			   se => $se });
  }
  return $hash;
}

1;

#!/usr/bin/env perl

=pod

=head1 NAME 

genbank2gff3.pl

=head1 SYNOPSIS

  genbank2gff3.pl [options] filename(s)


    Options:
        --filter   -x   genbank feature type(s) to ignore
        --ethresh  -e   error threshold for unflattener
                        set this high (>2) to ignore all unflattener errors
        --format   -f   Input format (SeqIO types): GenBank, Swiss or Uniprot, EMBL work
                        (GenBank is default)
        --quiet         don't talk about what is being processed 
        --typesource    SO sequence type for source (e.g. chromosome; region; contig)
        --help      -h  display this message
        --gffsource     String to use for GFF source (default: WormBase_imported)


=head1 DESCRIPTION

Locally hacked version of bp_genbank2gff3.pl

=cut

use strict;
use warnings;

#use Pod::Usage;
#use Bio::Root::RootI;
use Bio::SeqIO;
use File::Spec;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::IDHandler;
use Bio::Location::SplitLocationI;
use Bio::Location::Simple;

use Getopt::Long;
use List::Util qw(first);
use File::Basename; 

my (@FILTER, 
    @GFF_LINE_FEAT, 
    $DEBUG,
    $TYPE_MAP,
    $ETHRESH,
    $SOURCEID,
    $verbose, 
    $help,
    $source_type,
    $gene_id,
    $ncrna_id,
    $rna_id, 
    $tnum,
    $rnum,
    $outfile,
    $didheader,
    %method,
    %exonpar,

    );


use constant GM_NEW_TOPLEVEL => 2;
use constant GM_NEW_PART => 1;
use constant GM_DUP_PART => 0;
use constant GM_NOT_PART => -1;
  
$verbose = 1; # right default? -nov to turn off

my $CDSKEEP      = 1;
my $PROTEIN_TYPE = 'polypeptide';
my $FORMAT       = "GenBank";
    
my %TAG_MAP = (
  db_xref => 'Dbxref',
  name    => 'Name',
  note    => 'Note',
  synonym => 'Alias',
  symbol  => 'Alias',
);


my $quiet= !$verbose;
my $ok= GetOptions( 
  'h|help'        => \$help,
  'x|filter:s'    => \@FILTER,
  "ethresh|e=s"   => \$ETHRESH,
  'c|CDS!'        => \$CDSKEEP,
  'f|format=s'    => \$FORMAT,
  'typesource=s'  => \$source_type,
  'quiet!'        => \$quiet, # swap quiet to verbose
  'DEBUG!'        => \$DEBUG,
  'outfile=s'     => \$outfile,
  'gffsource=s'   => \$SOURCEID,
    );


$verbose= !$quiet;

$FORMAT   = "swiss" if $FORMAT =~/UniProt|trembl/;
$verbose  = 1 if($DEBUG);

my $tm  = Bio::SeqFeature::Tools::TypeMapper->new;
my $idh = Bio::SeqFeature::Tools::IDHandler->new;

$source_type ||= "region"; # should really parse from FT.source contents below
$SOURCEID ||= "WormBase_imported";

#my $FTSOmap = $tm->FT_SO_map();
my $FTSOmap;

my %terms = %{ $tm->FT_SO_map() };
while (my ($k,$v) = each %terms) {
  $FTSOmap->{$k} = ref($v) ? shift @$v : $v;
}


$TYPE_MAP = $FTSOmap;


# #convert $FTSOmap undefined to valid SO : moved to TypeMapper->map_types( -undefined => "region")

# stringify filter list if applicable

my $in = Bio::SeqIO->new(-fh => \*STDIN, -format => $FORMAT, -debug=>$DEBUG);
my $out;
if ($outfile) {
  open($out, ">$outfile") or die("Could not open $outfile for writing\n");
} else {
  $out = \*STDOUT;
}

while ( my $seq = $in->next_seq() ) {
  my $unflattener = Bio::SeqFeature::Tools::Unflattener->new; # for ensembl genomes (-trust_grouptag=>1);
  $unflattener->error_threshold($ETHRESH) if $ETHRESH;
  $unflattener->verbose(1) if($DEBUG);
  
  my $seq_name = $seq->accession_number;
  my $end = $seq->length;
  my @to_print;

  # filter out unwanted features
  my ($source_feat) = filter($seq);
  
  ($source_type,$source_feat)= 
      getSourceInfo( $seq, $source_type, $source_feat ) ;
  # always; here we build main prot $source_feat; # if @source;
  
  # abort if there are no features
  warn "$seq_name has no features, skipping\n" and next
      if !$seq->all_SeqFeatures;
  
  $FTSOmap->{source} = $source_type;

  # construct a GFF header
  # add: get source_type from attributes of source feature? chromosome=X tag
  # also combine 1st ft line here with source ft from $seq ..
  my($header,$info)= gff_header($seq_name, $end, $source_type, $source_feat);
  print $out $header;
  print "# working on $info\n" if($verbose);
  
  # unflatten gene graphs, apply SO types, etc; this also does TypeMapper ..
  unflatten_seq($unflattener, $seq);
  
  # Note that we use our own get_all_SeqFeatures function 
  # to rescue cloned exons
  
  @GFF_LINE_FEAT = ();
  
  my (@cds, @exons);
  for my $feature ( get_all_SeqFeatures($seq) ) {
    my $method = $feature->primary_tag;
    next if($SOURCEID =~/UniProt|swiss|trembl/i && $method ne $source_type);
    #next if $method eq 'chromosome';

    $feature->seq_id($seq->id) unless($feature->seq_id);
    $feature->source_tag($SOURCEID);
    
    # dgg; need to convert some Genbank to GFF tags: note->Note; db_xref->Dbxref;
    ## also, pull any GO:000 ids from /note tag and put into Ontology_term
    maptags2gff($feature);
    
    # current gene name.  The unflattened gene features should be in order so any
    # exons, CDSs, etc that follow will belong to this gene
    my $gene_name;
    if ( $method eq 'gene' || $method eq 'pseudogene' ) {
      @to_print= print_held($out, \@to_print);
      $gene_id = $gene_name= gene_name($feature); 
    } else {
      $gene_name= gene_name($feature);
    }
    
    #?? should gene_name from /locus_tag,/gene,/product,/transposon=xxx
    # be converted to or added as  Name=xxx (if not ID= or as well)
    ## problematic: convert_to_name ($feature); # drops /locus_tag,/gene, tags
    convert_to_name($feature); 
    
    ## dgg: extended to protein|polypeptide
    ## this test ($feature->has_tag('gene') ||) is not good: repeat_regions over genes
    ## in yeast have that genbank tag; why?
    ## these include pseudogene ...
    
    ## Note we also have mapped types to SO, so these RNA's are now transcripts:
    # pseudomRNA => "pseudogenic_transcript", 
    # pseudotranscript" => "pseudogenic_transcript", 
    # misc_RNA=>'processed_transcript',
    
    warn "#at: $method $gene_id/$gene_name\n" if $DEBUG;
    
    if ( $method =~ /(gene|RNA|CDS|exon|UTR|protein|polypeptide|transcript)/ 
         || ( $gene_id && $gene_name eq $gene_id ) ) {
      
      my $action = gene_features($feature, $gene_id, $gene_name);  # -1, 0, 1, 2 result
      if ($action == GM_DUP_PART) {
        # ignore, this is dupl. exon with new parent ...
        
      } elsif ($action == GM_NOT_PART) {
        add_generic_id( $feature, $gene_name, "nocount");
        my $gff = gff_string($feature);
        push @GFF_LINE_FEAT, $feature;
        
      } elsif ($action > 0) {
        # hold off print because exon etc. may get 2nd, 3rd parents
        @to_print= print_held($out, \@to_print) if ($action == GM_NEW_TOPLEVEL);
        push(@to_print, $feature);
      }
      
    }
    
    # otherwise handle as generic feats with IDHandler labels 
    else {
      add_generic_id( $feature, $gene_name, "");
      my $gff= gff_string($feature);
      push @GFF_LINE_FEAT, $feature;
    }
    
    
    if ($feature->primary_tag eq 'CDS')  {
      push @cds,  $feature;
    } elsif ($feature->primary_tag eq 'exon') {
      push @exons, $feature;
    }
    
  }
  
  @to_print= print_held($out, \@to_print);
  
  gff_validate(@GFF_LINE_FEAT);
  
  for my $feature (@GFF_LINE_FEAT) {
    my $gff= gff_string($feature);
    print $out "$gff\n" if $gff;
  }
  
}



##################################################

sub typeorder {
  return 1 if ($_[0] =~ /gene/);
  return 2 if ($_[0] =~ /RNA|transcript/);
  return 3 if ($_[0] =~ /protein|peptide/);
  return 4 if ($_[0] =~ /exon|CDS/);
  return 3; # default before exon (smallest part)
}

sub sort_by_feattype {
  my($at,$bt)= ($a->primary_tag, $b->primary_tag);
  return (typeorder($at) <=> typeorder($bt))  
      or ($at cmp $bt);
  ## or ($a->name() cmp $b->name());
}

sub print_held {  
  my($out,$to_print)= @_;
  return unless(@$to_print);
  @$to_print = sort sort_by_feattype  @$to_print; # put exons after mRNA, otherwise chado loader chokes
  while ( my $feature = shift @$to_print) {
    my $gff= gff_string($feature); # $gff =~ s/\'/./g; # dang bug in encode
    push @GFF_LINE_FEAT, $feature;
  }
  return (); # @to_print
}

sub maptags2gff {
  my $f = shift;
  ## should copy/move locus_tag to Alias, if not ID/Name/Alias already
  # but see below /gene /locus_tag usage
  foreach my $tag (keys %TAG_MAP) {
    if ($f->has_tag($tag)) {
      my $newtag= $TAG_MAP{$tag};
      my @v= $f->get_tag_values($tag);
      $f->remove_tag($tag);
      $f->add_tag_value($newtag,@v);
      
      ## also, pull any GO:000 ids from /note tag and put into Ontology_term
      ## ncbi syntax in CDS /note is now '[goid GO:0005886]' OR '[goid 0005624]'
      if ($tag eq 'note') {  
        map { s/\[goid (\d+)/\[goid GO:$1/g; } @v; 
        my @go= map { m/(GO:\d+)/g } @v; 
        $f->add_tag_value('Ontology_term',@go) if(@go);
      }  
    }
  }    
}


sub getSourceInfo {
  my ($seq, $source_type, $sf) = @_;  
  
  my $is_swiss = ($SOURCEID =~/UniProt|swiss|trembl/i);
  my $is_gene  = ($SOURCEID =~/entrezgene/i);
  my $is_rich  = (ref($seq) =~ /RichSeq/);
  my $seq_name = $seq->accession_number();
  
  unless($sf) { # make one
    $source_type=  $is_swiss ? $PROTEIN_TYPE 
        : $is_gene ? "eneg" # "gene"  # "region" # 
        : $is_rich ? $seq->molecule : $source_type;
    $sf = Bio::SeqFeature::Generic->direct_new();
    
    my $len = $seq->length(); $len=1 if($len<1); my $start = 1; ##$start= $len if ($len<1);
    my $loc= $seq->can('location') ? $seq->location() 
        :  new Bio::Location::Simple( -start => $start, -end => $len);
    $sf->location( $loc ); 
    $sf->primary_tag($source_type);
    $sf->source_tag($SOURCEID);
    $sf->seq_id( $seq_name);
    #? $sf->display_name($seq->id()); ## Name or Alias ?
    $sf->add_tag_value( Alias => $seq->id()); # unless id == accession
    $seq->add_SeqFeature($sf);
    ## $source_feat= $sf;
  }
  
  if ($sf->has_tag("chromosome")) {
    my ($chrname) = $sf->get_tag_values("chromosome");
    if ($chrname ne 'Unknown') {
      $source_type = "chromosome"; 
      $sf->add_tag_value( Alias => $chrname );
    }
  }
  
  # pull GB Comment, Description for source ft ...
  # add reference - can be long, not plain string...             
  warn "# $SOURCEID:$seq_name fields = ", join(",", $seq->annotation->get_all_annotation_keys()),"\n" if $DEBUG;
  # GenBank   fields: keyword,comment,reference,date_changed
  # Entrezgene fields 850293 =ALIAS_SYMBOL,RefSeq status,chromosome,SGD,dblink,Entrez Gene Status,OntologyTerm,LOCUS_SYNONYM
  
  # is this just for main $seq object or for all seqfeatures ?
  my %AnnotTagMap= ( 
    gene_name      => 'Alias',   
    ALIAS_SYMBOL   => 'Alias',
    LOCUS_SYNONYM  => 'Alias',
    symbol         => 'Alias',   
    synonym        => 'Alias',   
    dblink         => 'Dbxref',   
    product        => 'product',
    Reference      => 'reference',
    OntologyTerm   => 'Ontology_term',
    #comment        => 'Note',
    #comment1       => 'Note',
      );
  
  
  my ($desc)= $seq->annotation->get_Annotations("desc") || ( $seq->desc() );
  my ($date)= $seq->annotation->get_Annotations("dates") 
      || $seq->annotation->get_Annotations("update-date")
      || $is_rich ? $seq->get_dates() : ();
  my ($comment)= $seq->annotation->get_Annotations("comment");
  my ($species)= $seq->annotation->get_Annotations("species");
  if (!$species 
      && $seq->can('species') 
      && defined $seq->species() 
      && $seq->species()->can('binomial') ) {
    $species= $seq->species()->binomial();
  }
  
  # update source feature with main GB fields
  $sf->add_tag_value( ID => $seq_name ) unless $sf->has_tag('ID');
  $sf->add_tag_value( Note => $desc ) if($desc && ! $sf->has_tag('Note'));
  $sf->add_tag_value( organism => $species ) if($species && ! $sf->has_tag('organism'));
  $sf->add_tag_value( date => $date ) if($date && ! $sf->has_tag('date'));  
  $sf->add_tag_value( Dbxref => $SOURCEID.':'.$seq_name ) if $is_swiss || $is_gene;
  
  foreach my $atag (sort keys %AnnotTagMap) {
    my $gtag= $AnnotTagMap{$atag}; next unless($gtag);
    my @anno = map{ 
      if (ref $_ && $_->can('get_all_values')) { 
        split( /[,;] */, join ";", $_->get_all_values) 
      }
      elsif (ref $_ && $_->can('display_text')) { 
        split( /[,;] */, $_->display_text) 
      }
      elsif (ref $_ && $_->can('value')) { 
        split( /[,;] */, $_->value) 
      } 
      else {
        ();
      }
    } $seq->annotation->get_Annotations($atag);  
    foreach(@anno) { $sf->add_tag_value( $gtag => $_ ); }
  }
  
  return (wantarray)? ($source_type,$sf) : $source_type; #?
}


sub gene_features {
  my ($f, $gene_id, $genelinkID) = @_;
  local $_ = $f->primary_tag;
  $method{$_}++;
  
  if ( /gene/ ) {
    $f->add_tag_value( ID => $gene_id ) unless($f->has_tag('ID')); # check is same value!? 
    $tnum = $rnum= 0; $ncrna_id= $rna_id  = '';
    return GM_NEW_TOPLEVEL;
    
  } elsif ( /mRNA/ ) {  
    return GM_NOT_PART unless $gene_id;
    return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
    ($rna_id  = $gene_id ) =~ s/gene/mRNA/;
    $rna_id   .= '.t0' . ++$tnum;
    $f->add_tag_value( ID => $rna_id );
    $f->add_tag_value( Parent => $gene_id );
    
  } elsif ( /RNA|transcript/) { 
    ## misc_RNA here; missing exons ... flattener problem?
    #  all of {t,nc,sn}RNA can have gene models now
    ## but problem in Worm chr: mRNA > misc_RNA > CDS with same locus tag
    ## CDS needs to use mRNA, not misc_RNA, rna_id ...
    ## also need to fix cases where tRNA,... lack a 'gene' parent: make this one top-level
    
    if($gene_id) {
      return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
      ($ncrna_id = $gene_id) =~ s/gene/ncRNA/;
      $ncrna_id .= '.r0' . ++$rnum;
      $f->add_tag_value( Parent => $gene_id );
      $f->add_tag_value( ID => $ncrna_id );
    } else {
      unless ($f->has_tag('ID')) {
        if($genelinkID) {
          $f->add_tag_value( ID => $genelinkID  ) ;
        } else {
          $idh->generate_unique_persistent_id($f);
        }
      }
      ($ncrna_id)= $f->get_tag_values('ID'); 
      return GM_NEW_TOPLEVEL;
      # this feat now acts as gene-top-level; need to print @to_print to flush prior exons?
    }
    
  } elsif ( /exon/ ) { # can belong to any kind of RNA
    return GM_NOT_PART unless ($rna_id||$ncrna_id);
    return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
    ## we are getting duplicate Parents here, which chokes chado loader, with reason...
    ## problem is when mRNA and ncRNA have same exons, both ids are active, called twice
    ## check all Parents
    for my $expar ($rna_id, $ncrna_id) { 
      next unless($expar);
      if ( $exonpar{$expar} and $f->has_tag('Parent') ) {
        my @vals = $f->get_tag_values('Parent');
        next if (grep {$expar eq $_} @vals);
      }
      $exonpar{$expar}++;
      $f->add_tag_value( Parent => $expar); 
      # last; #? could be both
    }
    
  } elsif ( /CDS|protein|polypeptide/ ) {

    if( ! $CDSKEEP && /CDS/) {  
      $f->primary_tag($PROTEIN_TYPE); 
      
      ## duplicate problem is Location ..
      if ($f->location->isa("Bio::Location::SplitLocationI")) {
        # my($b,$e)=($f->start, $f->end); # is this all we need?
        my($b,$e)=(-1,0);
        foreach my $l ($f->location->each_Location) {
          $b = $l->start if($b<0 || $b > $l->start);
          $e = $l->end if($e < $l->end);
        }
        $f->location( Bio::Location::Simple->new(
                        -start => $b, -end => $e, -strand => $f->strand) );
      }
      
      $f->add_tag_value( Derives_from => $rna_id ); 
    }
    else {
      $f->add_tag_value( Parent => $rna_id );
      if ($f->primary_tag eq 'CDS') {
        my ($codon_start);
        eval {
          ($codon_start) = $f->get_tag_values('codon_start');
        };
        if (not $@) {
          $f->frame($codon_start - 1);
        } else {
          $f->frame(0);
        }
      }
    }
    
    return GM_NOT_PART unless $rna_id;
    return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id); #??

    (my $pro_id = $rna_id) =~ s/\.t/\.p/;
    
    $f->add_tag_value( ID => $pro_id );
    
    move_translation_fasta($f, $pro_id);
  } elsif ( /region/ ) {       
    $f->primary_tag('gene_component_region');
    $f->add_tag_value( Parent => $gene_id ); 
  } else {
    return GM_NOT_PART unless $gene_id;
    $f->add_tag_value( Parent => $gene_id );  
  }
  
  ## return GM_DUP_PART if /exon/ && ++$seen{$f} > 1;
  
  return GM_NEW_PART;
}

## was generic_features >  add_generic_id
sub add_generic_id {
  my ($f, $ft_name, $flags) = @_;
  my $method = $f->primary_tag;
  $method{$method}++ unless($flags =~ /nocount/); ## double counts GM_NOT_PART from above
  
  if ($f->has_tag('ID')) {
    
  }
  elsif ( $f->has_tag($method) ) {
    my ($name) = $f->get_tag_values($method);
    $f->add_tag_value( ID => "$method:$name" );
  }
  elsif($ft_name) { # is this unique ?
    $f->add_tag_value( ID => $ft_name ); 
  }
  else {
    $idh->generate_unique_persistent_id($f);
  }
  
  move_translation_fasta( $f, ($f->get_tag_values("ID"))[0] )
      if($method =~ /CDS/);
  
  # return $io->gff_string($f);
}

sub move_translation_fasta {
  my ($f, $ft_id) = @_;
  if( $ft_id && $f->has_tag('translation') ) {
    my ($aa) = $f->get_tag_values("translation");
    if($aa && $aa !~ /^length/) {
      $f->remove_tag("translation");
      $f->add_tag_value("translation","length.".length($aa)); # hack for odd chado gbl problem
    }
  }
}

sub gff_header {
  my ($name, $end, $source_type, $source_feat) = @_;
  $source_type ||= "region";
  
  my $info = "$source_type:$name";
  my $head = "##gff-version 3\n".
      "##sequence-region $name 1 $end\n".
      "# conversion-by bp_genbank2gff3.pl\n";
  if ($source_feat) {
    ## dgg: these header comment fields are not useful when have multi-records, diff organisms
    for my $key (qw(organism Note date)) {
      my $value;
      if ($source_feat->has_tag($key)) { 
        ($value) = $source_feat->get_tag_values($key);
      }
      if ($value) {
        $head .= "# $key $value\n";
        $info .= ", $value";
      }
    }
        $head = "" if $didheader;
  } else {
    $head .= "$name\t$SOURCEID\t$source_type\t1\t$end\t.\t.\t.\tID=$name\n";
  }
  $didheader++;
  return (wantarray) ? ($head,$info) : $head;
}

sub unflatten_seq {
  my ($unflattener, $seq) = @_;
  
  ## print "# working on $source_type:", $seq->accession, "\n"; 
  my $uh_oh = "Possible gene unflattening error with" .  $seq->accession_number .
      ": consult STDERR\n";
  
  eval {
    $unflattener->unflatten_seq( -seq => $seq, 
                                 -use_magic => 1 );
  };
  
  # deal with unflattening errors
  if ( $@ ) {
    warn $seq->accession_number . " Unflattening error:\n";
    warn "Details: $@\n";
    print "# ".$uh_oh;
  }
  
  return 0 if !$seq || !$seq->all_SeqFeatures;
  
  # map feature types to the sequence ontology
  ## $tm->map_types_to_SO( -seq => $seq );
  #$tm->map_types( -seq => $seq, -type_map => $FTSOmap, -undefined => "region" ); #dgg
  
  map_types( 
    $tm, 
    -seq => $seq, 
    -type_map  => $FTSOmap, 
    -undefined => "region" 
      ); #nml
  
}

sub filter {
  my $seq = shift;

  my (@feats, @sources);

  if (@FILTER) {
    for my $f ( $seq->remove_SeqFeatures ) {
      my $m = $f->primary_tag;
      push @sources, $f if ($m eq 'source'); # dgg? but leave in @feats ?
      push @feats, $f unless grep { $m eq $_ } @FILTER; 
    }
    $seq->add_SeqFeature($_) foreach @feats;
  } else {
    for my $f ( $seq->get_SeqFeatures ){
      my $m = $f->primary_tag;
      push @sources, $f if ($m eq 'source'); # dgg? but leave in @feats ?
    }
  }
  
  return @sources;
}


# The default behaviour of Bio::FeatureHolderI:get_all_SeqFeatures
# changed to filter out cloned features.  We have to implement the old
# method.  These two subroutines were adapted from the v1.4 Bio::FeatureHolderI
sub get_all_SeqFeatures  {
  my $seq = shift;
  my @flatarr;
  
  foreach my $feat ( $seq->get_SeqFeatures ){
    push(@flatarr,$feat);
    _add_flattened_SeqFeatures(\@flatarr,$feat);
  }
    return @flatarr;
}

sub gene_name {
  my $g = shift;
  my $gene_id = ''; # zero it;
  
  if ($g->has_tag('locus_tag')) {
    ($gene_id) = $g->get_tag_values('locus_tag');
  }
  
  elsif ($g->has_tag('gene')) {
    ($gene_id) = $g->get_tag_values('gene'); 
  }
  elsif ($g->has_tag('ID')) { # for non-Genbank > Entrezgene
    ($gene_id) = $g->get_tag_values('ID');
  }

  ## See Unflattener comment:
  # on rare occasions, records will have no /gene or /locus_tag
  # but it WILL have /product tags. These serve the same purpose
  # for grouping. For an example, see AY763288 (also in t/data)
  # eg. product=tRNA-Asp ;  product=similar to crooked neck protein
  elsif ($g->has_tag('product')) {
    my ($name)= $g->get_tag_values('product');
    ($gene_id) = $name unless($name =~ / /); # a description not name
  }
  
  ## dgg; also handle transposon=xxxx ID/name
  # ID=GenBank:repeat_region:NC_004353:1278337:1281302;transposon=HeT-A{}1685;Dbxref=FLYBASE:FBti0059746
  elsif ($g->has_tag('transposon')) {
    my ($name)= $g->get_tag_values('transposon');
    ($gene_id) = $name unless($name =~ / /); # a description not name
  }
  
  return $gene_id;
}

# same list as gene_name .. change tag to generic Name
sub convert_to_name {
  my $g = shift;
  my $gene_id = ''; # zero it;
  
  if ($g->has_tag('gene')) {
    ($gene_id) = $g->get_tag_values('gene'); 
    $g->remove_tag('gene');
    $g->add_tag_value('Name', $gene_id);
    }
  elsif ($g->has_tag('locus_tag')) {
    ($gene_id) = $g->get_tag_values('locus_tag');
    $g->remove_tag('locus_tag');
    $g->add_tag_value('Name', $gene_id);
  }
  
  elsif ($g->has_tag('product')) {
    my ($name)= $g->get_tag_values('product');
    ($gene_id) = $name unless($name =~ / /); # a description not name
    ## $g->remove_tag('product');
    $g->add_tag_value('Name', $gene_id);
  }
  
  elsif ($g->has_tag('transposon')) {
    my ($name)= $g->get_tag_values('transposon');
    ($gene_id) = $name unless($name =~ / /); # a description not name
    ## $g->remove_tag('transposon');
    $g->add_tag_value('Name', $gene_id);
  }
  elsif ($g->has_tag('ID')) { 
    my ($name)= $g->get_tag_values('ID');
    $g->add_tag_value('Name', $name);
  }    
  return $gene_id;
}


sub _add_flattened_SeqFeatures  {
  my ($arrayref,$feat) = @_;
  my @subs = ();
  
  if ($feat->isa("Bio::FeatureHolderI")) {
    @subs = $feat->get_SeqFeatures;
  } 
  elsif ($feat->isa("Bio::SeqFeatureI")) {
    @subs = $feat->sub_SeqFeature;
  }
  else {
    warn ref($feat)." is neither a FeatureHolderI nor a SeqFeatureI. ".
        "Don't know how to flatten.";
  }
  
  for my $sub (@subs) {
    push(@$arrayref,$sub);
    _add_flattened_SeqFeatures($arrayref,$sub);
  }
  
}

sub map_types {
  
  my ($self, @args) = @_;
  
  my($sf, $seq, $type_map, $syn_map, $undefmap) =
      $self->_rearrange([qw(FEATURE
                    SEQ
                    TYPE_MAP
                    UNDEFINED
                    )],
                        @args);
  
  if (!$sf && !$seq) {
    $self->throw("you need to pass in either -feature or -seq");
  }
  
  my @sfs = ($sf);
  if ($seq) {
    $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
    @sfs = $seq->get_all_SeqFeatures;
  }
  $type_map = $type_map || $self->typemap; # dgg: was type_map;
  foreach my $feat (@sfs) {
    
    $feat->isa("Bio::SeqFeatureI") || $self->throw("$feat NOT A SeqFeatureI");
    $feat->isa("Bio::FeatureHolderI") || $self->throw("$feat NOT A FeatureHolderI");
    
    my $primary_tag = $feat->primary_tag;
    
    
    my $mtype = $type_map->{$primary_tag};
    if ($mtype) {
      if (ref($mtype)) {
        if (ref($mtype) eq 'ARRAY') {
          my $soID;
          ($mtype, $soID) = @$mtype;
          
          if ($undefmap && $mtype eq 'undefined') { # dgg
            $mtype= $undefmap;
          }
          
          $type_map->{$primary_tag} = $mtype if $mtype;
        }
        elsif (ref($mtype) eq 'CODE') {
          $mtype = $mtype->($feat);
        }
        else {
          $self->throw('must be scalar or CODE ref');
        }
      }
      elsif ($undefmap && $mtype eq 'undefined') { # dgg
        $mtype= $undefmap;
      }
      $feat->primary_tag($mtype);
    } 
  }
}

sub gff_validate {
    warn "Validating GFF file\n" if $DEBUG;
    my @feat = @_;

    my (%parent2child, %all_ids, %descendants, %reserved);

    for my $f (@feat) {
        for my $aTags (['Parent', \%parent2child], ['ID', \%all_ids]) {
            map { push @{$$aTags[1]->{$_}}, $f } $f->get_tag_values($$aTags[0])
                if $f->has_tag($$aTags[0]); 
        }
    }

    id_validate(\%all_ids, \%reserved);        
}

sub id_validate {
    my ($hAll, $hReserved) = @_;


    for my $id (keys %$hAll) {

        #since 1 feature can have this id, 
        #let's just shift it off and uniquify
        #the rest (unless it's reserved)

        shift @{$hAll->{$id}} unless $hReserved->{$id};
        for my $feat (@{$hAll->{$id}}) {
            id_uniquify(0, $id, $feat, $hAll);
        }
    }

    for my $parentID (keys %$hReserved) {

        my @keys = keys %{$hReserved->{$parentID}};

        shift @keys;

        for my $k (@keys) {
            my $parent = $hReserved->{$parentID}{$k}{'parent'};
            my $aChildren= $hReserved->{$parentID}{$k}{'children'};

            my $value = id_uniquify(0, $parentID, $parent, $hAll);
            for my $child (@$aChildren) {

                my %parents;
                map { $parents{$_}++ } $child->get_tag_values('Parent');
                $child->remove_tag('Parent');
                delete $parents{$parentID};
                $parents{$value}++;
                my @parents = keys %parents;
                $child->add_tag_value('Parent', @parents);
            }

        }
    }
}

sub id_uniquify {
  my ($count, $value, $feat, $hAll) = @_;
  
  warn "UNIQUIFYING $value" if $DEBUG;
  
  if (! $count) {
    $feat->add_tag_value(Alias => $value);
    $value .= ('.' . $feat->primary_tag) 
  } elsif ($count == 1) {
    $value .= ".$count"; 
  } else { 
    chop $value;
    $value .= $count 
  }
  $count++;
  
  warn "ENDED UP WITH $value" if $DEBUG;
  if ( $hAll->{$value} ) { 
    warn "$value IS ALREADY TAKEN" if $DEBUG;
    id_uniquify($count, $value, $feat, $hAll);
  } else { 
    #warn "something's breaking ".$feat->primary_tag.' at '.$feat->start.' to '.$feat->end;
    $feat->remove_tag('ID');
    $feat->add_tag_value('ID', $value);
    push @{$hAll->{$value}}, $value;
  }
  
  $value;
}

my $i = 0;
my %GFF3_ID_Tags = map { $_ => $i++ } qw(ID Parent Target);
my %SKIPPED_TAGS = map { $_ => 1 } qw(score);

my $gff3_featureID = 0;

sub _incrementGFF3ID {
  return ++$gff3_featureID;
}

sub gff_string {
  my ($origfeat) = @_;
  my $feat;
  if ($origfeat->isa('Bio::SeqFeature::FeaturePair')){
    $feat = $origfeat->feature2;
  } else {
    $feat = $origfeat;
  }
  
  my $ID = _incrementGFF3ID();
  
  my ($score,$name,$strand);
  
  if( $feat->can('score') ) {
    $score = $feat->score();
  }
  $score = '.' unless defined $score;
  
  $strand = $feat->strand();
  
  if(! $strand) {
    $strand = ".";
  } elsif( $strand == 1 ) {
    $strand = '+';
  } elsif ( $feat->strand == -1 ) {
    $strand = '-';
  }
  
  if( $feat->can('seqname') ) {
    $name = $feat->seq_id();
    $name ||= 'SEQ';
  } else {
    $name = 'SEQ';
  }
  
  my @groups;
  
  # force leading ID and Parent tags
  my @all_tags =  grep { ! exists $GFF3_ID_Tags{$_} } $feat->all_tags;
  for my $t ( sort { $GFF3_ID_Tags{$b} <=> $GFF3_ID_Tags{$a} }
              keys %GFF3_ID_Tags ) {
    unshift @all_tags, $t if $feat->has_tag($t);
  }
  
  for my $tag ( @all_tags ) {
    next if exists $SKIPPED_TAGS{$tag};
    # next if $tag eq 'Target';
    if ($tag eq 'Target' && ! $origfeat->isa('Bio::SeqFeature::FeaturePair')){  
      # simple Target,start,stop
      my($target_id, $b,$e,$strand) = $feat->get_tag_values($tag); 
      next unless(defined($e) && defined($b) && $target_id);
      ($b,$e)= ($e,$b) if(defined $strand && $strand<0);
      $target_id =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;    
      push @groups, sprintf("Target=%s %d %d", $target_id,$b,$e);
      next;
    }
    
    my $valuestr;
    # a string which will hold one or more values 
    # for this tag, with quoted free text and 
    # space-separated individual values.
    my @v;
    for my $value ( $feat->get_tag_values($tag) ) {
      if(  defined $value && length($value) ) { 
        #$value =~ tr/ /+/;  #spaces are allowed now
        if ( ref $value eq 'Bio::Annotation::Comment') {
          $value = $value->text;
        }
        
        if ($value =~ /[^a-zA-Z0-9\,\;\=\.:\%\^\*\$\@\!\+\_\?\-]/) {
          $value =~ s/\t/\\t/g; # substitute tab and newline 
          # characters
          $value =~ s/\n/\\n/g; # to their UNIX equivalents
          
          # Unescaped quotes are not allowed in GFF3
          #                    $value = '"' . $value . '"';
        }
        $value =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;
      } else {
        # if it is completely empty, then just make empty double quotes
        $value = '""';
      }
      push @v, $value;
    }
    # can we figure out how to improve this?
    $tag = lcfirst($tag) unless ( $tag =~
                                  /^(ID|Name|Alias|Parent|Gap|Target|Derives_from|Note|Dbxref|Ontology_term)$/);
    
    push @groups, "$tag=".join(",",@v);
  }
  # Add Target information for Feature Pairs
  if( $feat->has_tag('Target') && 
      ! $feat->has_tag('Group') &&
      $origfeat->isa('Bio::SeqFeature::FeaturePair') ) {
    
    my $target_id = $origfeat->feature1->seq_id;
    $target_id =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;    
    
    push @groups, sprintf("Target=%s %d %d", 
                          $target_id,
                          ( $origfeat->feature1->strand < 0 ? 
                            ( $origfeat->feature1->end,
                              $origfeat->feature1->start) :
                            ( $origfeat->feature1->start,
                              $origfeat->feature1->end) 
                          ));
  }
  
  # unshift @groups, "ID=autogenerated$ID" unless ($feat->has_tag('ID'));
  if ( $feat->can('name') && defined($feat->name) ) {
    # such as might be for Bio::DB::SeqFeature
    unshift @groups, 'Name=' . $feat->name;
  }
  
  my $gff_string = "";
  my $prim_tag = $feat->primary_tag;
  if ($prim_tag eq 'pseudogene') {
    $prim_tag = 'gene';
  } elsif ($prim_tag =~ /^pseudo/) {
    $prim_tag = 'pseudogenic_transcript';
  }

  if ($feat->location->isa("Bio::Location::SplitLocationI")) {
    my @locs = $feat->location->each_Location;

    my %frames;
    if ($feat->primary_tag eq 'CDS') {
      my $prev_overhang = (3 - $feat->frame) % 3;
      my @these_locs = ($strand eq '-') ? sort { $b->start <=> $a->start } @locs : sort { $a->start <=> $b->start } @locs;
      foreach my $loc (@these_locs) {
        my $this_frame = (3 - $prev_overhang) % 3;
        $frames{$loc->start} = $this_frame;
        $prev_overhang = ($loc->length + 3 - $this_frame) % 3;
      }

    } else {
      map { $frames{$_->start} = "." } @locs; 
    }
    foreach my $loc (@locs) {
      $gff_string .= join("\t",
                          $name,
                          $feat->source_tag() || '.',
                          $prim_tag,
                          $loc->start(),
                          $loc->end(),
                          $score,
                          $strand,
                          $frames{$loc->start}, 
                          join(';', @groups)) . "\n";
    }
    chop $gff_string;
    return $gff_string;
  } else {
    my $frame = ($feat->primary_tag eq 'CDS') ? $feat->frame  : "."; 
    $gff_string = join("\t",
                       $name,
                       $feat->source_tag() || '.',
                       $prim_tag,
                       $feat->start(),
                       $feat->end(),
                       $score,
                       $strand,
                       $frame,
                       join(';', @groups));
  }
  return $gff_string;
  
}

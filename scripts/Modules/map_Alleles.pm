
=pod 

=head1 NAME

map_Alleles.pm

=head1 SYNOPSIS

use map_Alleles.pm

MapAlleles::setup(LogFile::,Wormbase::)

=head1 DESCRIPTION

modules containing functions for map_Alllele.pl
like for searching through GFF files and printing in Ace format

=head1 FUNCTIONS for MapAlleles::

=cut
package MapAlleles;

use strict;

# standard library
use IO::File;
#use Memoize; 

use Bio::Tools::CodonTable;
use Bio::Seq;

# Wormbase
use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/Modules";
use Feature_mapper;
use Sequence_extract;
use Coords_converter;
#memoize('Sequence_extract::Sub_sequence') unless $ENV{DEBUG}; 
use Modules::Remap_Sequence_Change;
use gff_model;

# needs to be attached in main code
my ($log,$wb,$errors,$extractor,$weak_checks, $db);

my ($index); # mem_index test

=head2 get_errors

    Title   : get_errors
    Usage   : $scalar=get_errors()
    Function: accessor for the $errors package variable
    Returns : (integer) number of errors occured
    Args    : none 

=cut
sub get_errors{return $errors}

=head2 setup

    Title   : setup
    Usage   : MapAlleles::setup($log,$wormbase)
    Function: sets up package wide variables
    Returns : 1 on success
    Args    : logfiles object , wormbase object

=cut

# used instead of an custom importer to keep the use statements in one place
sub setup{
  ($log,$wb,$db)=@_;

  $extractor=Sequence_extract->invoke($wb->autoace,undef,$wb);

  return 1;
}


=head2 get_all_alleles

    Title   : get_all_alleles
    Usage   : MapAlleles::get_all_alleles()
    Function: get a list of filtered alleles from acedb
    Returns : array_ref of Ace::Alleles
    Args    : none

=cut

# gets all alles and filters them
sub get_all_alleles {
  my $species=$wb->full_name;

  my @alleles = $db->fetch( -query =>"Find Variation WHERE Flanking_sequences AND Live AND species = \"$species\"");
  
  return \@alleles;
}        



=head2 get_all_allele_ids

    Title   : get_all_allele_ids
    Usage   : MapAlleles::get_all_allele_ids
    Function: get a list of filtered alleles from acedb
    Returns : array_ref of ids
    Args    : none

=cut

sub get_all_allele_ids {

  my $species=$wb->full_name;
  my $iter = $db->fetch_many( -query =>"Find Variation WHERE Flanking_sequences AND Live AND species = \"$species\"");

  my @all_names;
  
  while (my $var = $iter->next) {
    push @all_names, $var->name;
    $var->DESTROY;
  }

  return \@all_names;
}


=head2 get_all_allele_ids_table_maker

    Title   : get_all_allele_ids_table_maker
    Usage   : MapAlleles::get_all_allele_ids_table_maker
    Function: get a list of active alleles from acedb
    Returns : array_ref of ids
    Args    : none

=cut


=head2 get_all_allele

    Title   : get_all_allele
    Usage   : MapAlleles::get_all_allele('tm852')
    Function: get an allele from acedb
    Returns : array_ref containing ONE Ace::Allele
    Args    : (string) alllele name

=cut

# get allele for testing
sub get_allele {
  my ($allele)=@_;
  
  my @alleles = $db->fetch(Variation => "$allele");

  return \@alleles;
}

=head2 get_all_alleles_fromFile

    Title   : get_all_alleles_fromFile
    Usage   : MapAlleles::get_alleles_fromFile($filename)
    Function: get alleles from acedb
    Returns : array_ref containing Ace::Allele
    Args    : (string) file name

=cut

sub get_alleles_fromFile {
  my ($filename)=@_;
  my @alleles;
  
  my $fh = new IO::File $filename ,'r';
  while (<$fh>){
    chomp;
    push @alleles, $db->fetch(Variation => "$_");
  }
    
  return \@alleles;
}


=head2 filter_alleles

    Title   : _filter_alleles
    Usage   : _filter_alleles(array_ref of Ace::Alleles)
    Function: removes allleles failing some sanity checks from the array PRIVATE
    Returns : array_ref containing Ace::Alleles
    Args    : (array_ref) of Ace::Alleles 

=cut

# filter the alleles
# removes alleles without 2 flanking sequences, no Sequence or not Source attached to the parent sequence
sub filter_alleles {
  my ($alleles, $aggressive) = @_;
  
  my @good_alleles;
  
  foreach my $allele (@{$alleles}) {
    my $name = $allele->Public_name || "No_public_name";
    my $remark = ($allele->Remark||'none');

    my ($status1, $status2, $live1, $live2);
    $status1 = $allele->Status;
    ($status2) = $allele->Status;
    $live1 = $allele->Live;
    ($live2) = $allele->Live;

    if (not defined $allele->Status or $allele->Status ne 'Live') {
      $log->write_to("WARNING: $allele ($name) is not Live, so skipping\n");
      next;
    }

    # has no sequence connection
    if ( ! defined $allele->Mapping_target ) {
      $log->write_to("ERROR: $allele ($name) has missing Mapping_target tag (Remark: $remark)\n");
      $errors++;
      next;
    }
    
    # The bit below is commented out, as there are sadly some alleles connected to sequences without source that can't be moved :-( updated this commented out code to avoid iues in the future.
    
    #   # connected sequence has no source
    #        elsif ( !defined $allele->Mapping_target->name && !defined $weak_checks) { 
    #            $log->write_to("ERROR: $name connects to ${\$allele->Mapping_target} which has no Source tag (Remark: $remark)\n");$errors++
    #        }

    # no flanks at all
    elsif (!defined $allele->Flanking_sequences) {
      $log->write_to("ERROR: $allele ($name) has no/empty Flanking_sequence (Remark: $remark)\n");
      $errors++;
      next;
    }
    # no left flanking sequence
    elsif (!defined $allele->Flanking_sequences->name ) {
      $log->write_to("ERROR: $allele ($name) has no left Flanking_sequence (Remark: $remark)\n");
      $errors++;
      next;
    }
    
    # no right flanking sequence
    elsif (!defined $allele->Flanking_sequences->right || ! defined $allele->Flanking_sequences->right->name ) {
      $log->write_to("ERROR: $allele ($name) has no right Flanking_sequence (Remark: $remark)\n");
      $errors++;
      next;
    }       
    
    elsif (defined $allele->Type_of_mutation and $allele->Type_of_mutation eq 'Substitution') {
      if ($aggressive) {
        if (not defined $allele->Type_of_mutation->right){
          $log->write_to("WARNING: $allele ($name) has no FROM in the substitution tag (Remark: $remark)\n");
          $errors++;
          next;
        } elsif (not defined $allele->Type_of_mutation->right->right) {
          $log->write_to("ERROR: $allele ($name) has no TO in the substitution tag (Remark: $remark)\n");
          $errors++;
          next;
        } elsif ($allele->Type_of_mutation->right =~ /\d|\s/ or $allele->Type_of_mutation->right->right =~ /\d|\s/) {
          $log->write_to("ERROR: $allele ($name) has numbers an/or spaces in the FROM/TO tags (Remark: $remark)\n");
          $errors++;
          next;
        }
      }
    } elsif(! (defined $allele->Type_of_mutation)){
      $log->write_to("WARNING: $allele ($name) has no Type_of_mutation (Remark: $remark)\n");
    }

    push @good_alleles, $allele; 
  }
  return \@good_alleles;
}

=head2 map

    Title   : map
    Usage   : MapAlleles::map(array_ref of Ace::Alleles)
    Function: maps alleles to chromosome coordinates
    Returns : hash_ref containing {allele_name}->{alleles}->Ace::Alleles and ->{start},->{stop},->{chromosome},->{clone},->{orientation}
    Args    : (array_ref) of Ace::Alleles 
    
=cut

# get the chromosome coordinates
sub map {
  my ($alleles, $fix_unmapped) = @_;
  my %alleles;
  my $mapper = Feature_mapper->new( $wb->autoace, undef, $wb );

  my $count = scalar(@$alleles);
  $log->write_to("INFO: Mapping $count alleles...\n");

  my (%mapped_pos, @unmapped_inserts, @unmapped_others);

  foreach my $x (@$alleles) {
    $count--;
    if ($count % 1000 == 0) {
      $log->write_to("INFO:   ...$count to go...\n");
    }

    # $chromosome_name,$start,$stop
    my ($min_len, $max_len, %mut_types, @mut_types);

    map { $mut_types{$_->name} = $_ } $x->at('Sequence_details.Type_of_mutation');

    if ($mut_types{Substitution}) {
      my $from = $mut_types{Substitution}->right;
      if ($from) {
        $min_len = $max_len = length($from);
      }
    }
    

    my $target_seq = $x->Mapping_target->name;
    my $left_flank = $x->Flanking_sequences->name;
    my $right_flank = $x->Flanking_sequences->right->name;

    my @map;
    if ($left_flank eq '<') {
      $left_flank = "";
    }
    if ($right_flank eq '>') {
      $right_flank = "";
    }
    if (length($left_flank) < 10 and
        length($right_flank) < 10) {
      $log->write_to("WARNING: $x has at least both flanks < 10 residues, so will not attempt to map with flanks\n");
      @map = (0);
    } else {
      @map = $mapper->map_feature($target_seq,
                                  $left_flank,
                                  $right_flank,
                                  $min_len, 
                                  $max_len);
    }
    
    if ($map[0] eq '0'){
      if ($mut_types{Insertion} and scalar(keys %mut_types) == 1) {
        push @unmapped_inserts, $x;
      } else {
        push @unmapped_others, $x;
      }
      next;
    }

    if ($map[0] ne $x->Mapping_target->name and 
        $x->Mapping_target->name !~ /^CHROMOSOME/ and 
        $x->Mapping_target->name !~ /^chr/){
      $log->write_to("WARNING: moved $x (${\$x->Public_name}) from sequence ${\$x->Mapping_target->name} to $map[0]\n");
    }

    my ($cl_st, $cl_en, $cl_ori) = ($map[1], $map[2], '+');

    if ($cl_st > $cl_en) {
      ($cl_st, $cl_en) = ($cl_en, $cl_st);
      $cl_ori = '-';
    }

    $mapped_pos{$x->name} = {
      parent => $map[0],
      start  => $cl_st,
      end    => $cl_en,
      strand => $cl_ori,
    };
  }

  # decide what to do about the unmapped ones
  if (@unmapped_inserts or @unmapped_others) {

    if (not $fix_unmapped) {
      foreach my $x (@unmapped_inserts, @unmapped_others) {
        $log->write_to("ERROR: Failed to map $x and will not attempt to remap (${\$x->Public_name}) " .
                       "to sequence ${\$x->Mapping_target->name} " .
                       "with ${\$x->Flanking_sequences->name} and ${\$x->Flanking_sequences->right->name} " .
                       "(Remark: ${\$x->Remark})\n");
        $errors++;
      }
    } else {
      my %remapped;

      my $assembly_mapper = Remap_Sequence_Change->new($wb->get_wormbase_version - 1, 
                                                       $wb->get_wormbase_version, 
                                                       $wb->species, 
                                                       $wb->genome_diffs);

      foreach my $pair ([\@unmapped_inserts, 1], [\@unmapped_others,0]) {
        my ($list, $zero_length) = @$pair;

        my @ids = map { $_->name } @$list;
        next if not @ids;

        my $hash = $mapper->remap_and_generate_new_flanks_for_features($assembly_mapper, 
                                                                       'Variation',
                                                                       $zero_length,
                                                                       50,
                                                                       \@ids);
        foreach my $k (keys %$hash) {
          $remapped{$k} = $hash->{$k};
        }
      }

      foreach my $x (@unmapped_inserts, @unmapped_others) {
        if (exists $remapped{$x->name}) {
          my ($chr, $left_flank, $right_flank) = ($remapped{$x->name}->{chr},
                                                  $remapped{$x->name}->{left_flank},
                                                  $remapped{$x->name}->{right_flank});

          $log->write_to("ERROR: Failed to map $x but was able to generate new flanks using coord remapping (see below)\n\n");
          $log->write_to("ACE : \n");
          $log->write_to("ACE : Variation : \"$x\"\n");
          $log->write_to("ACE : Sequence $chr\n");
          $log->write_to("ACE : Flanking_sequences $left_flank $right_flank\n");
          $log->write_to("ACE : \n");
        } else {
          $log->write_to("ERROR: Failed to remap $x after attempted fix (${\$x->Public_name}) " .
                         "to sequence ${\$x->Mapping_target->name} " .
                         "with ${\$x->Flanking_sequences->name} and ${\$x->Flanking_sequences->right->name} " .
                         "(Remark: ${\$x->Remark})\n");
          $errors++;
        }
      }
    }
  }
  
  foreach my $x (@$alleles) {
    next if not exists $mapped_pos{$x->name};

    my $clone_seq    = $mapped_pos{$x->name}->{parent};
    my $clone_start  = $mapped_pos{$x->name}->{start};
    my $clone_end    = $mapped_pos{$x->name}->{end};
    my $clone_strand = $mapped_pos{$x->name}->{strand};

    if( abs($clone_start - $clone_end) > 100000) {
      # Note from Mary Ann 20 june 2008:
      # I have heard back from Mark Edgley and he confirms
      # that the few Massive deletions which map_alleles
      # is spitting out as an error are in fact massive and
      # OK.
      # Since this dataset is pretty much complete I think
      # it would be OK to exclude niDf* alleles from the "is it massive"
      # check.
      # 
      # Note from klh 2013-09:
      # Edgely data set are natural variants, so now have a WBVar public name;
      # will now instead check for "CGH_allele" method. Also, some new variants from
      # the MM project are also massive, so also exclude them
      my $len = abs($clone_end - $clone_start) + 1;

      if ($x->Method ne 'CGH_allele' and
          (not $x->Analysis or $x->Analysis ne 'Million_mutation_project_reanalysis')) {
        $log->write_to("ERROR: $x (${\$x->Public_name}) is massive ($len)\n");
        $errors++;
        next;
      }
    } 

    # get chromosome coords, just in case we're not already in them
    my ($chromosome, $chr_start, $chr_end) = ($clone_strand eq '-') 
        ? $mapper->LocateSpanUp($clone_seq, $clone_end, $clone_start)
        : $mapper->LocateSpanUp($clone_seq, $clone_start, $clone_end);

    my $chr_strand = '+';
    if ($chr_start > $chr_end) {
      ($chr_start, $chr_end) = ($chr_end, $chr_start);
      $chr_strand = '-';
    }

    # from flanks to variation. If the retured coords defined a 2bp span, then
    # the flanked feature is 0bp. In this case, we leave the coords as they are
    # (i.e. including 1bp from each flank), because this this the only way in 
    # which Acedb/GFF allows us to define 0-bp features
    if ($clone_end - $clone_start != 1) {
      $clone_start++;
      $clone_end--;

      $chr_start++;
      $chr_end--;
    }

    $alleles{$x->name} = {
      allele      => $x, 
      chromosome  => $chromosome, 
      start       => $chr_start, 
      stop        => $chr_end, 
      orientation => $chr_strand,
      clone       => $clone_seq,
      clone_start => ($clone_strand eq '-') ? $clone_end : $clone_start, 
      clone_stop  => ($clone_strand eq '-') ? $clone_start : $clone_end,
    };
    
    print "${\$x->name} ($chromosome ($chr_strand): $chr_start - $chr_end) clone: $clone_seq $clone_start-$clone_end ($clone_strand)\n" if $wb->debug;
    
  }
  return \%alleles;
}


sub remove_insanely_mapped_alleles {
  my ($alleles) = @_;

  while (my ($k, $v) = each %$alleles) {
    if ($v->{allele}->Type_of_mutation eq 'Substitution') {
      if ($v->{allele}->Type_of_mutation->right) {
        my ($from) = $v->{allele}->Type_of_mutation->right;
        my $len_of_from = length($from);
        my $mapped_region_len = $v->{stop} - $v->{start} + 1;

        if (length($from) != $mapped_region_len) {
          $log->write_to("ERROR: $k (${\$v->{allele}->Public_name}) maps to a region (length $mapped_region_len) that does not match length of FROM ($from) - non-unique flanks?\n");
          $errors++;
          delete $alleles->{$k};
        }
      }
    }
  }
}

=head2 print_genes

    Title   : print_genes
    Usage   : MapAlleles::print_genes(hash_ref {$gene_name}=[$allele_name,..])
    Function: print ACE format for Allele->gene connections
    Returns : hash_ref containing {$allele_name}=[$gene_name,..]
    Args    : (hash_ref) of {$gene_name}=[$allele_name,..] 

=cut

# create ace formatted gene links
sub print_genes {
    my ($genes,$fh)=@_;
    my %all2genes;
    while( my($key,$allele_list)=each %$genes){
        foreach my $allele(@$allele_list){
            $all2genes{$allele}||=[];
            push @{$all2genes{$allele}},$key;
        }
    }

    while( my($allele,$gene_list) = each %all2genes){
        print $fh "Variation : \"$allele\"\n";
        foreach my $gene(@$gene_list){
            print $fh "Gene $gene\n";
        }
        print $fh "\n";

	foreach my $gene(@$gene_list){
	    print $fh "Gene : $gene\n",
	          "Allele $allele Inferred_automatically map_Alleles.pl\n\n";
	}
    }
    return \%all2genes;
}

=head2 get_genes

    Title   : get_genes
    Usage   : MapAlleles::get_genes(array_ref of Ace::Alleles)
    Function: maps alleles to genes
    Returns : hash_ref containing {$gene_name}=[$allele_name,..]
    Args    : (array_ref) of mapped Ace::Alleles (from map() )

=cut

# map alleles to genes (test)
sub get_genes {
    my ($alleles)=@_;
    my %genes;
    while(my($k,$v)=each(%{$alleles})){
        my @hits;

        @hits=$index->search_genes($v->{'chromosome'},$v->{'start'},$v->{'stop'});
        foreach my $hit(@hits){
            $genes{$hit->{name}}||=[];
            push @{$genes{$hit->{name}}}, $k;
        }
    }
    if ($wb->debug){
        foreach my $y (keys %genes) {print "$y -> ",join " ",@{$genes{$y}},"\n"}
    }
    return \%genes;
}


=head2 get_cds

    Title   : get_cds
    Usage   : MapAlleles::get_cds(array_ref of Ace::Alleles)
    Function: maps alleles to cds
    Returns : hash_ref containing {$cds_name}{$type}{$allele_name}=1
    Args    : (array_ref) of mapped Ace::Alleles (from map())

=cut

# map the alleles to cdses 
sub get_cds {
  my ($alleles)=@_;
  my %cds;
  while(my($k,$v)=each(%{$alleles})){
    my $dont_calc_protein_effect = 0;
    my $potential_frameshift = 0;

    # Fetch allele data up front.
    my $allele        = $v->{allele};
    my $a_public_name = $allele->Public_name;
    my $a_remark      = $allele->Remark;
    my $a_chr         = $v->{chromosome};
    my $a_start       = $v->{start};
    my $a_stop        = $v->{stop};
    my $a_len         = $a_stop - $a_start + 1;

    my ($substitution, $deletion, $insertion, %a_type, $a_type);
    foreach my $tp ($allele->Type_of_mutation) {
      my @list;
      if ($tp->right) {
        push @list, lc($tp->right);
        if ($tp->right->right) {
          push @list, lc($tp->right->right);
        }
      }
      if ($tp eq 'Substitution') {
        $a_type{Substitution} = 1;
        $substitution = \@list;
      } elsif ($tp eq 'Insertion') {
        $a_type{Insertion} = 1;
        $insertion = $list[0];
      } elsif ($tp eq 'Deletion') {
        $a_type{Deletion} = 1;
        $deletion = 1;
      }
    }
    my $a_type = join("/", sort keys %a_type);

    # we can assert some things up front before comparing to annotaiotn
    if (defined $substitution and (defined $insertion or defined $deletion)) {
      $log->write_to(sprintf("ERROR: %s - Substitution defined together with Insertion/Deletion, ".
                             "so will not calculate prot effect (%s ; Remark:%s)\n", 
                             $k,
                             $a_public_name,
                             $a_remark));

        $errors++;
        $dont_calc_protein_effect = 1;

    } elsif (defined $substitution and scalar(@$substitution) != 2) {
      $log->write_to(sprintf("ERROR: %s - substitution (%d) has missing/empty FROM and/or TO, ".
                             "so will not calculate prot effect (%s; Remark:%s)\n", 
                             $k,
                             $a_len, 
                             $a_public_name,
                             $a_remark));
      $errors++;
      $dont_calc_protein_effect = 1;
    } elsif (defined $substitution and ($substitution->[0] =~ /\d|\s/ or $substitution->[1] =~ /\d|\s/)) {
      $log->write_to(sprintf("ERROR: %s - small substitution (%d) has numbers/spaces in FROM/TO, ".
                             "so will not calculate prot effect (%s; Remark:%s)\n", 
                             $k,
                             $a_len,
                             $a_public_name,
                             $a_remark));
      $errors++;
      $dont_calc_protein_effect = 1;
    } 

    my $reported_from_to_error = 0;
    my @hits = $index->search_cds( $a_chr, $a_start, $a_stop );
    foreach my $hit(@hits){
      print $hit->{name},"\n" if $wb->debug;
      my @exons   = grep {($a_stop >= $_->{start}) && ($a_start <= $_->{stop} ) } $hit->get_all_exons;
      my @introns = grep {($a_stop >= $_->{start}) && ($a_start <= $_->{stop} ) } $hit->get_all_introns;

      #
      # introns
      #
      if (@introns){

        $cds{$hit->{name}}{Intron}{$k}=1;

        # 10 bp allele overlapping a splice site
        if (abs($a_start - $a_stop) <= 10){

          my @types=('Acceptor','Donor');
          @types=('Donor','Acceptor') if $hit->{orientation} eq '+';
          my @intron_starts = grep { $a_start <= $_->{start}+1 && $a_stop >= $_->{start} } @introns;
          my @intron_stops  = grep { $a_start <= $_->{stop}    && $a_stop >= $_->{stop}-1 } @introns;
          
          # add the from to stuff:
          foreach my $intron (@intron_starts){
            my $site = &get_seq($intron->{chromosome},$intron->{start},$intron->{start}+1,$intron->{orientation});
            if ($a_start - $a_stop == 0 and defined $substitution) {
              my ($from_na, $to_na) = @$substitution;

              my $offset=$intron->{start} - $a_start;
              $offset=(1-$offset) if $intron->{orientation} eq '-';
              my $to_site   = $site;
              my $from_site = $site;
              substr($to_site,$offset,1,lc($to_na));
              substr($from_site,$offset,1,lc($from_na));
              my ($a,$b) = ($from_site eq $site) ? ($from_site,$to_site) : ($to_site,$from_site);
              $cds{$hit->{name}}{"$types[0] \"$a to $b\""}{$k}=1;
            }
            else {
              $cds{$hit->{name}}{"$types[0] \"$a_type disrupts $site\""}{$k}=1; 
            }
          }
          foreach my $intron(@intron_stops){
            my $site = &get_seq($intron->{chromosome},$intron->{stop}-1,$intron->{stop},$intron->{orientation});

            if ($a_start - $a_stop == 0 and defined $substitution) {
               my ($from_na, $to_na) = @$substitution;

               my $offset = 1-($intron->{stop} - $a_stop);
               $offset=(1-$offset) if $intron->{orientation} eq '-';
               my $to_site = $site;
               my $from_site = $site;
               substr($to_site,$offset,1,lc($to_na));
               substr($site,$offset,1,lc($from_na));
               my ($a,$b)=($from_site eq $site)?($from_site,$to_site):($to_site,$from_site);          
               $cds{$hit->{name}}{"$types[1] \"$a to $b\""}{$k}=1;
            }
            else {            
              $cds{$hit->{name}}{"$types[1] \"$a_type disrupts $site\""}{$k}=1;
            }
          }
        }
      }

      #
      # exons
      #
      if (@exons){
        $cds{$hit->{name}}{Coding_exon}{$k}=1;

        # if allele spans more than one exon, too complicated
        next if scalar(@exons) > 1;

        # if allele extends beyond the end of the exon, too complicated
        next if $a_start < $exons[0]->{start} or $a_stop > $exons[0]->{stop};

        next if $dont_calc_protein_effect;

        # get the cdna sequence of the CDS, and calc position of the allele within that CDS
        my @coding_exons;
        if ($hit->{orientation} eq '+'){ 
          @coding_exons = sort { $a->{start} <=> $b->{start} } $hit->get_all_exons;
          if ($coding_exons[0]->{phase} != 0) {
            $coding_exons[0]->{start} += $coding_exons[0]->{phase};
            $coding_exons[0]->{phase} = 0;
          }
        } else { 
          @coding_exons = sort { $b->{start} <=> $a->{start} } $hit->get_all_exons;
          if ($coding_exons[0]->{phase} != 0) {
            $coding_exons[0]->{stop} -= $coding_exons[0]->{phase};
            $coding_exons[0]->{phase} = 0;
          }
        }
        
        my $sequence = "";
        foreach my $ex (@coding_exons) {
          $sequence .= &get_seq($ex->{chromosome},$ex->{start},$ex->{stop},$ex->{orientation});
        }
        
        my $cds_start=0;
        my $exon=1;
        foreach my $c_exon(@coding_exons){
          if ($a_start <= $c_exon->{stop} and $a_stop >= $c_exon->{start}){
            # if allele extends beyond CDS limits, we will go no further
            if ($a_start < $c_exon->{start} or $a_stop > $c_exon->{stop}) {
              $cds_start = -1;
              last;
            }
            
            if ($c_exon->{orientation} eq '+') {
              $cds_start += ($a_start - $c_exon->{start} + 1);
            } else {
              $cds_start += ($c_exon->{stop} - $a_stop + 1); 
            }
            last;
          } else{
            $cds_start += ($c_exon->{stop} - $c_exon->{start} + 1);
          }
        }
        
        next if $cds_start < 0;
        my $cds_end = $cds_start + $a_len - 1;
        
        my $codon_start = $cds_start - (($cds_start-1) % 3);
        my $codon_end   = $cds_end + (3 - ($cds_end % 3)) % 3;
        
        #
        # Deal with frameshifts. 
        #
        if (defined $substitution and (length($substitution->[1]) - length($substitution->[0])) % 3 != 0) {
          $cds{$hit->{name}}{"Frameshift \"Unbalanced substitution at CDS position $cds_start-$cds_end\""}{$k} = 1;
          next;
        } elsif (defined $insertion and not defined $deletion and
                 (length($insertion) % 3 != 0)) {
          my $str = sprintf("%d bp insertion at CDS position %d-%d", length($insertion), $cds_start, $cds_end); 
          $cds{$hit->{name}}{"Frameshift \"$str\""}{$k} = 1;
        } elsif (defined $deletion and not defined $insertion and
                 $a_len % 3 != 0) {
          my $str = sprintf("%d bp deletion at CDS position %d-%d", $a_len, $cds_start, $cds_end);
          $cds{$hit->{name}}{"Frameshift \"$str\""}{$k} = 1;
          next;
        } elsif (defined $deletion and defined $insertion and 
                 (length($insertion) - $a_len) % 3 != 0) {
          my $str = sprintf("Unbalanced insertion/deletion at CDS position %d-%d", $cds_start, $cds_end);
          $cds{$hit->{name}}{"Frameshift \"$str\""}{$k} = 1;
          next;
        }

        if (defined $substitution or 
            defined $insertion and defined $deletion and $a_len == length($insertion)) {

          my $aa_start = ($codon_start - 1) / 3 + 1;
          my $aa_end   = $codon_end / 3;
          my $aa_range = ($aa_start == $aa_end) ? $aa_start : "$aa_start-$aa_end";

          my $table = Bio::Tools::CodonTable->new(-id=>1);
          
          my $from_na     = uc(substr($sequence, $cds_start - 1, $cds_end - $cds_start + 1));
          my $rev_from_na = Bio::Seq->new(-seq => $from_na, -alphabet => 'dna')->revcom->seq;

          my $to_na;

          if (defined $substitution) {
            my $sub_from_na = uc($substitution->[0]); 
            my $sub_to_na   = uc($substitution->[1]); 
            $to_na = $sub_to_na;

            my $rev_sub_from_na =  Bio::Seq->new(-seq => $sub_from_na, -alphabet => 'dna')->revcom->seq;
            my $rev_sub_to_na   =  Bio::Seq->new(-seq => $sub_to_na,   -alphabet => 'dna')->revcom->seq;
                      
            my $genome_from_na = ($hit->{orientation} eq '-') ? $rev_from_na : $from_na;
    
            if ($from_na ne $sub_from_na) {
              if ($hit->{orientation} eq '-' and
                $from_na eq $rev_sub_from_na) {
                # only accept this if the transcript is reverse strand, 
                # and the from matches the forward strand. In that case, 
                # we assert that FROM/TO have been defined in genomic sense 
                # (rather than transcript sense) and revcomp the from/to
                $to_na = $rev_sub_to_na;
              } else {
                if (not $reported_from_to_error) {
                  $log->write_to("ERROR: $k FROM/TO tags seem to be messed up; " . 
                                 "annotated FROM ($sub_from_na) does not match transcript FROM ($from_na) ".
                                 "or genome FROM ($genome_from_na)\n");
                  $reported_from_to_error = 1;
                  $errors++;
                }
                next;
              }
            }
          } elsif (defined $insertion) {
            # insertions should *always* be defined in forward genomic sense
            $to_na = $insertion;
            if ($hit->{orientation} eq '-') {
              $to_na =  Bio::Seq->new(-seq => $to_na, -alphabet => 'dna')->revcom->seq;
            }
          }

          my $to_sequence = $sequence;

          substr($to_sequence, $cds_start - 1, $cds_end - $cds_start + 1, $to_na);

          my $from_codons = substr($sequence, $codon_start -1, $codon_end - $codon_start + 1);
          my $to_codons   = substr($to_sequence, $codon_start -1, $codon_end - $codon_start + 1);

          my $from_aa = $table->translate($from_codons);
          my $to_aa =  $table->translate($to_codons);

          if ($to_aa eq $from_aa){
            # silent mutation
            $cds{$hit->{name}}{"Silent \"$to_aa ($aa_range)\""}{$k}=1;
            print "silent mutation: " if $wb->debug;
          } elsif (grep { /\*/ } split(//, $from_aa) and 
                   not grep { /\*/ } split(//, $to_aa)) {
            # readthrough; 
            $cds{$hit->{name}}{"Readthrough \"$from_aa to $to_aa\""}{$k} = 1;
          }
          # premature stop
          elsif (grep { /\*/ } split(//, $to_aa) and 
                 not grep { /\*/ } split(//, $from_aa)) {
            my $stop_aa_offset = index($to_aa, '*');

            my $stop_codon = substr($to_codons, $stop_aa_offset * 3, 3); 
            my $from_codon = substr($from_codons, $stop_aa_offset * 3, 3); 
            my $other_aa   = substr($from_aa, $stop_aa_offset, 1);

            my $stop_aa_pos = $aa_start + $stop_aa_offset;


            if (uc($stop_codon) =~ /[TWYKHDB]AG/ ){
              $cds{$hit->{name}}{"Nonsense Amber_UAG \"$other_aa to amber stop ($stop_aa_pos)\""}{$k}=1;
            }
            elsif (uc($stop_codon) =~ /[TWYKHDB]AA/){
              $cds{$hit->{name}}{"Nonsense Ochre_UAA \"$other_aa to ochre stop ($stop_aa_pos)\""}{$k}=1;
            }
            elsif (uc($stop_codon) =~ /[TWYKHDB]GA/){
              $cds{$hit->{name}}{"Nonsense Opal_UGA \"$other_aa to opal stop ($stop_aa_pos)\""}{$k}=1;
            }
            elsif (uc($stop_codon) eq 'TAR'){
              $cds{$hit->{name}}{"Nonsense Amber_UAG_or_Ochre_UAA \"$other_aa to amber or ochre stop ($stop_aa_pos)\""}{$k}=1;
            }
            elsif (uc($stop_codon) eq 'TRA') {
              $cds{$hit->{name}}{"Nonsense Ochre_UAA_or_Opal_UGA \"$other_aa to opal or ochre stop ($stop_aa_pos)\""}{$k}=1;
            }
            else {
              $log->write_to("ERROR: %k whatever stop $stop_codon is in $k ($a_public_name) it is not Amber/Opal/Ochre (Remark: $a_remark})\n");
              $errors++;
            }
          } else{
            $cds{$hit->{name}}{"Missense $aa_start \"$from_aa to $to_aa\""}{$k}=1;
            print "Missense: " if $wb->debug;
          }
        }
      }

    }
  }
  
  if ($wb->debug){
    foreach my $cds_name (keys %cds) {
      foreach my $type(keys %{$cds{$cds_name}}){
        print "$cds_name -> ",join " ",keys %{$cds{$cds_name}{$type}}," ($type)\n"
      }
    }
  }
  return \%cds;
}   

=head2 print_cds

    Title   : print_cds
    Usage   : MapAlleles::print_cds(hash_ref {$cds_name}{$type}{$allele_name}=1)
    Function: print ACE format for Allele->cds connections
    Returns : hash_ref containing {$allele_name}{$type}{$gene_name}=1
    Args    : (hash_ref) of {$cds _name}{$type}{$allele_name}=1 

=cut

# create ace formatted gene links
sub print_cds {
    my ($cds,$fh)=@_;
    my %inversecds;
    while(my($cds,$type_allele)=each %$cds){
         while(my($type,$allele)=each %$type_allele){
            foreach my $allele_name(keys %$allele){
                $inversecds{$allele_name}{$type}||=[];
                push @{$inversecds{$allele_name}{$type}},$cds;
            }
        }
    }

    while( my($allele_name, $type_cds) = each %inversecds){
        print $fh 'Variation : "',$allele_name,"\"\n";
        foreach my $type (keys %$type_cds){
            foreach my $cds(@{$type_cds->{$type}}){
                print $fh "Predicted_CDS $cds $type Inferred_automatically map_Alleles.pl\n";
                #print "Predicted_CDS $cds $type Inferred_automatically map_Alleles.pl\n" if $wb->debug;
            }
        }
        print $fh "\n";
    }
    return \%inversecds;
}

=head2 compare

    Title   : compare
    Usage   : MapAlleles::compare(hash_ref of mapped Ace::Allele, hash_ref of $gene->[$allele,..])
    Function: compare old and new gene connections
    Returns : 1 / writes errors to log
    Args    : hash_ref of mapped Ace::Allele, hash_ref of $gene->[$allele,..] 

=cut

# compare (old(ace objects) - new(hash))
sub compare { 
    my ($old,$new)=@_;
    my %check;
    foreach my $allele(keys %$old){
        foreach my $gene($old->{$allele}->{allele}->Gene){
            next if (qq/${\$old->{$allele}->{allele}->at("Affects.Gene.$gene")->col(1)}/ eq 'Genomic_neighbourhood');
	    next if (qq/${\$old->{$allele}->{allele}->at("Affects.Gene.$gene")->col(1)}/ eq 'Regulatory_feature');

	    $check{$allele}{$gene}=1; 
        }
    }
    foreach my $allele(keys %$new){
        foreach my $gene(@{$new->{$allele}}){
            $check{$allele}{$gene}+=2
        }
    }
    while(my($allele,$v)=each %check){
        while(my ($gene,$y)=each %$v){           
            if ($y==1) {
                my $remark=$old->{$allele}->{allele}->Remark;
                $log->write_to("ERROR: $allele (${\$old->{$allele}->{allele}->Public_name}) -> $gene connection is only in geneace (Remark: $remark)\n");$errors++}
            elsif($y==2){
                #$log->write_to("WARNING: $allele -> $gene connection created by script\n");
            }
            elsif($y==3){}
            else{die "comparison failed\n"}
        }
    }
    1;
}

=head2 get_seq

    Title   : get_seq
    Usage   : MapAlleles::get_seq(chromosome, start, stop, orientation)
    Function: get subsequences of the genome
    Returns : (string) sequence
    Args    : ('I','II',...) chromosome, (integer) start, (integer) stop, ('+','-') orientation 

=cut

# get sequence from chromosome through ACeDB
# get_seq('I',1,100,'+')
sub get_seq {
    my ($chromosome,$pos1,$pos2,$orientation)=@_; 
    
    my $seq=$extractor->Sub_sequence($chromosome);
    my $len=length $seq;
    
    $pos1--;
    $pos2--;
    
    my ($p1,$p2)=sort {$a<=>$b}($pos1,$pos2);
    
    my $raw_seq=substr($seq,$p1,$p2-$p1+1);
    $raw_seq= $extractor->DNA_revcomp($raw_seq) if ($orientation eq '-');
    return $raw_seq;
}

=head2 load_genes_and_cds

    Title   : _build_index
    Usage   : MapAlleles::build_index()
    Function: builds the Gene and Cds index
    Returns : 1
    Args    : none 

=cut

# iterates over the GFF_SPLIT and adds them to the index
sub load_genes_and_cds {

  my $gdir = $wb->gff_splits;

  my ($exon, $intron);

  foreach my $file (glob("$gdir/*_gene.gff"), 
                    glob("$gdir/*_curated.gff"), 
                    "$gdir/gene.gff", 
                    "$gdir/curated.gff") {
    next if not -e $file;

    my $fh = new IO::File $file, 'r';

    while(<$fh>) {
      chomp;
      s/\"//g;
      my @F = split;
      
      if ($F[1] eq 'gene' and $F[2] eq 'gene') {
        Store::store_gene( $F[0], $F[3], $F[4], $F[9], $F[6] );
      }
      if ( $F[1] eq 'curated' and  $F[2] eq 'CDS' ) {
        Store::store_cds( $F[0], $F[3], $F[4], $F[6], $F[9],undef );
      }        
      if ( $F[1] eq 'curated' and $F[2] eq 'coding_exon' ) {
        Store::store_exon( $F[0], $F[3], $F[4], $F[6], $exon++, $F[9], $F[7] );
      }
      if( $F[1] eq 'curated' and $F[2] eq 'intron' ) { 
        Store::store_intron( $F[0], $F[3], $F[4], $F[6], $intron++, $F[9], $F[7] );
      }
    }

    $fh->close;
  }

  $index=Mem_index::build;
}
                          
=head2 load_utr

        Title   : load_utr
        Usage   : MapAlleles::load_utr()
        Function: load the *UTR.gff files
        Returns : hash->chromosome->UTRs
        Args    : none

=cut

# load UTRs
sub load_utr{
  my $gff_dir = $wb->gff_splits;

  my @files = (glob("$gff_dir/*_UTR.gff"), "$gff_dir/UTR.gff");
  my %utrs;
  foreach my $file(@files){
    next if not -e $file;

    my $inf=new IO::File $file, 'r';
    #print "processing: $file\n" if $wb->debug;
    while (<$inf>) {
      next if /\#/;
      next if ! /UTR/;
      s/\"//g;
      my @fields=split;
      my ($chromosome,$start,$stop,$type,$transcript)=($fields[0],$fields[3],$fields[4],$fields[2],$fields[-1]);
      my $field={
        chromosome => $chromosome,
        start      => $start,
        stop       => $stop,
        type       => $type,
        transcript => $transcript,
      };
      $utrs{$chromosome}||=[];
      push @{$utrs{$chromosome}}, $field;
    }
  }
  return \%utrs;
}
             
=head2 search_utr

        Title   : search_utr
        Usage   : MapAlleles::search_utr(alleles,hashref from load_utr)
        Function: map alleles to UTRs
        Returns : hash->allele->transcript->type
        Args    : list of alleles,hashref of UTRs

=cut

# search UTRs
sub search_utr{                          
    my ($alleles,$utrs)=@_;
    my %allele_utr;
    while(my($k,$v)=each(%{$alleles})){
        my @hits = grep {$_->{start}<=$v->{stop} && $_->{stop}>=$v->{start}} @{$$utrs{$v->{chromosome}}};
        foreach my $hit(@hits){
            $allele_utr{$k}{$hit->{transcript}}{$hit->{type}}=1;
            print "$k -> ${\$hit->{transcript}} (${\$hit->{type}})\n" if $wb->debug;
        }
    }
    return \%allele_utr;
}

=head2 print_utr

        Title   : print_utr
        Usage   : MapAlleles::print_utr(hashref of allele->utr hits,filehandle)
        Function: print Allele->UTR connections into ace format
        Returns : nothing
        Args    : hashref of allele->utr hits,filehandle

=cut

# print UTRs
sub print_utr{
    my ($hits,$fh)=@_;
    while( my($allele_name, $type_cds) = each %$hits){
        print $fh 'Variation : "',$allele_name,"\"\n";
        foreach my $cds (keys %$type_cds){
            foreach my $type(keys %{$type_cds->{$cds}}){
                my $ace_type=$type eq 'three_prime_UTR'?'UTR_3':'UTR_5';
                print $fh "Transcript $cds $ace_type Inferred_automatically map_Alleles.pl\n";
            }
        }
        print $fh "\n";
    }
}

# load pseudogenes
sub load_pseudogenes{
  my $gff_dir = $wb->gff_splits;
  my @files = (glob("$gff_dir/*_Pseudogene.gff"), "$gff_dir/Pseudogene.gff");
  
  my %pgenes;
  foreach my $file(@files){
    next if not -e $file;
    my $inf=new IO::File $file, 'r';
    #print "processing: $file\n" if $wb->debug;
    while (<$inf>) {
      next if /\#/;
      next unless /Pseudogene/;
      s/\"//g;
      my @fields=split;
      my ($chromosome,$start,$stop,$type1,$type2,$gene)=($fields[0],$fields[3],$fields[4],$fields[1],$fields[2],$fields[-1]);
      next unless ($type1 eq 'Pseudogene' && $type2 eq 'Pseudogene');
      my $field={
        chromosome => $chromosome,
        start      => $start,
        stop       => $stop,
        pgene      => $gene,
      };
      $pgenes{$chromosome}||=[];
      push @{$pgenes{$chromosome}}, $field;
    }
  }
  return \%pgenes;
}
             
=head2 search_pseudogene

        Title   : search_pseudogene
        Usage   : MapAlleles::search_pseudogene(alleles,hashref from load_pseudogene)
        Function: map alleles to pseudogeness
        Returns : hash->allele->transcript->type
        Args    : list of alleles,hashref of pseudogenes

=cut

# search pseudogenes
sub search_pseudogenes{                          
    my ($alleles,$pgenes)=@_;
    my %allele_pgenes;
    while(my($k,$v)=each(%{$alleles})){
        my @hits = grep {$_->{start}<=$v->{stop} && $_->{stop}>=$v->{start}} @{$$pgenes{$v->{chromosome}}};
        foreach my $hit(@hits){
            $allele_pgenes{$k}{$hit->{pgene}}=1;
            print "$k -> ${\$hit->{pgene}}\n" if $wb->debug;
        }
    }
    return \%allele_pgenes;
}

=head2 print_pseudogene

        Title   : print_pseudogene
        Usage   : MapAlleles::print_pseudogene(hashref of allele->utr hits,filehandle)
        Function: print Allele->Pseudogene connections into ace format
        Returns : nothing
        Args    : hashref of allele->pseudogene hits,filehandle

=cut

# print Pseudogene
sub print_pseudogenes{
    my ($hits,$fh)=@_;
    while( my($allele_name, $pgenes) = each %$hits){
        print $fh 'Variation : "',$allele_name,"\"\n";
        foreach my $p_gene (keys %$pgenes){
                print $fh "Pseudogene $p_gene\n";
        }
        print $fh "\n";
    }
}


# search through non_coding Transcripts

# load ncrna
sub load_ncrnas{

  my @files;
  my $gff_dir = $wb->gff_splits;

  foreach my $type ('miRNA', 'pre_miRNA', 'ncRNA','rRNA','scRNA','snoRNA','snRNA','snlRNA','stRNA','tRNA','piRNA','asRNA','lincRNA') {
    push @files, (glob("$gff_dir/*_${type}.gff"), "$gff_dir/${type}.gff");
  }
  my %nc_rnas;
  foreach my $file(@files){
    next if not -e $file;

    my $inf=new IO::File $file, 'r';
    while (<$inf>) {
      next if /\#/;

      s/\"//g;
      my @fields=split;
      next if $fields[2] ne 'exon';

      my ($chromosome,$start,$stop,$transcript)=($fields[0],$fields[3],$fields[4],$fields[-1]);
      my $field= {
        chromosome => $chromosome,
        start      => $start,
        stop       => $stop,
        transcript => $transcript,
      };
      $nc_rnas{$chromosome}||=[];
      push @{$nc_rnas{$chromosome}}, $field;
    }
  }

  return \%nc_rnas;
}
             
=head2 search_ncrna

        Title   : search_ncrnas
        Usage   : MapAlleles::search_ncrnas(alleles,hashref from load_ncrna)
        Function: map alleles to ncrnas
        Returns : hash->allele->transcript->1
        Args    : list of alleles,hashref of ncrna

=cut

# search ncrnas
sub search_ncrnas{                          
    my ($alleles,$ncrnas)=@_;
    my %allele_ncrnas;
    while(my($k,$v)=each(%{$alleles})){
        my @hits = grep {$_->{start}<=$v->{stop} && $_->{stop}>=$v->{start}} @{$ncrnas->{$v->{chromosome}}};
        foreach my $hit(@hits){
            $allele_ncrnas{$k}{$hit->{transcript}}=1;
            print "$k -> ${\$hit->{transcript}}\n" if $wb->debug;
        }
    }
    return \%allele_ncrnas;
}

=head2 print_ncrna

        Title   : print_ncrnas
        Usage   : MapAlleles::print_ncrnas(hashref of allele->ncrna hits,filehandle)
        Function: print Allele->ncrna connections into ace format
        Returns : nothing
        Args    : hashref of allele->ncrna hits,filehandle

=cut

# print non coding RNAs
sub print_ncrnas{
  my ($hits,$fh)=@_;
  while( my($allele_name, $ncrnas) = each %$hits){
    print $fh 'Variation : "',$allele_name,"\"\n";
    foreach my $ncrna (keys %$ncrnas){
      print $fh "Transcript $ncrna\n";
    }
    print $fh "\n";
  }
}



1;

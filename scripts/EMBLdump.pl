#!/usr/bin/env perl
#
# EMBLdump.pl 
# 
#  Last updated on: $Date: 2015-06-17 11:45:03 $
#  Last updated by: $Author: klh $

use strict;
use Getopt::Long;
use File::Copy;
use Storable;
use Text::Wrap;

$Text::Wrap::columns = 60;
$Text::Wrap::unexpand = 0;

use lib $ENV{'CVS_DIR'};
use Bio::SeqIO;
use Wormbase;
use Coords_converter;

my %species_info = (

  strain => {
    elegans => 'Bristol N2',
    briggsae => 'AF16',
    brugia   => 'FR3',
    sratti   => 'ED321',
  },

);    

my ($test,
    $debug,
    $store,
    $species,
    $dump_raw,
    $dump_modified,
    $stage_to_repository,
    $full_species_name,
    $raw_dump_file,
    $mod_dump_file,
    $wormbase,
    $single,
    $database,
    $quicktest,
    $private,
    $sequencelevel,
    $agp_file,
    $clone2type, $gene2cgcname, $gene2gseqname, $trans2gene, $clone2dbid, $multi_gene_loci, $seleno_prots,
    $cds2proteinid_db, $cds2status_db, $trans2dbremark_db, $trans2briefid_db, $decorations_db,
    $cds2status_h, $cds2proteinid_h,  $trans2dbremark_h, $trans2type_h, $trans2briefid_h,
    $agp_segs,
    %additional_qualifiers,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "single=s"        => \$single,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "dumpraw"         => \$dump_raw,
  "dumpmodified"    => \$dump_modified,
  "stagetorepo"     => \$stage_to_repository,
  "decorationsdb=s" => \$decorations_db,
  "piddb=s"         => \$cds2proteinid_db,
  "dbremarkdb=s"    => \$trans2dbremark_db,
  "cdsstatusdb=s"   => \$cds2status_db,
  "briefiddb=s"     => \$trans2briefid_db,
  "modprivate=s"    => \$private,
  "moddumpfile:s"   => \$mod_dump_file,
  "rawdumpfile:s"   => \$raw_dump_file,
  "quicktest"       => \$quicktest,
  "sequencelevel"   => \$sequencelevel,
  "agpfile=s"       => \$agp_file,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
                             -autoace  => $database,
      );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
$full_species_name = $wormbase->full_name;
my $dbdir = ($database) ? $database : $wormbase->autoace;
my $tace = $wormbase->tace;
my $bioproject = $wormbase->ncbi_bioproject; 

$decorations_db =  $wormbase->database('current') if not defined $decorations_db;
$cds2proteinid_db = $decorations_db if not defined $cds2proteinid_db;
$cds2status_db = $decorations_db if not defined $cds2status_db;
$trans2dbremark_db = $decorations_db if not defined $trans2dbremark_db;
$trans2briefid_db = $decorations_db if not defined $trans2briefid_db;

###############################
# misc. variables             #
###############################


my ($delete_raw_dump_file, $delete_mod_dump_file) = (0,0);

if ($dump_raw) {
  #############################################
  # Use giface to dump EMBL files from camace #
  #############################################

  my $giface = $wormbase->giface;

  if (not defined $raw_dump_file) {    
    $raw_dump_file = "/tmp/EMBLdump.$$";
    $log->write_to("You are raw embl dumping from $dbdir to $raw_dump_file which will be delted\n\n");
    $delete_raw_dump_file = 1;
  } else {
    $log->write_to("You are raw embl dumping from $dbdir to $raw_dump_file which will be retained\n\n");
  }
  
  my $command;
  if ($single) {
    $command .= "query find Sequence $single\ngif EMBL $raw_dump_file\n";# find sequence and dump
    $command .= "quit\n";# say you don't want to save and exit
  } else {
    $command .= "query find Sequence EMBL_dump_info AND Species=\"$full_species_name\"";
    if ($species eq 'elegans') {
      if ($sequencelevel) {
        $command .= " AND Source";
      } else {
        $command .= " AND NOT Source"
      }
    }
    $command .= "\ngif EMBL $raw_dump_file\n";# find sequence and dump
    $command .= "quit\n";# say you don't want to save and exit
  }
  
  $log->write_to("$command\n");
  open (READ, "echo '$command' | $giface $dbdir |") or die ("Could not open $giface $dbdir\n");
  while (<READ>) {
    next if ($_ =~ /\/\//);
    next if ($_ =~ /acedb/);
  }
  close (READ);

  $raw_dump_file .= ".embl";
  $log->write_to("RAW dumped to $raw_dump_file\n");
} 


if ($dump_modified) {
  ######################################################################
  # cycle through the EMBL dump file, replacing info where appropriate #
  ######################################################################
  if (not $raw_dump_file) {
    $log->log_and_die("You cannot produce a modified dump file without producing a raw one first" .
                      " (either with -dumpraw or -rawdumpfile <file>\n");
  }

  if ($agp_file) {
    $agp_segs = &read_agp_files( $agp_file );
  }

  $trans2gene = &fetch_transcript2gene($dbdir);

  $clone2type = &fetch_clone2type($dbdir);
  $clone2dbid = &fetch_clone2dbid($dbdir);
  ($gene2cgcname, $gene2gseqname) = &fetch_gene2names($dbdir);

  if (not $sequencelevel) {
    $multi_gene_loci = &get_multi_gene_loci($gene2gseqname);
    
    $trans2type_h = &fetch_transcript2type($dbdir);
    
    # Some information is not yet available in autoace (too early in the build)
    # By default, pull cds2status, trans2dbremark and protein_ids from current_DB,
    # however, allow command-line override when necessary
    $cds2status_h = &fetch_cds2status($cds2status_db);
    $cds2proteinid_h = &fetch_cds2proteinid($cds2proteinid_db);
    $trans2dbremark_h = &fetch_trans2dbremark($trans2dbremark_db);
    $trans2briefid_h = &fetch_transcript2briefid($trans2briefid_db);
    
    $seleno_prots = &fetch_selenoproteins($dbdir);
  }

  if (not defined $mod_dump_file) {
    $mod_dump_file = "/tmp/EMBL_dump.$$.mod.embl";
    $log->write_to("You are MODIFIED embl dumping to $mod_dump_file which will be deleted\n");
    $delete_mod_dump_file = 1;
  } else {
    $log->write_to("You are MODIFIED embl dumping to $mod_dump_file which will be retained\n");
  }

  open(my $out_fh, ">$mod_dump_file") or $log->log_and_die("Could not open $mod_dump_file for writing\n");
  open(my $raw_fh, $raw_dump_file) or $log->log_and_die("Could not open $raw_dump_file for reading\n");
  
  my ($seqname, $chromosome, $seqlen, $idline_suffix, @accs, @features, $written_header, %written_version_tag);
  
  while (<$raw_fh>) {
    
    # Store the necessary default ID line elements ready for use in the new style EMBL ID lines.
    if(/^ID\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\d+)\s+\S+/){
      ($seqname, $seqlen) = ($1, $2);
      $idline_suffix = sprintf("SV XXX; linear; genomic DNA; %s; INV; $seqlen BP.", ($sequencelevel) ? "STD" : "CON");;
      @accs = ();
      @features = ();
      $written_header = 0;
      
      next;
    }
    
    if( /^AC\s+(\S+);/ ) {
      push @accs, $1;
      next;
    }
    
    if (/^DE/) {
      # should now have parsed everything necessary to write first block of entry
      #
      # ID
      #
      print $out_fh "ID   $accs[0]; $idline_suffix\n";
      print $out_fh "XX\n";
      
      #
      # ST * private
      #
      if ($private) {
        print $out_fh "ST * $private\n";
        print $out_fh "XX\n";
      }
      
      #
      # AC
      #
      print $out_fh "AC   $accs[0];";
          
      for(my $i=1; $i < @accs; $i++) {
        print $out_fh " $accs[$i];";
      }
      print $out_fh "\nXX\n";
      
      #
      # AC *
      #
      if ($clone2dbid->{$seqname}) {
        print $out_fh "AC * $clone2dbid->{$seqname}\n";
        print $out_fh "XX\n";
      }
      
      #
      # PR
      #
      print $out_fh "PR   Project:$bioproject\n";
      print $out_fh "XX\n";
      
      #
      # DE
      #
      my $de_line;
      
      if ($species eq 'elegans') {
        if ($seqname =~ /CHROMOSOME_(\S+)/) {
          $chromosome = $1;
          $de_line = "$full_species_name chromosome $chromosome";
        } elsif (!defined($clone2type->{$seqname})){
          $de_line =  "$full_species_name clone $seqname";
        } elsif (lc($clone2type->{$seqname}) eq "other") {
          $de_line = "$full_species_name clone $seqname";
        } elsif (lc($clone2type->{$seqname}) eq "yac") {
          $de_line = "$full_species_name YAC $seqname";
        } else {
          $de_line = "$full_species_name $clone2type->{$seqname} $seqname";
        }
      } elsif ($species eq 'briggsae') {
        $de_line = "$full_species_name AF16 supercontig from assembly CB4, $seqname";
      } elsif ($species eq 'brugia'){
        $de_line = "$full_species_name FR3 supercontig from assembly B_malayi-3.1 $seqname"
      }
      
      print $out_fh "DE   $de_line\n";
      $written_header = 1;
      next;
    }
    
    #
    # References
    #
    if (/^RN\s+\[1\]/) {
      my ($primary_RA, $primary_RL); 
      
      while(<$raw_fh>) {
        last if /^XX/;
        if (/^RA/) {
          $primary_RA = $_;
          next;
        }
        if (/^RL/) {
          $primary_RL = $_;
        }
      }

      my $ref_count = 1;
      
      if ($species eq 'elegans') {
        # elegans is a special case, because every clone has a different primary
        # author.
        printf $out_fh "RN   [%d]\n", $ref_count++;
        printf $out_fh "RP   1-%d\n", $seqlen;
        print $out_fh "RG   WormBase Consortium\n";
        if (defined $primary_RA) {
          print $out_fh $primary_RA;
        } else {
          print $out_fh "RA   Howe, K.L.;\n";
        }
        print $out_fh "RT   ;\n";
        print $out_fh $primary_RL if defined $primary_RL;
        print $out_fh "RL   Nematode Sequencing Project: Sanger Institute, Hinxton, Cambridge\n";
        print $out_fh "RL   CB10 1SA, UK and The Genome Institute at Washington University,\n"; 
        print $out_fh "RL   St. Louis, MO 63110, USA. E-mail: help\@wormbase.org\n";
        print $out_fh "XX\n";
      }

      my @refs = @{&get_references()};
      for(my $i=0; $i < @refs; $i++) {
        printf $out_fh "RN   [%d]\n", $ref_count++;
        print $out_fh "RP   1-$seqlen\n";
        map { print $out_fh "$_\n" } @{$refs[$i]};
        print $out_fh "XX\n";
      }

      next;
    }
    
    #
    # Comments
    #
    if (/^CC   / and not $sequencelevel) {
      if (not $written_version_tag{$seqname}) {
        printf $out_fh "CC   Annotated features correspond to WormBase release %s.\n", $wormbase->get_wormbase_version_name;
        print $out_fh "XX\n";
        $written_version_tag{$seqname} = 1;
      }
    }
    
    #
    # Feature table
    #
    if (/^FT   (\S+)\s+(.+)/) {
      my ($ftype, $content) = ($1, $2);
      
      push @features, {
        ftype     => $ftype,
        location  => [$content],
        quals     => [],
      };
      next;
    } elsif (/^FT\s+(.+)/) {
      my $content = $1;
      if ($content =~ /^\/\w+=/) {
        push @{$features[-1]->{quals}}, [$content];
      } else {
        if (not @{$features[-1]->{quals}}) {
          push @{$features[-1]->{location}}, $content;
        } else {        
          push @{$features[-1]->{quals}->[-1]}, $content;
        }
      }
      next;
    }
    
    if (/^SQ/) {
      &process_feature_table($out_fh, $seqname, @features);
      if ($sequencelevel) {
        print $out_fh $_;
      } elsif (defined $agp_segs) {
        $log->log_and_die("Could not find AGP segs for $seqname\n") 
            if not exists $agp_segs->{$seqname};

        my $con_string = &make_con_from_segs(@{$agp_segs->{$seqname}});

        my $sbreak = $Text::Wrap::break;
        my $scols  = $Text::Wrap::columns;

        $Text::Wrap::break = '(?<=[,])';
        $Text::Wrap::columns = 76;
        my $wrapped = wrap('', '', $con_string);

        $Text::Wrap::break   = $sbreak;
        $Text::Wrap::columns = $scols;

        foreach my $line (split(/\n/, $wrapped)) {
          print $out_fh "CO   $line\n";
        }
      }
      next;
    }
    
    if ($written_header and 
        (/^\S/ or $sequencelevel)) {
      print $out_fh $_;
    }
  }
  
  close($raw_fh);
  close($out_fh);
}


if ($stage_to_repository) {
  ######################################################################
  # stage the dump to the per-cosmid submissions repository
  ######################################################################

  if (not $mod_dump_file) {
    $log->log_and_die("You cannot stage without producing a modified dump file first " .
                      "(either with -dumpmodified or -moddumpfile <file>\n");
  }

  my $modified = &stage_dump_to_submissions_repository($mod_dump_file);
  $log->write_to("\nStaged EMBL dumps to submission repository - $modified files are now locally modified\n");
}

unlink $raw_dump_file if $delete_raw_dump_file;
unlink $mod_dump_file if $delete_mod_dump_file;

$log->mail();
exit(0); 

###################################################
#                 SUBROUTINES                     #
###################################################


############################
sub process_feature_table {
  my ($out_fh, $seqname, @feats) = @_;

  if ($sequencelevel) {
    @feats = grep { $_->{ftype} eq 'source' } @feats;
  }

  foreach my $feat (@feats) {

    if ($feat->{ftype} eq 'source') {
      printf $out_fh "FT   %-16s%s\n", $feat->{ftype}, $feat->{location}->[0];
      printf $out_fh "FT   %16s/db_xref=\"taxon:%d\"\n", " ", $wormbase->ncbi_tax_id;
      printf $out_fh "FT   %16s/strain=\"%s\"\n", " ", $species_info{strain}->{$species};
      printf $out_fh "FT   %16s/mol_type=\"genomic DNA\"\n", " ";
      if ($seqname =~ /CHROMOSOME_(\S+)/) {
        printf $out_fh "FT   %16s/chromosome=\"$1\"\n", " ";
      } elsif ($seqname =~ /chr(\S+)/) {
        printf $out_fh "FT   %16s/chromosome=\"$1\"\n", " ";
      }

      foreach my $tag (@{$feat->{quals}}) {
        foreach my $ln (@$tag) {
          printf $out_fh "FT   %16s%s\n", " ", $ln;
        }
      }
      if ($species eq 'briggsae' ||$species eq 'brugia') {
        printf $out_fh "FT   %16s/note=\"supercontig %s\"\n", " ", $seqname;
      }
      next;
    } elsif ($feat->{ftype} =~ /RNA$/) {
      #
      # 1st FT line should be one of 4
      # FT    ncRNA
      # FT    rRNA
      # FT    tRNA
      # FT    misc_RNA
      # Supported bio types for ncRNA
      #  /ncRNA_class="antisense_RNA"
      #  /ncRNA_class="miRNA"
      #  /ncRNA_class="piRNA"
      #  /ncRNA_class="scRNA"              
      #  /ncRNA_class="siRNA"
      #  /ncRNA_class="snoRNA"
      #  /ncRNA_class="snRNA"
      #  /ncRNA_class="other"
      # Nothing else counts

      my $mod_dir = $feat->{ftype};
      my ($rna_class, $rna_note);

      if ($mod_dir eq 'antisense_RNA' or
          $mod_dir eq 'miRNA' or
          $mod_dir eq 'piRNA' or
          $mod_dir eq 'scRNA' or
          $mod_dir eq 'siRNA' or 
          $mod_dir eq 'snoRNA' or 
          $mod_dir eq 'snRNA') {
        $rna_class = $mod_dir;
        $mod_dir = 'ncRNA';
      } elsif ($mod_dir eq 'ncRNA') {
        $rna_class = 'other';
      } elsif ($mod_dir eq 'tRNA' or  
               $mod_dir eq 'rRNA' or
               $mod_dir eq 'misc_RNA' or
               $mod_dir eq 'precursor_RNA' or
               $mod_dir eq 'prim_transcript') {
        # no class, do nothing
      } else {
        # for all other RNA types, pass them through as
        # ncRNA/other, but record the type in a note
        $rna_class = ($mod_dir eq 'lincRNA')? "lncRNA" : "other";
        $rna_note = "biotype:$mod_dir";
        $mod_dir = "ncRNA";
      }

      if (defined $rna_class) {
        $feat->{rna_class} = $rna_class;
        push @{$feat->{quals}}, ["/ncRNA_class=\"$rna_class\""];
      } 
      if (defined $rna_note) {
        $feat->{rna_note} = $rna_note;
      }
      $feat->{ftype} = $mod_dir;

    } elsif ($feat->{ftype} =~ /(\S+)_Pseudogene/) {
      my $new_dv = "$1";

      $feat->{ftype} = $new_dv;
      $feat->{is_pseudo} = 1;
      # we don't currently distinguish between the different
      # ENA pseudogene classes, so make them all "unknown" for now
      push @{$feat->{quals}}, ["/pseudogene=\"unknown\""];
    }

    #
    # Now do qualifiers. Pass through the list picking out information
    # we are interested in, and filter out the ones we don't want. 
    # The rest we keep as-are
    #
    
    my (%revised_quals,
        $wb_isoform_name,
        $wb_geneid,
        $gseq_name,
        $gene_qual, 
        $locus_tag_qual, 
        $standard_name_qual, 
        $product_qual);

    foreach my $qual (@{$feat->{quals}}) {
      if ($qual->[0] =~ /^\/standard_name=\"(\S+)\"/) {
        $wb_isoform_name = $1;
        $standard_name_qual = $qual;        
        $feat->{wb_isoform_name} = $wb_isoform_name;
        if (not exists $trans2gene->{$wb_isoform_name}) {
          $log->log_and_die("Could not find WBGene id for $wb_isoform_name\n");
        }
        $wb_geneid = $trans2gene->{$wb_isoform_name};

        if (not exists $gene2gseqname->{$wb_geneid}) {
          $log->log_and_die("Could not find gseqname for $wb_geneid\n");
        }
        $gseq_name =  $gene2gseqname->{$wb_geneid};
      } elsif ($qual->[0] =~ /\/product=/ or
               $qual->[0] =~ /\/gene=/ or
               $qual->[0] =~ /\/protein_id=/ or
               $qual->[0] =~ /\/note=/) {
        next;
      } elsif ($qual->[0] =~ /\/([^=]+)=?/) {
        push @{$revised_quals{$1}}, $qual;
      }
    }
    
    
    #
    # /product
    #
    my $prod_name = "";
    if ($feat->{ftype} eq 'CDS') {
      if (exists $trans2briefid_h->{$wb_isoform_name}) {
        $prod_name = $trans2briefid_h->{$wb_isoform_name};
      } else {
        if ($feat->{is_pseudo}) {
          $prod_name = "Pseudogenic transcript $wb_isoform_name";
        } else {     
          # new guidelines from ENA/UniProt; do not put gene name
          # in product qualifier. Default to Uncharacterized protein
          $prod_name = "Uncharacterized protein";

          #if (exists $gene2cgcname->{$wb_geneid}) {
          #  $prod_name .= " " . uc($gene2cgcname->{$wb_geneid});
          #}

          # The following is odd, but it is what UniProt wanted...
          #my $isoform_suf;
          #if ($wb_isoform_name =~ /^$gseq_name(.+)$/) {
          #  $isoform_suf = $1;
          #}          
          #if (defined $isoform_suf) {
          #  $prod_name .= " (isoform $isoform_suf)";
          #}
        }
      }
    } elsif ($feat->{ftype} eq 'mRNA') {
      my $isoform_suf;
      if ($wb_isoform_name =~ /^$gseq_name(.+)$/) {
        $isoform_suf = $1;
        $isoform_suf =~ s/^\.//;
      }
      if ($gene2cgcname->{$wb_geneid}) {
        $prod_name = uc($gene2cgcname->{$wb_geneid});
      }
      $prod_name .= " primary transcript";
      if ($isoform_suf) {
        $prod_name .= " $isoform_suf";
      }
      
    } elsif ($feat->{ftype} eq 'misc_RNA') {
      if ($gene2cgcname->{$wb_geneid}) {
        $prod_name = $gene2cgcname->{$wb_geneid};
      } else {
        $prod_name = $gseq_name;
      }
      $prod_name = "Non-coding transcript of $prod_name";
    } elsif ($feat->{ftype} =~ /RNA/) {
      # note that prim_transcript will not match this, but this is fine
      # because prim_transcript features are not allowed to have a product
      my $pname = ($gene2cgcname->{$wb_geneid}) ? $gene2cgcname->{$wb_geneid} : $wb_isoform_name;

      if ($feat->{ftype} eq 'precursor_RNA') {
        $prod_name = "microRNA $pname precursor";
      } elsif (($feat->{ftype} eq 'tRNA' or $feat->{ftype} eq 'rRNA') and not $feat->{is_pseudo}) {
        if ($trans2briefid_h->{$wb_isoform_name}) {
          $prod_name = $trans2briefid_h->{$wb_isoform_name};
        } else {
          $prod_name = $wb_isoform_name;
        }
      } elsif (exists $feat->{rna_class} and $feat->{rna_class} eq 'miRNA') {
        $prod_name = "microRNA $pname";
        if ($wb_isoform_name =~ /([ab])$/) {
          $prod_name .= " ($1)";
        }
      } else {
        $prod_name = "RNA transcript $pname";
      }
    }
    $product_qual = ["/product=\"$prod_name\""] if $prod_name;

    #
    # /protein_id and prediction_status notes
    #
    if ($feat->{ftype} eq 'CDS') {
      if (exists $cds2proteinid_h->{$wb_isoform_name} and
          exists $cds2proteinid_h->{$wb_isoform_name}->{$seqname} 
          
          # dont include protein-ids for multi-clone objects; they cause trouble
          #not exists $spans->{$wb_isoform_name}) {
          ) {
        my $pid = $cds2proteinid_h->{$wb_isoform_name}->{$seqname};
        push @{$revised_quals{protein_id}}, ["/protein_id=\"$pid\""];
      }

      my $status_note;
      if (defined $cds2status_h->{$wb_isoform_name}) {
        if ($cds2status_h->{$wb_isoform_name} eq 'Confirmed') {
          $status_note = "Confirmed by transcript evidence";
        } elsif ($cds2status_h->{$wb_isoform_name} eq 'Partially_confirmed') {
          $status_note = "Partially confirmed by transcript evidence";
        } else {
          $status_note = "Predicted";
        }
      }
      if (defined $status_note) {
        push @{$revised_quals{note}}, ["/note=\"$status_note\""];
      }
    }

    #
    # /gene
    #
    my $gene_qual_val;
    if ($gene2cgcname->{$wb_geneid}) {
      $gene_qual_val = $gene2cgcname->{$wb_geneid};
    } else {
      $gene_qual_val = $gseq_name;
    }
    $gene_qual = ["/gene=\"$gene_qual_val\""];

    #
    # locus_tag
    #
    # note that for multi-gene loci (i.e. "isoforms" of a single locus that
    # been attched to distinct gene; there are a handful of these in C. elegans)
    # we cannot reliably use the gene sequence name as the locus tag, because
    # the /gene and /locus_tag qualifiers must be consistent across features
    # (i.e. two features with the same /locus_tag must also have the same /gene)
    # For these then, we just treat them as two separate gene, and use the thing
    # we put in the gene qualifier as the locus_tag suffix

    my $lt = $gseq_name;
    if (exists $multi_gene_loci->{$lt}) {
      $lt = $gene_qual_val;
    }
    if ($species eq 'briggsae' and $lt =~ /CBG(\S+)/) {
      $lt = "CBG_$1";
    } elsif ($species eq 'elegans') {
      $lt = "CELE_${lt}";
    } elsif ($species eq 'brugia'){
      $lt = "BM_${lt}";
    }
    $locus_tag_qual = [sprintf("/locus_tag=\"%s\"", $lt)];

   
    #
    # Other qualifiers
    #
    &add_additional_qualifiers( $feat );

    if ($additional_qualifiers{$wb_isoform_name}) {
      foreach my $qk (keys %{$additional_qualifiers{$wb_isoform_name}}) {
        foreach my $qual (@{$additional_qualifiers{$wb_isoform_name}->{$qk}}) {
          my $qual_string = "/$qk=$qual";
          my @wqs = split(/\n/, wrap('','',$qual_string));
          push @{$revised_quals{$qk}}, \@wqs;
        }
      }
    }

    if ($feat->{ftype} eq 'CDS') {      
      if ($trans2dbremark_h->{$wb_isoform_name}) {
        my $rem = sprintf("%s", $trans2dbremark_h->{$wb_isoform_name});
        $rem =~ s/\s+/ /g;
        
        my $rem_line = "/note=\"$rem\"";
        my @wl = split(/\n/, wrap('','',$rem_line));
        
        push @{$revised_quals{note}}, \@wl;
      }
    } elsif ($feat->{ftype} eq 'ncRNA' and 
             $feat->{rna_class} eq 'other') { 
      
      # We can only have one note for ncRNA features
      if ($feat->{rna_note}) {
        push @{$revised_quals{note}}, ["/note=\"$feat->{rna_note}\""];
      } elsif (exists $trans2type_h->{$wb_isoform_name} and 
          exists $trans2type_h->{$wb_isoform_name}->{note} and
          $trans2type_h->{$wb_isoform_name}->{note} ) {
        my $type_rna_note = $trans2type_h->{$wb_isoform_name}->{note};
        push @{$revised_quals{note}}, ["/note=\"$type_rna_note\""];
      } elsif ($trans2briefid_h->{$wb_isoform_name} and 
          $trans2briefid_h->{$wb_isoform_name} ne "ncRNA") {
        my $bid_rna_note = $trans2briefid_h->{$wb_isoform_name};
        push @{$revised_quals{note}}, ["/note=\"$bid_rna_note\""];
      }
    }
    

    #
    # Finally, print them all out in a consistent sensible order
    #

    printf $out_fh "FT   %-16s%s\n", $feat->{ftype}, $feat->{location}->[0];
    for(my $i=1; $i < @{$feat->{location}}; $i++) {
      printf $out_fh "FT   %16s%s\n", " ", $feat->{location}->[$i];
    }
    
    foreach my $qual ($gene_qual, $locus_tag_qual, $standard_name_qual, $product_qual) {
      foreach my $line (@$qual) {
        printf $out_fh "FT   %16s%s\n", " ", $line;
      }
    }
    foreach my $k (sort keys %revised_quals) {
      next if $k eq 'note';
      foreach my $qual (@{$revised_quals{$k}}) {
        foreach my $line (@$qual) {
          printf $out_fh "FT   %16s%s\n", " ", $line;
        }
      }
    }
    if (exists ($revised_quals{note})) {
      foreach my $qual (@{$revised_quals{note}}) {
        foreach my $line (@$qual) {
          printf $out_fh "FT   %16s%s\n", " ", $line;
        }
      }
    }
  }
  print $out_fh "XX\n";
}

###########################
sub add_additional_qualifiers {
  my ($feat) = @_;

  my $wb_isoform_name = $feat->{wb_isoform_name};

  #
  # Selenoproteins
  #
  if (exists $seleno_prots->{$wb_isoform_name}) {
    
    my ($strand, @exons) = &parse_location($feat->{location});
    
    my $cds_cdna_len = $seleno_prots->{$wb_isoform_name}->{cds_cdna_length};
    foreach my $seleno_pos_start (@{$seleno_prots->{$wb_isoform_name}->{selenos}}) {
      if ($strand < 0) {
        $seleno_pos_start = $cds_cdna_len - $seleno_pos_start - 1;
      }
      my $seleno_pos_end = $seleno_pos_start + 2;

      my $cur_cdna_pos = 0;
      my @segs;
      foreach my $ex (@exons) {
        for(my $gen_pos = $ex->{start}; $gen_pos <= $ex->{end}; $gen_pos++) {
          $cur_cdna_pos++;
          
          if ($cur_cdna_pos >= $seleno_pos_start and $cur_cdna_pos <= $seleno_pos_end) {
            if (not @segs or $segs[-1]->{end} + 1 < $cur_cdna_pos) {
              push @segs, {
                start => $gen_pos,
                end   => $gen_pos,
              };
            } else {
              $segs[-1]->{end} = $gen_pos,
            }
          }
        }
      }
      @segs = map { sprintf("%d..%d", $_->{start}, $_->{end}) }@segs;
      my $pos_str = (scalar(@segs) > 1) ? "join(" . join(",", @segs) . ")" : $segs[0];

      # transl_except => ["(pos:4060..4062,aa:Sec)"],      
      push @{$additional_qualifiers{$wb_isoform_name}->{transl_except}},  "(pos:$pos_str,aa:Sec)";
    }
  }
}



##########################
sub parse_location {
  my ($loc_line_arr) = @_;

  my $loc_string = "";
  foreach my $line (@$loc_line_arr) {
    if ($line =~ /(\S+)/)  {
      $loc_string .= $1;
    }
  }

  my $strand = 1;
  if ($loc_string =~ /complement/) {
    $strand = -1;
  }

  while(1) {
    if ($loc_string =~ /^[^\(]*(complement|join)\(.+\)[^\)]*$/) {
      $loc_string =~ s/^([^\(]*)(complement|join)\((.+)\)([^\)]*)$/$1$3/;
    } else {
      last;
    }
  }

  my @exons;
  foreach my $comp (split(/,/, $loc_string)) {
    my ($st, $en) = $comp =~ /^(\d+)\.\.(\d+)$/;

    push @exons, {
      start => $st,
      end   => $en,
    }
  }

  return ($strand, @exons);
}

##########################
sub get_multi_gene_loci {
  my ($gene2gseqname) = @_;

  # multi gene loci: when two or more transcripts have
  # the same sequence name, but multiple genes

  my (%multi_gene_loci, %gseqnames);

  foreach my $wbg (sort keys %{$gene2gseqname}) {
    my $gseqname = $gene2gseqname->{$wbg};

    $gseqnames{$gseqname}->{$wbg} = 1;
  }

  foreach my $gseqname (keys %gseqnames) {
    if (scalar(keys %{$gseqnames{$gseqname}}) > 1) {
      $multi_gene_loci{$gseqname} = 1;
    }
  }

  return \%multi_gene_loci;
}

##########################
sub get_references {

  my %primary_references = (
    elegans =>  [
      [ 
        "RX   PUBMED; 9851916.",
        "RG   Caenorhabditis elegans Sequencing Consortium",
        "RA   Sulson J.E., Waterston R.;",
        "RT   \"Genome sequence of the nematode C. elegans: a platform for investigating",
        "RT   biology\";",
        "RL   Science 282(5396):2012-2018(1998).",
      ],
    ],
    
    briggsae => [
      ["RG   WormBase Consortium",
       "RA   Howe K.L.;",
       "RT   ;",
       "RL   Submitted (03-OCT-2011) to the INSDC.",
       "RL   WormBase group, European Bioinformatics Institute,",
       "RL   Cambridge, CB10 1SA, UNITED KINGDOM.",
       ],
      [
       "RX   PUBMED; 14624247.",
       "RA   Stein L.D., Bao Z., Blasiar D., Blumenthal T., Brent M.R., Chen N.,",
       "RA   Chinwalla A., Clarke L., Clee C., Coghlan A., Coulson A., D'Eustachio P.,",
       "RA   Fitch D.H., Fulton L.A., Fulton R.E., Griffiths-Jones S., Harris T.W.,",
       "RA   Hillier L.W., Kamath R., Kuwabara P.E., Mardis E.R., Marra M.A.,",
       "RA   Miner T.L., Minx P., Mullikin J.C., Plumb R.W., Rogers J., Schein J.E.,",
       "RA   Sohrmann M., Spieth J., Stajich J.E., Wei C., Willey D., Wilson R.K.,",
       "RA   Durbin R., Waterston R.H.;",
       "RT   \"The genome sequence of Caenorhabditis briggsae: a platform for comparative",
       "RT   genomics\";",
       "RL   PLoS Biol. 1(2):E45-E45(2003).",
      ],
      [ 
        "RX   PUBMED; 21779179.",
        "RA   Ross J.A., Koboldt D.C., Staisch J.E., Chamberlin H.M., Gupta B.P.,",
        "RA   Miller R.D., Baird S.E., Haag E.S.;",
        "RT   \"Caenorhabditis briggsae recombinant inbred line genotypes reveal",
        "RT   inter-strain incompatibility and the evolution of recombination\";",
        "RL   PLoS Gen. 7(7):E1002174-E1002174(2011).",
      ],
    ],
    brugia => [
      ["RG   WormBase Consortium",
       "RA   Ghedin E., Paulini M.;",
       "RT   ;",
       "RL   Submitted (12-Dec-2012) to the INSDC.",
       "RL   WormBase Group, European Bioinformatics Institute,",
       "RL   Cambridge, CB10 1SA, UNITED KINGDOM.",
       ],
      ["RX   PUBMED; 17885136.",
       "RA   Ghedin E., Wang S., Spiro D., Caler E., Zhao Q., Crabtree J., Allen J.E.,",
       "RA   Delcher A.L., Guiliano D.B., Miranda-Saavedra D., Angiuoli S.V., Creasy T.,",
       "RA   Amedeo P., Haas B., El-Sayed N.M., Wortman J.R., Feldblyum T., Tallon L.,",
       "RA   Schatz M., Shumway M., Koo H., Salzberg S.L., Schobel S., Pertea M., Pop M.,",
       "RA   White O., Barton G.J., Carlow C.K., Crawford M.J., Daub J., Dimmic M.W.,",
       "RA   Estes C.F., Foster J.M., Ganatra M., Gregory W.F., Johnson N.M., Jin J.,",
       "RA   Komuniecki R., Korf I., Kumar S., Laney S., Li B.W., Li W., Lindblom T.H.,",
       "RA   Lustigman S., Ma D., Maina C.V., Martin D.M., McCarter J.P., McReynolds L.,",
       "RA   Mitreva M., Nutman T.B., Parkinson J., Peregrin-Alvarez J.M., Poole C., Ren Q.,",
       "RA   Saunders L., Sluder A.E., Smith K., Stanke M., Unnasch T.R., Ware J., Wei A.D.,",
       "RA   Weil G., Williams D.J., Zhang Y., Williams S.A., Fraser-Liggett C., Slatko B.,",
       "RA   Blaxter M.L., Scott A.L.;",
       "RT   \"Draft genome of the filarial nematode parasite Brugia malayi.\";",
       "RL   Science. 2007 Sep 21;317(5845):1756-60.",
      ],
    ],
    sratti => [
      [
       "RA   Aslett M.;",
       "RT   ;",
       "RL   Submitted (04-NOV-2013) to the INSDC.",
       "RL   Pathogen Sequencing Unit, Wellcome Trust Sanger Institute, Wellcome Trust",
       "RL   Genome Campus, Hinxton, Cambridge, Cambridgeshire CB10 1SA, UNITED KINGDOM.",
      ],
    ],
    );

  return $primary_references{$species};
}


###################################
sub stage_dump_to_submissions_repository {
  my ($dump_file) = @_;

  my (%records, $current_lines, $cosmid);

  my $submit_repo = $wormbase->submit_repos;

  open(my $dfg, $dump_file) or $log->log_and_die("Could not open $dump_file for reading\n");

  while(<$dfg>) {
    /^ID/ and do {
      $cosmid = undef;
      $current_lines = [];
    };

    /^DE\s+.+\s+(\S+)$/ and do {
      if (not defined $cosmid ) {
        $cosmid = $1;
        foreach my $line (@$current_lines) {
          push @{$records{$cosmid}->{embl}}, $line;
        }
        $current_lines = $records{$cosmid}->{embl};
      }
    };
    
    if (/^\s+(.+)$/) {
      my $seq = $1;
      $seq =~ s/\s+//g;
      $records{$cosmid}->{seq} .= uc($seq);
    }
    
    push @$current_lines, $_;
  }
  
  foreach my $cosmid (sort keys %records) {    
    next if not exists $records{$cosmid}->{seq};
    
    my $hash = $wormbase->submit_repos_hash_from_sequence_name($cosmid); 

    my $seq_file  = "$submit_repo/$hash/$cosmid/$cosmid.fasta";      
    if (not -e $seq_file) {
      $log->log_and_die("Could not find the FASTA file for $cosmid ($seq_file)");
    }
    
    my $sio = Bio::SeqIO->new(-format => 'fasta',
                              -file   => $seq_file);
    my $seqobj = $sio->next_seq;
    
    if (uc($seqobj->seq) ne $records{$cosmid}->{seq}) {
      # sequence has changed
      my $newseqobj = Bio::PrimarySeq->new(-id => $cosmid,
                                           -seq => $records{$cosmid}->{seq});
      my $sioout = Bio::SeqIO->new(-format => 'fasta',
                                   -file   => ">$seq_file");
      $sioout->write_seq($newseqobj);
      $sioout->close();
    }
  }
  
  foreach my $cosmid (sort keys %records) {
     my $hash = ($species eq 'elegans' and not $sequencelevel) 
         ? "0"
         : $wormbase->submit_repos_hash_from_sequence_name($cosmid);
      
     my $embl_file = "$submit_repo/$hash/$cosmid/$cosmid.embl";
     if (not -e $embl_file) {
       $log->log_and_die("Could not find the EMBL file for $cosmid ($embl_file)");
     }
       
     open my $emblfh, ">$embl_file" or $log->log_and_die("Could not open $embl_file for writing\n");
     foreach my $ln (@{$records{$cosmid}->{embl}}) {
       print $emblfh $ln;
     }
     close($emblfh);
  }


  my $modified_count = 0;

  open(my $gitcmd, "cd $submit_repo && git ls-files -m |")
      or $log->write_to("Warning: could not run git command to get list of locally modified files\n");
  while(<$gitcmd>) {
    /^(\S+)/ and do {
      $log->write_to("   $1 locally modified\n");
      $modified_count++;
    };
  }

  return $modified_count;
}


##############################################
sub read_agp_files {
  my @f = @_;

  my %segs;

  foreach my $f (@f) {
    my $fh;
    if ($f =~ /\.gz$/) {
      open($fh, "gunzip -c $f |") or $log->log_and_die("Could not open zipped agp file $f for reading\n");
    } else {
      open($fh, $f) or $log->log_and_die("Could not open agp file $f for reading\n");
    }
    while(<$fh>) {
      my @l = split(/\t/, $_);

      # Grr, having to deal with that CHROMOSOME prefix again...
      if ($species eq 'elegans' and grep { $l[0] eq $_ } ('I','II','III','IV','V','X')) {
        $l[0] = "CHROMOSOME_$l[0]";
      }

      my $seg =  {
        chr       => $l[0],
        start     => $l[1], 
        end       => $l[2],
        type      => ($l[4] eq 'N') ? 'gap' : 'seq',
      };

      if ($seg->{type} eq 'seq') {
        $seg->{ctg_acc} = $l[5];
        $seg->{ctg_start} = $l[6];
        $seg->{ctg_end} = $l[7];
        $seg->{ctg_ori} = $l[8];
      };

      push @{$segs{$l[0]}}, $seg;
    }
  }

  return \%segs;
}

##############################################
sub make_con_from_segs {
  my @segs = @_;

  my @components;
  foreach my $seg (@segs) {
    if ($seg->{type} eq 'gap') {
      push @components, sprintf("gap(%d)", $seg->{end} - $seg->{start} + 1);
    } else {
      my $reg = sprintf("%s:%d..%d", $seg->{ctg_acc}, $seg->{ctg_start}, $seg->{ctg_end});
      if ($seg->{ctg_ori} eq '-') {
        $reg = "complement($reg)";
      }
      push @components, $reg;
    }
  }
  return "join(" . join(",", @components) . ")";
}

##############################################
## Database fetching methods.
##############################################

################################
sub fetch_clone2type {
  my ($db) = @_;

  my %clone2type;
  $wormbase->FetchData('clone2type', \%clone2type, "$db/COMMON_DATA");
  
  return \%clone2type;
}

################################
sub fetch_clone2dbid {
  my ($db) = @_;

  my %clone2dbid;
  $wormbase->FetchData('clone2dbid', \%clone2dbid, "$db/COMMON_DATA");

  return \%clone2dbid;
}

################################
sub fetch_transcript2gene {
  my ($db) = @_;

  my (%transcript2gene);

  $wormbase->FetchData('worm_gene2geneID_name', \%transcript2gene, "$db/COMMON_DATA");
  return \%transcript2gene;
}


#####################
sub fetch_selenoproteins {
  my ($ref_db) = @_;
  
  my %selenos;
  
  if ($quicktest) {
    return \%selenos;
  }

  my $db = Ace->connect(-path => $ref_db);
  my $iter = $db->fetch_many(-query => "Find CDS WHERE Genetic_code = \"Selenocysteine\"");
  while(my $cds = $iter->next) {
    my $pep_seq = $cds->asPeptide;

    $pep_seq =~ s/>\S+\n//;
    $pep_seq =~ s/\n//g; 
    # EMBL feature includes implicit STOP feature
    $pep_seq .= "*";
    my $pep_len = length($pep_seq);
    my $cds_len = $pep_len * 3;

    $selenos{$cds->name} = {
      cds_cdna_length => $cds_len,
      selenos    => [],
    };
    
    while ($pep_seq =~ /U/g) {
      my $match_pos = pos($pep_seq);
      
      # assumption: selenoproteins will always have a complete CDS
      # (so that we can trivially convert to CDS coords)
      my $cds_cdna_end = $match_pos * 3;
      my $cds_cdna_start = $cds_cdna_end - 2;

      push @{$selenos{$cds->name}->{selenos}}, $cds_cdna_start;
    }
  }

  $db->close();

  return \%selenos;
}


################################
sub fetch_gene2names {
  my ($db) = @_;

  my (%gene2cgc, %gene2gseqname);
  
  if (not $quicktest and -d $db) {

    $log->write_to("You are using $db to get gene names\n");

    my $query = &get_gene_name_query();
    my $command = "Table-maker -p $query\nquit\n";
    open(TACE, "echo '$command' | $tace $db |");
    while(<TACE>) {
      chomp;
      s/\"//g;
      next if /^acedb/;
      next if /^\/\//;
      next if not /\S/;

      my ($gene, $cgc_name, $gseqname) = split(/\t/, $_);

      $log->log_and_die("Could not get gseqname for $gene\n") if not $gseqname;
      
      if (defined $gene) {
        $gene2cgc{$gene} = $cgc_name if $cgc_name;
        $gene2gseqname{$gene} = $gseqname if $gseqname;
      }
    }
    unlink $query;
  }

  return (\%gene2cgc, \%gene2gseqname);
}


#######################
sub fetch_cds2status {
  my ($ref_db)= @_;

  my %cds2status;

  if ($quicktest) {
    return \%cds2status;
  }

  if (-d $ref_db) {

    $log->write_to("You are using $ref_db to get CDS status\n");
    
    my $query = &get_status_query();
    my $command = "Table-maker -p $query\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ref_db |");
    while (<TACE>) {
      chomp;
      s/\"//g;
      next unless (/^([A-Z,0-9,.]+?\w)\s+(\w+)/) ;
      my ($cds, $status) = ($1, $2);
      $cds2status{$cds} = $status;
    }
    close TACE;

    unlink $query;
  } elsif (lc($ref_db) eq 'none') {
    $log->write_to("You have specifed 'none' for the cds2status db, so will not be added\n");
  } else {
    $log->log_and_die("Could not find cds2status database '$ref_db'");
  }

  return \%cds2status;
}

#######################
sub fetch_cds2proteinid {
  my ($ref_db) = @_;

  my %cds2proteinid;

  if ($quicktest) {
    return \%cds2proteinid;
  }

  if (-d $ref_db) {
    $log->write_to("You are using $ref_db to get CDS protein_id\n");
    
    my $query = &get_protein_id_query();
    my $command = "Table-maker -p $query\nquit\n";
    open(TACE, "echo '$command' | $tace $ref_db |");
    while(<TACE>) {
      chomp;
      s/\"//g;
      if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)$/) {
        my ($cds, $seqname, $protein_id, $version) = ($1, $2, $3, $4);
        $cds2proteinid{$cds}->{$seqname} = "${protein_id}.${version}";
      }   
    }
    close(TACE);
    
    unlink $query;
  } elsif (lc($ref_db) eq 'none') {
    $log->write_to("You have specifed 'none' for the cds2protein_id db, so will not be added\n");
  } else {
    $log->log_and_die("Could not find cds2protein_id database '$ref_db'");
  }

  return \%cds2proteinid;
}

#######################
sub fetch_trans2dbremark {
  my ($ref_db) = @_;

  my %trans2dbremark;

  if ($quicktest) {
    return \%trans2dbremark;
  }

  if (-d $ref_db) {
    $log->write_to("You are using $ref_db to get CDS/Transcript db_remarks\n");
    
    my @qfiles;
    foreach my $class ("CDS", "Transcript") {
      my $query = &get_db_remark_query($class);
      push @qfiles, $query;
      my $command = "Table-maker -p $query\nquit\n";
      open(TACE, "echo '$command' | $tace $ref_db |");
      while(<TACE>) {
        chomp;
        if (/^\"(\S+)\"\s+\"(.+)\"$/) {
          my ($obj, $dbremark) = ($1, $2);
          $dbremark =~ s/\\//g;
          $trans2dbremark{$obj} = $dbremark;
        }
      }
    }
    
    unlink @qfiles;
  } elsif (lc($ref_db) eq 'none') {
    $log->write_to("You have specifed 'none' for the trans2dbremark db, so will not be added\n");
  } else {
    $log->log_and_die("Could not find trans2dbremark database '$ref_db'");
  }

  return \%trans2dbremark;

}

########################
sub fetch_transcript2type {
  my ($ref_db) = @_;

  my %tran2type;

  if (not $quicktest and -d $ref_db) {
    $log->write_to("You are using $ref_db to get Transcript types\n");
    
    my $query = &get_ncrna_type_query();
    my $command = "Table-maker -p $query\nquit\n";
    open(TACE, "echo '$command' | $tace $ref_db |");
    while(<TACE>) {
      chomp;
      s/\"//g;
      my ($obj, $type, $other_comment) = split(/\t/, $_);

      if ($type) {
        $tran2type{$obj}->{type} = $type;
        if ($other_comment) {
          $other_comment =~ s/\\//g;
          $tran2type{$obj}->{note} = $other_comment;
        }
      }
    }
    unlink $query;
  }

  return \%tran2type;

}

########################
sub fetch_transcript2briefid {
  my ($ref_db) = @_;

  my (%res, %tran2briefid);

  if (not $quicktest and -d $ref_db) {
    $log->write_to("You are using $ref_db to get CDS/Transcript Brief_identification\n");
    
    foreach my $class ("CDS", "Transcript", "Pseudogene") {
      my $query = &get_brief_id_query($class);
      my $command = "Table-maker -p $query\nquit\n";
      open(TACE, "echo '$command' | $tace $ref_db |");
      while(<TACE>) {
        chomp;
        next if not /^\"\S+\"/;
        s/\"//g; 
        my ($obj, $brief_id, $evi_database) = split(/\t/, $_);

        if (not $evi_database) {
          $evi_database = "no_evidence";
        }
        
        $brief_id =~ s/C\. elegans //; 
        $brief_id =~ s/\\//g; 

        next if $brief_id =~ /^Uncharacterized protein/;

        $res{$class}->{$obj}->{lc($evi_database)} = $brief_id;
      }
      unlink $query;
    }
  }

  foreach my $class (keys %res) {    
    foreach my $obj (keys %{$res{$class}}) {
      if ($class eq 'CDS') {
        # only include the UniProt ones for CDS for time being;
        # InterPro and no-evidence ones have questionable source
        if (exists $res{$class}->{$obj}->{uniprot}) {
          $tran2briefid{$obj} = $res{$class}->{$obj}->{uniprot};
        }
      } else {
        my ($evi) = values %{$res{$class}->{$obj}};
        $tran2briefid{$obj} = $evi;
      }
    }
  }

  return \%tran2briefid;
}



########################
sub get_db_remark_query {
  my ($class) = @_;

  my $tmdef = "/tmp/dbremark_tmquery.$class.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");

  my $condition = "";
  if ($single) {
    $condition = "Condition Sequence = \"$single\""
  }

  my $db_remark_tablemaker_query = <<"EOF";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class $class
From 1 
$condition

Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Text 
From 1 
Tag DB_remark 

EOF

  print $qfh $db_remark_tablemaker_query;
  close($qfh);

  return $tmdef;
}


########################
sub get_brief_id_query {
  my ($class) = @_;

  my $tmdef = "/tmp/briefid_tmquery.$class.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");

  my $condition = "";
  if ($single) {
    $condition = "Condition Sequence = \"$single\""
  }

  my $briefid_tablemaker_query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class $class 
From 1 
$condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Text 
From 1 
Tag Brief_identification 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Database 
Right_of 2 
Tag  HERE  # Accession_evidence 
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  

EOF

  print $qfh $briefid_tablemaker_query;
  close($qfh);

  return $tmdef;
}


##########################
sub get_protein_id_query {

  my $tmdef = "/tmp/pid_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";
  if ($single) {
    $condition = "Condition Sequence = \"$single\""
  }

  my $protein_id_tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class CDS
From 1 
$condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Sequence 
From 1 
Tag Protein_id 
 
Colonne 3 
Width 12 
Optional 
Visible 
Text 
Right_of 2 
Tag  HERE  
 
Colonne 4 
Width 12 
Optional 
Visible 
Integer 
Right_of 3 
Tag  HERE  
 
EOF

  print $qfh $protein_id_tablemaker_template;
  close($qfh);

  return $tmdef;

}

##########################
sub get_status_query {

  my $tmdef = "/tmp/status_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  


  my $condition = "Condition Live AND Species = \"$full_species_name\"";
  if ($single) {
    $condition .= " AND Sequence = \"$single\""
  }

  my $status_query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 30 
Optional 
Hidden 
Class 
Class Gene 
From 1 
$condition
 
Colonne 2 
Width 30 
Optional 
Visible 
Class 
Class CDS 
From 1 
Tag Corresponding_CDS 
 
Colonne 3 
Width 30 
Optional 
Visible 
Next_Tag 
From 2 
Tag Prediction_status 

EOF

  print $qfh $status_query;
  close($qfh);

  return $tmdef;
}

##########################
sub get_ncrna_type_query {

  my $tmdef = "/tmp/rnatype_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  


  my $condition = "Condition Live AND Species = \"$full_species_name\"";
  if ($single) {
    $condition .= " AND Sequence = \"$single\""
  }

  my $rnatype_query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Transcript 
From 1 
$condition
 
Colonne 2 
Width 12 
Optional 
Hidden 
Show_Tag 
From 1 
Tag Transcript 
 
Colonne 3 
Width 12 
Optional 
Visible 
Next_Tag 
Right_of 2 
 
Colonne 4 
Width 12 
Optional 
Visible 
Text 
Right_of 3 
Tag ncRNA 

EOF
 
  print $qfh $rnatype_query;
  close($qfh);

  return $tmdef;
  
}


##########################
sub get_gene_name_query {

  my $tmdef = "/tmp/gname_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  


  my $condition = "Condition Live AND Species = \"$full_species_name\"";
  if ($single) {
    $condition .= " AND Sequence = \"$single\"";
  } else {
    $condition .= " AND Sequence";
  }

  my $gname_query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
$condition
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag CGC_name 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag Sequence_name 

EOF
 
  print $qfh $gname_query;
  close($qfh);

  return $tmdef;
  
}



__END__

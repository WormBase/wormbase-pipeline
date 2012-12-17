#!/usr/local/bin/perl5.8.0 -w
#
# EMBLdump.pl :  makes modified EMBL dumps from camace.
# 
#  Last updated on: $Date: 2012-12-17 10:34:16 $
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
  genome_project_id => {
    elegans  => 13758,
    briggsae => 10731,
    brugia   => 10729,
  },

  taxon_id => {
    elegans  => 6239,
    briggsae => 473542,
    brugia   => 6279,
  },

  strain => {
    elegans => 'Bristol N2',
    briggsae => 'AF16',
    brugia   => 'FR3',
  },

  genome_project_accs => { },
);    

#
# Yuk. The below is a selenoprotein, the only one
# currently in elegans or briggsae. Since we do not store
# in the database the information required to populate the
# the transl_exep qualifier, we need to hard-code it here. 
# Like I said: yuk.
#

my %additional_qualifiers = (

  "C06G3.7" => { 
    transl_except => ["(pos:4060..4062,aa:Sec)"],
  },
  'Bm2025a' => {
    transl_except => ['(pos:823..825,aa:Sec)']
  },
  'Bm2025b' => {
    transl_except => ['(pos:1989..1991,aa:Sec)']
  },
  'Bm2025c' => {
    transl_except => ['(pos:2134..2136,aa:Sec)']
  },
  'Bm2025d' => {
    transl_except => ['(pos:1123..1125,aa:Sec)']
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
    $keep_redundant,
    $exclude_multis,
    %clone2type, %clone2acc, %cds2cgc, %rnagenes, %clone2dbid, %pseudo2cgc, $multi_gene_loci,
    $cds2proteinid_db, $cds2status_db, $cds2dbremark_db,
    $cds2status_h, $cds2proteinid_h,  $cds2dbremark_h
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "single=s@"       => \$single,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "dumpraw"         => \$dump_raw,
  "dumpmodified"    => \$dump_modified,
  "stagetorepo"     => \$stage_to_repository,
  "piddb=s"         => \$cds2proteinid_db,
  "dbremarkdb=s"    => \$cds2dbremark_db,
  "cdsstatusdb=s"   => \$cds2status_db,
  "modprivate=s"    => \$private,
  "moddumpfile:s"   => \$mod_dump_file,
  "rawdumpfile:s"   => \$raw_dump_file,
  "quicktest"       => \$quicktest,
  "keepredundant"   => \$keep_redundant,
  "excludemultis"   => \$exclude_multis,

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
    $log->write_to("You are raw embl dumping from $dbdir to $raw_dump_file wihch will be delted\n\n");
    $delete_raw_dump_file = 1;
  } else {
    $log->write_to("You are raw embl dumping from $dbdir to $raw_dump_file which will be retained\n\n");
  }
  
  my $command;
  if ($single) {
    $command .= "query find Sequence $single\ngif EMBL $raw_dump_file\n";# find sequence and dump
    $command .= "quit\n";# say you don't want to save and exit
  } else {
    $command .= "query find Sequence Genomic_canonical AND (From_laboratory = \"HX\" OR From_Laboratory = \"RW\")\n";
    $command .= "gif EMBL $raw_dump_file\n";# find sequence and dump
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

  %clone2type      = $wormbase->FetchData('clone2type');
  %cds2cgc         = $wormbase->FetchData('cds2cgc');
  %rnagenes        = $wormbase->FetchData('rna2cgc');
  %clone2dbid      = $wormbase->FetchData('clone2dbid');
  %pseudo2cgc      = $wormbase->FetchData('pseudo2cgc');

  $multi_gene_loci = &get_multi_gene_loci(\%cds2cgc, \%rnagenes, \%pseudo2cgc);

  # Some information is not yet available in autoace (too early in the build)
  # By default, pull cds2status, cds2dbremark and protein_ids from current_DB,
  # however, allow command-line override when necessary
  $cds2status_h = &fetch_cds2status( defined $cds2status_db ? $cds2status_db : $wormbase->database('current') );
  $cds2dbremark_h = &fetch_cds2dbremark( defined $cds2dbremark_db ? $cds2dbremark_db : $wormbase->database('current') );
  $cds2proteinid_h = &fetch_cds2proteinid( defined $cds2proteinid_db ? $cds2proteinid_db : $wormbase->database('current') );
  
  if (not defined $mod_dump_file) {
    $mod_dump_file = "/tmp/EMBL_dump.$$.mod.embl";
    $log->write_to("You are MODIFIED embl dumping to $mod_dump_file which will be deleted\n");
    $delete_mod_dump_file = 1;
  } else {
    $log->write_to("You are MODIFIED embl dumping to $mod_dump_file which will be retained\n");
  }

  my $span_hash = &get_locs_for_multi_span_objects($raw_dump_file);
  
  open(my $out_fh, ">$mod_dump_file") or $log->log_and_die("Could not open $mod_dump_file for writing\n");
  open(my $raw_fh, $raw_dump_file) or $log->log_and_die("Could not open $raw_dump_file for reading\n");
  
  my ($clone, $seqlen, $idline_suffix, @accs, @features, $written_header);
  
  while (<$raw_fh>) {
    
    # Store the necessary default ID line elements ready for use in the new style EMBL ID lines.
    if(/^ID\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\d+)\s+\S+/){
      ($clone, $seqlen) = ($1, $2);
      $idline_suffix = sprintf("SV XXX; linear; genomic DNA; %s; INV; $seqlen BP.", 
                               ($species eq 'briggsae') ? "CON" : "STD");
      
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
      if (exists $species_info{genome_project_accs}->{$species}) {
        foreach my $acc (@{$species_info{genome_project_accs}->{$species}}) {
          print $out_fh " $acc;";
        }
      }
      print $out_fh "\nXX\n";
      
      #
      # AC *
      #
      if ($clone2dbid{$clone}) {
        print $out_fh "AC * $clone2dbid{$clone}\n";
        print $out_fh "XX\n";
      }
      
      #
      # PR
      #
      print $out_fh "PR   Project:$species_info{genome_project_id}->{$species};\n";
      print $out_fh "XX\n";
      
      #
      # DE
      #
      my $de_line;
      
      if ($species eq 'elegans') {
        if (!defined($clone2type{$clone})){
          $de_line =  "$full_species_name clone $clone";
          $log->write_to("WARNING: no clone type for $_");
        } elsif (lc($clone2type{$clone}) eq "other") {
          $de_line = "$full_species_name clone $clone";
        } elsif (lc($clone2type{$clone}) eq "yac") {
          $de_line = "$full_species_name YAC $clone";
        } else {
          $de_line = "$full_species_name $clone2type{$clone} $clone";
        }
      } elsif ($species eq 'briggsae') {
        $de_line = "$full_species_name AF16 supercontig from assembly CB4, $clone";
      } elsif ($species eq 'brugia'){
        $de_line = "$full_species_name FR3 supercontig from assembly B_malayi-3.0 $clone"
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

      my @refs = @{&get_references()};
      for(my $i=0; $i < @refs; $i++) {
        printf $out_fh "RN   [%d]\n", $ref_count++;
        print $out_fh "RP   1-$seqlen\n";
        map { print $out_fh "$_\n" } @{$refs[$i]};
        print $out_fh "XX\n";
      }
      
      if ($species eq 'elegans') {
        printf $out_fh "RN   [%d]\n", $ref_count;
        printf $out_fh "RP   1-%d\n", $seqlen;
        print $out_fh "RG   WormBase Consortium\n";
        if (defined $primary_RA) {
          print $out_fh $primary_RA;
        } else {
          print $out_fh "RA   ;\n";
        }
        print $out_fh "RT   ;\n";
        print $out_fh $primary_RL if defined $primary_RL;
        print $out_fh "RL   Nematode Sequencing Project: Sanger Institute, Hinxton, Cambridge\n";
        print $out_fh "RL   CB10 1SA, UK and The Genome Institute at Washington University,\n"; 
        print $out_fh "RL   St. Louis, MO 63110, USA. E-mail: help\@wormbase.org\n";
        print $out_fh "XX\n";
      }
      next;
    }
    
    #
    # Comments
    #
    if (/^CC   For a graphical/) {
      while(<$raw_fh>) {
        last if /^XX/;
        
        if ($species eq 'elegans') {
          printf $out_fh "CC   Annotated features correspond to WormBase release %s.\n", $wormbase->get_wormbase_version_name;
          if ($exclude_multis) {
            print $out_fh "CC   Features spanning 2 or more clones excluded in this round.\n";
          }
          print $out_fh "XX\n";
          print $out_fh  "CC   For a graphical representation of this sequence and its analysis\n";
          print $out_fh  "CC   see:- http://www.wormbase.org/species/c_elegans/clone/$clone\n";
          print $out_fh  "XX\n";
        }elsif($species eq 'brugia'){
          print $out_fh "XX\n",
                        "CC   For a graphical representation of this sequence and its analysis\n",
                        "CC   see:- http://www.wormbase.org/species/c_elegans/clone/$clone\n",
                        "XX\n";
        }
      }
      next;
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
      &process_feature_table($out_fh, $clone, $span_hash, @features);
      if ($species eq 'elegans') {
        print $out_fh $_;
      }
      next;
    }
    
    if ($written_header and 
        (/^\S/ or $species eq 'elegans')) {
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
  my ($out_fh, $clone, $spans, @feats) = @_;

  foreach my $feat (@feats) {

    if ($feat->{ftype} eq 'source') {
      printf $out_fh "FT   %-16s%s\n", $feat->{ftype}, $feat->{location}->[0];
      printf $out_fh "FT   %16s/db_xref=\"taxon:%d\"\n", " ", $species_info{taxon_id}->{$species};
      printf $out_fh "FT   %16s/strain=\"%s\"\n", " ", $species_info{strain}->{$species};
      printf $out_fh "FT   %16s/mol_type=\"genomic DNA\"\n", " ";
      foreach my $tag (@{$feat->{quals}}) {
        foreach my $ln (@$tag) {
          printf $out_fh "FT   %16s%s\n", " ", $ln;
        }
      }
      if ($species eq 'briggsae' ||$species eq 'brugia') {
        printf $out_fh "FT   %16s/note=\"supercontig %s\"\n", " ", $clone;
      }
      next;
    } elsif ($feat->{ftype} =~ /RNA$/) {
      #
      # 1st FT line should be one of 3
      # FT    ncRNA
      # FT    rRNA
      # FT    tRNA
      # Supported bio types for ncRNA
      #  /ncRNA_class="miRNA"
      #  /ncRNA_class="siRNA"
      #  /ncRNA_class="scRNA"              
      #  /ncRNA_class="other"
      #  /ncRNA_class="snoRNA"
      #  /ncRNA_class="snRNA"
      # Nothing else counts

      my $mod_dir = $feat->{ftype};
      my $rna_class;

      if ($mod_dir eq 'snoRNA' or 
          $mod_dir eq 'miRNA' or
          $mod_dir eq 'siRNA' or 
          $mod_dir eq 'scRNA' or
          $mod_dir eq 'snRNA') {
        $rna_class = $mod_dir;
        $mod_dir = 'ncRNA';
      } elsif ($mod_dir eq 'misc_RNA' or
               $mod_dir eq 'ncRNA') {
        $rna_class = 'other';
        $mod_dir = 'ncRNA';
      } elsif ($mod_dir eq 'tRNA' or 
               $mod_dir eq 'rRNA') {
        # no class, do nothing
      } else {
        # for all other RNA types, pass them through as
        # ncRNA/other, but record the type in a note
        push @{$feat->{quals}}, ["/note=\"$mod_dir\""];
        $mod_dir = "ncRNA";
        $rna_class = "other";
      }

      if (defined $rna_class) {
        push @{$feat->{quals}}, ["/ncRNA_class=\"$rna_class\""];
        $feat->{class} = $rna_class;
      }
      $feat->{ftype} = $mod_dir;
    } elsif ($feat->{ftype} =~ /Pseudogene/) {
      # skip psudogenes for now; currently in briggsae, the only pseudogenes
      # are tRNAs, but we cannot assume this; therefore no way of determining
      # reliably what sort of psseudogene we are looking at
      next if $species eq 'briggsae';
        
      my $new_dv = "CDS";

      # hack: in C.elegans, all Pseudogenes with .t\d+ suffices are
      # tRNA pseudogenes; 
      foreach my $tg (@{$feat->{quals}}) {
        if ($tg->[0] =~ /\/gene=\"\S+\.t\d+\"/) {
          $new_dv = "tRNA";
        } 
      }

      $feat->{ftype} = $new_dv;
      # we don't currently distinguish between the different
      # ENA psedogene classes, so make them all "unknown" for now
      push @{$feat->{quals}}, ["/pseudogene=\"unknown\""];
    }

    #
    # Now do qualifiers. Pass through the list picking out information
    # we are interested in, and filter out the ones we don't want. 
    # The rest we keep as-are
    #
    
    my (%revised_quals,
        $wb_isoform_name,
        $gene_qual, 
        $locus_tag_qual, 
        $standard_name_qual, 
        $product_qual, 
        $db_remark, 
        $brief_ident);

    foreach my $qual (@{$feat->{quals}}) {
      if ($qual->[0] =~ /^\/standard_name=\"(\S+)\"/) {
        $wb_isoform_name = $1;
        $standard_name_qual = $qual;

        if ($cds2dbremark_h->{$wb_isoform_name}) {
          my $rem = sprintf("%s", $cds2dbremark_h->{$wb_isoform_name});
          $rem =~ s/\s+/ /g;

          my $rem_line = "/note=\"$rem\"";
          my @wl = split(/\n/, wrap('','',$rem_line));
          
          $db_remark = \@wl;
        }        
      } elsif ($qual->[0] =~ /\/note=\"similar to (.+\S)\s*$/) {
        my $bi = "\"$1";
        for(my $i=1; $i < @$qual; $i++) {
          $qual->[$i] =~ /^\s*(.+\S)\s*$/ and $bi .= " $1";
        }
        $bi =~ s/^\"//;
        $bi =~ s/\"$//;

        my @wl = split(/\n/, wrap('','',"/note=\"$bi\""));

        $brief_ident = \@wl;
        $feat->{brief_identification} = $bi;

      } elsif ($qual->[0] =~ /\/note=\"(.+)\-RNA\"/) {
        $feat->{rna_identification} = $1;
      } elsif ($qual->[0] =~ /\/note=\"preliminary prediction\"/ or
               $qual->[0] =~ /\/note=\"cDNA EST/ or 
               $qual->[0] =~ /\/product=/ or
               $qual->[0] =~ /\/gene=/ or
               $qual->[0] =~ /\/protein_id=/) {
        next;
      } elsif ($qual->[0] =~ /\/([^=]+)=?/) {
        push @{$revised_quals{$1}}, $qual;
      }
    }
    
    
    #
    # /product
    #
    if ($feat->{ftype} eq 'CDS') {
      my ($isoform_suf, $prod_name);
      if ($wb_isoform_name =~ /^(\S+)([a-z])$/) {
        $prod_name = $1;
        $isoform_suf = $2;
      } else {
        $prod_name = $wb_isoform_name;
      }
      if (exists $cds2cgc{$wb_isoform_name}) {
        $prod_name = uc($cds2cgc{$wb_isoform_name});
      }

      $product_qual = sprintf("/product=\"Protein %s%s\"",
                              $prod_name,
                              defined $isoform_suf ? ", isoform $isoform_suf" : "");
    } elsif ($feat->{ftype} =~ /RNA/) {
      if ($feat->{ftype} eq 'tRNA' and 
          exists $feat->{brief_identification}) {
        $product_qual = "/product=\"$feat->{brief_identification}\"";
      } elsif ($feat->{ftype} eq 'rRNA' and 
               exists $feat->{rna_identification}) {
        $product_qual = "/product=\"$feat->{rna_identification}\"";
      } else {
        if ($rnagenes{$wb_isoform_name}) {
          $product_qual = "/product=\"RNA transcript $rnagenes{$wb_isoform_name}\"";
        } else {
          $product_qual = "/product=\"RNA transcript $wb_isoform_name\"";
        }
      }
    }
    $product_qual = [$product_qual];

    #
    # /protein_id and prediction_status note
    #
    if ($feat->{ftype} eq 'CDS') {
      if (exists $cds2proteinid_h->{$wb_isoform_name} and
          exists $cds2proteinid_h->{$wb_isoform_name}->{$clone} and 
          # dont include protein-ids for multi-clone objects; they cause trouble
          not exists $spans->{$wb_isoform_name}) {
        my $pid = $cds2proteinid_h->{$wb_isoform_name}->{$clone};
        push @{$revised_quals{protein_id}}, ["/protein_id=\"$pid\""];
      }
      
      my $status_note;
      if (defined $cds2status_h->{$wb_isoform_name}) {
        if ($cds2status_h->{$wb_isoform_name} eq 'Confirmed') {
          $status_note = "/note=\"Confirmed by transcript evidence\"";
        } elsif ($cds2status_h->{$wb_isoform_name} eq 'Partially_confirmed') {
          $status_note = "/note=\"Partially confirmed by transcript evidence\"";
        } else {
          $status_note = "/note=\"Predicted\"";
        }
      }
      if (defined $status_note) {
        push @{$revised_quals{note}}, [$status_note];
      }
    }

    #
    # /gene
    #
    if ($cds2cgc{$wb_isoform_name}) {
      $gene_qual = ["/gene=\"$cds2cgc{$wb_isoform_name}\""];
    } elsif ($rnagenes{$wb_isoform_name}) {
      $gene_qual = ["/gene=\"$rnagenes{$wb_isoform_name}\""];
    } elsif ($pseudo2cgc{$wb_isoform_name}) {
      $gene_qual = ["/gene=\"$pseudo2cgc{$wb_isoform_name}\""];
    } else {
      if ($wb_isoform_name =~ /^(\S+\.\d+)[a-z]?/) {
        $gene_qual = ["/gene=\"$1\""];
      } elsif (($species eq 'briggsae' || $species eq 'brugia') and $wb_isoform_name =~ /^(\S+\d+)[a-z]$/) {
        $gene_qual = ["/gene=\"$1\""];
      } else {
        $gene_qual = ["/gene=\"$wb_isoform_name\""];
      }
    }
      
    #
    # locus_tag
    #
    my $lt = $wb_isoform_name;
    if ($lt =~ /^(\S+\.\d+)[a-z]?/) {
      $lt = $1;
    }
    $locus_tag_qual = [];
    if (not exists $multi_gene_loci->{$lt}) {
      if ($species eq 'briggsae' and $lt =~ /CBG(\d+)/) {
        $lt = "CBG_$1";
      } elsif ($species eq 'elegans') {
        $lt = "CELE_${lt}";
      } elsif ($species eq 'brugia'){
        $lt = "BM_${lt}";
      }
      push @{$locus_tag_qual}, sprintf("/locus_tag=\"%s\"", $lt);
    }
   
    #
    # Other qualifiers
    #
    
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
      if (defined $brief_ident) {
        push @{$revised_quals{note}}, $brief_ident;
      } elsif (defined $db_remark) {
        push @{$revised_quals{note}}, $db_remark;
      }
    } elsif ($feat->{ftype} eq 'ncRNA' and $feat->{class} eq 'other') {
      # only attempt to add notes to ncRNAs defined with "other" class
      if (exists $feat->{rna_identification}) {
        push @{$revised_quals{note}}, ["/note=\"$feat->{rna_identification}\""];
      } elsif (exists $feat->{brief_identification} and 
               $feat->{brief_identification} ne 'ncRNA' and
               $feat->{brief_identification} !~ /non-coding RNA gene$/) {
        push @{$revised_quals{note}}, $brief_ident;
      } elsif (defined $db_remark) {
        push @{$revised_quals{note}}, $db_remark;
      }
    }


    #
    # Finally, print them all out in a consistent sensible order
    #
    #

    if (exists $spans->{$wb_isoform_name} and 
        not exists $spans->{$wb_isoform_name}->{$clone}) {
      # this feature is annotated on multiple clones, 
      # and this is not the clone we have chosen to 
      # place it on
      if ($debug) {
        $log->write_to("Discarding duplicate feat: $wb_isoform_name will not be placed on $clone\n");
      }
      next;
    }

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
# this method obtains a list of multi-span objects, i.e. objects
# that span multiple clones. Since ENA no longer accept these
# annotations on multiple clones, we need to choose one clone
# to attach them to, and do this by selecting the clone where
# most of the object lies, using a series of rules
#
# Additionally, when sequence has ben updayted and there are mutual 
# dependencies between clones (i.e. clone X refers to clone Y and vice 
# versa), this can cause problems because clone X will be referring to 
# the "new" version of clone Y which has not been loaded yet. We therefore 
# provide the option of excluding all multi-clone features for one round
# of submissions. A second submission round will therefore be required to
# process the clones with foreign references

sub get_locs_for_multi_span_objects {
  my $file = shift;

  my ($clone, $acc, $type, $loc, $obj_name, %locs, %best_obj_locs);

  open(my $fh, $file) or $log->log_and_die("Could not open $file for reading\n");
  while(<$fh>) {
    /^ID\s+(\S+)/ and do {
      $clone = $1;
      $acc = "";
      next;
    };

    /^AC\s+(\S+);/ and do {
      $acc = $1 if not $acc;
    };

    /^FT   (\S+)\s+(\S+)$/ and do {
      ($type, $loc) = ($1, $2);
      while(<$fh>) {
        last if /^FT\s+\//;
        /FT\s+(\S+)/ and $loc .= $1;
      }
    };
      
    /^FT\s+\/standard_name=\"(\S+)\"/ and do {
      $obj_name = $1;
      
      push @{$locs{"$type:$obj_name"}}, {
        clone => $clone,
        acc   => $acc,
        location => $loc,
      };

      next;
    }
  }

  foreach my $obj (keys %locs) {
    my @locs = @{$locs{$obj}};
    my ($obj_name) = $obj =~ /:(\S+)/;

    next if scalar(@locs) == 1;

    $best_obj_locs{$obj_name} = {};

    next if $exclude_multis;
    
    foreach my $loc (@locs) {
      my ($local_extent, $foreign_extents) = &parse_location($loc->{location});
      $loc->{local_extent} = $local_extent;
      $loc->{foreign_extents} = $foreign_extents;
    }
    
    @locs = grep { $_->{local_extent} > 0 } @locs;
    @locs = sort { $b->{local_extent} <=> $a->{local_extent} } @locs;

    my $best_loc = 0;
    for(my $i=0; $i < @locs; $i++) {
      my $loc = $locs[$i];
      my $clone = $loc->{clone};

      if ($obj_name =~ /$clone/) {
        # sequence name of object matches clone name, 
        # so this is the best choice for this obj
        $best_loc = $i;
        last;
      }
    }


    my @best_locs = ($keep_redundant) ? @locs : ($locs[$best_loc]);

    foreach my $loc (@best_locs) {    
      my $clone = $loc->{clone};
      $best_obj_locs{$obj_name}->{$clone} = 1;
    }

  }

  return \%best_obj_locs;
}
  

##################################
# read_gff_annotation
#
# This method scans the GFF files for CDS features
# spanning clones in different orientations. AceDB
# will not include these in the dumps, so we need
# to add them in by hand (which of course begs the
# the question of why we need Acedb at all for the
# EMBL dumping...)
##################################
sub read_gff_annotation {
  my @types = ('curated');

  my (%additional_feats_by_clone);

  my $cc = Coords_converter->invoke($wormbase->autoace, 0, $wormbase);

  foreach my $chr ($wormbase->get_chromosome_names(-prefix => 1)) {
    my $gff_fh = $wormbase->open_GFF_file($chr, "curated", $log);
      
    my %trans;

    while(<$gff_fh>) {
      next if (/^\#/);
      my @data = split(/\t+/, $_);
      
      if ($data[2] eq 'exon') {
        my ($class, $obj_name) = ($data[8] =~ /(\S+)\s+\"(\S+)\"/);

        push @{$trans{$obj_name}->{exons}}, [$data[3], $data[4]];
        $trans{$obj_name}->{chr} = $chr;
        $trans{$obj_name}->{strand} = $data[6];

        $trans{$obj_name}->{start} = $data[3]
            if not exists $trans{$obj_name}->{start} or $trans{$obj_name}->{start} > $data[3];
        
        $trans{$obj_name}->{end} = $data[4]
            if not exists $trans{$obj_name}->{end} or $trans{$obj_name}->{end} < $data[4];
        
      }      
    }

    CDS: foreach my $cds (keys %trans) {
      my ($seq) = $cc->LocateSpan($trans{$cds}->{chr},
                                  $trans{$cds}->{start},
                                  $trans{$cds}->{end});
      
      if (not $cc->isa_clone($seq)) {
        # this CDS crosses boundaries, do break it down to constituent exons
        my (@exons, @mapped_exons, %strands);

        @exons = @{$trans{$cds}->{exons}};
        if ($trans{$cds}->{strand} eq '+') {
          @exons = sort { $a->[0] <=> $b->[0] } @exons;
        } else {
          @exons = sort { $b->[0] <=> $a->[0] } @exons;
        }

        foreach my $ex (@exons) {
          my ($ex_seq, $ex_st, $ex_en) = $cc->LocateSpan($trans{$cds}->{chr},
                                                         @$ex);
          if ($ex_en > $ex_st) {
            push @mapped_exons, {
              clone => $ex_seq,
              start => $ex_st,
              end   => $ex_en,
              strand => 1,
            }; 
          } else {
            push @mapped_exons, {
              clone => $ex_seq,
              start => $ex_en,
              end   => $ex_st,
              strand => -1,
            }; 
          }
          
          if (not $cc->isa_clone($ex_seq)) {
            # This gene has an exon crossing a supercontig boundary. Complicated. Messy. Ignore
            $log->write_to("Gene $cds has an exon crossing a supercontig boundary; ignoring it\n");
            next CDS;
          }
        }

        map { $strands{$_->{strand}} = 1 } @mapped_exons;
        
        if (scalar(keys %strands) > 1) {
          # this gene spans across neighbouring sctgs that are in different orientations. 
          # Acedb does (not) currently include these in the raw dump, so we will have to
          # insert them ourselves. 
          my $parent_clone = $mapped_exons[0]->{clone};

          my @join_els; 
          foreach my $ex (@mapped_exons) {
            if ($ex->{strand} < 0) {
              if ($ex->{clone} eq $parent_clone) {
                push @join_els, sprintf("complement(%d..%d)", $ex->{start}, $ex->{end});
              } else {
                push @join_els, sprintf("complement(%s:%d..%d)", $clone2acc{$ex->{clone}}, $ex->{start}, $ex->{end});
              }
            } else {
              if ($ex->{clone} eq $parent_clone) {
                push @join_els, sprintf("%d..%d", $ex->{start}, $ex->{end});
              } else {
                push @join_els, sprintf("%s:%d..%d", $ex->{clone}, $ex->{start}, $ex->{end});
              }
            }
          }
          my $join_str = sprintf("join(%s)", join(", ", @join_els));
          my @loc_lines =split(/\n/, wrap('','',$join_str));
          map { s/\s+//g } @loc_lines;

          my $feat = {
            location => \@loc_lines,
            ftype    => 'CDS',
            quals    => [["/standard_name=\"$cds\""]],
            clone    => $parent_clone,            
          };

          push @{$additional_feats_by_clone{$parent_clone}}, $feat;
        }
      }
    }
  }

  return \%additional_feats_by_clone;
}


##########################
sub parse_location {
  my ($loc_string) = @_;

  while(1) {
    if ($loc_string =~ /^[^\(]*(complement|join)\(.+\)[^\)]*$/) {
      $loc_string =~ s/^([^\(]*)(complement|join)\((.+)\)([^\)]*)$/$1$3/;
    } else {
      last;
    }
  }
  
  my @comps = split(/,/, $loc_string);
  my @local_components = grep { $_ =~ /^\d+\.\.\d+$/ } @comps;

  my (%remote_len);
  my $locallen = 0;
  foreach (@comps) {
    /^(\d+)\.\.(\d+)/ and do {
      $locallen += ($2 - $1 + 1);
      next;
    };
    /^(\S+):(\d+)\.\.(\d+)/ and do {
      $remote_len{$1} += ($3 - $2 + 1);
      next;
    };
  }
  return ($locallen, \%remote_len);
}

##########################
sub get_multi_gene_loci {
  my @hashes = @_;

  my (%locus2gene, %multi_gene_loci);
  foreach my $h (@hashes) {
    foreach my $iso (keys %$h) {
      my $locus_id = $iso;
      if ($iso =~ /^(\S+\d)[a-z]$/) {
        $locus_id = $1;
      }

      $locus2gene{$locus_id}->{$h->{$iso}} = 1;
    }
  }

  foreach my $locus_id (keys %locus2gene) {
    if (scalar(keys %{$locus2gene{$locus_id}}) > 1) {
      # pathological case: multiple gene at same
      # "locus". In these cases, we could choose to
      # (a) omit the locus_tag, or 
      # (b) use the locus_id as the gene_name (so that it
      # matches /locus_tag) and pass through the gene 
      # name as a note. Pass back a structure that allows
      # both of these possibilities. 
      $multi_gene_loci{$locus_id} = 1;
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
        "RA   ;",
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
       "RA   Paulini M;",
       "RT   ;",
       "RL   Submitted (12-Dec-2012) to the INSDC.",
       "RL   WormBase Group, European Bioinformatics Institute,",
       "RL   Cambridge, CB10 1SA, UNITED KINGDOM.",
       ],
    ]

    );

  return $primary_references{$species};
}

#######################
sub fetch_cds2status {
  my ($ref_db)= @_;

  my %cds2status;

  if ($quicktest) {
    return \%cds2status;
  }

  if (-d $ref_db) {

    $log->write_to("You are using $ref_db to get CDS status\n\n");
    
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
    $log->write_to("You are using $ref_db to get CDS protein_id\n\n");
    
    my $query = &get_protein_id_query();
    my $command = "Table-maker -p $query\nquit\n";
    open(TACE, "echo '$command' | $tace $ref_db |");
    while(<TACE>) {
      chomp;
    s/\"//g;
      if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)$/) {
        my ($cds, $clone, $protein_id, $version) = ($1, $2, $3, $4);
        $cds2proteinid{$cds}->{$clone} = "${protein_id}.${version}";
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
sub fetch_cds2dbremark {
  my ($ref_db) = @_;

  my %cds2dbremark;

  if ($quicktest) {
    return \%cds2dbremark;
  }

  if (-d $ref_db) {
    $log->write_to("You are using $ref_db to get CDS/RNA db_remarks\n\n");
    
    my @qfiles;
    foreach my $class ("CDS", "RNA") {
      my $query = &get_db_remark_query($class);
      push @qfiles, $query;
      my $command = "Table-maker -p $query\nquit\n";
      open(TACE, "echo '$command' | $tace $ref_db |");
      while(<TACE>) {
        chomp;
        if (/^\"(\S+)\"\s+\"(.+)\"$/) {
          my ($obj, $dbremark) = ($1, $2);
          $dbremark =~ s/\\//g;
          
          $cds2dbremark{$obj} = $dbremark;
        }
      }
    }
    
    unlink @qfiles;
  } elsif (lc($ref_db) eq 'none') {
    $log->write_to("You have specifed 'none' for the cds2dbremark db, so will not be added\n");
  } else {
    $log->log_and_die("Could not find cds2dbremark database '$ref_db'");
  }

  return \%cds2dbremark;

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
Class ${species}_${class}
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
Class ${species}_CDS
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
  return $tmdef;
}
 
###################################
sub stage_dump_to_submissions_repository {
  my ($dump_file) = @_;

  my (%records, $current_lines, $cosmid);

  my $submit_repo = $wormbase->submit_repos;

  open(my $dfg, $dump_file) or $log->log_and_die("Could not open $dump_file for reading\n");

  while(<$dfg>) {
    /^ID/ and do {
      $current_lines = [];
    };

    /^DE\s+/ and do {
      if (/^DE\s+.+\s+(\S+)$/) {
        $cosmid = $1;
      } else {
        $log->log_and_die("Could not parse contig name from DE line: $_\n");
      }

      foreach my $line (@$current_lines) {
        push @{$records{$cosmid}->{embl}}, $line;
      }
      $current_lines = $records{$cosmid}->{embl};
    };
    
    if (/^\s+(.+)$/) {
      my $seq = $1;
      $seq =~ s/\s+//g;
      $records{$cosmid}->{seq} .= uc($seq);
    }
    
    push @$current_lines, $_;
  }
  
  # only do sequences for C.elegans, because we do not do the
  # submission of primary sequence records for other species
  if ($species eq 'elegans') {
    foreach my $cosmid (sort keys %records) {
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
  }
  
  foreach my $cosmid (sort keys %records) {
     my $hash = $wormbase->submit_repos_hash_from_sequence_name($cosmid);
      
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



__END__

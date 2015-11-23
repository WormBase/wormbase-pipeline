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
  #
  # STRAIN
  #
  strain => {
    elegans  => 'Bristol N2',
    briggsae => 'AF16',
    brugia   => 'FR3',
    sratti   => 'ED321',
  },
  
  #
  # REFERENCES
  #
  references => {
    elegans =>  [
      [
        "RA   Sulson J.E., Waterston R.;",
        "RL   Nematode Sequencing Project: Sanger Institute, Hinxton, Cambridge",
        "RL   CB10 1SA, UK and The Genome Institute at Washington University,", 
        "RL   St. Louis, MO 63110, USA.",
      ],
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
      [
        "RX   PUBMED; 17885136.",
        "RA   Ghedin E., Wang S., Spiro D., Caler E., Zhao Q., Crabtree J., Allen J.E.,",
        "RA   Delcher A.L., Guiliano D.B., Miranda-Saavedra D., Angiuoli S.V., Creasy T.,",
        "RA   Amedeo P., Haas B., El-Sayed N.M., Wortman J.R., Feldblyum T., Tallon L.,",
        "RA   Schatz M., Shumway M., Koo H., Salzberg S.L., Schobel S., Pertea M.,", 
        "RA   Pop M., White O., Barton G.J., Carlow C.K., Crawford M.J., Daub J.,",
        "RA   Dimmic M.W., Estes C.F., Foster J.M., Ganatra M., Gregory W.F.,",
        "RA   Johnson N.M., Jin J., Komuniecki R., Korf I., Kumar S., Laney S., Li B.W.,",
        "RA   Li W., Lindblom T.H., Lustigman S., Ma D., Maina C.V., Martin D.M.,",
        "RA   McCarter J.P., McReynolds L., Mitreva M., Nutman T.B., Parkinson J.,",
        "RA   Peregrin-Alvarez J.M., Poole C., Ren Q., Saunders L., Sluder A.E.,",
        "RA   Smith K., Stanke M., Unnasch T.R., Ware J., Wei A.D., Weil G.,",
        "RA   Williams D.J., Zhang Y., Williams S.A., Fraser-Liggett C., Slatko B.,",
        "RA   Blaxter M.L., Scott A.L.;",
        "RT   \"Draft genome of the filarial nematode parasite Brugia malayi.\";",
        "RL   Science. 2007 Sep 21;317(5845):1756-60.",
      ],
    ],

    sratti => [
      [
        "RA   Aslett M.;",
        "RT   ;",
        "RL   Submitted (16-SEP-2014) to the INSDC.",
        "RL   Pathogen Sequencing Unit, Wellcome Trust Sanger Institute, Wellcome Trust",
        "RL   Genome Campus, Hinxton, Cambridge, Cambridgeshire CB10 1SA, UNITED KINGDOM.",
      ],
      [
        "RA   Coghlan A., Holroyd N., Foth B.J., Tsai I.J., Kikiuchi T., Tracey A.,",
        "RA   Nichol S., Brooks K., Beasley H., Stanley E.J., Sanchez-Flores A.,",
        "RA   Cotton J.A., Bennett H.M., Gordon D., Ribeiro D., Protasio A.V., Hunt V.,",
        "RA   Randle N., Kulkarni A., Wastling J.A., Lok J.B., Stoltzfus J.B., Doyle S.R.,",
        "RA   Grant W.N., Streit A., Viney M., Berriman M.;", 
        "RT   \"Sequencing and analysis of the genomes of the parasitic nematode",
        "RT   Strongyloides ratti and its close relatives for drug and vaccine discovery\"",
        "RL   ",
      ],
    ],

    ovolvulus => [
      [
        "RA   Aslett M.;",
        "RT   ;",
        "RL   Submitted (04-NOV-2013) to the INSDC.",
        "RL   Pathogen Sequencing Unit, Wellcome Trust Sanger Institute, Wellcome Trust",
        "RL   Genome Campus, Hinxton, Cambridge, Cambridgeshire CB10 1SA, UNITED KINGDOM.",
      ],
      [

        "RA   Cotton J., Tsai J., Stanley E., Tracey A., Holroyd N., Lustigman S.,",
        "RA   Berriman M.;",
        "RT   \"Genome sequencing of Onchocerca volvulus\"",
        "RL   Unpublished.",
      ],
    ]
  }    
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
    $cds2proteinid_db, $cds2status_db, $trans2dbremark_db, $decorations_db,
    $cds2status_h, $cds2proteinid_h,  $trans2dbremark_h, $trans2briefid_h, $trans2type_h, $gclass2desc_h,
    $agp_segs,
    $no_annotation,
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
  "modprivate=s"    => \$private,
  "moddumpfile:s"   => \$mod_dump_file,
  "rawdumpfile:s"   => \$raw_dump_file,
  "quicktest"       => \$quicktest,
  "sequencelevel"   => \$sequencelevel,
  "noannotation"    => \$no_annotation,
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

if ($species eq 'sratti') {
  $sequencelevel = 1;
} elsif ($species ne "elegans") {
  $sequencelevel = 0;
}

$decorations_db =  $wormbase->database('current') if not defined $decorations_db;
$cds2proteinid_db = $decorations_db if not defined $cds2proteinid_db;
$cds2status_db = $decorations_db if not defined $cds2status_db;
$trans2dbremark_db = $decorations_db if not defined $trans2dbremark_db;

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

  if (not $no_annotation) {
    $multi_gene_loci = &get_multi_gene_loci($gene2gseqname);
    
    $trans2type_h = &fetch_transcript2type($dbdir);
    $gclass2desc_h = &fetch_geneclass2desc($dbdir);
    $trans2briefid_h = &fetch_transcript2briefid($dbdir);
    $seleno_prots = &fetch_selenoproteins($dbdir);

    # Some information is not yet available in autoace (too early in the build)
    # By default, pull cds2status, trans2dbremark and protein_ids from current_DB,
    # however, allow command-line override when necessary
    $cds2status_h = &fetch_cds2status($cds2status_db);
    $cds2proteinid_h = &fetch_cds2proteinid($cds2proteinid_db);
    $trans2dbremark_h = &fetch_trans2dbremark($trans2dbremark_db);
   
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
  
  my ($seqname, $seqlen, $idline_suffix, @accs, @features, %source_quals, $written_header, %written_version_tag);
  
  while (<$raw_fh>) {
    
    # Store the necessary default ID line elements ready for use in the new style EMBL ID lines.
    if(/^ID\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\d+)\s+\S+/){
      ($seqname, $seqlen) = ($1, $2);
      $idline_suffix = sprintf("SV XXX; linear; genomic DNA; %s; INV; $seqlen BP.", ($sequencelevel) ? "STD" : "CON");
      @accs = ();
      @features = ();
      %source_quals = ();
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
      } else {
        print $out_fh "AC * $accs[0]\n";
      }
      print $out_fh "XX\n";
      
      #
      # PR
      #
      print $out_fh "PR   Project:$bioproject\n";
      print $out_fh "XX\n";
      
      #
      # DE
      #
      my (@de_line);
      
      if ($species eq 'elegans') {
        if ($seqname =~ /CHROMOSOME_(\S+)/) {
          @de_line = ("$full_species_name chromosome $1");
          $source_quals{chromosome} = $1;
        } elsif (!defined($clone2type->{$seqname})){
          @de_line = ("$full_species_name clone $seqname");
          $source_quals{clone} = $seqname;
        } elsif (lc($clone2type->{$seqname}) eq "other") {
          @de_line = ("$full_species_name clone $seqname");
          $source_quals{clone} = $seqname;
        } elsif (lc($clone2type->{$seqname}) eq "yac") {
          @de_line = ("$full_species_name YAC $seqname");
          $source_quals{clone} = $seqname;
        } else {
          @de_line = ("$full_species_name $clone2type->{$seqname} $seqname");
          $source_quals{clone} = $seqname;
        }
      } elsif ($species eq 'briggsae') {
        @de_line = ("$full_species_name AF16 supercontig from assembly CB4, $seqname");
        $source_quals{note} = "supercontig $seqname";
      } elsif ($species eq 'brugia'){
        @de_line = ("$full_species_name FR3 supercontig from assembly B_malayi-3.1 $seqname");
        $source_quals{note} = "supercontig $seqname";
      } elsif ($species eq 'sratti') {
        my $de_prefix = "$full_species_name assembly S_ratti_ED321";
        if ($seqname =~ /_chr(\S)_scaffold(\d+)/) {
          @de_line = ("$de_prefix, chr $1 scaffold $2");
          $source_quals{chromosome} = $1;
        } elsif ($seqname =~ /chr(\S)/) {
          @de_line = ("$de_prefix, chromosome $1");
          $source_quals{chromosome} = $1;
        } elsif ($seqname =~ /scaffold(\d+)/) {
          @de_line = ("$de_prefix, scaffold $1"); 
          $source_quals{note} = "supercontig $seqname";
        }
      }
    
      foreach my $de_line (@de_line) {
        print $out_fh "DE   $de_line\n";
      }
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
      
      printf $out_fh "RN   [%d]\n", $ref_count++;
      printf $out_fh "RP   1-%d\n", $seqlen;
      print $out_fh "RG   WormBase Consortium\n";
      print $out_fh "RA   WormBase;\n";
      print $out_fh "RT   ;\n";
      print $out_fh "RL   WormBase Group, European Bioinformatics Institute,\n";
      print $out_fh "RL   Cambridge, CB10 1SA, UK. Email: help\@wormbase.org\n";
      print $out_fh "XX\n";

      my @refs = @{$species_info{references}->{$species}};
      for(my $i=0; $i < @refs; $i++) {
        printf $out_fh "RN   [%d]\n", $ref_count++;
        printf $out_fh "RP   1-%d\n", $seqlen;
        map { print $out_fh "$_\n" } @{$refs[$i]};
        print $out_fh "XX\n";
      }

      next;
    }
    
    #
    # Comments
    #
    if (/^CC   / and not $no_annotation) {
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
      if ($ftype eq 'source') {
        foreach my $k (sort keys %source_quals) {
          push @{$features[-1]->{quals}}, ["/$k=\"$source_quals{$k}\""];
        }
      }

      next;
    } elsif (/^FT\s+(.+)/) {
      my $content = $1;
      next if $content eq '/pseudo';
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

  if ($no_annotation) {
    @feats = grep { $_->{ftype} eq 'source' } @feats;
  }

  foreach my $feat (@feats) {

    if ($feat->{ftype} eq 'source') {
      printf $out_fh "FT   %-16s%s\n", $feat->{ftype}, $feat->{location}->[0];
      printf $out_fh "FT   %16s/db_xref=\"taxon:%d\"\n", " ", $wormbase->ncbi_tax_id;
      printf $out_fh "FT   %16s/strain=\"%s\"\n", " ", $species_info{strain}->{$species};
      printf $out_fh "FT   %16s/mol_type=\"genomic DNA\"\n", " ";
      foreach my $tag (@{$feat->{quals}}) {
        foreach my $ln (@$tag) {
          printf $out_fh "FT   %16s%s\n", " ", $ln;
        }
      }
      next;
    } elsif ($feat->{ftype} =~ /RNA$/) {
      my $mod_dir = $feat->{ftype};
      my ($rna_class);

      if ($mod_dir eq 'antisense_RNA' or
          $mod_dir eq 'miRNA' or
          $mod_dir eq 'piRNA' or
          $mod_dir eq 'scRNA' or
          $mod_dir eq 'siRNA' or
          $mod_dir eq 'snoRNA' or
          $mod_dir eq 'snRNA') {
        $rna_class = $mod_dir;
        $mod_dir = 'ncRNA';
      } elsif ($mod_dir eq 'lincRNA') {
        $rna_class = "lncRNA";
        $mod_dir = "ncRNA";
      } elsif ($mod_dir eq 'ncRNA') {
        $rna_class = 'other';        
      } elsif ($mod_dir eq 'tRNA' or  
               $mod_dir eq 'rRNA' or
               $mod_dir eq 'misc_RNA' or
               $mod_dir eq 'precursor_RNA' or
               $mod_dir eq 'prim_transcript') {
        # no class
      } else {
        $log->log_and_die("Unexpected feature type " . $feat->{ftype});
      }

      if (defined $rna_class) {
        $feat->{rna_class} = $rna_class;
        push @{$feat->{quals}}, ["/ncRNA_class=\"$rna_class\""];
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
    my ($prod_name, $class_based_name, @prod_notes);
    if ($feat->{ftype} eq 'CDS') {
      if (exists $trans2briefid_h->{$wb_isoform_name}) {
        $prod_name = $trans2briefid_h->{$wb_isoform_name};
      } elsif (exists $gene2cgcname->{$wb_geneid}) {
        my $gclass = $gene2cgcname->{$wb_geneid};
        $gclass =~ s/-[^-]+$//;
        if (exists $gclass2desc_h->{$gclass}) {
          $prod_name = $gclass2desc_h->{$gclass};
          $prod_name =~ s/\\//g; 
          if ($feat->{is_pseudo}) {
            $prod_name .= " pseudogene";
          }
          $class_based_name = $gclass;
        }
      } 
      if ($feat->{is_pseudo}) {
        if (defined $prod_name) {
        # no product allowed for pseudogenes; add them as a note
          push @prod_notes, $prod_name;
          undef $prod_name;
        }
      } else {
        if (defined $prod_name) {
          if ($class_based_name) {
            push @prod_notes, "Product from WormBase gene class $class_based_name";
          }
        } else {
          $prod_name = "Uncharacterized protein";
        }
      }
      #if ($trans2dbremark_h->{$wb_isoform_name}) {
      #  my $rem = sprintf("%s", $trans2dbremark_h->{$wb_isoform_name});
      #  $rem =~ s/\s+/ /g;
      #  push @prod_notes, $rem;
      #}
    } elsif ($feat->{ftype} eq 'mRNA') {
      # unclear what product should be for mRNAs. Defer for now
    } elsif ($feat->{ftype} =~ /RNA/) {      
      # note that prim_transcript will not match this, but this is fine
      # because prim_transcript features are not allowed to have a product
      my $pname = ($gene2cgcname->{$wb_geneid}) ? $gene2cgcname->{$wb_geneid} : $wb_isoform_name;
      my $unclassified = 0;
      
      if ($feat->{ftype} eq 'misc_RNA') {
        $prod_name = "Non-coding transcript of protein-coding gene $pname";
      } elsif ($feat->{ftype} eq 'rRNA' or $feat->{ftype} eq 'tRNA') {
        if ($trans2briefid_h->{$wb_isoform_name}) {
          $prod_name = $trans2briefid_h->{$wb_isoform_name};
        } else {
          $prod_name = $feat->{ftype};
        }
        if ($feat->{is_pseudo}) {
          $prod_name .= " pseudogene";
          push @prod_notes, $prod_name;
          undef $prod_name;
        }        
      } elsif ($feat->{ftype} eq 'precursor_RNA') {
        $prod_name = "pre-microRNA $pname";
      } elsif (exists $feat->{rna_class}) {
        my $rclass = $feat->{rna_class};
        if ($rclass eq 'miRNA') {
          if ($wb_isoform_name =~ /([ab])$/) {
            $pname .= ",$1";
          }
          $prod_name = "microRNA $pname";
        } elsif ($rclass  eq 'antisense_RNA') {
          $prod_name = "antisense RNA $pname";
        } elsif ($rclass eq 'piRNA') {
          $prod_name = "piwi-interacting RNA $pname";
        } elsif ($rclass eq 'scRNA') {
          $prod_name = "small cytoplasmic RNA $pname";
        } elsif ($rclass eq 'siRNA') {
          $prod_name = "small interfering RNA $pname";
        } elsif ($rclass eq 'snoRNA') {
          $prod_name = "small nucleolar RNA $pname";
        } elsif ($rclass eq 'snRNA') {
          $prod_name = "small nuclear RNA $pname";
        } elsif ($rclass eq 'lincRNA') {
          $prod_name = "long non-coding RNA $pname";
        } else {
          $prod_name = "Unclassified non-coding RNA $pname";
          $unclassified = 1;
        }
      } else {
        $prod_name = "Unclassified non-coding RNA $pname";
        $unclassified = 1;
      }
      
      # The following will work once Brief_identification and
      # RNA type remarks have been cleaned up. Until then, disabled. 
      if ($unclassified and 0) {
        if (exists $trans2type_h->{$wb_isoform_name} and 
            exists $trans2type_h->{$wb_isoform_name}->{note} and
            $trans2type_h->{$wb_isoform_name}->{note} ) {
          my $type_rna_note = $trans2type_h->{$wb_isoform_name}->{note};
          push @{$revised_quals{note}}, ["/note=\"$type_rna_note\""];
          push @prod_notes, $type_rna_note;
        } elsif ($trans2briefid_h->{$wb_isoform_name} and 
                 $trans2briefid_h->{$wb_isoform_name} ne "ncRNA") {
          my $bid_rna_note = $trans2briefid_h->{$wb_isoform_name};
          push @prod_notes, $bid_rna_note;
        }
      }

    }
    $product_qual = [];
    if ($prod_name) {
      my $prod_line = "/product=\"$prod_name\"";
      my @wl =  split(/\n/, wrap('','',$prod_line));
      push @$product_qual, @wl;
    }
    if (@prod_notes) {
      foreach my $note (@prod_notes) {
        my $note_line = "/note=\"$note\"";
        my @wl2 = split(/\n/, wrap('','',$note_line));
        
        push @$product_qual, @wl2;
      }
    }

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
    
    push @$current_lines, $_;
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
  my ($ref_db) = @_;

  my (%gene2cgc, %gene2gseqname);
  
  if (not $quicktest and -d $ref_db) {
    foreach my $aref (&_run_table_maker($ref_db, &get_gene_name_query() )) {
      my ($gene, $cgc_name, $gseqname) = @$aref;

      $log->log_and_die("Could not get gseqname for $gene\n") if not $gseqname;
      
      if (defined $gene) {
        $gene2cgc{$gene} = $cgc_name if $cgc_name;
        $gene2gseqname{$gene} = $gseqname if $gseqname;
      }
    }
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
    foreach my $aref (&_run_table_maker($ref_db, &get_status_query())) {
      my ($cds, $status) = @$aref;
      $cds2status{$cds} = $status;
    }
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
    foreach my $aref (&_run_table_maker($ref_db,&get_protein_id_query() )) {
      my ($cds, $seqname, $protein_id, $version) = @$aref;
      $cds2proteinid{$cds}->{$seqname} = "${protein_id}.${version}";
    }
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
    foreach my $class ("CDS", "Transcript") {
      foreach my $aref (&_run_table_maker($ref_db, &get_db_remark_query($class))) {
        my ($obj, $dbremark) = @$aref;
        $dbremark =~ s/\\//g;
        $trans2dbremark{$obj} = $dbremark;
      }
    }
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
    foreach my $aref (&_run_table_maker($ref_db, &get_ncrna_type_query())) {
      my ($obj, $type, $other_comment) = @$aref;

      if ($type) {
        $tran2type{$obj}->{type} = $type;
        if ($other_comment) {
          $other_comment =~ s/\\//g;
          $tran2type{$obj}->{note} = $other_comment;
        }
      }
    }
  } 

  return \%tran2type;

}

########################
sub fetch_transcript2briefid {
  my ($ref_db) = @_;

  my (%res, %tran2briefid);

  if (not $quicktest and -d $ref_db) {
    foreach my $class ("CDS", "Transcript", "Pseudogene") {
      foreach my $aref (&_run_table_maker($ref_db, &get_brief_id_query($class) )) {
        my ($obj, $brief_id, $evi_database) = @$aref;

        if (not $evi_database) {
          $evi_database = "no_evidence";
        }
        
        $brief_id =~ s/C\. elegans //; 
        $brief_id =~ s/\\//g; 
        
        next if $brief_id =~ /^Uncharacterized protein/;
        
        $res{$class}->{$obj}->{lc($evi_database)} = $brief_id;
      }
    }
  }

  foreach my $class (keys %res) {    
    foreach my $obj (keys %{$res{$class}}) {
      if ($class eq 'CDS') {
        # only include the UniProt ones for CDS for time being;
        # InterPro and no-evidence ones have questionable source
        if ($species eq 'elegans') {
          if (exists $res{$class}->{$obj}->{uniprot}) {
            $tran2briefid{$obj} = $res{$class}->{$obj}->{uniprot};
          }
        } else {
          foreach my $evi (keys %{$res{$class}->{$obj}}) {
            $tran2briefid{$obj} = $res{$class}->{$obj}->{$evi};
          }
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
sub fetch_geneclass2desc {
  my ($ref_db) = @_;

  #
  # These are ones that we do not want to use as product names, 
  # because they are uninformative (e.g. tag) or based mutant-phenotype
  # (e.g. unc). There will be others, but these are the big hitters in
  # terms of gen numbers. 
  #
  my @blacklist = (
    "age",
    "tag",
    "unc",
    "eat",
    "egl",
    "lin",
    "pfk",
    "pfkb",
      );

  my %class2desc;

  if (not $quicktest and -d $ref_db) {
    foreach my $aref (&_run_table_maker($ref_db, &get_gene_class_query())) {
      my ($gclass, $desc) = @$aref;

      next if grep { $gclass eq $_ } @blacklist;
      next if $desc =~ /abnormal/i;
      next if $desc =~ /defective/i;
      next if $desc =~ /sensitiv/i;
      next if $desc =~ /lethal/i;
      next if $desc =~ /resist/i;

      $class2desc{$gclass} = $desc;
    }
  }

  return \%class2desc;

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

  my $condition = "Condition Method";
  if ($single) {
    $condition .= " AND Sequence = \"$single\""
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

  my $tmdef = "/tmp/cdsstatus_tmquery.$$.def";
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


  my $condition = "Condition Species = \"$full_species_name\"";
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


  my $condition = "Condition Live AND Species = \"$full_species_name\" AND Sequence_name";
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

####################################
sub get_gene_class_query {

  my $tmdef = "/tmp/gclass_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  


  my $gclass_query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene_class
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Text 
From 1 
Tag Description
 
EOF
 
  print $qfh $gclass_query;
  close($qfh);

  return $tmdef;

}


##############################################
sub _run_table_maker {
  my ($ref_db, $def_file) = @_;

  $log->write_to("Running table-maker on $ref_db with $def_file\n");

  my $command = "Table-maker -p $def_file\nquit\n";
  open(my $tacefh, "echo '$command' | $tace $ref_db |");

  my @res;
  while(<$tacefh>) {
    chomp;
    s/\"//g;
    next if /^acedb/;
    next if /^\/\//;
    next if not /\S/;
    
    my @line = split(/\t/, $_);

    push @res, \@line;
  }

  unlink $def_file;
  return @res;
}


__END__

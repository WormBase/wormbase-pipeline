#!/usr/bin/env perl
#
# map_Oligo_sets
#
#
#############################################################################################

use lib $ENV{'CVS_DIR'};
use strict;
use Wormbase;
use Getopt::Long;
use warnings;
use Modules::Map_Helper;


my $config = {
  CDS => [ ['curated', 'curated', 'exon'] ],

  Transcript => [ ['Coding_transcript',     'Coding_transcript',     'exon'],
                  ['Non_coding_transcript', 'non_coding_transcript', 'exon'],
                  ['ncRNA',                 'ncRNA',                 'exon'],
                  ['Pseudogene',            'Pseudogene',            'exon'] ],
  
  Pseudogene => [ ['Pseudogene', 'Pseudogene', 'exon'] ],
};


my $help;       # Help perldoc
my $test;       # Test mode
my $debug;      # Debug mode, output only goes to one user
my $acefile;    # specify a custom acefile
my $store;      # specify a frozen configuration file
my $noload;   # don't parse the acefile
my $species;

GetOptions(
  'debug=s'   => \$debug,
  'test'      => \$test,
  'help'      => \$help,
  'acefile=s' => \$acefile,
  'store=s'   => \$store,
  'noload'    => \$noload,
  'species=s' => \$species,
);

my @oligo_set_file_prefixes = ('Oligo_set', 'Oligo_set_mapping');

############################
# recreate configuration   #
# ##########################
my $wb;
if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
} else { 
  $wb = Wormbase->new( -debug    => $debug, 
                       -test     => $test,
                       -organism => $species);
}

my $gffdir  = $wb->gff_splits;
if (not defined $acefile) {
  $acefile = $wb->acefiles . "/Oligo_sets_to_genes.ace";
}

my $log = Log_files->make_build_log($wb);


my @oligo_set_files;
if ($wb->assembly_type eq 'contig') {
  foreach my $oligo_set_file_prefix (@oligo_set_file_prefixes) {
    my $file = "$gffdir/${oligo_set_file_prefix}.gff";
    if (not -e $file) {
      $log->write_to("Could not find file $file - assuming not present for species\n");
    } else {
      push @oligo_set_files, $file;
    }
  }
} else {
  foreach my $oligo_set_file_prefix (@oligo_set_file_prefixes) {
    foreach my $chr_name ($wb->get_chromosome_names(-mito => 1, -prefix => 1)) {
      my $file = "$gffdir/${chr_name}_${oligo_set_file_prefix}.gff";
      if (not -e $file) {
        $log->write_to("Could not find file $file - assuming not present for species\n");
      } else {
        push @oligo_set_files, $file;
      }
    }
  }
}

if (@oligo_set_files) {
  open(my $acefh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

  foreach my $class (keys %$config) {
    $log->write_to("Mapping to $class...\n");
    
    my $fm = Map_Helper->new();
    
    foreach my $stanza (@{$config->{$class}}) {
      my ($file_prefix, $method, $type) = @$stanza;
      
      my @files;
      if ($wb->assembly_type eq 'contig') {
        my $file = $gffdir . "/${file_prefix}.gff";
        push @files, $file;
      } else {
        foreach my $chr_name ($wb->get_chromosome_names(-mito => 1, -prefix => 1)) {
          my $file = $gffdir . "/${chr_name}_${file_prefix}.gff";
          push @files, $file;
        }
      }
      
      foreach my $file (@files) {
        if (not -e $file) {
          warn("Could not find file $file - assuming this datatype does not exist for $species\n");
          next;
        }
        $log->write_to(" Reading $file...\n");
        $fm->populate_from_GFF($file, $method, $type, sprintf('%s \"(\S+)\"', $class));
      }
    }
    
    $log->write_to(" Building index...\n");
    $fm->build_index;
    
    my %object_map;
    
    foreach my $file (@oligo_set_files) {
      $log->write_to(" Mapping $file to $class...\n");
      open( my $fh, $file) or $log->log_and_die("Could not open $file for reading\n");
      while(<$fh>) {
        next if /^\#/;
        my @l = split(/\t+/, $_);
        
        next if $l[1] ne 'Oligo_set' or $l[2] ne 'reagent';
        my ($name) = $l[8] =~ /Oligo_set \"(\S+)\"/;
        if (not defined $name) {
          ($name) = $l[8] =~ /Target \"Oligo_set:(\S+)\"/;
        }
        
        next if not defined $name;
        
        my @matches = @{$fm->search_feature_segments($l[0], $l[3], $l[4], $l[6])};
        map { $object_map{$_}->{$name} = 1 } @matches;
      }
    }
    
    foreach my $obj_name (sort keys %object_map) {
      print $acefh "\n$class : \"$obj_name\"\n";
      foreach my $os (sort keys %{$object_map{$obj_name}}) {
        print $acefh "Corresponding_oligo_set $os\n";
      }
    }
  }
  
  close($acefh) or $log->log_and_die("Did not cleanly close output file $acefile\n");
  
  if (not $noload) {
    $wb->load_to_database( $wb->autoace, $acefile, 'map_Oligo_set', $log );
  }
} else {
  $log->write_to("Could not find any oligo_set files to map - doing nothing\n");
}

###############
# hasta luego #
###############

$log->mail();

exit(0);


__END__


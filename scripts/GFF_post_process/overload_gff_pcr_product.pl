#!/usr/bin/env perl
#
# overload_gff_pcr_product.pl
#
# Overloads PCR_product mapping lines with extra info (Lab of clone etc)
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-07-22 11:37:35 $      

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;
use Ace;

my ($debug,$test,$species,$store,$wormbase);
my ($infile,$outfile,$gff3, $changed_lines);

GetOptions(
  'debug=s'   => \$debug,
  'test'      => \$test,
  'species:s' => \$species,
  'store:s'   => \$store,
  'infile:s'  => \$infile,
  'outfile:s' => \$outfile,
  'gff3'      => \$gff3,
)||die(@!);

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

my $clone_info = &get_clone_info();
my $amplified = &get_amplified_info();

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while ( <$gff_in_fh> ) {
  unless(/(Orfeome|GenePair_STS)\s+PCR_product/){
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @l = split(/\t/, $_);

  my ($name);
  if ($l[8] =~ /PCR_product \"(\S+)\"/) {
    $name = $1;
  } elsif ($l[8] =~ /Name=PCR_product:(\S+)/) {
    $name = $1;
  } else {
    $log->log_and_die("Could not extract product name from line @l\n");
  }

  my $old_attr = $l[8];
  if ($gff3) {
    $l[8] .=
    $l[8] .= ";VendorID=$clone_info->{$name}" if exists $clone_info->{$name};
  } else {
    $l[8] .= " ; Amplified \"$amplified->{$name}\"" if exists $amplified->{$name};
    $l[8] .= " ; VendorID \"$clone_info->{$name}\"" if exists $clone_info->{$name};
  }
  $changed_lines++ if $l[8] ne $old_attr;
  print $gff_out_fh join("\t", @l), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);

#################################
sub get_amplified_info {
  my %amplified;

  my $def_file = &write_amplified_def_file();
  my $table = $wormbase->table_maker_query($wormbase->autoace, $def_file);
  while(<$table>) {
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    chomp;
    my ($pcr_prod, $amplified) = split(/\t/,$_);
    if (defined $pcr_prod) {
      $amplified{$pcr_prod} = $amplified;
    }
  }
  unlink $def_file;

  return \%amplified;
}

#################################
sub get_clone_info {
  my %clone_acc;

  my $def_file = &write_clone_def_file();
  my $table = $wormbase->table_maker_query($wormbase->autoace, $def_file);
  while(<$table>) {
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    chomp;
    my ($pcr_prod, $db_key, $db_val) = split(/\t/,$_);
    if (defined $pcr_prod and $db_key eq 'HGMP_location') {
      $clone_acc{$pcr_prod} = $db_val;
    }
  }
  unlink $def_file;

  return \%clone_acc;
}


#################################
sub write_amplified_def_file {
  my $def = '/tmp/overload_gff_pcr_product.def';
  open my $tmpfh,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class PCR_product 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Integer 
From 1 
Tag Amplified 
  
// End of these definitions
END

  print $tmpfh $txt;
  close($tmpfh);
  return $def;
}



#####################
sub write_clone_def_file {
  my $def = '/tmp/overload_gff_pcr_product.def';
  open my $tmpfh,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class PCR_product 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class Clone 
From 1 
Tag Clone 
 
Colonne 3 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 2 
Tag Database 
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Database_field 
Right_of 3 
Tag  HERE  
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 4 
Tag  HERE  
 
// End of these definitions
END

  print $tmpfh $txt;
  close($tmpfh);
  return $def;
}



1;

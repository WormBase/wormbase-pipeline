#!/usr/bin/env perl
#
# overload_gff_tf.pl
# 
# overloads GFF TF_binding_site lines with info about the TF itself
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-08-27 21:50:10 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ( $debug, $test, $verbose, $store, $wormbase, $species, $database );
my ( $infile, $outfile, $gff3, %TF, $changed_lines);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "gff3"       => \$gff3,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "database:s" => \$database,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$database = $wormbase->autoace if not defined $database;
$species = $wormbase->species;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

my $table = $wormbase->table_maker_query($database, &write_def_file);
while(<$table>) {
  s/\"//g; #"
  next if (/acedb/ or /\/\//);
  chomp;
  my ($feature, $tf_id, $tf_name) = split(/\t/,$_);
  if (defined $feature) {
    $TF{$feature}->{id} = $tf_id;
    $TF{$feature}->{name} = $tf_name;
  }
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while( <$gff_in_fh> ) {
  /^\#/ and do {
    print $gff_out_fh $_;
    next;
  };
  chomp;
  my @f = split(/\t+/, $_);
  
  my ($fid);
  if ($f[8] =~ /Feature\s+\"(\S+)\"/) {
    $fid = $1;
  } elsif ($f[8] =~ /Feature:([^;]+)/) {
    $fid = $1;
  }
  
  #CHROMOSOME_V    TF_binding_site  TF_binding_site     155950  155951  .       +       .       Feature "WBsf046537"
  if (defined $fid and exists $TF{$fid}) {
    my $id = $TF{$fid}{id};
    my $name = $TF{$fid}{name};
    if ($gff3) {
      $f[8] .= ";tf_id=$id;tf_name=$name";
    } else {
      $f[8] .= " ; TF_ID \"$id\" ; TF_name \"$name\"";
    }
    $changed_lines++;
  }
  print $gff_out_fh join("\t", @f), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);


##########################################

sub write_def_file {
	my $def = '/tmp/overload_TF_GFF.def';
	open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
	my $txt = <<END;
// Spread sheet definition for the ACeDB software 
// User: gw3
// Date: 2012-07-13_14:29:33

// %n (%%n in the graphic) are parameter to be given on the command line in tace
// or by default by the Parameters given in this file
// \%n (%n in the graphic) are substituted by the value of column n at run time
// Line starting with // are ignored, starting with # are comments

Sortcolumn 1

Colonne 1 
Subtitle Feature 
Width 12 
Optional 
Visible 
Class 
Class Feature 
From 1 
Condition Associated_with_transcription_factor
 
Colonne 2 
Subtitle TF 
Width 40 
Optional 
Visible 
Class 
Class Transcription_factor 
From 1 
Tag Associated_with_transcription_factor 
 
Colonne 3 
Subtitle Name 
Width 40 
Optional 
Visible 
Text 
From 2 
Tag Name 
 
 

// End of these definitions
END

	print TMP $txt;
	close TMP;
	return $def;
}

1;

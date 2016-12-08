#!/usr/bin/env perl
#
# map_microarray.pl
#
# Add information to Microarray_results objects based on overlaps in GFF files 
#
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2013-10-14 10:16:24 $        

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

######################################
# variables and command-line options #
######################################
my ($test, 
    $debug,
    $noload,     
    $store,
    $species,
    $outfile,
    $wb);

GetOptions ("debug=s"   => \$debug,
	    "test"      => \$test,
	    "noload"    => \$noload,
            "species=s" => \$species,
	    "acefile=s"	=> \$outfile,
	    'store=s'	=> \$store
    );


############################
# recreate configuration   #
############################
if ($store){
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n");
} else {
  $wb = Wormbase->new(-debug    => $debug,
                      -test     => $test,
                      -organism => $species,
      );
}

my $log = Log_files->make_build_log($wb);

my $dbdir= $wb->autoace; 
$outfile = $outfile ? $outfile : $wb->acefiles."/microarray_mappings.ace";

open (my $out_fh, ">$outfile") or $log->log_and_die("Can't open the output file $outfile\n");

###########################################
# PCR_products and SMD_microarray results #
###########################################

foreach my $class ("PCR_product", "Oligo_set") {
  my %data;

  foreach my $pair (["CDS", "Overlaps_CDS"], 
                    ["Pseudogene", "Overlaps_pseudogene"],
                    ["Transcript", "Overlaps_transcript"]) {
    my ($tclass, $overlap_q) = @$pair;

    $log->write_to("Making $tclass connections to Microarray_results objects based on $class\n");
                      
    my $tm_def = &get_tm_def( $class, $tclass, $overlap_q );

    my $tm_query = $wb->table_maker_query( $dbdir, $tm_def );
    while(<$tm_query>) {
      chomp;
      s/\"//g; 
      next if (/acedb/ or /\/\//);
      next if /^\s*$/;
      
      my ($obj1, $mic_obj, $t_obj, $g_obj) = split(/\t/, $_);
        
      $data{$mic_obj}->{$tclass}->{$t_obj} = 1;
      $data{$mic_obj}->{Gene}->{$g_obj} = 1;
    }

    unlink $tm_def;
  }

  foreach my $obj (sort keys %data) {
    print $out_fh "\nMicroarray_results : \"$obj\"\n";
    foreach my $tag (sort keys %{$data{$obj}}) {
      foreach my $oobj (sort keys %{$data{$obj}->{$tag}}) {
        print $out_fh "$tag $oobj\n";
      }
    }
  }

}

close($out_fh) or $log->log_and_die("Problem when closing output file (probably out of space)\n");


unless ($noload) {
  $log->write_to("Loading file to autoace\n");
  $wb->load_to_database($wb->autoace, $outfile, 'microarray_mappings', $log);
}

$log->mail();
exit(0);


############################################

sub get_tm_def {
  my ($qclass, $tclass, $overlapst) = @_;

  my $txt = <<"END";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class $qclass
From 1 

Colonne 2
Width 12
Mandatory
Visible
Class
Class Microarray_results
From 1
Tag Microarray_results
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class $tclass 
From 1 
Tag $overlapst 
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 3 
Tag Gene 

END
    
  my $def_file = "/tmp/tm_def.$qclass.$tclass.$overlapst.$$.def";
  open(my $fh, ">$def_file") or $log->log_and_die("Could not open $def_file for writing\n");
  print $fh $txt;
  close($fh) or $log->log_and_die("Could not close $def_file after writing\n");

  return $def_file;
}


__END__

=pod

=head2 NAME - map_microarray.pl

=head1 USAGE

=over 4

=item map_microarray.pl [-options]

=back

map_microarray.pl connects PCR_products and micro_array_results  objects to CDS/transcript / Pseudogene /gene objects
via the defined_by tag in the database.

mandatory arguments:

=over 4

=item none

=back

optional arguments:

=over 4

=item -debug username, Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=item -outace specify location for the output Ace (if you want to use it outside of the build)

=item -store specifiy a configuration file

=back

=cut


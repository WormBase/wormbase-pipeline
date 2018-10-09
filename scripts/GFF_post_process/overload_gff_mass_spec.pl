#!/usr/bin/env perl
#
# overload_gff_mass_spec.pl
# 
# Overloads mass_spec_genome lines with extra info (peptide match etc)
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-08-27 21:50:10 $      

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($help, $debug, $test, $store, $wormbase, $database, $species);
my ( $gff3,$infile,$outfile,$changed_lines);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "gff3"       => \$gff3,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "database:s" => \$database,
  "species:s"  => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -organism => $species,
			     );
}

$database = $wormbase->autoace if not defined $database;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

my $tace = $wormbase->tace;

# get the set of CDS, history versions and resulting protein IDs
my $matches = &get_protein_match_data();
my $protein_history_aref = &get_protein_history();   

my $db = Ace->connect (-path => $database, 
                       -program => $tace) || die "cannot connect to database at $wormbase->database('current')\n";

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  
  
while (<$gff_in_fh>) {
  if (/^\#/) {
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @f = split /\t/, $_;
  my $id;

  if ($f[1] eq 'mass_spec_genome') {
    # get the ID name
    if ($gff3) {
      ($id) = ($f[8] =~ /Target=(\S+)/);
    } else {
      ($id) = ($f[8] =~ /Target \"Mass_spec_peptide:(\S+)\"/);
    }
    
    if (exists $matches->{$id}) { # for this peptide id
      my %times_observed;
      my $proteins = "";
      my $cdss = "";
      
      foreach my $prot (keys %{ $matches->{$id} }) {
        if ($proteins ne "") {$proteins .= " ";}
        $proteins .= $prot;
        
        # get the CDS name for this protein
        my $prot_obj = $db->fetch(Protein => $prot);
        if (! defined $prot_obj) {
          $log->log_and_die( "Can't fetch Protein object for $prot\n");
        }
        my $cds = $prot_obj->Corresponding_CDS;
        if (! defined $cds) { 
          $cds = get_previous_wormpep_ids($prot, $protein_history_aref);
        }
        $prot_obj->DESTROY();
        if ($cdss ne "") {$cdss .= " ";}
        $cdss .= $cds;
        
        foreach my $experiment (keys %{ $matches->{$id}->{$prot} }) {
          $times_observed{$experiment} = 1;	# count the number of unique experiments that this peptide has been seen in
          }
      }
      
      my $times_observed = (keys %times_observed); # count the unique experiments
      
      if ($gff3) {
        $f[8] .= ";Note=$id";
        $f[8] .= ";protein_matches=$proteins";
        $f[8] .= ";cds_matches=$cdss";
        $f[8] .= ";times_observed=$times_observed";
      } else {
        $f[8] .= " ; Note \"$id\"";
        $f[8] .= " ; Protein_matches \"$proteins\"";
        $f[8] .= " ; CDS_matches \"$cdss\"";
        $f[8] .= " ; Times_observed \"$times_observed\"";
      }
      $changed_lines++;
    }
  }
  
  print $gff_out_fh join("\t", @f), "\n";
}
$db->close;
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing - $changed_lines lines modified\n");
$log->mail();
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub get_protein_match_data {
  my $ace_dir = $database;
  my $tace    = $wormbase->tace;

  my $cmd1 = "Query Find Mass_spec_peptide\nshow -a\nquit";
  
  my (%match, $id, $count);
  
  $log->write_to("Finding Mass_spec data\n");
  open (TACE, "echo '$cmd1' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Mass_spec_peptide\s+:\s+\"(\S+)\"/) {
      $id = $1;
    } elsif (/\"(.+)\"\s+Protein\s+\"(.+)\"/) {
      $match{$id}->{$2}->{$1} = 1;
    }
  }
  close TACE;
  
  return \%match;
}



##########################################
sub get_protein_history {
  my @protein_history;
                            
  my $release = $wormbase->get_wormbase_version;
  my $pepdir_base = $wormbase->peproot;
  my $pepdir = $wormbase->pepdir_prefix . "pep${release}";
  my $pephistory =  $wormbase->pepdir_prefix . "pep.history${release}";
  my $data_file = "$pepdir_base/$pepdir/$pephistory";
 
  open (HIST, "< $data_file") || die "Can't open $data_file\n";
  while (my $line = <HIST>) {
    chomp $line;
    my @f = split /\s+/, $line;
    # $cds_id, $wormpep_id, $version1, $version2 <- ($version2 is undef if current version)
    push @protein_history, [@f];
  }
  close(HIST);
 
  return \@protein_history;
}

##########################################
sub get_previous_wormpep_ids {
 
  my ($protein, $protein_history_aref) = @_;
 
  my @protein_history = @{$protein_history_aref};
 
  my @wormpep_ids;
  my @versions;

  foreach my $line (@protein_history) {
    my ($cds_id, $wormpep_id, $version1, $version2) = @{$line};
    if ($wormpep_id eq $protein) {
      return $cds_id;
    }
  }
 
  print STDERR "Still can't find the CDS for protein $protein_name\n";
  return "unknown";
}
 
1;

#!/usr/bin/env perl
#
# map_Interaction.pl
#
# Add information to Interaction objects via aceperl follows:
# checks genes connected to linked variations and adds them to the interaction if they don't exist yet.
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-03-18 10:16:06 $

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;

#####################################
# variables and command-line options #
######################################

my ($output, $test, $debug, $noload, $store, $acefile, $wb, $dbdir,$testInteraction);

GetOptions(
    'debug=s'   => \$debug,
    'test'      => \$test,
    'noload'    => \$noload,
    'acefile=s' => \$output,
    'store=s'   => \$store,
    'database=s' => \$dbdir,
    'interaction=s' = \$testInteraction,
)||die("invalid commandline optioon\n");

if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n");
}
else { 
  $wb = Wormbase->new( -debug => $debug, -test => $test ); 
}

$dbdir ||= $wb->autoace; # Database path
$output ||= $wb->acefiles.'/Interaction_connections.ace';

my $log = Log_files->make_build_log($wb);

# open a connection
my $db = Ace->connect( -path    => "$dbdir", -program => $wb->$tace)
    or $log->log_and_die("Could not successfully connect to Acedb $dbdir\n");

open(my $outfh, ">$output" ) || $log->log_and_die("Can't open output file $output\n");

my @interactions = $testInteraction ? $db->fetch(Interaction => $testInteraction) 
   : $db->fetch( -query => 'Find Interaction WHERE Variation_Interactor OR Unaffiliated_variation');

foreach my $interaction (@interactions) {
  my (%results, %current_genes, %delete_unaffiliated);

  # Interactor Interactor_overlapping_gene WBGeneXYZ Interactor_type Effector
  if ($interaction->Interactor_overlapping_Gene) {
    foreach my $g ($interaction->get('Interactor_overlapping_Gene')) {
      my @evis;
      foreach my $col ($g->col) { # foreach type tag
        my @evi_cmps = ($col->name);
        if ($col->right) {
          push @evi_cmps, $col->right->name; # push the type
        }
      
        my $evidence = join(" ", @evi_cmps);
        push @evis, $evidence;
      }
      $current_genes{$g} = \@evis;
    }
  }

  # Interactor Interactor_overlapping_gene WBGeneXYZ Interactor_type Effector
  if ($interaction->Variation_interactor ) {
    foreach my $v ($interaction->get('Variation_interactor')) {
      my @evis;
      foreach my $col ($v->col) { # foreach type
        my @evi_cmps = ($col->name);
        if ($col->right) {
          push @evi_cmps, $col->right->name; # push the type
        }
        
        my $evidence = join(" ", @evi_cmps);
        push @evis, $evidence;
      }
      
      my $var = $db->fetch( Variation => $v );
      if (defined $var->Gene) {
        my @gene = $var->Gene;
        my (%found, %not_found);

        foreach my $g (@gene) {
          if (exists $current_genes{$g->name}) {
            $found{$g} = 1;
          } else {
            $not_found{$g} = 1;
          }
        }
        if (keys %found) {
          # transfer evidence if the gene does not already have some
          foreach my $g (keys %found) {
            if (@evis and not @{$current_genes{$g}}) {
              $log->write_to("$interaction : Transferring Info only from $var to existing gene connection $g\n") if $debug;
              $results{Interactor_overlapping_gene}->{$g->name} = \@evis;
            } elsif (@{$current_genes{$g}} and not @evis) {
              $log->write_to("$interaction : Transferring Info only from $g to existing Variation_interactor $v\n") if $debug;
              $results{Variation_interactor}->{$v} = $current_genes{$g};
            }
          } 
        } else {
          $log->write_to(sprintf("%s : Adding new gene connections (%s) with Info via %s\n", 
                                 $interaction, 
                                 join(',', keys %not_found), 
                                 $v)) if $debug;
          foreach my $g (keys %not_found) {
            $results{Interactor_overlapping_gene}->{"$g"} = \@evis;
          }
        }
      }
    }
  }

  if ($interaction->Unaffiliated_variation) {
    foreach my $v ($interaction->get('Unaffiliated_variation')) {
      my $var = $db->fetch( Variation => $v );
      if (defined $var->Gene) {
        foreach my $g ($var->Gene) {
          if (exists $current_genes{$g->name}) {
            $log->write_to("$interaction : Promoting unaffilated variation $v to Variation_interactor\n") if $debug;
            $results{Variation_Interactor}->{$v} = $current_genes{$g->name};
            $delete_unaffiliated{$v} = 1;
          }
        }
      }
    }
  }

  if (keys %results) {
    print $outfh "\nInteraction : \"$interaction\"\n";
    foreach my $tag (keys %results) {
      foreach my $obj (keys %{$results{$tag}}) {
        my @evi = @{$results{$tag}->{$obj}};
        if (@evi) {
          foreach my $evi (@evi) {
            print $outfh "$tag \"$obj\" $evi\n";
          }
        } else {
            print $outfh "$tag \"$obj\"\n";
        }
      }
    }
  }
  foreach my $v (sort keys %delete_unaffiliated) {
    print $outfh "\nInteraction : \"$interaction\"\n";
    print $outfh "-D Unaffiliated_variation \"$v\"\n";
  }
}

close($outfh) or $log->log_and_die("Could not succesfully close $acefile\n");
$db->close;

unless ($noload) {
  $log->write_to("Loading file to autoace\n") if $debug;
  $wb->load_to_database( $wb->autoace, $output, 'map_interaction_script', $log );
}

$log->mail();
exit(0);

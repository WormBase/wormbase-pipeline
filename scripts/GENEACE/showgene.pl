#!/software/bin/perl -w
#
# showgene.pl
# 
# by Gary Williams                  
#
# Small utility for showing the current information of a gene.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/perl';
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use FileHandle;          # Or `IO::Handle' or `IO::'-anything-else used for autoflush

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;


=pod

=head batch_pname_update.pl

=item Parameters


    -gene genename 
     show information about a gene
     genename can be the WBGene ID, the sequence name or the CGC name

    -variation varname
     show information about a variation
     varname can be WBVarID or variation name

    -seqfeature featurename
     show information about a sequence feature
     featurename is WBsfID

    -variation strainname
     show information about a strain
     strainname can be WBStrainID or strain name

    -entity
     lists all entities in the nameserver


    -all
     Prints all of the change events, not just the first and last


=cut


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $verbose, $store, $wormbase);
my ($species, $gene, $variation, $seqfeature, $strain, $entity, $all);
my $BATCH_SIZE = 500; # maximum entries to put into any one batch API call

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "gene:s"     => \$gene,
	    "variation:s"=> \$variation,
	    "sequencefeature:s" => \$seqfeature,
	    "strain:s"   => \$strain,
	    "entity"     => \$entity,
	    "all"        => \$all,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);


my $db = NameDB_handler->new($wormbase);


if (defined $gene) {
  
  find_gene($gene, $all);

} elsif (defined $variation) {
  
  find_variation($variation);

} elsif (defined $seqfeature) {
  
  find_seqfeature($seqfeature);

} elsif (defined $strain) {
  
  find_strain($strain);

} elsif ($entity) {
  
  find_entity();
  
} else {
  die "Please specify the gene name with: -gene name\n";
}

$db->close;


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


##########################################

sub find_gene {
  my ($gene, $all) = @_;
  my $info = $db->info_gene("$gene");
  
  print "// Genes matching '$gene'\n";
  my $id = $info->{'id'};
  my $seqname = $info->{'sequence-name'};
  if (!defined $seqname) {$seqname = '.'}
  my $cgcname = $info->{'cgc-name'};
  if (!defined $cgcname) {$cgcname = '.'}
  my $species = $info->{'species'};
  my $status = $info->{'status'};
  my $biotype = $info->{'biotype'};
  if (!defined $biotype) {$biotype = '.'}
  print "Gene: $id\nSequence-name $seqname\nCGC name $cgcname\n";
  print "Species: $species\nStatus $status\nBiotype $biotype\n";
  # order the history events by time, so element [0] is the creation and element [-1] is the last event
  my $no_changes = scalar @{$info->{'history'}};
  my @events = sort {$a->{'t'} cmp $b->{'t'}} @{$info->{'history'}};
  my $created_when = $events[0]{'when'};
  if (!defined $created_when) {$created_when = $events[0]{'t'}}  # {'t'} holds the date the data was put in the database
  my $created_by = $events[0]{'who'}{'name'};
  my $created_id = $events[0]{'who'}{'id'};
  my $creation_thing = $events[0]{'what'};
  print "Created: $creation_thing $created_when by $created_by ($created_id)\n";
  my $change = '';
  foreach my $change_type (@{$events[0]{'changes'}}) {
    $change .= $change_type->{'attr'} . " : " . $change_type->{'value'} . "\t"
  }
  print "\t$change\n";

  if ($no_changes > 1) {
    if ($all) {
      foreach my $element (1 .. $no_changes-1) {
	my $last_changed = $events[$element];
	my $last_changed_when = $last_changed->{'when'};
	if (!defined $last_changed_when) {$last_changed_when = $last_changed->{'t'}} # {'t'} holds the date the data was put in the database
	my $last_changed_by = $last_changed->{'who'}{'name'};
	my $last_changed_id = $last_changed->{'who'}{'id'};
	my $last_changed_thing = $last_changed->{'what'};
	my $last_changed_changes = $last_changed->{'changes'};
	print "Change: $last_changed_thing $last_changed_when by $last_changed_by ($last_changed_id)\n";
	$change = '';
	foreach my $change_type (@{$last_changed_changes}) {
	  $change .= $change_type->{'attr'} . " : " . $change_type->{'value'} . "\t"
	}
	print "\t$change\n";
      }

    } else {
      
      my $last_changed = $events[$no_changes-1];
      my $last_changed_when = $last_changed->{'when'};
      if (!defined $last_changed_when) {$last_changed_when = $last_changed->{'t'}} # {'t'} holds the date the data was put in the database
      my $last_changed_by = $last_changed->{'who'}{'name'};
      my $last_changed_id = $last_changed->{'who'}{'id'};
      my $last_changed_thing = $last_changed->{'what'};
      my $last_changed_changes = $last_changed->{'changes'};
      print "Last change: $last_changed_thing $last_changed_when by $last_changed_by ($last_changed_id)\n";
      $change = '';
      foreach my $change_type (@{$last_changed_changes}) {
	$change .= $change_type->{'attr'} . " : " . $change_type->{'value'} . "\t"
      }
      print "\t$change\n";
    }
  }

}

##########################################
sub find_variation {

  my ($variation) = @_;

  my $info = $db->{'db'}->curl('GET', "entity/variation/$variation");

  if (exists $info->{'message'} && $info->{'message'} eq 'Entity lookup failed') {
    print "$variation Not found\n";
  } else {
    print "ID: ".$info->{'id'}."\n";
    print "Name: ".$info->{'name'}."\n";
    print "Status: ".$info->{'status'}."\n";
  }
}
##########################################
sub find_seqfeature {

  my ($seqfeature) = @_;

  my $info = $db->{'db'}->curl('GET', "entity/sequence-feature/$seqfeature", undef, 1);
  if (exists $info->{'message'} && $info->{'message'} eq 'Entity lookup failed') {
    print "$seqfeature Not found\n";
  } else {
    print "ID: ".$info->{'id'}."\n";
    print "Name: ".$info->{'name'}."\n";
    print "Status: ".$info->{'status'}."\n";
  }
}
##########################################
sub find_strain {

  my ($strain) = @_;

  my $info = $db->{'db'}->curl('GET', "entity/strain/$strain", undef, 1);

  if (exists $info->{'message'} && $info->{'message'} eq 'Entity lookup failed') {
    print "$strain Not found\n";
  } else {
    print "ID: ".$info->{'id'}."\n";
    print "Name: ".$info->{'name'}."\n";
    print "Status: ".$info->{'status'}."\n";
  }
}
##########################################
sub find_entity {

  my $info = $db->{'db'}->curl('GET', 'entity', undef, 1);
  foreach my $entity (@{$info->{'entity-types'}}) {
    my $name = $entity->{'entity-type'};
    print "$name\n";
  }

}
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut

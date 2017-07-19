#!/usr/bin/perl
#
# dumps all interactions into a flatfile, except the no_interaction ones
#

use strict;

use Getopt::Long;
use Dumper;
use IO::File;
use Storable;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;


my ($store,$debug,$test,$database,$species,$outfile);
GetOptions(
       'store=s' => \$store,
       'debug=s' => \$debug,
       'test'    => \$test,
       'species=s'  => \$species,
       'database=s' => \$database,
       'outfile=s'  => \$outfile,
)||die(@!);

my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, 
                                 -test => $test,
                                 -organism => $species)}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;

$log->write_to("connecting to $database\n");
my $dbh = Ace->connect(-path => $database ) or $log->log_and_die("Could not connect to $database\n");

$outfile = $wormbase->reports . '/interactions.txt' if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

print $of "# WormBase gene interactions\n";
print $of "# WormBase version: " . $dbh->version . "\n";
print $of '# Generated: ',&get_date,"\n";
print $of '# ' .  join("\t",qw/WBInteractionID Interaction_type Interaction_subtype Summary Citation Interactor1 Common-name Role1 Interactor2 Common-name Role2 .../),"\n";

# ignore objects with invalid name and Predicted
my @interactions = $dbh->fetch(-query=>'find Interaction WBInteraction????????? ! Interaction_type = Predicted'); 

foreach my $interaction (@interactions) {
    my ($brief_citation,$db_field,$db_acc) = eval { 
          $interaction->Paper->Brief_citation,$interaction->Paper->Database(2),$interaction->Paper->Database(3) 
    };

    #my $reference = "[$db_field:$db_acc] $brief_citation";
    my $interaction_type = $interaction->Interaction_type;
    my $subtype = ($interaction_type=~/Genetic|Regulatory/) ? &right_tip($interaction_type) : 'N/A';

    # exclude negative results
    next if $subtype=~/No_interaction/;
    next if $interaction->Regulation_result=~/Does_not_regulate/;
    
    my @cols = ($interaction,$interaction_type,$subtype);

    my $summary = eval{$interaction->Interaction_summary};
    push @cols,$summary;

    my $reference = eval{$interaction->Paper->Brief_citation};
    push @cols,$reference;

    foreach my $interactor_type ($interaction->Interactor) {# e.g. PCR_interactor, Interactor_overlapping_gene
	my $role = '';
	my $count = 0;
	foreach my $interactor ($interactor_type->col) {
	    my @interactors = $interactor_type->col;
	    my @tags = eval{ $interactors[$count++]->col }; # Interactor_info
	    my %info;
	    $info{obj} = $interactor;
	    if (@tags) {
		map { $info{"$_"} = $_->at; } @tags;

                # Exclude those that have been translated to Gene
		if ($interactor_type =~ /Other_regulator|Interactor_overlapping_gene|Molecule_regulator|Other_regulated|Rearrangement/) {
                    $role = $info{Interactor_type};
		    my @interactor_names;
		    my $public_name = eval{$interactor->Public_name} ||'N/A';
		    push (@interactor_names,$interactor,$public_name);
		    push (@cols,@interactor_names,$role);
		}
	    }
	}
    }
    print $of join("\t",@cols) . "\n";
}

$log->mail;
$of->close;
exit(0);

# as right() doesn't take -1
sub right_tip{
  my ($node) = (@_);
  my $i=0;
  my $last;
  while (my $value = $node->right($i)){
     $last = $value;
     $i++;
  }
  return $last;
}

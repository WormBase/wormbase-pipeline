#!/usr/bin/env perl 
#
# transcriptmasker.pl
#
# masks out ?Feature data spans in mRNA/ESTs prior to the BLAT analysis
# (essentially to remove TSL and polyA sequences)

# 031023 dl1

# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2013-07-22 09:26:26 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Bio::SeqIO;
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Carp;
use File::Path;
use Storable;

$|=1;

##############################
# command-line options       #
##############################

my $debug;    # debug mode
my $database; #
my $test;     # Test mode
my $store;    # Wormbase.pm storable
my $help;     # Help/Usage page
my $mrna;     # mRNA data
my $ncrna;    # ncRNA data
my $est;      # EST data
my $ost;      # OST data
my $nematode; # non-washu or nembase ests.  Just to get them dumped from database - not masked.
my $all;      # all of the above
my $file;

my ($species, $mol_type, $qspecies);

GetOptions (
	    "all"            => \$all,#not a valid option
	    "mrna"           => \$mrna,#not a valid option
	    "ncrna"          => \$ncrna, #not a valid option
	    "est"            => \$est, #not a valid option
	    "ost"            => \$ost,#not a valid option
	    "nematode"       => \$nematode,#not a valid option
	    "debug:s"        => \$debug,
	    "help"           => \$help,
	    "database:s"     => \$database,
	    "test"           => \$test,
	    "store:s"        => \$store,
	    "file:s"         => \$file,
	    "species:s"      => \$species,
	    "qspecies:s"     => \$qspecies,
	    "mol_type:s"     => \$mol_type
	   );

# Help pod if needed
&usage("Help") if ($help);

my $wormbase;

#if running from BLAT_controller and masking a different species some jiggery-pokery is needed to 
# make sure the correct paths get used and debug/test info is not lost.

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test     => $test,
                             -organism => $species
                             );
}

if( $qspecies and ($qspecies ne $wormbase->species) ){
  $wormbase = Wormbase->new( -debug   => $wormbase->debug,
                             -test     => $wormbase->test,
                             -organism => $qspecies
                             );
}						

my $log = Log_files->make_build_log($wormbase);
$log->write_to("mol_type must be specifies\n") unless $mol_type;

$species = $wormbase->species;

# valid Feature_data methods
my %valid_methods = (
		     "TSL" 		 => 1, # TSL data to flag data that is of unknown SL type.
		     "SL1" 		 => 1, # SL1 TSL data
		     "SL2" 		 => 1, # SL2 TSL data
		     "polyA"		 => 1, # polyA tail sequences
		     "poly_nucleotide"   => 1, # poly nucleotide sequences (e.g. 5' AAAAAAA)
		     "BLAT_discrepancy"  => 1, # General method for BLAT problems
		     "vector"  		 => 1, # Vector sequences
		     "chimeric"          => 1, # Chimeric EST reads (1 part of the EST is masked out).
		     "rRNA_contamination"=> 1, # rRNA_contamination of library.
		     "Oligo_cap"         => 1, # Oligo cap sequence masking.
		     "low_complexity"    => 1, #low-complexity or poor quality/unclipped sequence.
		     'low'               => 1, #low-complexity or poor quality/unclipped sequence.
		    );

# which database?
if (-e $wormbase->orgdb."/database/block1.wrm") {
  $database = $wormbase->orgdb unless $database;
}
elsif ($species eq "elegans"){
  $database = $wormbase->database('camace');
}
else {
  $database = $wormbase->database($species);
}

my $blat_dir =  $wormbase->blat;
my $tace  = $wormbase->tace;                          # tace executable path
my $id;                                               # id for the entry
my $seq;                                              # raw sequence for the entry
my %sequence;                                         # 
my $feature;                                          #
my $seqmasked;                                        #
my $seqlength;                                        #
my $masked;                                           # No of entries masked
my $ignore;
my %seq2feature;                                      #stores feature data info

#remove all old masked data
&remove_masked_files($species, $mol_type);

#get all of the Feature_data via table maker
&fetch_features;	

# which data file to parse
$masked = &MaskSequence($species, $mol_type);
$log->write_to("\n=============================================\n\n");

#########################################
# hasta luego                           #
#########################################

$log->mail;
exit(0);


############################################################
######################## Subroutines #######################
############################################################



sub remove_masked_files {
  my $species  = shift;
  my $mol_type = shift;
  if ((-e $wormbase->maskedcdna."/$mol_type.masked") or (-e $wormbase->maskedcdna."/$mol_type.masked_1")) {
    $log->write_to("Removing old $mol_type masked files for $species in ".$wormbase->maskedcdna."\n");
    $wormbase->run_command ("rm ".$wormbase->maskedcdna."/$mol_type.masked*", $log);
  }
  else {
    $log->write_to("No old $mol_type masked files for $species to be removed from ".$wormbase->maskedcdna."\n");
  }
}

#_ MaskSequence -#
# 
# pass type of transcript data to be masked (e.g. mRNA, EST etc)

sub MaskSequence {
  my $species  = shift;
  my $mol_type = shift;
  my $masked = 0;
  #left ignored and ignore after merge as they seem to be meant for different things - ignored seems to be doing nothing though!
  my $ignored = 0;
  my $ignore ;

  #$log->write_to("masking $mol_type for $species\n");

  # set input record seperator
  $/ = ">";
  
  # assign output file
  my $output_file = $file ? $file."masked" : $wormbase->maskedcdna."/$mol_type.masked";
  mkpath ($wormbase->maskedcdna) unless -e $wormbase->maskedcdna;
  open (OUTPUT, ">$output_file") || $log->log_and_die("ERROR: Can't open output file: $output_file");
  
  # input file loop structure
  my $file2mask = $file ? $file : $wormbase->cdna_dir."/$mol_type";
  $log->write_to("masking $mol_type for $species input from $file2mask output to $output_file\n");
  my $seq_in = Bio::SeqIO->new(-file => "$file2mask" , '-format' => 'Fasta');
 SEQ: while (my $seqobj = $seq_in->next_seq) {
    my $seqmasked = $seqobj->seq;		# copy sequence to masked file and
    my $seqlength = length ($seqmasked);	# calculate the length of the sequence
	my $id = $seqobj->id;
    if( $seq2feature{$id} ) {
      foreach my $feature (keys %{$seq2feature{$id}}) {
      	my $coords = $seq2feature{$id}->{$feature};
		if (defined($feature)) {
	  	next unless ($valid_methods{$feature}); # only mask valid Feature_data methods
	  	my ($method,$start,$stop,$length,$remark);
	  	my ($cut_to,$cut_from,$cut_length);
	  	$start = $coords->[0]; # start coord
	  	$stop  = $coords->[1]; # stop coord
		if ($stop > $seqlength) {
		  $log->write_to("$id: $feature feature is past end of sequence\n"); 
		  $stop = $seqlength; # correct the end position of the feature
		}
	  	$cut_to     = $start - 1; # manipulations for clipping 
	  	$cut_from   = $stop;
	  	$cut_length = $stop - $start + 1;
	  	if ($cut_to < 0 ) {
	  	  $cut_to = 0;
	  	}			# fudge to ensure non-negative clipping coords
	  
	  	print "$id: '$feature' $start -> $stop [$cut_to : $cut_from ($cut_length)]\n" if ($wormbase->debug);
	  	print "// # [$id] $feature:" . (substr($seqmasked,$cut_to,$cut_length)) . " [$start - $stop]\n\n" if ($debug);
	 	my $newseq = (substr($seqmasked,0,$cut_to)) . ('n' x $cut_length)  . (substr($seqmasked,$cut_from));
	 	$seqmasked = $newseq;
		}
	# increment count of sequences masked
		$masked++;
      }		
    }
    # output masked sequence
    print OUTPUT ">$id\n$seqmasked\n";
  }
  $/ = "\n";
  close OUTPUT;
  return ($masked);
}
#_ end MaskSequence _#



###############################################################

sub usage {
  my $error = shift;
    
  if ($error eq "Help") {
    # Help menu
    exec ('perldoc',$0);
  }
  elsif ($error == 0) {
    # Normal help menu
    exec ('perldoc',$0);
  }
}


sub fetch_features {
  my $tm_data = $wormbase->table_maker_query($database,$wormbase->autoace."/wquery/SCRIPT:transcriptmasker.def");
  while( <$tm_data> ) {
    s/\"//g;  #"
    next if (/acedb/ or /\/\// or /^$/);
    my ($seq, $type, $start, $end) = split;
    if ($seq and $type and $start and $end) {
      $seq2feature{$seq}->{$type} = [($start, $end)];
    }
    else {
      $log->error("something wrong with $seq:$type\n");
    }
  }
  $log->write_to( scalar(keys %seq2feature) ." sequences to be masked\n");
}

__END__

=pod

=head2   NAME - transcriptmasker.pl

=head1 USAGE

=over 4

=item transcriptmasker.pl [-options]

=back

transcriptmasker.pl takes the fasta file generated using SRS and queries autoace
for ?Feature_data classes. It will mask out any annotated features (currently TSL
acceptor sites, and polyA+ tails) with N's. These masked files are then used by 
the BLAT scripts to map back to the genome.

transcriptmasker.pl mandatory arguments:

=over 4

=item none,

=back

transcriptmasker.pl OPTIONAL arguments:

=over 4

=item -all, process EST, mRNA & OST datafiles 

=item -est, process EST datafile 

=item -ost, process OST datafile 

=item -mrna, process mRNA datafile 

=item -debug <user>, Debug mode - log file only goes to named user

=item -help, these help pages

=back

=cut

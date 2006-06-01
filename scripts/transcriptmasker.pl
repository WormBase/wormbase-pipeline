#!/usr/local/bin/perl5.8.0 -w
#
# transcriptmasker.pl
#
# masks out ?Feature data spans in mRNA/ESTs prior to the BLAT analysis
# (essentially to remove TSL and polyA sequences)

# 031023 dl1

# Last edited by: $Author: pad $
# Last edited on: $Date: 2006-06-01 10:40:04 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Ace;
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
my $all;      # all of the above
my $file;

GetOptions (
	    "all"            => \$all,
	    "mrna"           => \$mrna,
	    "ncrna"          => \$ncrna,
	    "est"            => \$est,
	    "ost"            => \$ost,
	    "debug:s"        => \$debug,
	    "help"           => \$help,
	    "database:s"     => \$database,
	    "test"           => \$test,
	    "store:s"        => \$store,
	    "file:s"         => \$file,
	   );

# Help pod if needed
&usage("Help") if ($help);

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);


# datafiles for input
my $EST_dir = $wormbase->wormpub."/analysis/ESTs";

our %datafiles = (
		  "mrna"  => "elegans_mRNAs",
		  "ncrna" => "elegans_ncRNAs",
		  "est"   => "elegans_ESTs",
		  "ost"   => "elegans_OSTs"
		 );

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
		     "Oligo_cap"         => 1 # Oligo cap sequence masking.
		    );
# transcript accessions to names from a hash in common data

print "// Reading EST_names.dat hash\n\n" if ($debug);
our %EST_name = $wormbase->FetchData('NDBaccession2est');
print "// Finished reading EST_names.dat hash\n\n" if ($debug);

# which database
$database = $wormbase->autoace unless $database;
my $blat_dir =  $wormbase->blat;
my $tace  = $wormbase->tace;                          # tace executable path
my $acc;                                              # accession for the entry
my $id;                                               # id for the entry
my $seq;                                              # raw sequence for the entry
my %sequence;                                         # 
my $feature;                                          #
my $seqmasked;                                        #
my $seqlength;                                        #
my $masked;                                           # No of entries masked
my $ignore;
my %seq2feature;                                      #stores feature data info

#get all of the Feature_data via table maker
&fetch_features;	

# connect to database
print  "\nOpening $database for masking ..\n" if ($wormbase->debug);
my $db = Ace->connect(	-path=>$database,
			-program =>$tace) || $log->log_and_die("Connection failure: ".Ace->error."\n");

# which data file to parse
$masked = &MaskSequence($datafiles{mrna}) if ($mrna || $all);
$log->write_to($wormbase->runtime." : masked $masked mRNA sequences\n") if ($mrna || $all);

$masked = &MaskSequence($datafiles{ncrna}) if ($ncrna || $all);
$log->write_to($wormbase->runtime." : masked $masked ncRNA sequences\n") if ($ncrna || $all);

$masked = &MaskSequence($datafiles{est})  if ($est || $all);
$log->write_to($wormbase->runtime." : masked $masked EST sequences\n") if ($est || $all);

$masked = &MaskSequence($datafiles{ost})   if ($ost || $all);
$log->write_to($wormbase->runtime." : masked $masked OST sequences\n") if ($ost || $all);

$log->write_to("\n=============================================\n\n");

#########################################
# hasta luego                           #
#########################################

$log->mail;
$db->close;
exit(0);


############################################################
######################## Subroutines #######################
############################################################

#_ MaskSequence -#
# 
# pass type of transcript data to be masked (e.g. mRNA, EST etc)

sub MaskSequence {
  my $data   = shift;#gets populated by value from system call.
  my $masked = 0;
  #left ignored and ignore after merge as they seem to be meant for different things - ignored seems to be doing nothing though!
  my $ignored = 0;
  my $ignore ;
  
  # set input record seperator
  $/ = ">";
  
  # assign output file
  my $output_file = $file ? $file."masked" : "$blat_dir/$data.masked";
  open (OUTPUT, ">$output_file") || $log->log_and_die("ERROR: Can't open output file: $output_file");
  
  # input file loop structure
  my $skip = 1;
  
  my $file2mask = $file ? $file : "$EST_dir/$data";  
  open (INPUT, "<$file2mask")     || $log->log_and_die("ERROR: Can't open input file: $file2mask");
 SEQ: while (<INPUT>) {
    chomp;
    next if ($_ eq "");		# catch empty lines
    if ($skip == 1) {		# skip first bad entry
      $skip--;
      next;
    }
    next unless /\w/; #skip blank lines
    if (/^(\S+)\s+\S+.+\n/) {	# deal with accessions {$acc} and WormBase internal names {$id}
      $acc = $1;
      if (defined $EST_name{$acc}) {
	$id  = $EST_name{$acc};
      } 
      else {
	#		print "// ERROR: No accession-id connection for this sequence [$acc]\n" if ($wormbase->debug);
	$EST_name{$acc} = $acc;
	$id = $acc;
      }
    }
    
    my $seq = "$'";		# assign the rest of the string to $seq
    $seq =~ s/[^gatcn]//g;	# remove non-gatcn characters (i.e. newlines)
    my $seqmasked = $seq;		# copy sequence to masked file and
    my $seqlength = length ($seq);	# calculate the length of the sequence
	 
    # fetch the sequence object from the database. push the feature_data objects to memory
    my $obj = $db->fetch(Sequence=>$id);
    if (!defined ($obj)) {
      $log->write_to("ERROR: Could not fetch sequence $id \n\n");
      next;
    }

    # Is the Ignore tag set?
    if (defined $obj->Ignore) {
      print "\n// Ignore tag set for $acc $id\n\n" if ($debug);
      next SEQ;
    }	
	
    if( $seq2feature{$id} ) {
      foreach my $feature (keys %{$seq2feature{$id}}) {
      	my $coords = $seq2feature{$id}->{$feature};
	if (defined($feature)) {
	  next unless ($valid_methods{$feature}); # only mask valid Feature_data methods
	  my ($method,$start,$stop,$length,$remark);
	  my ($cut_to,$cut_from,$cut_length);
	  $start = $coords->[0]; # start coord
	  $stop  = $coords->[1]; # stop coord
	  $cut_to     = $start - 1; # manipulations for clipping 
	  $cut_from   = $stop;
	  $cut_length = $stop - $start + 1;
	  if ($cut_to < 0 ) {
	    $cut_to = 0;
	  }			# fudge to ensure non-negative clipping coords
	  
	  print "$acc [$id]: '$feature' $start -> $stop [$cut_to : $cut_from ($cut_length)]\n" if ($wormbase->debug);
	  print "// # $acc [$id] $feature:" . (substr($seq,$cut_to,$cut_length)) . " [$start - $stop]\n\n" if ($debug);
	  my $newseq = (substr($seqmasked,0,$cut_to)) . ('n' x $cut_length)  . (substr($seqmasked,$cut_from));
	  $seqmasked = $newseq;
	}
	# increment count of sequences masked
	$masked++;
      }		
    }
    # output masked sequence
    print OUTPUT ">$acc $id\n$seqmasked\n";
  }
  close INPUT;
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
  my $tm_data = $wormbase->table_maker_query($database,"$database/wquery/SCRIPT:transcriptmasker.def");
  while( <$tm_data> ) {
    s/\"//g;
    my ($seq, $type, $start, $end) = split;
    if ($seq and $type and $start and $end) {
      $seq2feature{$seq}->{$type} = [($start, $end)];
    }
    else {
      $log->error("something wrong with $seq:$type\n");
    }
  }
  $log->write_to(scalar(keys %seq2feature), " sequences to be masked\n");
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

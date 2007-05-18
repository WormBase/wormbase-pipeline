#!/usr/local/bin/perl5.8.0 
#
# makesuperlinks.pl
#
# constructs superlink objects for camace based on Overlap_right tags
# dl
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2007-05-18 13:45:00 $
 
$!=1;
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace ;


##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);
 
# Database name for databases other than ~wormpub/DATABASES/camace
my $db;         
# output file for ace
my $acefile;		      

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "db=s"       => \$db,
	    "acefile=s"  => \$acefile,
            );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}
 
# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


##############################
# superlink data for camace  #
##################Structure.Subsequence############

our %clone2super = ("H06O01"  => "I",
		    "Y95D11A" => "IR",
		    "T05A6"   => "II",
		    "F45H7"   => "IIIL",
		    "K01F9"   => "IIIR",
		    "F21D5"   => "IV",
		    "F07D3"   => "V",
		    "B0272"   => "X",
		    ) ;

##############################
# open output file
##############################

if (! $acefile) {die "-acefile is not specified\n";}
open (ACE, ">$acefile") || die "Can't open file $acefile\n";

##############################
# get data from camace       #
##############################

my $tace = $wormbase->tace;


my $campath = "/nfs/disk100/wormpub/DATABASES/camace";
($campath = $db) if ($db);
# check database path
&usage("Bad_database_name") unless (-e "$campath/database/database.map");
warn "Acessing database $campath\n" if ($verbose);

# connect to database
my $camdb   = Ace->connect(-path=>"$campath",
			   -program => $tace) || die "failed to connect to database\n";


# first get lots of info about Genome_sequence objects

my (%isGenomeSequence,%length,%right,%rightOffset,%left,%isLinkCandidate,%isExternal);
my (%currSource,%currStart,%currEnd);
my (%CDSSource,%CDSStart,%CDSEnd);
my (%PseudoSource,%PseudoStart,%PseudoEnd);
my (%TransposonSource,%TransposonStart,%TransposonEnd);
my (%TranscriptSource,%TranscriptStart,%TranscriptEnd);
my ($it,$obj,$seq,$start,$end);

my $error = 0;			# error status to return from program

$it = $camdb->fetch_many(Genome_sequence => '*') ;
while ($obj = $it->next) {
  if (! defined $obj) {
    $log->write_to("Genome_sequence not defined\n");
    $error = 1;
  }
    $isGenomeSequence{$obj} = 1 ;
    $length{$obj}           = $obj->DNA(2) ;
    if (!$length{$obj})     { 
      $log->write_to("No length for $obj\n") ; 
      $error = 1;
    }
    $right{$obj}            = $obj->Overlap_right ;
    $rightOffset{$obj}      = $obj->Overlap_right(2) ;
    $left{$obj}             = $obj->Overlap_left ;

    foreach ($obj->Confidential_remark) {
	if (/^not in Cambridge LINK$/) { $isExternal{$obj} = 1; }
    }
    $isLinkCandidate{$obj}  = (!$isExternal{$obj} && $length{$obj} > 0);
}

# then some info about current links

$it = $camdb->fetch_many(Sequence => 'SUPERLINK*') ;
while ($obj = $it->next) {
    foreach $a ($obj->at('Structure.Subsequence')) {
      if (! defined $a) {
	$log->write_to("superlink Structure.Subsequence not defined\n"); 
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("Structure.Subsequence row not defined\n"); 
	$error = 1;
      }
      $currSource{$seq}    = $obj;
      $currStart{$seq}     = $start;
      $currEnd{$seq}       = $end;
	
#	print "// push $seq to subsequence hash\n";

    }

    foreach $a ($obj->at('SMap.S_child.CDS_child')) {
      if (! defined $a) {
	$log->write_to("superlink SMap.S_child.CDS_child not defined\n");
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("SMap.S_child.CDS_child row not defined\n");
	$log->write_to($obj->name."\n");
	if (defined$seq) {$log->write_to("$seq\n")}
	$error = 1;
      }
      $CDSSource{$seq}    = $obj;
      $CDSStart{$seq}     = $start;
      $CDSEnd{$seq}       = $end;
#	print "// push $seq to CDS hash\n";

    }

    foreach $a ($obj->at('SMap.S_child.Pseudogene')) {
      if (! defined $a) {
	$log->write_to("superlink SMap.S_child.Pseudogene not defined\n");
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("SMap.S_child.Pseudogene row not defined\n");
	$error = 1;
      }
      $PseudoSource{$seq}    = $obj;
      $PseudoStart{$seq}     = $start;
      $PseudoEnd{$seq}       = $end;
#	print "// push $seq to Pseudo hash\n";

    }

    foreach $a ($obj->at('SMap.S_child.Transcript')) {
      if (! defined $a) {
	$log->write_to("superlink SMap.S_child.Transcript not defined\n");
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("SMap.S_child.Transcript row not defined\n");
	$error = 1;
      }
      $TranscriptSource{$seq}    = $obj;
      $TranscriptStart{$seq}     = $start;
      $TranscriptEnd{$seq}       = $end;
#	print "// push $seq to Transcript hash\n";

    }

    foreach $a ($obj->at('SMap.S_child.Transposon')) {
      if (! defined $a) {
	$log->write_to("superlink SMap.S_child.Transposon not defined\n");
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("SMap.S_child.Transposon row not defined\n");
	$error = 1;
      }
      $TransposonSource{$seq}    = $obj;
      $TransposonStart{$seq}     = $start;
      $TransposonEnd{$seq}       = $end;
#	print "// push $seq to Transposon hash\n";
    }

}

warn "Stored data to hash\n" if ($verbose);

###########################################
# make links
###########################################

my ($lk,%start,%end,%link,$parent,$startright);

foreach $seq (keys %isGenomeSequence) {
    
    # only keep seeds
    next if (!$isLinkCandidate{$seq} ||
	     ($left{$seq} && $isLinkCandidate{$left{$seq}} && $rightOffset{$left{$seq}}) ||
	     !$rightOffset{$seq} ||
	     !$isLinkCandidate{$right{$seq}});

    # don't make link objects for these clones
    # these should be in St Louis SUPERLINKS
    next if (($seq eq "6R55") || ($seq eq "F38A1") || ($seq eq "cTel52S") || ($seq eq "cTel7X"));

    # print LINK header
    $lk = "SUPERLINK_CB_$clone2super{$seq}";
    print ACE "\nSequence $lk\n";
    print ACE "From_laboratory HX\n";

    # loop over subsequences
    $startright = 1;
    while ($isLinkCandidate{$seq}) {
	$start{$seq}  = $startright; 
	$end{$seq}    = $startright + $length{$seq} - 1;
	print ACE "Subsequence $seq $start{$seq} $end{$seq}\n";
	$link{$seq}   = $lk;
	if (!$rightOffset{$seq}) {
	    warn "ending loop here because rightOffset{$seq} is not set\n" if ($verbose);
	    last;
	}     # POSS EXIT FROM LOOP
	$startright   = $startright + $rightOffset{$seq} - 1;
	$seq          = $right{$seq};
    }		
}

###########################################
# Re-map subsequences back onto the new links 
###########################################

foreach $seq (keys %CDSStart) {
    next if ($isGenomeSequence{$seq});
    
    if ($seq =~ /(\S+)\.\d+/) {
	$parent = $1; 
    } elsif ($seq =~ /(\S+)\.gc\d+/) {	# allow genefinder-style IDs to be read
	$parent = $1; 
    } else { 
	$log->write_to("no dot in subsequence name $seq\n");
	$error = 1;
	next; 
    }
    
    # assign parent for problem child (it spans clones Y41C4A to Y66A7AR)
    if ($seq eq "Y66A7A.8") {$parent = "Y66A7AR";}
    
    # next if no coordinate for parent clone
    if (!$currStart{$parent}) { 
      if ($seq !~ /.tw$/) {	# don't complain about the twinscan CDSs - just ignore them
	$log->write_to("no coord in link for parent $parent of $seq\n");
	$error = 1;
      }
      next;
    }
    # next if parent and child map to different links
    if (!($CDSSource{$seq} eq $currSource{$parent})) { 
	$log->write_to("parent $parent and child $seq in different links\n"); 
	$error = 1;
	next;
    }
    # next if parent has no home to go to
    if (!$link{$parent}) {
	$log->write_to("no new link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if sequence is a superlink
    next if ($seq =~ /^SUPERLINK/);

    $start{$seq} = $start{$parent} - $currStart{$parent} + $CDSStart{$seq};
    $end{$seq}   = $start{$parent} - $currStart{$parent} + $CDSEnd{$seq};
    $link{$seq}  = $link{$parent};

    print ACE "\nSequence $link{$parent}\n";
    print ACE "CDS_child $seq $start{$seq} $end{$seq}\n";
}

foreach $seq (keys %PseudoStart) {
    next if ($isGenomeSequence{$seq});
    
    if ($seq =~ /(\S+)\.\d+/) {
	$parent = $1; 
    }
    else { 
	$log->write_to("no dot in subsequence name $seq\n");
	$error = 1;
	next; 
    }
    
    # next if no coordinate for parent clone
    if (!$currStart{$parent}) { 
	$log->write_to("no coord in link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if parent and child map to different links
    if (!($PseudoSource{$seq} eq $currSource{$parent})) { 
	$log->write_to("parent $parent and child $seq in different links\n"); 
	$error = 1;
	next;
    }
    # next if parent has no home to go to
    if (!$link{$parent}) {
	$log->write_to("no new link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if sequence is a superlink
    next if ($seq =~ /^SUPERLINK/);

    $start{$seq} = $start{$parent} - $currStart{$parent} + $PseudoStart{$seq};
    $end{$seq}   = $start{$parent} - $currStart{$parent} + $PseudoEnd{$seq};
    $link{$seq}  = $link{$parent};

    print ACE "\nSequence $link{$parent}\n";
    print ACE "Pseudogene $seq $start{$seq} $end{$seq}\n";
}

# Transposons

foreach $seq (keys %TransposonStart) {
    next if ($isGenomeSequence{$seq});
    
    if ($seq =~ /(\S+)\.\d+/) {
	$parent = $1; 
    }
    else { 
	$log->write_to("no dot in subsequence name $seq\n");
	$error = 1;
	next; 
    }
    
    # next if no coordinate for parent clone
    if (!$currStart{$parent}) { 
	$log->write_to("no coord in link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if parent and child map to different links
    if (!($TransposonSource{$seq} eq $currSource{$parent})) { 
	$log->write_to("parent $parent and child $seq in different links\n"); 
	$error = 1;
	next;
    }
    # next if parent has no home to go to
    if (!$link{$parent}) {
	$log->write_to("no new link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if sequence is a superlink
    next if ($seq =~ /^SUPERLINK/);

    $start{$seq} = $start{$parent} - $currStart{$parent} + $TransposonStart{$seq};
    $end{$seq}   = $start{$parent} - $currStart{$parent} + $TransposonEnd{$seq};
    $link{$seq}  = $link{$parent};

    print ACE "\nSequence $link{$parent}\n";
    print ACE "Transposon $seq $start{$seq} $end{$seq}\n";
}


foreach $seq (keys %TranscriptStart) {
    next if ($isGenomeSequence{$seq});
    
    if ($seq =~ /(\S+)\.\d+/) {
	$parent = $1; 
    }
    else { 
	$log->write_to("no dot in subsequence name $seq\n");
	$error = 1;
	next; 
    }
    
    # next if no coordinate for parent clone
    if (!$currStart{$parent}) { 
	$log->write_to("no coord in link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if parent and child map to different links
    if (!($TranscriptSource{$seq} eq $currSource{$parent})) { 
	$log->write_to("parent $parent and child $seq in different links\n"); 
	$error = 1;
	next;
    }
    # next if parent has no home to go to
    if (!$link{$parent}) {
	$log->write_to("no new link for parent $parent of $seq\n");
	$error = 1;
	next;
    }
    # next if sequence is a superlink
    next if ($seq =~ /^SUPERLINK/);

    $start{$seq} = $start{$parent} - $currStart{$parent} + $TranscriptStart{$seq};
    $end{$seq}   = $start{$parent} - $currStart{$parent} + $TranscriptEnd{$seq};
    $link{$seq}  = $link{$parent};

    print ACE "\nSequence $link{$parent}\n";
    print ACE "Transcript $seq $start{$seq} $end{$seq}\n";
}



###########################################
# Do a quick check to make sure that everything 
# that was in a link has been put back in one
###########################################

foreach $seq (keys %currSource) {
    if (!$link{$seq}) { 
      $log->write_to("$seq not put back into a link\n"); 
      $error = 1;
    }
}


$camdb->close;

close (ACE);


# Close log files and exit
$log->mail();
print "Finished.\n" if ($verbose);
exit($error);			# return the error status


############# end of file ################



#######################################################################
# Help and error trap outputs                                         #
#######################################################################
 
sub usage {
    my $error = shift;
    
    if ($error eq "Bad_database_name") {
        # Named database mode: database name not recognised
        print "\nNo database of this name exists ('$db')\n";
	exit(0);
    }
    elsif ($error eq "Help") {
        # Normal help menu
        system ('perldoc',$0);
        exit (0);
    }
}






__END__

=pod

=head2   NAME - makesuperlinks.pl

=head1 USAGE

=over 4

=item makesuperlinks.pl [-options]

=back

makesuperlinks queries an ACEDB database and generates the SUPERLINK
objects to an .acefile. This uses the Overlap_right tags within the
database and a hard-coded hash of starting clones within this
script

makesuperlinks mandatory arguments:

=over 4

=item B<-acefile>, file to output ACE to

=back

makesuperlinks OPTIONAL arguments:

=over 4

=item B<-db text>, database mode. Only dumps acefiles for the named database. The default is /nfs/disk100/wormpub/DATABASES/camace

=item B<-debug>, send output to specified user only 

=item B<-verbose>, verbose report

=item B<-help>, this help page

=back

=head1 AUTHOR (& person to blame)

=over 4

=item Dan Lawson dl1@sanger.ac.uk

=back

=cut

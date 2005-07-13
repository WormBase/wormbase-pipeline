#!/usr/local/bin/perl5.8.0 
#
# fetch_seqs_for_blatting.pl
#
# by Dan Lawson
#
# Attempt to unify all of the diverse scripts to fetch ESTs, OSTs, mRNAs etc. used by blat
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2005-07-13 14:25:23 $
 
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Data::Dumper;
 
##############################
# command-line options       #
##############################
 
my $help;                # Help/Usage page
my $ace;                 # Dump acefile
my $verbose;             # turn on extra output
my $blastdb;             # make blast database using pressdb?
my $ftp;                 # also copy to ftp site
my $debug;               # For sending output to just one person
my $maintainers = "All"; # log file recipients
my ($est, $mrna, $ncrna, $ost, $nematode, $embl, $all); # the main options

GetOptions (
            "est"      => \$est,
            "mrna"     => \$mrna,
            "ncrna"    => \$ncrna,
            "ost"      => \$ost,
            "nematode" => \$nematode,
            "embl"     => \$embl,
	    "all"      => \$all,
            "verbose"  => \$verbose,
            "blastdb"  => \$blastdb,
            "ace"      => \$ace,
            "ftp"      => \$ftp,
            "debug=s"  => \$debug,
            "help"     => \$help
            );

# Help pod if needed
&usage(0) if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

##############################
# Other script variables     #
##############################

my $acc;                 # EMBL accession
my $id;                  # EMBL ID
my $sv;                  # EMBL sequence version
my $def = "";            # EMBL description, needs to be initialised to ""
my $protid;              # EMBL Protein_ID
my $protver;             # EMBL Protein_ID version
my $org;                 # EMBL species
my $ost_seq;             # tag for OST/EST split

my %ISNB2wormbase = &FetchData('NDBaccession2est');      # EST accession => name
my %wormbase2ISNB = reverse %ISNB2wormbase;


my $log;                                                       # for log file
my $dir      = "/nfs/disk100/wormpub/analysis/ESTs";           # path for output files
my $ftpdir   = "/nfs/disk69/ftp/pub/wormbase/sequences/ESTS";  # path for ftp site
my $getz     = "/usr/local/pubseq/bin/getzc";                  # getz binary
our $tace    = &tace;                                          # tace binary
our $ace_dir = "/wormsrv2/autoace";                            # database to extract data from

if ($debug) {
    $tace = "/nfs/team71/acedb/edgrif/TEST/DAN/tace";       # tace binary
}

# get masking data.....

our %mask;
our %bmasked;

&get_masking_data if ($est || $ost || $mrna || $all);


#########################################
# MAIN   PART   OF   SCRIPT
#########################################

#&create_log_files;
&make_ests          if ($est || $ost || $all);
&make_mrnas         if ($mrna || $all);
&make_ncrnas        if ($ncrna || $all);

# still to deal with....
&make_embl_cds      if ($embl);
&make_nematode_ests if ($nematode);
 
exit(0);

########################################################################


#################
## Subroutines ##
#################

#################################################

sub make_ests {
    my $count = 0;

    print "// Fetching EST sequences\n";
    
    my $est_file = "elegans_ESTs.masked";
    my $ost_file = "elegans_OSTs.masked";
    my $name;
    
    my ($name,$cut_to,$cut_from,%sequence,$cut_length,$newseq,$seq);
    my ($type,$seqtag);

    open (OUT_EST, ">$dir/$est_file");
    open (OUT_OST, ">$dir/$ost_file");
    # assign tace query for this data set 
    my $command = "quiet -on\nquery find sequence where cDNA_EST AND NOT Ignore\ndna\nquit\n";
    # set carriage return char to >
    $/=">";
    # tace call to the database to dump dna
    open (SEQUENCES, "echo '$command' | $tace $ace_dir |");
    while (<SEQUENCES>) {	
	chop;

	s/acedb//g;
	# increment count
	$count++;        
	# print sequence to file
	($name) = (/^(\S+)/);
	$seq = $_;
	$seq =~ s/$name//g;
 	$seq =~ s/\W//g;
	$seq =~ s/[^gatcn]/n/g;
	
	$sequence{$name} = $seq;
    }
    # close file handles
    close (SEQUENCES);

    # Now do the masking

    foreach $name (sort keys %sequence) {

	next if ($sequence{$name} eq "Abientot");
	next if ($name eq "");

	foreach (@{$bmasked{$name}}) {
	    chomp;
	    $type = $_;
	    
	    if ($type eq "SL1") {
		$seqtag = $name . ":TSL";
	    }
	    elsif ($type eq "SL2") {
		$seqtag = $name . ":TSL";
	    }
	    else {
		$seqtag = $name . ":" . $type;
	    }
	    
	    $cut_to     = $mask{$seqtag}{START} - 1;                        # manipulations for clipping 
	    $cut_from   = $mask{$seqtag}{STOP};
	    $cut_length = $cut_from -  $cut_to;
	    
	    if ($cut_to < 0 ) {$cut_to = 0;}                 # fudge to ensure non-negative clipping coords
	    
	    $seq = $sequence{$name};

            print "// Sequence: $name ($seqtag) [$cut_to - $cut_from :" . substr($seq,$cut_to,$cut_length) ."\n";
	    $newseq = substr($seq,0,$cut_to) . ('n' x $cut_length)  . substr($seq,$cut_from);
	    
	    $sequence{$name} = $newseq;
	    
	}
	
	print  OUT_EST ">$name \n$sequence{$name}\n\n" unless ($name =~ /^OST/);
	print  OUT_OST ">$name \n$sequence{$name}\n\n" if     ($name =~ /^OST/);

    }
    
    # reset carriage return char to normal
    $/="\n";



    close (OUT_EST);
    close (OUT_OST);
}

#################################################

sub make_mrnas {
    my $count = 0;

    print "// Fetching mRNA sequences\n";
    
    my $mRNA_file = "elegans_mRNAs.masked";

    my ($name,$cut_to,$cut_from,%sequence,$cut_length,$newseq,$seq);

    my ($type,$seqtag);

    open (OUT_mRNA, ">$dir/$mRNA_file");
    # assign tace query for this data set 
    my $command = "quiet -on\nquery find sequence where mRNA AND NOT Ignore\ndna\nquit\n";
    # set carriage return char to >
    $/=">";
    # tace call to the database to dump dna
    open (SEQUENCES, "echo '$command' | $tace $ace_dir |");
    while (<SEQUENCES>) {
	chop;

	s/acedb//g;
	# increment count
	$count++;        
	# print sequence to file
	($name) = /^(\S+)/;
	$seq = $_;
	$seq =~ s/$name//g;
 	$seq =~ s/\W//g;
	$seq =~ s/[^gatcn]/n/g;
	
	$sequence{$name} = $seq;
    }
    close SEQUENCES;
    
    
    # Now do the masking

    foreach $name (sort keys %sequence) {

	next if ($sequence{$name} eq "Abientot");
	next if ($name eq "");

	foreach (@{$bmasked{$name}}) {
	    chomp;
	    $type = $_;
	    
	    if ($type eq "SL1") {
		$seqtag = $name . ":TSL";
	    }
	    elsif ($type eq "SL2") {
		$seqtag = $name . ":TSL";
	    }
	    else {
		$seqtag = $name . ":" . $type;
	    }
	    
	    $cut_to     = $mask{$seqtag}{START} - 1;                        # manipulations for clipping 
	    $cut_from   = $mask{$seqtag}{STOP};
	    $cut_length = $cut_from -  $cut_to;
	    
	    if ($cut_to < 0 ) {$cut_to = 0;}                 # fudge to ensure non-negative clipping coords
	    
	    $seq = $sequence{$name};

            print "// Sequence: $name ($seqtag) [$cut_to - $cut_from :" . substr($seq,$cut_to,$cut_length) ."\n";
	    $newseq = substr($seq,0,$cut_to) . ('n' x $cut_length)  . substr($seq,$cut_from);
	    
	    $sequence{$name} = $newseq;
	    
	}
	
	print  OUT_mRNA ">$name $name\n$sequence{$name}\n\n";
	
    }
    
    # reset carriage return char to normal
    $/="\n";
    # close file handles
    close (OUT_mRNA);
}

#################################################


sub make_ncrnas {

    my $count = 0;
    
    print "// Fetching ncRNA sequences\n";
    
    my $ncRNA_file = "elegans_ncRNAs.masked";
    my $name;
    my $seq;

    open (OUT_ncRNA, ">$dir/$ncRNA_file");

    # assign tace query for this data set 
    my $command = "quiet -on\nquery find sequence where RNA AND NOT mRNA\ndna\nquit\n";
    # set carriage return char to >
    $/=">";
    # tace call to the database to dump dna
    open (SEQUENCES, "echo '$command' | $tace -silent $ace_dir |");
    while (<SEQUENCES>) {	
	chop;
	
	# remove acedb lines...
	s/acedb//g;

	# end loop if last message
	last if (/bientot/);

	# increment count
	$count++;        
	# short-circuit header info
	next if ($count < 5);

	# print sequence to file
	($name) = (/^(\S+)/);
	$seq = $_;
	$seq =~ s/$name//g;
 	$seq =~ s/\W//g;
	

	print OUT_ncRNA ">$name $name\n$seq\n";
    }

    # close file handles
    close (SEQUENCES);

    # reset carriage return char to normal
    $/="\n";

    # close file handles
    close (OUT_ncRNA);
}

#############################################################

sub make_nematode_ests {

    my $count = 0;
    
    print "// Fetching nematode sequences\n";
    
    my $nematode_file = "other_nematode_ESTs";
    my $name;
    my $seq;

    open (OUT_nematode, ">$dir/$nematode_file");

    # assign tace query for this data set 
    my $command = "quiet -on\nquery find sequence where method = EST_nematode\ndna\nquit\n";
    # set carriage return char to >
    $/=">";
    # tace call to the database to dump dna
    open (SEQUENCES, "echo '$command' | $tace -silent $ace_dir |");
    while (<SEQUENCES>) {	

	print;
	chomp;
	
	# remove acedb lines...
	s/acedb//g;

	# end loop if last message
	last if (/bientot/);

	# increment count
	$count++;        
	# short-circuit header info
	next if ($count < 5);

	# print sequence to file
	($name) = (/^(\S+)/);
	$seq = $_;
	$seq =~ s/$name//g;
 	$seq =~ s/\W//g;
	

	print OUT_nematode ">$name $name\n$seq\n";
    }

    # close file handles
    close (SEQUENCES);

    # reset carriage return char to normal
    $/="\n";

    # close file handles
    close (OUT_nematode);
}

#############################################################

sub make_embl_cds {

    my $count = 0;
    
    print "// Fetching EMBL CDS sequences\n";
    
    my $embl_file = "elegans_embl_cds";
    my $name;
    my $seq;

    open (OUT_nematode, ">$dir/$embl_file");

    # assign tace query for this data set 
    my $command = "quiet -on\nquery find sequence where method = NDB_CDS\ndna\nquit\n";
    # set carriage return char to >
    $/=">";
    # tace call to the database to dump dna
    open (SEQUENCES, "echo '$command' | $tace -silent $ace_dir |");
    while (<SEQUENCES>) {	

	print;
	chomp;
	
	# remove acedb lines...
	s/acedb//g;

	# end loop if last message
	last if (/bientot/);

	# increment count
	$count++;        
	# short-circuit header info
	next if ($count < 5);

	# print sequence to file
	($name) = (/^(\S+)/);
	$seq = $_;
	$seq =~ s/$name//g;
 	$seq =~ s/\W//g;
	

	print OUT_nematode ">$name\n$seq\n";
    }

    # close file handles
    close (SEQUENCES);

    # reset carriage return char to normal
    $/="\n";

    # close file handles
    close (OUT_embl);
}



#################################################

sub get_masking_data {

    print "// Fetching masking data\n";
    
    my $count = 0;
    my ($seqtag,$target,$start,$stop,$span,$desc,$seq,$type);

    # assign tace query for this data set 
    my $command = "Table-maker -p /nfs/disk100/wormpub/DATABASES/current_DB/wquery/FeatureMasker.def\nquit\n";
    # tace call to the database to dump dna
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {

	# Exclusions from the TableMaker output
	next if (/inverted/);
	next if (/tandem/);
	next if (/intron/);
	next if (/polyA_site/);
	next if (/polyA_signal_sequence/);

	s/\"//g;

#	print "// $_\n" if ($debug);

	if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+.+)/) {
	    ($seqtag,$target,$type,$start,$stop,$span,$desc) = ($1,$2,$3,$4,$5,$6,$7);

	    # output data to array and hash
	    push (@{$bmasked{$target}},$type);
	    $mask{$seqtag}{SEQ}   = $target;
	    $mask{$seqtag}{START} = $start;
	    $mask{$seqtag}{STOP}  = $stop;
	    $mask{$seqtag}{DESC}  = $desc;

	    print "$bmasked{$target} [$target = '$seqtag'] : $mask{$seqtag}{START} $mask{$seqtag}{STOP} $mask{$seqtag}{DESC} \n"  if ($debug);
	    
	}
    }
    # close file handles
    close (TACE);
    
}


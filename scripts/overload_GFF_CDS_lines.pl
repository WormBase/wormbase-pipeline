#!/usr/local/bin/perl5.8.0
#
# overload_GFF_CDS_lines.pl
# 
# by Dan Lawson
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2008-03-28 10:56:32 $
#

#
#    1. Brief_identification
#    2. Locus names when available
#    3. WormPep IDs
#    4. Confirmed_by status
#    5. WBGene IDs
#    6. Genbank accesions (for region:Genomic_canonical)
#
#	        [1]                                                                      [3]                    [2]             [4]                   [5]
# CDS "JC8.10a" ; Note "Inositol polyphosphate phosphatase, catalytic domain homologues" ; WormPep "WP:CE28239" ; Note "unc-26" ; Confirmed_by "cDNA" ; Gene "WBGene00006763"
#

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($chrom, $splitdir, $gffdir, $release);

GetOptions (
            "help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "chrom:s"    => \$chrom,
	    "release:s"  => \$release,
	    "splits:s"   => \$splitdir,
	    "gff:s"      => \$gffdir,
	    "species:s"  => \$species	    
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
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
$release = $wormbase->get_wormbase_version;

#################################
# Set up some useful paths      #
#################################


$splitdir           = $wormbase->gff_splits unless $splitdir;
$gffdir             = $wormbase->gff unless $gffdir;
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP

# check that the supplied $release is numeric
unless ($release =~ /\d+/) {
    exit(1);
}

print "// Working with wormpep release $release\n" if ($verbose);


# get data from wormpep release
my $CDS;
my %wormpep;
my %geneID;
my %status;
my $line;
my %locus;
my %briefID;
my %RNAgenes;
my %seqname2geneid;

#&get_RNA_info;
&get_wormpep_info;

$wormbase->FetchData("worm_gene2geneID_name",\%seqname2geneid);
# test output from wormpep
#foreach $i (sort keys %wormpep) {
#    print "CDS: $i\t$wormpep{$i}\t$geneID{$i}\t$locus{$i}\t$status{$i}\t\'$briefID{$i}\'\n";
#}

# parse GFF lines
my @gff_files;

if (defined($chrom)){
    push(@gff_files,$chrom);
}
else {
    @gff_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
}

foreach my $file (@gff_files) {
    
    # Check for existance of GFFfile
	$log->log_and_die("$file doesnt exist\n") unless (-e "$gffdir/$file.gff") ;
 
    open (OUT, ">$gffdir/${file}.CSHL.gff");

    open (GFF, "<$gffdir/${file}.gff")  || die "Cannot open $file.gff\n";
    while (<GFF>) {
	chomp;
	
		#skip header lines of file
		unless  (/^\S+\s+(curated|miRNA|ncRNA|snRNA|snlRNA|snoRNA|tRNAscan-SE-1\.23|snl)\s+(\w+primary_transcript|CDS)/) {
		    print OUT "$_\n";
		    next;
		}
	
		my ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$name) = split /\t/;
	
		if( $source eq 'curated') {
			my ($i) = $name =~ (/CDS \"(\S+)\"/);
	
			print OUT "$chromosome\t$source\t$feature\t$start\t$stop\t$score\t$strand\t$other\tCDS \"$i\" ;";
			print OUT " Note \"$briefID{$i}\" ;"        if ($briefID{$i} ne "");
			print OUT " WormPep \"".$wormbase->pep_prefix.":$wormpep{$i}\" ; " if ($wormpep{$i} ne "");
			print OUT " Note \"$locus{$i}\" ; "         if ($locus{$i} ne "");
			print OUT " Status \"$status{$i}\" ; "      if ($status{$i} ne "");
			print OUT " Gene \"$geneID{$i}\" ; "        if ($geneID{$i} ne "");
	    }
	    else {
	    	#non-coding genes
	    	my ($transcript) = $name =~ (/Transcript \"(\S+)\"/);
	    	print OUT "$chromosome\t$source\t$feature\t$start\t$stop\t$score\t$strand\t$other\tTranscript \"$transcript\" ;";
	    	print OUT " Note \"".$RNAgenes{$transcript}->{'remark'}."\" ; " if $RNAgenes{$transcript}->{'remark'} ;
	    	print OUT " Note \"".$RNAgenes{$transcript}->{'locus'}."\" ; "  if $RNAgenes{$transcript}->{'locus'} ;;
	    	print OUT " Gene \"".$seqname2geneid{$transcript}."\" ; "  if $seqname2geneid{$transcript} ;
	    }
		print OUT "\n";
	}
    close GFF; #_ end of input GFF file
    close OUT; #_ end of output GFF file
}


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub get_RNA_info {
	my $rna_file = $wormbase->wormrna."/wormrna".$wormbase->get_wormbase_version.".rna";
	open (RNA,"<$rna_file") or $log->log_and_die("cant open $rna_file\t$!\n");
	while(<RNA>) {
		#I couldnt think of a way to do this in one regex!
		my ($locus, $remark, $cds);
		if(/locus:(\S+)/){
			$locus = $1;
			/>(\S+)\s+(.*)\s+locus/;
			$cds = $1;
			$remark = $2;
		}
		elsif(/>(\S+)\s+(.*)$/) {
			$cds = $1;
			$remark = $2;
		}
		$RNAgenes{$cds}->{'remark'} = $remark if $remark;
		$RNAgenes{$cds}->{'locus'} = $locus if $locus;
	}
	close RNA;
}


sub get_wormpep_info {
	my $file = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep".$wormbase->get_wormbase_version;
	open (WORMPEP, "<$file") or $log->log_and_die("cant open $file $!\n");
	while (<WORMPEP>) {
	    if (/^>(\S+)\s+(\S+)\s+(WBGene\d+)\s+status\:(\S+)/) {
			$CDS           = $1;
			$wormpep{$CDS} = $2;
			$geneID{$CDS}  = $3;
			$status{$CDS}  = $5;
		
			$line = $4;
		
			if ($line =~ /locus:(\S+)\s+(\S+.+)/) {
			    $locus{$CDS}   = $1;
			    $briefID{$CDS} = $2;
			}
			elsif ($line =~ /locus:(\S+)$/) {
			    $locus{$CDS}   = $1;
			}
			else {
			    $briefID{$CDS} = $line;
			}
	    }
	    elsif (/^>(\S+)\s+(\S+)\s+(WBGene\d+)\s+status\:(\S+)/) {
			$CDS           = $1;
			$wormpep{$CDS} = $2;
			$geneID{$CDS}  = $3;
			$status{$CDS}  = $4;
	    }
    
	}
	close WORMPEP;
}




sub usage {
 my $error = shift;
 
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }

}

=pod

=head2 NAME - overload_GFF_CDS_lines.pl

=head1 USAGE

=over 4

=item  overload_GFF_CDS_lines.pl [-options]

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

=item -test, Test mode, run the script, but don't change anything.

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

=item Unknown

=back

=cut

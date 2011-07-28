#!/usr/local/bin/perl5.8.0 -w
#
# GFFmunger.pl
# 
# by Dan Lawson
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2011-07-28 16:01:23 $
#
# Usage GFFmunger.pl [-options]


#################################################################################
# variables                                                                     #
#################################################################################

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
use Ace;

##################################################
# Script variables and command-line options      #
##################################################
my ($help, $debug, $test, $verbose, $store, $wormbase,$gffdir,$datadir);

my $all;       # All of the following:
my $landmark;  #   Landmark genes
my $motifs;    #   Protein motifs mapped down to genome
my $gmap2pmap; #   Physical positions for genes based on interpolation using GMap
my $UTR;       #   UTRs 
my $WBGene;    #   WBGene spans
my $CDS;       #   CDS overload
my $rnai;      #   RNAi
my $chrom;     # single chromosome mode
my $version;

GetOptions (
  "all"        => \$all,
  "landmark"   => \$landmark,
  "motifs"     => \$motifs,
  "gmap2pmap"  => \$gmap2pmap,
  "UTR"        => \$UTR,
  "CDS"        => \$CDS,
  "chrom:s"    => \$chrom,
  "gff:s"      => \$gffdir,
  "splits:s"   => \$datadir,
  "release:s"  => \$version,
  "help"       => \$help,
  "debug=s"    => \$debug,
  "test"       => \$test,
  "verbose"    => \$verbose,
  "store:s"    => \$store,
  'rnai'       => \$rnai,
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


# get version number
my $WS_version;

if ($version) {
  $WS_version = $version;
} else {
  $WS_version = $wormbase->get_wormbase_version;
}



##############################
# Paths etc                  #
##############################

$datadir ||= $wormbase->gff_splits;
$gffdir  ||= $wormbase->gff;
my @file_prefices;

# prepare array of file names and sort names

if (defined($chrom)){
  push(@file_prefices,$chrom);
}
else {
  if ($wormbase->assembly_type eq 'contig'){
    @file_prefices = ($wormbase->species);
  } else {
    @file_prefices = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
  }
}

# check to see if full chromosome gff dump files exist
foreach my $filep (@file_prefices) {
  unless (-e "$gffdir/$filep.gff") {
    &usage("No GFF file");
  }
  if (-e -z "$gffdir/$filep.gff") {
    &usage("Zero length GFF file");
  }
}


my $gffpath;


#################################################################################
# Main Loop                                                                     #
#################################################################################


if ($CDS || $all) {
  $log->write_to("# Overloading CDS lines\n");
  if (defined($chrom)){
    $log->write_to("overload_GFF_CDS_lines.pl -chrom $chrom -gff $gffdir\n");
    $wormbase->run_script("overload_GFF_CDS_lines.pl -chrom $chrom -gff $gffdir", $log); # generate *.CSHL.gff files
    
  } else {
    $log->write_to("overload_GFF_CDS_lines.pl -gff $gffdir\n");
    $wormbase->run_script("overload_GFF_CDS_lines.pl -gff $gffdir", $log);
  }

  foreach my $filep (@file_prefices) {
    next if not $filep;
    system ("mv -f $gffdir/$filep.CSHL.gff $gffdir/$filep.gff");        # copy *.CSHL.gff files back to *.gff name
  }
}

############################################################
# remove unwanted lines from the GFF files
############################################################

$log->write_to("# Removing unwanted lines\n");
if (defined($chrom)){
  $log->write_to("remove_unwanted_GFF_lines.pl -chrom $chrom -gff $gffdir\n");
  $wormbase->run_script("remove_unwanted_GFF_lines.pl -chrom $chrom -gff $gffdir", $log);
  
} else {
  $log->write_to("remove_unwanted_GFF_lines.pl -gff $gffdir\n");
  $wormbase->run_script("remove_unwanted_GFF_lines.pl -gff $gffdir", $log);
}


############################################################
# loop through each GFF file                               #
############################################################

foreach my $filep (@file_prefices) {

  next if ($filep eq "");               # end loop if no filename

  $gffpath = "$gffdir/${filep}.gff";
  $log->write_to("# File $filep\n");
    
  my @addfiles;
  
  if (($landmark || $all) && ($filep ne "CHROMOSOME_MtDNA") && ($wormbase->species eq 'elegans')) {
    push @addfiles, "$datadir/${filep}_landmarks.gff";
  }

  if (($gmap2pmap or $all) and $wormbase->species eq 'elegans') {
    push @addfiles, "$datadir/${filep}_gmap2pmap.gff";
  }

  if ($motifs or $all) {
    push @addfiles, $wormbase->assembly_type eq 'contig' ? "$datadir/proteinmotifs.gff" : "$datadir/${filep}_proteinmotifs.gff";
  }

  if ($UTR || $all) {
    push @addfiles, $wormbase->assembly_type eq 'contig' ? "$datadir/UTR.gff" : "$datadir/${filep}_UTR.gff";
  }

  foreach my $addfile (@addfiles) {
    if (-e $addfile) {
      $log->write_to("# Adding $addfile file\n");
      &addtoGFF($addfile,$gffpath);
    } else {
      $log->log_and_die("Could not find $addfile - badness\n");
    }
  }

  unless ($filep =~ (/MtDNA/)) {
    &check_its_worked($gffpath) if $all;
  }
}

##################
# Check the files
##################


if ($wormbase->assembly_type eq 'contig') {
  my ($file) = @file_prefices;
  my $prefix = $wormbase->chromosome_prefix;
  $wormbase->check_file("$gffdir/${file}.gff", $log,
			minsize => 1500000,
			lines => ['^##',
				  "^$prefix\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			);

} else {
  foreach  my $file (@file_prefices) {
    my $minsize = ($file=~/random|un/)?350000:1500000;
    $wormbase->check_file("$gffdir/${file}.gff", $log,
			  minsize => $minsize,
			  lines => ['^##',
				    "^$file\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			  );
  }
}

# Tidy up
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);



###############################
# subroutines                 #
###############################

sub check_its_worked {
	my $file = shift;
        my $fiveprime = qx{grep five_prime  $file | wc -l };
	my $partially = qx{grep Partially  $file | wc -l };
	
	my $msg;

        if (($landmark || $all) && ($file ne "CHROMOSOME_MtDNA") && ($wormbase->species eq 'elegans')) {
           my $landmark_genes = qx{grep landmark $file | wc -l }; #qx{} captures system command output.

	   if ($landmark_genes < 10) {
		$msg .= "landmark genes are not present\n";
	   }
        }
	if ( $fiveprime  < 10 ) {
		$msg .= "UTRs are not present\n";
	}
	if ( $partially  < 10 ) {
		$msg .= "CDS overloading not present\n";
	}	
	
	if( defined ($msg) ) {
		$log->write_to("GFFmunging failed : $msg");
		$log->error;
	}
	else {
		$log->write_to("GFFmunging appears to have worked\n");
	}
}

sub addtoGFF {
    
    my $addfile = shift;
    my $GFFfile = shift;
    
    system ("cat $addfile >> $GFFfile") && warn "ERROR: Failed to add $addfile to the main GFF file $GFFfile\n";

}


##########################################
sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  elsif ($error eq "No GFF file") {
      # No GFF file to work from
      print "One (or more) GFF files are absent from $gffdir\n\n";
      exit(0);
  }
  elsif ($error eq "Zero length GFF file") {
      # Zero length GFF file
      print "One (or more) GFF files are zero length. The GFF dump may not have worked\n\n";
      exit(0);
  }
  elsif ($error eq "Debug") {
    # No debug person named
    print "You haven't supplied your name\nI won't run in debug mode until I know who you are\n\n";
    exit (0);
  }
}

################################################
#
# Post-processing GFF routines
#
################################################



__DATA__
clone_path
CDS
WBGene
pseudogenes
transposon
rna
worm_genes
Coding_transcript
coding_exon
exon
exon_tRNA
exon_pseudogene
exon_noncoding
intron
intron_all
Genefinder
history
repeats
assembly_tags
TSL_site
polyA
oligos
RNAi
TEC_RED
SAGE
allele
clone_ends
PCR_products
cDNA_for_RNAi
BLAT_EST_BEST
BLAT_EST_OTHER
BLAT_OST_BEST
BLAT_OST_OTHER
BLAT_TRANSCRIPT_BEST
BLAT_mRNA_BEST
BLAT_mRNA_OTHER
BLAT_EMBL_BEST
BLAT_EMBL_OTHER
BLAT_NEMATODE
BLAT_NEMBASE
BLAT_TC1_BEST
BLAT_TC1_OTHER
BLAT_ncRNA_BEST
BLAT_ncRNA_OTHER
BLAT_WASHU
Expr_profile
BLASTX
WABA_BRIGGSAE
operon
Oligo_set
rest
__END__



=pod

=head2 NAME - GFFsplitter.pl

=back 

=head1 USAGE

=over 4

=item GFFsplitter.pl <options>

=back

This script splits the large GFF files produced during the build process into
smaller files based on a named set of database classes to be split into.
Output written to autoace/GFF_SPLITS/WSxx

=over 4

=item MANDATORY arguments: 

None.

=back

=over 4

=item OPTIONAL arguments: -help, this help page.

= item -debug <user>, only email report/logs to <user>

= item -archive, archives (gzips) older versions of GFF_SPLITS directory

= item -verbose, turn on extra output to screen to help track progress


=back


=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut

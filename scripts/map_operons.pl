#!/usr/local/bin/perl5.8.0 -w
#
# map_operons.pl

# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2003-12-16 13:49:44 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Data::Dumper;
use Getopt::Long;
use File::Copy;
use Ace;

##############################
# Script variables (run)     #
##############################

# who will receive log file?
my $maintainers = "dl1";

our $dir    = "/wormsrv2/autoace/OPERONS";
our $gff    = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
our $WS_version =  &get_wormbase_version_name;

our $output1 = "/wormsrv2/autoace/OPERONS/operons_${WS_version}.ace";
our $output2 = "/wormsrv2/wormbase/misc/misc_operons.ace";

# do everything in /wormsrv2/autoace/OPERONS
chdir("$dir");

our %operon;
our %cds;
our %est;
our $errors = 0;

my $verbose;           # verbose mode
my $help;              # Help/Usage page
my $debug;             # debug mode


GetOptions (
            "verbose"        => \$verbose,
            "help"           => \$help,
            "debug:s"        => \$debug
	    );
 
# Help pod if needed
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
    print "// DEBUG = \"$debug\"\n\n";
    ($maintainers = $debug . '\@sanger.ac.uk');
}

# Set up logfile
my $rundate    = &rundate;
our $log = "/wormsrv2/logs/map_operons.pl.WS${WS_version}.${rundate}.$$";
open (LOG,">$log");
LOG->autoflush();

# print logfile header
print LOG "# map_operons.pl\n";

&recreate_hashes;

print "// create output file\n\n" if ($verbose);

open (OUTPUT, ">$output1") || die "Can't open file for output\n";

&acedump_operons;

close OUTPUT;

print "// end output file\n\n" if ($verbose);

# copy this file to correct place

my $status = copy($output1, $output2);

if ($status == 0) {
    print "// Failed to copy file to $output2\n"; 
    print LOG "Failed to copy file: $!\n";   
    $errors++;
}

my $subject_line = "map_operons.pl";

# warn about errors in subject line if there were any
if ($errors == 1) {
  $subject_line .= " : $errors ERROR";
}
elsif ($errors > 1) {
  $subject_line .= " : $errors ERRORS";
}

print "// mail log file\n\n" if ($verbose);

&mail_maintainer("$subject_line",$maintainers,$log);

# hasta luego

exit(0);





##################################################################################################################

sub recreate_hashes {

    open (FH, "<$dir/operon.dat") or die "operon.dat : $!";
    undef $/;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}

sub acedump_operons {
    
    my $operon_start;
    my $operon_stop;
    my $reset;
    my $gene_count;
    my @f;


    for my $operon_lookup (sort keys %operon) {
	
	print OUTPUT "\n";
	print OUTPUT "Operon : \"$operon_lookup\"\n";
	print OUTPUT "Species \"Caenorhabditis elegans\"\n";

	print "// Dump operon $operon_lookup\n";
	
	$gene_count = 1;
	foreach my $gene_lookup (@{$operon{$operon_lookup}{CDS}}) {
	    
#	    print "// $gene_lookup\n";

	    if ($gene_count == 1) {
		$operon_start = 0;
		open (GFF_1, "grep -w '$gene_lookup' $gff/CHROMOSOME_${operon{$operon_lookup}->{CHROMOSOME}}.genes.gff |");
		while (<GFF_1>) {
		    @f = split /\t/;
		    if ($f[6] eq "+") {$operon_start = $f[3];}
		    if ($f[6] eq "-") {$operon_start = $f[4];}
		}
		close GFF_1;
#		print "OP_start = $operon_start\n";
	    }
	    
	    if ($gene_count == $operon{$operon_lookup}->{NO_GENES}) {
		$operon_stop = 0;
		open (GFF_2, "grep -w '$gene_lookup' $gff/CHROMOSOME_${operon{$operon_lookup}->{CHROMOSOME}}.genes.gff |");
		while (<GFF_2>) {
		    @f = split /\t/;
		    ($operon_stop = $f[4]) if ($f[6] eq "+");
		    ($operon_stop = $f[3]) if ($f[6] eq "-");
		}
		close GFF_2;
#		print "OP_stop = $operon_stop\n";
	    }
	    
	    if (scalar (@{$cds{$gene_lookup}{SL1_EST}}) > 1) {
		print OUTPUT "Contains_CDS \"$gene_lookup\" SL1 \"EST clones @{$cds{$gene_lookup}{SL1_EST}}\"\n";
		$reset = 1;
	    }
	    if (scalar (@{$cds{$gene_lookup}{SL2_EST}}) > 1) {
		print OUTPUT "Contains_CDS \"$gene_lookup\" SL2 \"EST clones @{$cds{$gene_lookup}{SL2_EST}}\"\n";
		$reset = 1;
	    }
	    if ($cds{$gene_lookup}{SL2_MICROARRAY} eq "++") {
		print OUTPUT "Contains_CDS \"$gene_lookup\" SL2 \"Microarray strong\"\n";
		$reset = 1;
	    }
	    elsif ($cds{$gene_lookup}{SL2_MICROARRAY} eq "+") {
		print OUTPUT "Contains_CDS \"$gene_lookup\" SL2 \"Microarray weak\"\n";
		$reset = 1;
	    }
	    
	    unless ($reset == 1) {
		print OUTPUT "Contains_CDS \"$gene_lookup\"\n";
	    }
	    
	    $reset = 0;
	    $gene_count++;
	}
	print OUTPUT "Method operon\n";
	print OUTPUT "\n";
	print OUTPUT "Sequence : \"CHROMOSOME_$operon{$operon_lookup}->{CHROMOSOME}\"\n";
	print OUTPUT "S_child Operon $operon{$operon_lookup}->{ACC} $operon_start $operon_stop\n";
	print OUTPUT "\n";
	
	if (($operon_start == 0) || ($operon_stop == 0)) {
	    print LOG "\n//ERROR: MISSING COORDINATE for operon $operon{$operon_lookup}->{ACC}\n\n";
	    $errors++;
	}
	
	$gene_count = 1;
	undef ($operon_start);
	undef ($operon_stop);
    }
    
}

sub usage {
    my $error = shift;
    if ($error eq "Help") {
	exec ('perldoc',$0);
    }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
    }


}




__END__

=pod

=head2   NAME - map_operons.pl

=head1 USAGE

=over 4

=item autoace_minder.pl 

=back

map_operons.pl will generate an .acefile for the operon data using the 
operon.dat hash in /wormsrv2/autoace/OPERONS/.


map_operons.pl mandatory arguments:

=over 4

=item none,

=back

map_operons.pl optional arguments:

=over 4

=item none,

=back

=cut


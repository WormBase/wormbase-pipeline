#!/nfs/team71/worm/mh6/bin/perl
#===============================================================================
#
#         FILE:  fix_gff.pl
#
#        USAGE:  ./fix_gff.pl
#
#  DESCRIPTION:  to fix gffs with interpolated map positions
#
#      OPTIONS:  -help, -test, -debug USER, -snps, -genes, -clones, -all, -store FILE
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  $Author: mh6 $ (), <>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  13/02/06 09:37:00 GMT
#     REVISION:  $Revision: 1.10 $
#===============================================================================

# BACS / SNPS / GENEs

use strict;
use lib $ENV{'CVS_DIR'};
use Modules::Physical_map;
use Wormbase;
use Getopt::Long;
use IO::File;

my ( $store, $test, $prep, $debug, $alleles, $genes, $clones, $all, $help, $wormbase,$chromosome );    #options

GetOptions(
    'help'    => \$help,
    'test'    => \$test,
    'debug:s' => \$debug,
    'alleles'    => \$alleles,
    'genes'   => \$genes,
    'clones'  => \$clones,
    'all'     => \$all,
    'store:s' => \$store,
    'chromosome:s'=> \$chromosome,
    'prepare' => \$prep
);

die `perldoc $0` if $help;

if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test ) }

my $log         = Log_files->make_build_log($wormbase);    # prewarning will be misused in a global way
$log->{SCRIPT}="$0 : [@ARGV]";

my $acedb    = $wormbase->{'autoace'};
my $chromdir = $wormbase->{'gff_splits'};
my $outdir   = "$acedb/acefiles/";
##############################################################################################
#generate a new mapper based on the files (also needs to point to better gene/gene files)
my $mapper = Physical_mapper->new( $acedb, glob("$chromdir/CHROMOSOME_*_gene.gff") );

# check the mappings
if ($prep) {
	$mapper->check_mapping($log,$acedb);
	$log->mail( );
	exit;
}

###############################################################################################
my $rev_genes = Map_func::get_phys($acedb);                # hash of all rev_map genes

# snp/gene/whatever loop over gff
$log->write_to("\n\ngenerating acefiles:\n");
$log->make_line;

my @chromosomes=('CHROMOSOME_I', 'CHROMOSOME_II', 'CHROMOSOME_III', 'CHROMOSOME_IV', 'CHROMOSOME_V', 'CHROMOSOME_X');
@chromosomes=("CHROMOSOME_$chromosome",) if $chromosome;

foreach my $chrom ( @chromosomes) {
###################################################################################
# $chrom,$chromdir,$snps,$genes,$clones,$rev_genes
#
# hmm .... log?
# 
# SHOULD really go into a method/function
#

	
    # Input files
    my @data;
    push( @data, "$chromdir/${chrom}_allele.gff" )          if ( $alleles|| $all ); # GFF_method_dump.pl -method Allele
    push( @data, "$chromdir/${chrom}_gene.gff" )            if ( $genes  || $all );
    push( @data, "$chromdir/${chrom}_clone_acc.gff" ) if ( $clones || $all );
    foreach my $file (@data) {
        $chrom=~/_(.*)$/;
	
	&dump_alleles($wormbase,$1) if ($alleles && (! -e $file));
	
   	my $fh = IO::File->new( $file, "r" ) || ( $log->write_to("cannot find: $file\n") && next);
	$file =~ /${chrom}_(allele|gene|clone)/; # hmpf
        my $of = IO::File->new( "> $outdir/interpolated_$1_$chrom.ace" );

	$log->write_to("writing to: interpolated_$1_$chrom.ace\n");

        while (<$fh>) {
            next if /\#/;
            s/\"//g;
            my @fields = split;

            # dumb assumption that f[9] is always the id
            my ( $chr, $source, $feature, $id , $ctag) = ( $fields[0], $fields[1], $fields[2], $fields[9],$fields[8]);

            my $class;
            if    ( $source eq 'Genomic_canonical' && $feature eq 'region' ) { $class = 'Sequence' }
            elsif ( $source eq 'Allele'            && $ctag eq 'Variation' ) { $class = 'Variation' }
            elsif ( $source eq 'gene'              && $feature eq 'gene' )   { $class = 'Gene';
                next if $rev_genes->{$id} # need to check for existing reverse maps for genes
            }
            else { next }
	    
	    my $pos = ( $fields[3] + $fields[4] ) / 2; # average map position
            my $aceline = $mapper->x_to_ace( $id, $pos, $chr, $class );
	    
            print $of $aceline if $aceline; # mapper returns undef if it cannot be mapped (like on the telomers)
            $log->write_to("cannot map $class : $id (might be on a telomer) - phys.pos $chr : $pos\n") if ( !$aceline ); #--
        }
	
        close $of;
        close $fh;
    }
}
###########################################################################################
$log->mail();

exit 0;

###############################
sub dump_alleles {
	my ($wormbase,$chromosome)=@_;
#	my $cmd = "GFF_method_dump.pl -database ".$wormbase->autoace." -method Allele -dump_dir ".$wormbase->autoace."/GFF_SPLITS -chromosome $chromosome";
	my $cmd = "grep Allele ".$wormbase->autoace."/CHROMOSOMES/CHROMOSOME_$chromosome.gff >".$wormbase->autoace."/GFF_SPLITS/CHROMOSOME_${chromosome}_allele.gff";

	print `$cmd`
}

package Log_files;

sub make_line {
	my ($self)=@_;
	$self->write_to("-----------------------------------\n");
}

package Physical_mapper;

sub x_to_ace {
    my ( $self, $id, $map, $chr, $x ) = @_;
    $chr =~ s/CHROMOSOME_//;
    my $mpos = $self->map( $map, $chr );
    if ($mpos) { return "$x : \"$id\"\nInterpolated_map_position \"$chr\" $mpos\n\n" }
}

sub check_mapping {
	my ($self,$logger,$acedir)=@_;

	# print it to some dark place
	my $revh = IO::File->new("$acedir/logs/rev_phys.log","w");	
	$logger->write_to("have a look at $acedir/logs/rev_phys.log to resolve:\n");
	$logger->make_line;

	foreach my $key ( keys %{$self->{pmap}} ) { 
		my $last;
		foreach my $i(@{ $self->{smap}->{$key}}){
 			if ($last&&($self->{pmap}->{$key}->{$i}->[0] < $self->{pmap}->{$key}->{$last}->[0])){
				print $revh "----------------------------------\n";
				print $revh "$key: $i\t",$self->{pmap}->{$key}->{$i}->[0],"\t",$self->{pmap}->{$key}->{$i}->[1]," ERROR (conflict with last line) \n";
				$logger->write_to("$key: $i\t".$self->{pmap}->{$key}->{$i}->[0]."\t".$self->{pmap}->{$key}->{$i}->[1]." ERROR (conflict with last line) \n");
				print $revh "----------------------------------\n"
			}
			else {	print $revh "$key: $i\t",$self->{pmap}->{$key}->{$i}->[0],"\t",$self->{pmap}->{$key}->{$i}->[1],"\n"}

			$last=$i
		}
	}
}
__END__

=head1 NAME

=head1 DESCRIPTION

creates mapping files in AceDB format for alleles, clones and genes by using 
the physical and genetic position of reference genes.

=head1 USAGE 

=over 2

=item * -prep checks and prepares the reverse physical mapping

=item * -help

=item * -test

=item * -debug NAME

=item * -alleles

=item * -genes

=item * -clones

=item * -all	all of the above snp+genes+clones

=item * -store FILENAME	specifies a stored Wormbase configuration

=item * -chromosome [I II III IV V X] specify ONE chromosome to process

=back 

=head1 DEPENDENCIES

=over 2

=item * Modules::Physical_map

=item * Wormbase

=back 





#!/software/bin/perl
#===============================================================================
#
#         FILE:  make_UTR_GFF.pl
#
#        USAGE:  ./make_UTR_GFF.pl
#
#  DESCRIPTION:
#
#      OPTIONS:  [-help, -verbose, -store FILE, -debug USER, -test, -noload, -chromosome CHROMOSOME_NUMBER]
# REQUIREMENTS:  Wormbase.pm, Modules::GFF_sql.pm
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  $Author: ar2 $
#      COMPANY:
#      VERSION:  2 
#      CREATED:  21/02/06 14:11:30 GMT
#     REVISION:  $Revision: 1.18 $ 
#===============================================================================

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Modules::GFF_sql;
use Wormbase;
use Getopt::Long;
use IO::File;

my ( $help, $debug, $test, $store, $wormbase, $chrom, $load , $verbose);

GetOptions(
    "help"         => \$help,
    "debug=s"      => \$debug,
    "test"         => \$test,
    "store:s"      => \$store,
    "chromosome:s" => \$chrom,
    "noload"       => \$load,
    'verbose'	   => \$verbose
);

die `perldoc $0` if $help;
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test ) }

my $log = Log_files->make_build_log($wormbase);    # prewarning will be misused in a global way

# do a single chromosome if prompted else do the lot......
my @chromosome;
if ($chrom) {
	@chromosome = split(/,/,join(',',$chrom));
}
else {
    @chromosome = $wormbase->get_chromosome_names(-mito => 1);
}

# global setup setup

my $gffdir = $wormbase->gff_splits;
my $db     = GFF_sql->new( { -build => 1 } );
my $outdir = $gffdir;
my %cds_cache;  # crude hack to speed up the cds lookup
my $contig_assembly = ($wormbase->assembly_type eq 'contig') ? 1 : undef;
my $suffix = $contig_assembly ? $$ : '';
# main
CHROM:foreach my $chr (@chromosome) {
    %cds_cache=(); # clean out old data
    my $pref = $wormbase->chromosome_prefix;
    $chr =~ s/$pref//; #seems silly but is safest way to ensure that its not missing or duplicated!
    my $file_prefix = $contig_assembly ? $wormbase->species : $wormbase->chromosome_prefix.$chr;
   	my $long_name = $wormbase->chromosome_prefix.$chr;
   		
   	&write_tmp_gff($file_prefix,$long_name) if $contig_assembly;
   	
    # load    
    $db->load_gff( "$gffdir/${file_prefix}_curated.gff$suffix", $long_name, 1 ) if !$load;
    $db->load_gff( "$gffdir/${file_prefix}_Coding_transcript.gff$suffix", $long_name ) if !$load;

    my $outfile = IO::File->new( "$outdir/${file_prefix}_UTR.gff$suffix",          ">>" ); # append or create
    my $infile  = IO::File->new( "$gffdir/${file_prefix}_Coding_transcript.gff$suffix", "r" );

    my $n_exons=0;
    
    # iterate over exons GFF
    while (<$infile>) {

        next if /\#/;
        s/\"//g;#"
        my @f = split;
        my ( $chrom, $start, $stop, $ori, $name ) = ( $f[0], $f[3], $f[4], $f[6], $f[9] );
        next if ( $f[1] ne 'Coding_transcript' || $f[2] ne 'exon' );
	print "processing transcript: $name \n" if $verbose;
        #make the gff
        print $outfile &make_gff( $db, $chrom, $start, $stop, $ori, $name );
	$n_exons++;
    }

    $log->write_to("processed $n_exons exons\n");
    $wormbase->run_command("cat $outdir/${file_prefix}_UTR.gff$$ >> $outdir/UTR.gff$$", 'no_log');
    
    #clean up
    $wormbase->run_command("rm $gffdir/".$wormbase->species."*$$",'no_log');
}

$log->mail();

##########################

#
sub write_tmp_gff {
	my $prefix = shift;
	my $contig = shift;
	open(OUT,">$gffdir/${prefix}_curated.gff$suffix") or $log->log_and_die("cant make curated tmp GFF for $prefix $contig: $!\n");
	my $handle = $wormbase->open_GFF_file($contig,'curated', $log);
	while (<$handle>) {
		print OUT;
	}
	close $handle;
	close OUT;
	open(OUT,">$gffdir/${prefix}_Coding_transcript.gff$suffix") or $log->log_and_die("cant make Coding tmp GFF for $prefix $contig: $!\n");
	$handle = $wormbase->open_GFF_file($contig,'Coding_transcript', $log);
	while (<$handle>) {
		print OUT;
	}
	close $handle;
	close OUT;
}


# get CDS name (lifted from original version)
sub short_name {
    my ($name) = @_;
    my $cds_regex = $wormbase->cds_regex;
    my ($cdsname) = $name =~ /($cds_regex)/;
    return $cdsname;
}

# UTR GFF line generator
sub utr {
    my ( $chrom, $name, $start, $stop, $type, $ori ) = @_;
    return "$chrom\tCoding_transcript\t${type}_prime_UTR\t$start\t$stop\t.\t$ori\t.\tTranscript \"$name\"\n";
}

# CDS GFF line generator
sub cds {
    my ( $chrom, $start, $stop, $orientation, $frame, $name, $fluff ) = @_;
    return "$chrom\tCoding_transcript\tcoding_exon\t$start\t$stop\t.\t$orientation\t$frame\tTranscript \"$name\" ; $fluff\n";
}

# get CDS start/stop coordinates
sub get_cds {
    my ( $db, $chrom, $id ) = @_;

    if ( $cds_cache{ &short_name($id) } ) { return @{ $cds_cache{ &short_name($id) } } } # hmpf
    else {
        my @hits = $db->get_chr( $chrom, { feature => 'curated', source => 'coding_exon', fluff => '"' . &short_name($id) . '"' } );
        my @sorted_hits = sort { $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} } @hits;
	$cds_cache{ &short_name($id)}=[$sorted_hits[0]->{start}, $sorted_hits[-1]->{stop}];
        return $sorted_hits[0]->{start}, $sorted_hits[-1]->{stop};
    }
}

# creates the gff-lines and returns them
sub make_gff {
    my ( $db, $chrom, $start, $stop, $orientation, $name ) = @_;
    
    # get_hits , get cds
    my @hits =
      $db->get_chr( $chrom, { start => $start, stop => $stop, feature => 'curated', source => 'coding_exon', 'fluff' => '"' . &short_name($name) . '"' } );
    my @cds = get_cds( $db, $chrom, &short_name($name) );

    if (@hits) {
        if ( $start < $cds[0] || $stop > $cds[1] ) {    # need to create both utrs and cds exons
            my ( $type, $utr, $cds );
	    if ( $start < $cds[0] && $stop > $cds[1]){
		my @types = ( $orientation eq '+' ) ? ('five','three') : ('three','five');
		#utr1
		my $utr1=&utr( $chrom, $name, $start, $cds[0] - 1, $types[0], $orientation );
		#cds
		$cds = &cds( $chrom, $cds[0], $cds[1], $orientation, $hits[0]->{'frame'}, $name, $hits[0]->{'fluff'} );
		#utr2
		my $utr2 = &utr( $chrom, $name, $cds[1] + 1, $stop, $types[1], $orientation );
		return $utr1,$cds,$utr2;
	    }
            elsif ( $start < $cds[0] ) {
                $type = ( $orientation eq '+' ) ? 'five' : 'three';
                $utr = &utr( $chrom, $name, $start, $cds[0] - 1, $type, $orientation );
                $cds = &cds( $chrom, $cds[0], $stop, $orientation, $hits[0]->{'frame'}, $name, $hits[0]->{'fluff'} );
            }
            else { # has to be stop > cds[1]
                $type = ( $orientation eq '+' ) ? 'three' : 'five';
                $utr = &utr( $chrom, $name, $cds[1] + 1, $stop, $type, $orientation );
                $cds = &cds( $chrom, $start, $cds[1], $orientation, $hits[0]->{'frame'}, $name, $hits[0]->{'fluff'} );
            }
            return $utr, $cds;
        }
        else {    # is cds only => only cds exon
            return &cds( $chrom, $start, $stop, $orientation, $hits[0]->{'frame'}, $name, $hits[0]->{'fluff'} );
        }
    }
    else {        # is utr only => only utr exon
        my $type;
        if    ( $start < $cds[0] ) { $type = ( $orientation eq '+' ) ? 'five'  : 'three' }
        elsif ( $stop > $cds[1] )  { $type = ( $orientation eq '+' ) ? 'three' : 'five' }
        else { $log->write_to("urghs: $name start/stop CDS coordinates seem to be messed up\n") }
        return &utr( $chrom, $name, $start, $stop, $type, $orientation );
    }
}

__END__


=head1 NAME 

make_UTR_GFF.pl

=head1 USAGE

make_UTR_GFF.pl [-test -debug NAME -noload -chromosome CHROMOSOME_NUMBER -store STOREFILE]

=head1 DESCRIPTION

creates a GFF file per chromosome containing UTRs and coding exons (connected to CDS and Transcripts) from   _curated.gff and _Coding_transcript.gff (in GFF_split).
It does that by creating an SQL table for the GFF files in the build mySQL database (can theoretically overwrite each other as there are no locks on the tables, but only if they run on the same chromosome).
Mixed exons (partially coding) will be split into 2 exons for the GFF (one coding/one noncoding).

=head1 Arguments

=over 

=item -store : loads from a stored configuration file

=item -chromosome C_NAME: limits the processing to one chromosome (like I)

=item -test : will use the TEST_BUILD directories

=item -debug NAME : mail logs to NAME as email

=item -noload : does not reload the GFF-database (faster startup)

=item -verbose : prints out which transcript it is processing

=back

=head1 Output

GFF_SPLITS/CHROMOSOME_nn_UTR.gff

=head1 Database

=over

=item mySQL 4.x

=item host mcs2a

=item name mh6_build

=item port 3316

=back

=head1 Dependencies

=over

=item Wormbase.pm

=item Modules/GFF_sql.pm

=back

=cut

#!/usr/bin/perl -w
#===============================================================================
#         FILE:  build_webdb.pl
#
#        USAGE:  ./build_webdb.pl 
#
#  DESCRIPTION:  build Todd's GFF databases
#
#      $AUTHOR:$
#      COMPANY:  WormBase
#      CREATED:  11/27/09 10:10:08 GMT
#      CHANGED: $Date: 2011-02-09 10:44:32 $
#    $Revision: 1.8 $
#===============================================================================

# need to set the PERl5LIb  to pick up bioperl for the load_gff thing ... maybe have to hardcode it

use Getopt::Long;
use lib $ENV{CVS_DIR};
use Wormbase;
use strict;

my @species;
my %tmpfiles;
my ($debug,$currentDB,$ftp,$version);

GetOptions(
      'species:s' => \@species,
      'debug:s'   => \$debug,
      'currentDB' => \$currentDB, # use current_DB instead of the build databases
      'ftp'       => \$ftp, # grab it from the ftp-site
      'version:s' => \$version, # grab a specific version from the FTP site
)||die(@!);

# postprocess the collected GFF files per species

@species = qw(elegans) unless $species[0];

## for each species
foreach my $worm(@species){
    my $wormbase = Wormbase->new(
        -organism => $worm,
        -debug    => $debug,
    );

    if ($currentDB){
       $wormbase = Wormbase->new(
        -organism => $worm,
        -debug    => $debug,
        -autoace  => glob('~wormpub/DATABASES/current_DB'),
       );
    }

    $wormbase->version($version) if $version; # override the default one

    my $log = Log_files->make_build_log($wormbase);

    $$wormbase{log}=$log; # crude, but works
    $$wormbase{ftp}=1 if $ftp; # even cruder

    # munge the GFF collection
    $log->write_to("... creating munged GFFs\n");
    &process_worm($wormbase);

    # deletes the tmp directory
    &clean_tmpfiles;

    $log->mail();
}


print "Hasta Luego\n";

# tars up the database directory and copies it to the ftp server
sub tar_and_feather {
    my ($wb)=@_;
    my $sp = $wb->species;
    
    my $from = "/tmp/$sp.gff";
    my $to = $wb->ftp_site .'/WS'. $wb->version .'/genomes/'.$wb->full_name(-g_species => 1)
             .'/genome_feature_tables/GFF2/'.$wb->full_name(-g_species => 1).'.WS'.$wb->version.'GBrowse.gff.gz';
    
    system("gzip -9 $from") && die(@!);
    system("cp /$from.gz $to") && die(@!);
}

# takes GFFs from $wormbase and creates a big one in /tmp/$species.gff
sub process_worm {
    my ($wb)=@_;

    my @raw_gffs;
    my @gz_gffs;
    if ($$wb{ftp}){
        # ~ftp/pub2/wormbase/WS209/genomes/c_elegans/genome_feature_tables/SUPPLEMENTARY_GFF/*.gff
        
        # commented out, as we probably don't need them
        
        # @raw_gffs = glob($wb->ftp_site .'/WS'. $wb->version .'/genomes/'.$wb->full_name(-g_species => 1)
        # .'/genome_feature_tables/SUPPLEMENTARY_GFF/*.gff');

        # ~ftp/pub2/wormbase/WS209/genomes/c_elegans/genome_feature_tables/GFF2/c_elegans.WS209.gff.gz
        @gz_gffs = glob($wb->ftp_site .'/WS'. $wb->version .'/genomes/'.$wb->full_name(-g_species => 1)
        .'/genome_feature_tables/GFF2/'.$wb->full_name(-g_species => 1).'.WS'.$wb->version.'.gff.gz');
    }
    else {
         my $buildDataDir = '/nfs/wormpub/BUILD_DATA/SUPPLEMENTARY_GFF';
         @raw_gffs = glob($wb->chromosomes.'/*.gff');
         @gz_gffs  = glob($wb->chromosomes.'/*.gff.gz');
         push @gz_gffs, glob("$buildDataDir/*.gff2.gz");
         push @raw_gffs,glob("$buildDataDir/*.gff");
    }

    $$wb{log}->write_to("... processing @gz_gffs\n");
    $$wb{log}->write_to("... processing @raw_gffs\n");

    my $gffile = concatenate_gff($wb,\@gz_gffs,\@raw_gffs);

    process_gff($wb->species,$gffile);
}

# concatenate the files
sub concatenate_gff {
    my ($wb,$gz,$gff)=@_;
    my $outfile = "/tmp/${\$wb->species}.gff_inf";

    &tmpfile_hook($outfile);

    # gzipped ones
    if ($$gz[0]){
      $gz = join(' ',@$gz);
      system("zcat $gz>$outfile") && die(@!);
    }        

    # normal ones
    if ($$gff[0]){
      $gff = join(' ',@$gff);
      system("cat $gff>>$outfile") && die(@!);
    }

    return $outfile;
}

# just the processing bit
sub process_gff {
    my($species,$file)=@_;

    &tmpfile_hook("/tmp/$species.gff");

    print STDERR "/software/bin/perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl ".
           ">/tmp/$species.gff <$file\n" if $debug;
    system("/software/bin/perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl ".
           ">/tmp/$species.gff <$file") && die(@!);
}

# unlink the file and add it to the final cleanup
sub tmpfile_hook {
    my ($tmpfile)=@_;
    unlink $tmpfile if -e $tmpfile;
    $tmpfiles{$tmpfile}=1;
}

# pack it up for Todd

sub clean_tmpfiles {
    foreach my $key (keys %tmpfiles){
         unlink($key) if (-e $key);
    }
}

# just in case, you never know
END {
    &clean_tmpfiles();
}

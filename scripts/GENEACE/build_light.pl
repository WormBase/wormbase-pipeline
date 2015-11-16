#!/software/bin/perl -w
#
# build_light
# 
#  Last updated on: $Date: 2010-07-14 14:37:58 $
#  Last updated by: $Author: pad $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Storable;
use Log_files;

my ($wormbase,$debug,$test,$store,$all,$quick,$sourcedb,$noinit);

GetOptions (
	    "store"      => \$store,
	    "test"       => \$test,
	    "debug:s"    => \$debug,
	    "all"        => \$all,
	    "quick"      => \$quick,
	    "source:s"   => \$sourcedb,
	    "noinit"     => \$noinit,
	   );

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);
my $dbdir = "/nfs/wormpub/DATABASES/build_light_full";
unless (defined $sourcedb) {
  $sourcedb = "/nfs/wormpub/BUILD/elegans";
}
unless ($all || $quick){ $log->log_and_die("You must specify all or quick command line options or else you dont have any data to load\n");}

# create the necessary dir structure
unless ($noinit) {
  if (-e "$dbdir/wspec/models.wrm") {
    $wormbase->run_command("rm -rf $dbdir", $log) && die "Failed to remove the database $dbdir\n";
  }
  unless (-e "$dbdir/wspec/models.wrm") {
    $wormbase->run_command("mkdir $dbdir", $log) && die "Failed to create directory $dbdir/wspec/models.wrm\n";
    $wormbase->run_command("cp -rf /nfs/wormpub/wormbase-pipeline/wspec $dbdir/", $log) && die "Failed to copy the wspec directory\n";
    $wormbase->run_command("mkdir $dbdir/database", $log) && die "Failed to create directory database\n";
  }
}

####################################
# What data should we load         #
####################################

my @filenames;
if (($all) || ($quick)){
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Sequence.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Pseudogene.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Transcript.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_DNA.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_CDS.ace";      
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Genetic_code.ace";
  push @filenames, "$sourcedb/acefiles/primaries/geneace/geneace_Gene.ace";
  push @filenames, "$sourcedb/acefiles/primaries/geneace/geneace_Gene_class.ace";
  push @filenames, "/nfs/wormpub/wormbase/autoace_config/misc_autoace_methods.ace";

}

if ($all){
  push @filenames, "$sourcedb/acefiles/primaries/citace/caltech_Person.ace";
  push @filenames, "$sourcedb/acefiles/primaries/geneace/geneace_Feature.ace";
  push @filenames, "$sourcedb/acefiles/primaries/geneace/geneace_Variation.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_EST_features.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Feature.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Transposon.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Clone.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Transposon_fam.ace"; 
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_NDB_features.ace";
  push @filenames, "$sourcedb/acefiles/primaries/camace/camace_Motif.ace";
  push @filenames, "$sourcedb/acefiles/map_alleles_output.ace";
  push @filenames, "$sourcedb/acefiles/WBgene_spans.ace";
  push @filenames, "$sourcedb/acefiles/repeat_homologies.ace";
  push @filenames, "$sourcedb/acefiles/inverted_repeats.ace";
  push @filenames, "$sourcedb/acefiles/feature_binding_site.ace";                      
  push @filenames, "$sourcedb/acefiles/feature_histone_binding_site.ace";              
  push @filenames, "$sourcedb/acefiles/feature_promoter.ace";                          
  push @filenames, "$sourcedb/acefiles/feature_TF_binding_site.ace"; 
  push @filenames, "$sourcedb/acefiles/feature_binding_site_region.ace";               
  push @filenames, "$sourcedb/acefiles/feature_histone_binding_site_region.ace";       
  push @filenames, "$sourcedb/acefiles/feature_regulatory_region.ace";                 
  push @filenames, "$sourcedb/acefiles/feature_TF_binding_site_region.ace"; 
  push @filenames, "$sourcedb/acefiles/feature_Corrected_genome_sequence_error.ace";   
  push @filenames, "$sourcedb/acefiles/feature_micro_ORF.ace";                         
  push @filenames, "$sourcedb/acefiles/feature_segmental_duplication.ace";             
  push @filenames, "$sourcedb/acefiles/feature_three_prime_UTR.ace"; 
  push @filenames, "$sourcedb/acefiles/feature_DNAseI_hypersensitive_site.ace";        
  push @filenames, "$sourcedb/acefiles/feature_polyA_signal_sequence.ace";             
  push @filenames, "$sourcedb/acefiles/feature_SL1.ace";                               
  push @filenames, "$sourcedb/acefiles/feature_transcription_end_site.ace"; 
  push @filenames, "$sourcedb/acefiles/feature_Genome_sequence_error.ace";             
  push @filenames, "$sourcedb/acefiles/feature_polyA_site.ace";                        
  push @filenames, "$sourcedb/acefiles/feature_SL2.ace";                               
  push @filenames, "$sourcedb/acefiles/feature_transcription_start_site.ace"; 
  push @filenames, "$sourcedb/acefiles/elegans_blastx.ace";
  push @filenames, "$sourcedb/acefiles/elegans_blastp.ace";
}


####################################
# Re-initialise the ACEDB database #
####################################

unless ($noinit) {
  $log->log_and_die("*Reinitdb error - lock.wrm file present..\n") if (-e "$dbdir/database/lock.wrm");
  
  if( -e "$dbdir/database/log.wrm" ) {
    unlink glob("$dbdir/database/new/*") or $log->write_to("ERROR: Couldn't unlink file $dbdir/database/new/ : $!\n");
    unlink glob("$dbdir/database/touched/*") or $log->write_to( "ERROR: Couldn't unlink file $dbdir/database/touched/ : $!\n");
    
    if( -e "$dbdir/database/log.wrm") {
      my $status = move("$dbdir/database/log.wrm", "$dbdir/database/log.old");
      print "ERROR: Couldn't move file: $!\n" if ($status == 0);
    }
    unlink glob("$dbdir/database/*.wrm") or $log->write_to("ERROR: Couldn't run rm command $dbdir/database/*.wrm: $!\n");
  }
  else {
    eval{
      mkdir("$dbdir/database");
    };
    $log->log_and_die(@!) if @!;
  }
  my $command="y\n";
  $log->write_to("* Reinitdb: reinitializing the database ..\n");
  &DbWrite($command,$wormbase->tace,$dbdir,"ReInitDB");
}

##############################
# Upload the .ace files      #
##############################

my $filename;
foreach $filename (@filenames) { 

  # .ace_files can have URLs in the .ace, so quote '//'
  #my $tmp = "$filename.tmp";
  #system("/bin/sed 's#http:\\/\\/#http:\\\\/\\\\/#g' < $filename > $tmp");
  #$wormbase->run_command("mv -f $tmp $filename", $log);
  
  
  my $command=<<END;
pparse $filename
save 
quit
END
  if (-e $filename) {
    $log->write_to( "Loading $filename into $dbdir\n");
    &DbWrite($command,$wormbase->tace,$dbdir,"ParseFile");
  } else {
    $log->write_to("$filename is not existent - skipping ..\n");
    next;
  }
}


sub DbWrite {
  my ($command,$exec,$dir,$name,$logfile)=@_;
  open (WRITEDB,"| $exec $dir ") or $log->log_and_die("$name DbWrite failed\n");
  print WRITEDB $command;
  close WRITEDB;
}


#close things down and end
print "Diaskeda same Poli\n"; #we had alot of fun#
$log->mail();
exit(0);
__END__

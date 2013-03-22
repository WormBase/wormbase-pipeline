#!/software/bin/perl -w
#
# fasta_dumper.pl
# 
# by Paul Davis
#
# Script to dump a fasta file for a specified species/class/method
# Option to add WBGeneIDs to the header of the file.
# >2L52.1 gene=WBGene00007063
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2013-03-22 15:11:05 $
use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Wormbase;
use File::Path;
use Storable;

my ($debug, $store, $verbose, $database, $test, $wormbase, $species, $class, $method, $out_file, $target_dir, $geneid, $process, $version);

GetOptions ( "debug:s"      => \$debug,      #email a specified user only, so you don't flood all your colleagues with test run spam.
	     "verbose"      => \$verbose,    #verbose quces a little more info to screen
	     "class:s"      => \$class,      #specify an acedb class
	     "method:s"     => \$method,     #specify an acedb method
	     "database:s"   => \$database,   #specify a database to query
	     "output:s"     => \$out_file,   #specify an output file name
	     "outdir:s"     => \$target_dir, #specify an output location or else the file will be written to BUILD/species/SEQUENCES/
	     "test"         => \$test,       #invoke test env
	     "store:s"      => \$store,      #supply a storable
	     "species:s"    => \$species,    #needed to work out what species is being processed
	     "gene"         => \$geneid,     #adds the WBGeneID to the header line (retrieved via a wormpep mapping)
	     "processonly"  => \$process,    #process the specified output file without all the delay of dumping
	     "version:s"    => \$version,    #adds flexibility 
	   ) ;



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species,
                           );
}
# Create log file
my $log = Log_files->make_build_log($wormbase);

# Define database to be queried
my $db;
if (not defined $database or $database eq "autoace") {
  $db = $wormbase->autoace;
}else{$db = $database}

my $runtime = $wormbase->runtime;
my $tace = $wormbase->tace;
my (%gene_ids);

unless (defined $class) {$log->write_to("!!!!Warning - No class specified so defaulting to Transcript\n");
			 $class = "Transcript";}
unless (defined $method) {$log->write_to("!!!!Warning - No method specified so defaulting to Coding_transcript\n");
			  $method = "Coding_transcript";}

unless (defined $target_dir) {$target_dir = $wormbase->sequences;}
unless (defined $out_file){$out_file     = "coding_transcripts.dna";}
if ($geneid) {
  unless (defined $version) {$log->write_to("The current build wormpep will be used for WBGene data\n")}
}
my $out = "$target_dir/$out_file";
my $target= "${out}.mod";
my $gspecies = $wormbase->full_name('-g_species' => 1);
my $full_name = $wormbase->full_name;

############
# Main body
############


if ($process) {
  $log->write_to("\n$runtime: Only processing $out as -processonly option used on command line\n\n");
}

unless ($process) {
  $log->write_to("\n$runtime: Started making $class($method) fasta file:\n$out\n");
  if (-e "$out") {
    $log->write_to("\nOld $out exists in this location, removing.......\n");
    $wormbase->run_command("rm -f $out", $log);
  }
  open FAOUT,">$out" ||die($!);

  my $object;
  my $connection = Ace->connect(-path => "$db") || die (Ace->error);
  my $object_it = $connection->fetch_many(-query => "Find $class; Species=\"${full_name}\"; Method=\"$method\";");
  while(my $object=$object_it->next){
    if ($verbose){print $object->name."\n";}
    my $dna =$object->asDNA();
    print FAOUT "$dna";
  }
  close FAOUT;
}

if ($geneid){
  &fetch_gene_ids;
  my $cds_regex = $wormbase->cds_regex_noend;
  if (-e $out) {
    open(out_fh, "<$out") || die "Failed to open $out\n" ;
    open(target_fh, ">$target") || die "Failed to open $target\n";
    while(<out_fh>) {
      chomp;
      s/>//;
      if (/($cds_regex)/) {
	if (exists $gene_ids{$1}) {
	  print target_fh "\n>$_ gene=$gene_ids{$1}\n";
	} else {
	  print target_fh "\n$_\n";
	}
      } elsif (/(\S+)/) {
	print target_fh "$_\n";
      }
    }
    close(target_fh) or $log->log_and_die("Could not successfully close target_fh\n");
    close(out_fh)or $log->log_and_die("Could not successfully close out_fh\n");
  }
  else {
    $log->error("ERROR: Could not find transcript file for $gspecies ($out)\n");
    die;
  }
  $wormbase->run_command("mv $out ${out}_bk", $log);
  $log->write_to("\nMoving $out -> ${out}_bk\n");
  $wormbase->run_command("mv ${out}.mod $out", $log);
  $log->write_to("\nMoving ${out}.mod -> $out\n");
}


if (-e "${out}.gz") {
  $log->write_to("\nOld zipped data exists ${out}.gz exists in this location, removing.......\n");
  $wormbase->run_command("rm -f ${out}.gz", $log);
}


$wormbase->run_command("gzip -9 -f $out", $log);
$runtime = $wormbase->runtime;
$log->write_to("$runtime: Finished making $class($method) fasta file\n\nFile location:$out\n\n");

$log->mail();
exit(0);

sub fetch_gene_ids {
  my $wsid;
  unless ($version) {
    $wsid = $wormbase->get_wormbase_version;
  }
  else {
    $wsid = $version;
  }
  my $peppre = $wormbase->pepdir_prefix;
  my $pep_dir = $wormbase->peproot;
  my $source_pepfile = "$pep_dir/${peppre}pep${wsid}/${peppre}pep${wsid}";
  if (-e $source_pepfile) {
    open(my $source_pep_fh, $source_pepfile);
    while(<$source_pep_fh>) {
      /^\>(\S+)\s+\S+\s+(\S+)/ and do {
	$gene_ids{$1} = $2;
      };
    }
    close($source_pep_fh);
  }
  else {
    $log->log_and_die("Failed to open $source_pepfile\n");
  }
}

__END__

=pod

=head2 NAME - fasta_dumper.pl

=head3 USAGE

=head2 DESCRIPTION

This script "does exactly what it says on the tin". ie it dumps fasta sequences for a given method class combination.

=head2 OPTIONS:

=over 4

B<-verbose:> 

Terminal output

B<-databases:> 

Either use autoace if used in build process or give the full database path. Basically retrieves paired read info for ESTs from that database.

B<-gene:> 

This option relies upon wormpep to get a mapping between CDS ID and WBGene ID.

B<-processonly:> 

This will take the output file and add the "gene=WBGene00012345" to the headers.

B<-class:> 

Specify an acedb class for data to be dumped

B<-method:> 

Narrows the keyset based on method

B<-output:> 

Specify a filename for the results to be saved in

B<-debug:> 

Sends the script log to the specifies user.

B<-species:> 

this option is necessary to establish the source of information

=back

=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@ebi.ac.uk)

=back

=cut

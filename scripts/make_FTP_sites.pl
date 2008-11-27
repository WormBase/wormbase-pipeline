#!/usr/local/bin/perl5.8.0 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase & wormpep FTP sites
# 
# Originally written by Dan Lawson
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2008-11-27 15:58:28 $
#
# see pod documentation (i.e. 'perldoc make_FTP_sites.pl') for more information.
#
##########################################################################################

=pod

=head1 NAME - make_FTP_sites.pl

=back


=head1 USAGE

=over 4

=item make_FTP_sites.pl

=back


This script does :

 [01] - make a new directory for the WS release
 [02] - copy the WS release files to the target directory
 [03] - make a new directory for the chromosome DNA/GFF/AGP files
 [04] - copy the chromosome DNA/GFF/AGP files to the target directory
 [05] - copy the models.wrm file across (also present in the database.*.tar.gz files)
 [06] - copy the relevant dbcomp file across
 [07] - copy across latest wormpep release
 [08] - make wormpep FTP site
 [09] - copy WormRNA release across
 [10] - extract confirmed genes from autoace and make a file on FTP site
 [11] - delete the old symbolic link and make the new one
 [12] - delete the old WS release version directory
 [13] - makes a file of cDNA2orf connections
 [14] - makes a file of all gene IDs with CGC names and Sequence names (where present)
 [15] - exit gracefully


=over 4

=item MANDATORY arguments:

none

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)


=head1 AUTHOR - Originally written by Dan Lawson (dl1@sanger.ac.uk) with extensive rewrite
and upgrade by Keith Bradnam (krb@sanger.ac.uk)

=cut
use strict;
use lib $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'}."/Modules";
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use IO::Handle;
use File::Path;



#################################################################################
# Command-line options and variables                                            #
#################################################################################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $release; # only copy across release files
my $chroms;  # only copy across chromosome files
my $ont;     # only copy across ontology files
my $misc;    # only copy misc files
my $wormpep; # only copy wormpep files
my $genes;   # only copy confirmed genes
my $cDNA;    # only copy cDNA2orf file
my $geneIDs; # only copy file of gene IDs
my $pcr;     # only copy file of PCR products
my $homols;  # only copy best blast hits 
my $manifest;# check everything has been copied.
my $all;     # copy everything across
my $dump_ko; # create the knockout consortium file
my $dna;
my $rna;
my $gff;
my $supplementary;
my $clustal;

GetOptions ("help"     => \$help,
	    "debug=s"  => \$debug,
	    "test"     => \$test,
	    "verbose"  => \$verbose,
	    "store:s"    => \$store,
	    "release"  => \$release,
	    "dna"      => \$dna,
	    "rna"      => \$rna,
	    "gff"      => \$gff,
	    "supplementary"      => \$supplementary,
	    "ont"      => \$ont,
	    "misc"     => \$misc,
	    "wormpep"  => \$wormpep,
	    "genes"    => \$genes,
	    "cDNAlist" => \$cDNA,
	    "geneIDs"  => \$geneIDs,
	    "clustal"  => \$clustal,
	    "pcr"      => \$pcr,
	    "homols"   => \$homols,
	    "manifest"=> \$manifest,
	    "knockout" => \$dump_ko,
	    "all"      => \$all);


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


# using -all option?
($clustal=$release=$dump_ko=$dna=$gff=$supplementary=$rna=$misc=$wormpep=$genes=$cDNA=$geneIDs=$pcr=$homols=$manifest=$ont = 1 ) if ($all);

my $base_dir = $wormbase->basedir;    # The BUILD directory
my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $citace_dir = $wormbase->primaries."/citace";
my $targetdir = $wormbase->ftp_site;  # default directory, can be overidden
my $WS              = $wormbase->get_wormbase_version();      # e.g.   132
my $WS_name         = $wormbase->get_wormbase_version_name(); # e.g. WS132
my $annotation_dir = "$targetdir/$WS_name/genomes/".$wormbase->full_name('-g_species' => 1)."/annotation";
my $maintainers     = "All";
my $runtime;



#################################################################################
# Main                                                                          #
#################################################################################

#if this fails its really important that the next step doesnt run - create lockat start and remove at end
# will be checked by finish_build.pl

my $lockfile = $wormbase->autoace."/FTP_LOCK";
$log->write_to("writing lock file\n");

open (FTP_LOCK,">$lockfile") or $log->log_and_die("cant write lockfile $!\n");

print FTP_LOCK "If this file exists something has failed in make_ftp_site.pl\n DO NOT continue until you know what happend and have fixed it\n\n";
close FTP_LOCK;

&copy_release_files if ($release);    # make a new directory for the WS release and copy across release files

&copy_dna_files if ($dna);  		  
&copy_rna_files if ($rna);
&copy_gff_files if ($gff);  		  
&copy_supplementary_gff_files if ($supplementary);  		  

&copy_ontology_files if ($ont);       # make a new /ONTOLOGY directory and copy files across

&copy_misc_files if ($misc);          # copy across models.wrm and other misc. files, e.g. wormRNA

&copy_wormpep_files if ($wormpep);    # copied from ~wormpub/WORMPEP

&extract_confirmed_genes if ($genes); # make file of confirmed genes from autoace and copy across

&make_cDNA2ORF_list if ($cDNA);       # make file of cDNA -> ORF connections and add to FTP site

&make_geneID_list if ($geneIDs);      # make file of WBGene IDs -> CGC name & Sequence name and add to FTP site

&make_pcr_list if ($pcr);             # make file of PCR products -> WBGene IDs, CDS, CGC name

&copy_homol_data if ($homols);        # copies best blast hits files across

&check_manifest if ($manifest);       # compares whats on the FTP site with what should be

&extract_ko if ($dump_ko);            # dumps KO-consortium data to FTP site

&copy_clustal if ($clustal);          # copies the clustal sql-dump to the FTP site

################################
#
# Tidy up and exit
#
################################


# warn about errors in subject line if there were any
$wormbase->run_command("rm -f $lockfile",$log) if ($log->report_errors == 0);
$log->mail;
print "Finished.\n" if ($verbose);
exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################


##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying release files\n");

  my $ftp_acedb_dir = "$targetdir/$WS_name/acedb";
  mkpath("$ftp_acedb_dir",1,0775);

  my $filename;

  opendir (RELEASE,"$ace_dir/release") or croak ("Could not open directory $ace_dir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    $wormbase->run_command("scp $ace_dir/release/$filename $ftp_acedb_dir/$filename", $log);

    my $O_SIZE = (stat("$ace_dir/release/$filename"))[7];
    my $N_SIZE = (stat("$ftp_acedb_dir/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir RELEASE;
  
  # Copy across the models.wrm file
  $wormbase->run_command("scp $ace_dir/wspec/models.wrm $ftp_acedb_dir/models.wrm.$WS_name", $log);

  # copy some miscellaneous files across
  my $old_release = $WS -1;
  $wormbase->run_command("scp ".	$wormbase->compare."/WS$old_release-$WS_name.dbcomp $ftp_acedb_dir", $log);
  $wormbase->run_command("scp $base_dir/autoace_config/INSTALL $ftp_acedb_dir", $log);  
  
  # change group ownership
  $wormbase->run_command("chgrp -R  worm $ftp_acedb_dir", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished\n\n");
  
}

##################################################
# copy the DNA, GFF, and agp files across
##################################################

sub copy_dna_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying dna and agp files\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $wb (values %accessors) {
    my $gspecies = $wb->full_name('-g_species'=>1);
    my $chromdir = $wb->chromosomes;

    if (-e "$chromdir") {
      my $dna_dir = "$targetdir/$WS_name/genomes/$gspecies/sequences/dna";
      mkpath($dna_dir,1,0775);
      #todd wants all species to have whole genome in one file
      if ($wb->assembly_type eq 'contig') {
	my $dna_file = "$chromdir/supercontigs.fa";
	my $species = $wb->species;
	my $masked_file = "$chromdir/".$species."_masked.dna.gz";
	my $soft_file = "$chromdir/".$species."_softmasked.dna.gz";
		
	$wormbase->run_command("cp -f $dna_file $dna_dir/".$gspecies.".$WS_name.dna.fa",$log);
	$wormbase->run_command("/bin/gzip -f $dna_dir/".$gspecies.".$WS_name.dna.fa",$log);
	$wormbase->run_command("cp -f $soft_file $dna_dir/".$gspecies."_softmasked.$WS_name.dna.fa.gz", $log);
	$wormbase->run_command("cp -f $masked_file $dna_dir/".$gspecies."_masked.$WS_name.dna.fa.gz", $log);
		
      } elsif ($wb->assembly_type eq 'chromosome') {
	$wormbase->run_command("cp -R $chromdir/*.dna* $dna_dir/", $log);
      } else {$log->error("$gspecies : unknown assembly_type\n")}
      $wb->run_command("cp -R $chromdir/*.agp $dna_dir/", $log) if (scalar glob("$chromdir/*.agp"));

      # change group ownership
      $wb->run_command("chgrp -R  worm $dna_dir", $log);
    }
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

sub copy_gff_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying gff files\n");
  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $wb (values %accessors) {
    my $species = $wb->species;
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $chromdir = $wb->chromosomes;

    if (-e "$chromdir") {
      my $gff_dir = "$targetdir/$WS_name/genomes/$gspecies/genome_feature_tables/GFF2";
      mkpath($gff_dir,1,0775);

      # we don't want to end up with thousands of files for species with DNA still in contigs, so make one file
      my @contigs = $wb->get_chromosome_names(-prefix => 1, -mito => 1);
	my $whole_filename = "$gspecies.$WS_name.gff"; #concatenated whole genome file require for all species
	$wormbase->run_command("rm -f $gff_dir/$whole_filename", $log);
	if($wb->assembly_type eq 'contig') {
	      if (-e "$chromdir/$species.gff") { # tierII does it this way
			$wormbase->run_command("cp -R $chromdir/$species.gff $gff_dir/$whole_filename", $log);
      	}
      	else {
      		$log->error("$chromdir/$species.gff missing\n");
      	}
      }
      else {
           $wormbase->run_command("cat $chromdir/*.gff* > $gff_dir/$whole_filename", $log);
           $wormbase->run_command("cp $chromdir/*.gff* $gff_dir/", $log); #individual files too
      }

	#add supplementary and nGASP GFF
	my $ngaspdir;
	if($species eq 'elegans') {
		my $supdir = $wb->build_data."/SUPPLEMENTARY_GFF";
		my @gfffiles = glob("$supdir/*.gff");
		foreach my $sup (@gfffiles){
			$wb->run_command("cat $sup >> $gff_dir/$whole_filename", $log);
		}
		$ngaspdir = $supdir;
	}
	$ngaspdir = $wb->database("$species")."/nGASP" unless $ngaspdir;
	
	#nGASP - zcat files stored under primaries (or BUILD_DATA for C_ele) on to FTP full GFF file.
	if(-e $ngaspdir){
		my @ngasp_methods = qw(augustus fgenesh jigsaw mgene);
		foreach my $method(@ngasp_methods){
			my $file = "$ngaspdir/$species.$method.gff2.gz";
			if(-e $file){
				$wb->run_command("zcat $file >> $gff_dir/$whole_filename", $log);
			}else {
				$log->error("$file missing\n");
			}
		}
	}
	else{
		$log->write_to("no ngasp for $gspecies\n");
	}
      
	$wormbase->run_command("/bin/gzip -f $gff_dir/$whole_filename",$log);
	$wormbase->run_command("cp -R $chromdir/composition.all $gff_dir/", $log) if (-e "$chromdir/composition.all");
      $wormbase->run_command("cp -R $chromdir/totals $gff_dir/", $log) if (-e "$chromdir/totals");

      # change group ownership
      $wormbase->run_command("chgrp -R  worm $gff_dir", $log);  
    }
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

sub copy_supplementary_gff_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying supplementary gff files\n");
  if($wormbase->species eq'elegans') {
	my $chromdir = $wormbase->chromosomes;
	if (-e "$chromdir/SUPPLEMENTARY_GFF") {
		my $sgff_dir = "$targetdir/$WS_name/genomes/c_elegans/genome_feature_tables/SUPPLEMENTARY_GFF";
		mkpath($sgff_dir,1,0775);
		$wormbase->run_command("cp -R $chromdir/SUPPLEMENTARY_GFF/*.gff $sgff_dir/", $log);

		# change group ownership
		$wormbase->run_command("chgrp -R  worm $sgff_dir", $log);
	}
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

sub copy_clustal{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying clustal sql dump\n");
  if($wormbase->species eq 'elegans') {
	my $chromdir = $wormbase->autoace;
	if (-e "$chromdir/wormpep_clw.sql.bz2") {
		my $sgff_dir = "$targetdir/$WS_name/";
		mkpath($sgff_dir,1,0775);
		$wormbase->run_command("cp -R $chromdir/wormpep_clw.sql.bz2 $sgff_dir/", $log);

		# change group ownership
		$wormbase->run_command("chgrp  worm $sgff_dir/wormpep_clw.sql.bz2", $log);
	}
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}
sub copy_rna_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying rna files\n");

  # run through all possible organisms
  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $wb (values %accessors) {

    my $rnadir = $wb->wormrna;

    if( -e "$rnadir") {
      my $ftprna_dir = "$targetdir/$WS_name/genomes/".$wb->full_name('-g_species' => 1)."/sequences/rna";
      mkpath($ftprna_dir,1,0775);
      $wormbase->run_command("cp -R $rnadir/* $ftprna_dir/", $log);
      chdir "$ftprna_dir" or $log->write_to("Couldn't cd $ftprna_dir\n");
      my $prefix = $wb->pepdir_prefix;
      $wormbase->run_command("/bin/tar -cf ${prefix}rna$WS.tar README ${prefix}rna$WS.rna", $log);
      $wormbase->run_command("/bin/gzip -f ${prefix}rna$WS.tar",$log);

      # change group ownership
      $wormbase->run_command("chgrp -R  worm $ftprna_dir", $log);  

    }
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}
############################################
# copy across ontology files
############################################

sub copy_ontology_files {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying ontology files\n");

  my $obo_dir = "$citace_dir/temp_unpack_dir/home/citace/Data_for_${WS_name}/Data_for_Ontology/";
  mkpath("$ace_dir/ONTOLOGY",1,0775);
  mkpath("$targetdir/$WS_name/ONTOLOGY",1,0775);

  $wormbase->run_command("cp $obo_dir/*.obo $ace_dir/ONTOLOGY", $log);
  $wormbase->run_command("cp -R $ace_dir/ONTOLOGY $targetdir/$WS_name", $log);

  # change group ownership
  $wormbase->run_command("chgrp -R  worm $ace_dir/ONTOLOGY", $log);  
  $wormbase->run_command("chgrp -R  worm $targetdir/$WS_name/ONTOLOGY", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");
}

############################################
# copy across annotation files
#############################################

sub copy_misc_files{
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying misc files\n");
  # zip and copy the microarray oligo mapping files.
  chdir "$ace_dir";
  $wormbase->run_command("/bin/gzip -f *oligo_mapping", $log);
  my $annotation_dir = "$targetdir/$WS_name/genomes/".$wormbase->full_name('-g_species' => 1)."/annotation";
  mkpath($annotation_dir,1,0775);
  $wormbase->run_command("cp $ace_dir/*oligo_mapping.gz $annotation_dir/", $log);

  # copy the compara tarball
  $wormbase->run_command("cp $base_dir/autoace/compara.tar.bz2 $targetdir/$WS_name/",$log);

  # change group ownership
  $wormbase->run_command("chgrp -R  worm $targetdir/$WS_name", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying\n\n");

}


############################################
# copy across wormpep files
#############################################

sub copy_wormpep_files {

	$runtime = $wormbase->runtime;
	$log->write_to("$runtime: copying wormpep files\n");

	my $wp_source_dir = $wormbase->wormpep;
	my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
	my $wp_ftp_dir = "$wormpep_ftp_root/wormpep$WS";
	mkpath($wp_ftp_dir,1,0775);

	&_copy_pep_files($wormbase);#elegans

	# copy wormpep file if one exists for that species.
	$log->write_to("zip and copy other species\n");
	my %accessors = ($wormbase->species_accessors);
	foreach my $species (keys %accessors){
		$log->write_to("copying $species protein data to FTP site\n");
		&_copy_pep_files($accessors{$species})
	}

	#need to update ftp/databases/wormpep for website links
	if($wormbase->species eq 'elegans') {
		$log->write_to("updating ftp/database/wormpep//\n");
		my @wormpep_files = $wormbase->wormpep_files;
		my $WS = $wormbase->get_wormbase_version;
		my $source = $wormbase->basedir . "/WORMPEP/".$wormbase->pepdir_prefix."pep$WS";
		foreach my $file ( @wormpep_files ){
			my $sourcefile = "$source/$file$WS";
			$wormbase->run_command("cp $sourcefile $wp_ftp_dir/$file$WS", $log);
			&CheckSize("$sourcefile","$wp_ftp_dir/$file$WS");
			$wormbase->run_command("ln -sf $wp_ftp_dir/$file$WS $wormpep_ftp_root/$file", $log);
		}
	}

	
	$runtime = $wormbase->runtime;
	$log->write_to("$runtime: Finished copying\n\n");

	# change group ownership
	$wormbase->run_command("chgrp -R  worm $wp_ftp_dir", $log);
}

sub _copy_pep_files {
	my $wb = shift;
	my $source = $wb->basedir . "/WORMPEP/".$wb->pepdir_prefix."pep$WS";
	my $target = "$targetdir/$WS_name/genomes/".$wb->full_name('-g_species' => 1)."/sequences/protein";
	mkpath($target,1,0775);

	# if wwormpep has not been rebuilt it may need to be carried from previous build.  Possibly done earlier but if not do it here.
	unless (-e $source) {
		require CarryOver;
		my $carrier = CarryOver->new($wb, $log);
		$carrier->carry_wormpep($WS,$wb->version);
	}
	# tar up the latest wormpep release and copy across (files added in next loop)
	my $tgz_file = "$source/".$wb->pepdir_prefix."pep$WS.tar.gz";
	my $command = "tar -c -z -h -C \"$base_dir/WORMPEP/\" -f $tgz_file";
	
	# copy em over
	my @wormpep_files = $wb->wormpep_files;
	foreach my $file ( @wormpep_files ){
		my $sourcefile = "$source/$file$WS";
		$wb->run_command("cp $sourcefile $target/$file$WS", $log);
		&CheckSize("$sourcefile","$target/$file$WS");
		$command .= " $sourcefile";
	}
	$wb->run_command("$command", $log);
	$wb->run_command("cp $tgz_file $target", $log);

	#single gzipped fasta file for Todd
	my $WS_name = $wormbase->get_wormbase_version_name;
	my $pepfile = "$target/".$wb->pepdir_prefix."pep$WS";
	my $toddfile ="$target/".$wb->pepdir_prefix."pep.$WS_name.fa.gz";
	$wb->run_command("gzip -c $pepfile > $toddfile", $log);

	# change group ownership
	$wormbase->run_command("chgrp -R  worm $target", $log);
}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting confirmed genes\n");

  my $db = Ace->connect(-path => "$ace_dir/") || die print "ACE connection failure: ", Ace->error;
  my $query = "Find elegans_CDS; Confirmed";
  my @confirmed_genes   = $db->fetch(-query=>$query);
  mkpath("$annotation_dir",1,0775);

  open(OUT,">$annotation_dir/confirmed_genes.$WS_name") or $log->write_to("Couldn't write to $annotation_dir/confirmed_genes.$WS_name\n");

  foreach my $seq (@confirmed_genes){
    my $dna = $seq->asDNA();
    print OUT "$dna";
  }

  close(OUT);
  $wormbase->run_command("/bin/gzip -f $annotation_dir/confirmed_genes.$WS_name", $log);

  $db->close;
  $wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished extracting\n\n");
  return(0);
}


#######################################
# extracting knockout consortium data
#######################################

sub extract_ko {

	$runtime = $wormbase->runtime;
  	$log->write_to("$runtime: Extracting Knockout Consortium Data\n");
	mkpath("$annotation_dir",1,0775);

	my $outfile= "$annotation_dir/knockout_consortium_alleles.$WS_name.xml.bz2";

	$wormbase->run_script("dump_ko.pl -file $outfile",$log);
	$wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  

	$runtime = $wormbase->runtime;
  	$log->write_to("$runtime: Finished dumping\n\n");
}


################################################################################
# make list of cDNA -> orf connections
################################################################################

sub make_cDNA2ORF_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making cDNA2ORF files\n");
  # simple routine to just get cDNA est names and their correct ORFs and make an FTP site file
  # two columns, second column supports multiple ORF names

  my $tace = $wormbase->tace;
  my $command=<<EOF;
Table-maker -p "$ace_dir/wquery/cDNA2CDS.def" quit
EOF

  my $dir = "$ace_dir";

  my %cDNA2orf;
  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    chomp;
    if (m/^\"/){
      s/\"//g;
      m/(.+)\s+(.+)/;     
      $cDNA2orf{$1} .= "$2 ";
    }
  }
  close(TACE);
  
  # output to ftp site

  
  my $out = "$annotation_dir/cDNA2orf.$WS_name";
  mkpath("$annotation_dir",1,0775);
  open(OUT, ">$out") or $log->write_to("Couldn't open $out\n");

  foreach my $key (keys %cDNA2orf){
    print OUT "$key,$cDNA2orf{$key}\n";
  }
  close(OUT);

  $wormbase->run_command("/bin/gzip -f $out", $log);
  $wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making files\n\n");  
  
}

################################################################################
# make list of WBGene IDs to CGC name and Sequence name
################################################################################

sub make_geneID_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making Gene ID list\n");
  # For each 'live' Gene object, extract 'CGC_name' and 'Sequence_name' fields (if present)

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/gene2cgc_name_and_sequence_name.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$annotation_dir/geneIDs.$WS_name";
  mkpath("$annotation_dir",1,0775);

  open (OUT, ">$out") || croak "Couldn't open $out\n";
  open (TACE, "echo '$command' | $tace $dir | ") or $log->write_to("cant access $dir\n");  
  while (<TACE>){
      if (m/^\"/){
	  s/\"//g;
	  tr/\t/,/;
	  print OUT "$_";
      }
  }
  close(TACE);
  close(OUT);

  $wormbase->run_command("/bin/gzip -f $out", $log);
  $wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making list\n\n");
}


################################################################################
# make list of PCR_product connections to CDS and Gene ID plus CGC name
################################################################################

sub make_pcr_list {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making PCR product 2 gene list list\n");

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/pcr_product2gene.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$annotation_dir/pcr_product2gene.$WS_name";
  mkpath("$annotation_dir",1,0775);

  # hashes needed because one pcr product may hit two or more genes
  my %pcr2gene;
  my %gene2sequence;
  my %gene2cgc;

  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    if (m/^\"/){
      s/\"//g;
      chomp;
      # split the line into various fields
      my ($pcr,$gene,$cgc,$sequence) = split(/\t/, $_) ;

      # fill hashes and arrays, allow for multiple genes per PCR product
      $pcr2gene{$pcr}      .= "$gene,";
      ($gene2cgc{$gene} = $cgc) if($cgc);
      $gene2sequence{$gene} = "$sequence";
    }
  }
  close(TACE);



# Now write output, cycling through list of pcr products in %pcr2gene

  open (OUT, ">$out") or $log->log_and_die("Couldn't open $out\n");

  foreach my $pcr (keys %pcr2gene){
    my @genes = split(/,/,$pcr2gene{$pcr});
    my $counter =0;
    print OUT "$pcr";
    foreach my $gene (@genes){

      # remove next element if it is the same gene ID and start loop again
      if (defined($genes[$counter+1]) && $genes[$counter+1] eq $genes[$counter]){
	splice(@genes, $counter+1,1);
	redo;
      }
      $counter++;
      print OUT "\t$gene";
      print OUT "($gene2cgc{$gene})" if (exists($gene2cgc{$gene}));
      
      # now print sequence name
      print OUT ",$gene2sequence{$gene}";
      
    }
    print OUT "\n";
  }

  close(OUT);

  $wormbase->run_command("/bin/gzip -f $out", $log);
  $wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making list\n\n");
}

############################################################
# copy best blast hits file to ftp site
############################################################

sub copy_homol_data {

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying homol files\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $wb (values %accessors) {
    my $blast_dir = $wb->acefiles;
    my $species = $wb->species;
    my $source_file = "$blast_dir/${species}_best_blastp_hits";
    if(-e $source_file || -e "$source_file.gz") { # this script might be run more than once if there are problems
      my $protein_dir = "$targetdir/$WS_name/genomes/".$wb->full_name('-g_species'=>1)."/sequences/protein";
      mkpath($protein_dir,1,0775);
      my $target_file = "$protein_dir/best_blastp_hits_${species}.$WS_name.gz";
      $wormbase->run_command("/bin/gzip -f $source_file",$log) if (! -e "$source_file.gz"); # this script might be run more than once if there are problems
      $wormbase->run_command("scp $source_file.gz $target_file", $log);
    }
    $runtime = $wormbase->runtime;
    $log->write_to("$runtime: Finished copying\n\n");
  } 

  # change group ownership
  $wormbase->run_command("chgrp -R  worm $targetdir/$WS_name", $log); 

}

########################

sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    $log->write_to("\tERROR: $first SRC: $F_SIZE TGT: $S_SIZE - not same size, please check\n");
  } 
}

##################################################################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##################################################################################


sub check_manifest {
	my $rel = $wormbase->get_wormbase_version;
	my $ftp_dir = $wormbase->ftp_site."/WS$rel";
	my $path = $ftp_dir;
	my $count;
	my %otherspecies = $wormbase->species_accessors;
	$otherspecies{$wormbase->species} = $wormbase;
	my $gspecies =0;
	while(<DATA>) {
		next unless /\w/;
		next if/^#/;
		chomp;
		s/REL/$rel/g;

		if(/\.\/(\S+)/) {
			$path = "$ftp_dir/$1";
			if($path =~ /gspecies/){
				$gspecies = 1;
			}
			else {
				$gspecies = 0;
			}
			next;
		}

		my $file = $_;
		if($gspecies ==1) {
			foreach my $worm (keys %otherspecies){
				my $gspecies = $otherspecies{$worm}->full_name('-g_species'=>1);
				my $gWORM = $otherspecies{$worm}->pepdir_prefix;
				my $species = $otherspecies{$worm}->species;
				my $wormpath = $path;
				my $wormfile = $file;
				$wormpath =~ s/gspecies/$gspecies/;
				$wormfile =~ s/gspecies/$gspecies/;
				$wormfile =~ s/WORM/$gWORM/;
				$wormfile =~ s/SPECIES/$species/;
				$count += &checkfile("$wormpath/$wormfile");
			}
		}
		else {
			$count += &checkfile("$path/$file");
		}
	}
	$log->write_to("$count files in place on FTP site\n");
}

sub checkfile {
	my $file = shift;
	if( -s "$file") {
      	$log->write_to("OK: $file\n");
      	return 1;
	} else {
      	$log->error("ERROR: $file has a problem\n");
      	return 0;
      }
}

#for manifest checks
# gspecies = $wb->full_name('-gspecies' => 1)
# WORM = $wb->pepdir_prefix

__DATA__
./acedb
files_in_tar
md5sum.WSREL
models.wrm.WSREL

#elegans specific stuff
./genomes/c_elegans/annotation
affy_oligo_mapping.gz
agil_oligo_mapping.gz
pcr_product2gene.WSREL.gz
geneIDs.WSREL.gz
gsc_oligo_mapping.gz
cDNA2orf.WSREL.gz
confirmed_genes.WSREL.gz

./genomes/c_elegans/sequences/dna
CHROMOSOME_I.agp
CHROMOSOME_I.dna.gz
CHROMOSOME_II.agp
CHROMOSOME_II.dna.gz
CHROMOSOME_III.agp
CHROMOSOME_III.dna.gz
CHROMOSOME_III_masked.dna.gz
CHROMOSOME_II_masked.dna.gz
CHROMOSOME_IV.agp
CHROMOSOME_IV.dna.gz
CHROMOSOME_IV_masked.dna.gz
CHROMOSOME_I_masked.dna.gz
CHROMOSOME_MtDNA.dna.gz
CHROMOSOME_V.agp
CHROMOSOME_V.dna.gz
CHROMOSOME_V_masked.dna.gz
CHROMOSOME_X.agp
CHROMOSOME_X.dna.gz
CHROMOSOME_X_masked.dna.gz
intergenic_sequences.dna.gz

./genomes/gspecies/sequences/dna
gspecies.WSREL.dna.fa.gz

./genomes/gspecies/sequences/protein
WORMpepREL.tar.gz
WORMpep.WSREL.fa.gz
best_blastp_hits_SPECIES.WSREL.gz

./genomes/gspecies/sequences/rna
WORMrnaREL.tar.gz

./genomes/c_elegans/genome_feature_tables/GFF2
CHROMOSOME_I.gff
CHROMOSOME_II.gff
CHROMOSOME_III.gff
CHROMOSOME_IV.gff
CHROMOSOME_MtDNA.gff
CHROMOSOME_V.gff
CHROMOSOME_X.gff
composition.all
totals

./genomes/gspecies/genome_feature_tables/GFF2
gspecies.WSREL.gff.gz

./genomes/c_elegans/genome_feature_tables/SUPPLEMENTARY_GFF
RNAz.gff
genemark.gff
mSplicer_orf.gff
mSplicer_transcript.gff
miranda.gff
pictar.gff
elegans_curation_anomalies.gff

./genomes/c_briggsae/genome_feature_tables/GFF2
chrI.gff
chrII.gff
chrIII.gff
chrIII_random.gff
chrII_random.gff
chrIV.gff
chrIV_random.gff
chrI_random.gff
chrUn.gff
chrV.gff
chrV_random.gff
chrX.gff

./ONTOLOGY
anatomy_association.WSREL.wb
anatomy_ontology.WSREL.obo
gene_association.WSREL.wb.cb
gene_association.WSREL.wb.ce
gene_association.WSREL.wb
gene_ontology.WSREL.obo
phenotype_association.WSREL.wb
phenotype_ontology.WSREL.obo




#!/usr/bin/env perl
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2015-03-09 16:31:34 $

use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Log_files;
use Storable;
use Bio::SeqIO;
use Net::FTP;

my ($test, $debug);
my ($fly, $yeast, $human, $uniprot, $swissprot, $trembl, $cleanup, $all, $default);
my $store;
my ($species, $qspecies, $nematode);

GetOptions (
	    'debug:s'     => \$debug,
	    'test'        => \$test,
	    'store:s'     => \$store,
	    'species:s'   => \$species,
	    'fly'	  => \$fly,
	    'yeast' 	  => \$yeast,
	    'human'	  => \$human,
	    'uniprot'	  => \$uniprot,
            'swissprot'   => \$swissprot,
            'trembl'      => \$trembl,
	    'cleanup'     => \$cleanup,
	    'all'         => \$all,
            'default'     => \$default,
	    );

my $wormbase;
if( $store ) {
    $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test     => $test,
			       -organism => $species
			       );
}

$species = $wormbase->species;   #for load
my $log = Log_files->make_build_log($wormbase);

my $farm_loc = $ENV{'PIPELINE'};
my $blastdir    = "$farm_loc/BlastDB";
my $acedir      = "$farm_loc/ace_files";
my $swalldir    = "$farm_loc/swall_data";

if ($default) {
  $human=$fly=$yeast=$swissprot=$cleanup=1;
} elsif ($all) {
  $human=$fly=$yeast=$uniprot=$cleanup=1;
}

if( $human ) { &process_human; } 

if($uniprot or $swissprot or $trembl) {
  #get current ver.
  my $cver = determine_last_vers('slimswissprot');
  
  #find latest ver
  my $reldate_file = "/ebi/ftp/pub/databases/uniprot/knowledgebase/reldate.txt";
  if (-e $reldate_file) { 
    open (WG, "< $reldate_file") or $log->log_and_die("cant get Uniprot page\n");
  } else {
    open (WG,"wget -O - -q ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/reldate.txt |") or $log->log_and_die("cant get Uniprot page\n");
  }
  my $lver;
  while(<WG>) { #to make processing easier we use the uniprot release no.rather than separate SwissProt and Trembl
    if (/UniProt\s+Knowledgebase\s+Release\s+(\d+)_(\d+)/){
      my $newver = sprintf("%d%d", $1, $2);
      if (!-e '/software/worm' && ($newver != $cver)) { # only run this is the EBI - Supported/uniprot is provided for us at the Sanger
        &process_uniprot($newver, $swissprot, $trembl);
      }
      if($newver != $cver and ($uniprot or $swissprot)) {
        &process_swissprot($newver);
      } else {
        $log->write_to("Not updating swissprot ($newver)");
      }
      if ($newver != $cver and ($uniprot or $trembl)) {
        &process_trembl($newver);
      } else { 
        $log->write_to("\tNot updating trembl ($newver)\n"); 
      }
    }
  }
  close WG;
}


if ($yeast) {
    $log->write_to("Updating yeast\n");

    my $update = 0;
    my $target = '/tmp/download.yeast.gz';
    my $source='http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz';
    
    eval {
      die "Could not fetch file" if $wormbase->run_command("wget -O $target $source",$log);
      die "Could not unzip file" if $wormbase->run_command("gunzip -f $target",$log);

      $target =~ s/\.gz//; #no longer gzipped
      
      my $ver = &determine_last_vers('yeast');
      #check if number of proteins has changed
      $log->write_to("\tcomparing\n");
      my $old_file = "$blastdir/yeast$ver.pep";
      my $old_cnt = qx{grep -c '>' $old_file};
      my $new_cnt = qx{grep -c '>' $target};
      
      die "Could not work out protein counts\n" unless( $old_cnt =~ /^\d+$/ and $new_cnt =~ /^\d+$/);
      
      $update = 1 if $old_cnt != $new_cnt;

    };
    if ($@) {
      $log->write_to("Could not successfully fetch yeast file ($@) so defaulting to previous version\n");
    } elsif (not $update) {
      $log->write_to("yeast is up to date\n");
    } else {
      $log->write_to("\tupdating yeast . . .\n");
      &process_yeast($target);
      $log->write_to("\tremoving download\n");
      $wormbase->run_command("rm -f $target", $log);
    }
}	

if ($fly) {
    $log->write_to("Updating fly . . \n");
    my $update = 0;
    my $fly_download = '/tmp/flybase.gz';
    my $ver;
    eval {
       # find the release version
	my $page_download = '/tmp/page_download';
	my $fly_version;
	$log->write_to("\tdownloading flybase listing\n");
	die "Could not get FlyBase listing" if $wormbase->run_command("wget -O $page_download ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/md5sum.txt", $log);
	open (PAGE, "<$page_download") || $log->log_and_die("Can't open $page_download\n");
	while (my $line = <PAGE>) {
	    if ($line =~ /dmel-all-translation-r(\d+)\.(\d+)\.fasta.gz/) {
		$fly_version = "$1.$2";
		last;
	    }
	}
	close(PAGE);
	$wormbase->run_command("rm -f $page_download", $log);

	#get the file 
	$log->write_to("\tdownloading flybase file\n");
	die "Could not fetch file" if $wormbase->run_command("wget -O $fly_download ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-translation-r${fly_version}.fasta.gz", $log);
	die "Could not unzip file" if $wormbase->run_command("gunzip -f $fly_download", $log);
	$fly_download = '/tmp/flybase';
    
	#check if number of proteins has changed
	$log->write_to("\tcomparing\n");
	$ver = &determine_last_vers('gadfly'); 
	my $old_file = "$blastdir/gadfly$ver.pep";
	my $old_cnt = qx{grep -c '>' $old_file};
	my $new_cnt = qx{grep -c '>' $fly_download};
	die "Could not work out protein counts\n" unless( $old_cnt =~ /^\d+$/ and $new_cnt =~ /^\d+$/);

	$update = 1 if $old_cnt != $new_cnt;
    };

    if ($@) {
	$log->write_to("Could not successfully fetch flybase file ($@) so defaulting to previous version\n");
    }
    elsif (!$update) {
	$log->write_to("flybase is already up to date\n");
	$wormbase->run_command("rm -f $fly_download", $log);
    }
    else {
	#update the file
	$log->write_to("\tupdating flybase . . .\n");
	$ver++;
	my $pepfile = "$blastdir/gadfly${ver}.pep";

	my $record_count = 0;	
	my $problem_count =0;
	my ($gadID, $FBname, $FBgn);
	my $count;
	my %genes;
	my $seqs = Bio::SeqIO->new('-file'  => $fly_download, 
                                   '-format'=>'fasta') or $log->log_and_die("cant open SeqIO from $fly_download:$!\n");
        SEQ:while (my $seq = $seqs->next_seq){
	  $count++;
	  $record_count++;
	  
	  my %fields = $seq->primary_seq->desc =~ /(\w+)=(\S+)/g;
          
	  if ($fields{name}){
            $FBname = $fields{name};
            $FBname =~ s/;//g;
	  }
	  ($FBgn) = $fields{'parent'} =~ /(FBgn\d+)/;
	  foreach ( split(/,/,$fields{'dbxref'}) ) {
            my($key, $value) = split(/:/);
            if( $key eq 'FlyBase_Annotation_IDs') {
              ($gadID) = $value =~ /(\w+)-/;
            }
	  }
	  
	  # some old style names still exist eg pp-CT*****.  In these cases
	  # we need to use the 1st field of the "from_gene" fields.
          
	  if( $gadID ){
            if($genes{$gadID}) {
              next SEQ if(length($genes{$gadID}->{'pep'}) > $seq->length);
            }
            $genes{$gadID}->{'fbname'} = $FBname if $FBname;
            $genes{$gadID}->{'fbgn'} = $FBgn if ($FBgn);
            $genes{$gadID}->{'pep'} = $seq->seq;
	  }
	  else {
            # report problems?
            $log->write_to("PROBLEM : $_\n");
	  }
        }
	
	my $acefile  = "$acedir/flybase.ace";
	open (ACE,">$acefile") or die "cant open $acefile\n"; 

        my $seqs_out = Bio::SeqIO->new('-file'   => ">$pepfile", 
                                       '-format' => 'fasta');
	
	foreach my $gadID (keys %genes){
	    #print ace file
	    my $FBname = $genes{$gadID}->{'fbname'};
	    my $FBgn   = $genes{$gadID}->{'fbgn'};
	    print ACE "\n\nProtein : \"FLYBASE:$gadID\"\n";
	    print ACE "Peptide \"FLYBASE:$gadID\"\n";
	    print ACE "Species \"Drosophila melanogaster\"\n";
	    print ACE "Gene_name \"$FBname\"\n" if $FBname;
	    print ACE "Database \"FlyBase\" FlyBase_gn \"$FBgn\"\n" if ($FBgn);
	    print ACE "Database \"FlyBase\" FlyBase_ID \"$gadID\"\n" if $gadID;
	    print ACE "Description \"Flybase gene name is $FBname\"\n" if $FBname;

	    print ACE "\nPeptide : \"FLYBASE:$gadID\"\n";
	    print ACE $genes{$gadID}->{'pep'}."\n";

            my $outseq = Bio::PrimarySeq->new( -id => $gadID, 
                                               -seq => $genes{$gadID}->{'pep'} );
            $seqs_out->write_seq($outseq);
	}
	
	#write database file
	close ACE;
	
	my $redundant = scalar keys %genes;
	$log->write_to("\t$record_count ($redundant) proteins\nflybase updated\n\n");
	$wormbase->run_command("rm -f $fly_download", $log);
    }
}


#
# And finally, clean up all the old database files
#
# If you have updated any of the blast databases ie not Interpro, which
# are elsewhere, in the steps above you must ensure that you remove the
# old versions from /lustre/scratch101/ensembl/wormpipe/BlastDB/. E.g. if
# gadfly3.2 is a new database, then you should end up with gadfly3.2.pep
# and remove gadfly3.1.pep.
#
# This will ensure that only the latest databases get copied across
# the ensembl compute farm. Also remove the old blast database index
# files for any old database (*.ahd, *.atb, *bsq, *.xpt, *.xps, *.xpd,
# *.psq, *.pin, *.phr).
#
if ($cleanup) {

  $log->write_to("  Removing old blast databases . . .\n");
  
  # root name regular expressions of the databases to check
  my @roots = (
	       'wormpep\d+.pep',
	       'brugpep\d+.pep',
	       'brepep\d+.pep',
	       'brigpep\d+.pep',
               'jappep\d+.pep',
               'ppapep\d+.pep',
	       'remapep\d+.pep',
	       'ovolpep\d+.pep',
               'gadfly\d+.pep',
               'ipi_human_\d+_\d+.pep',
	       'slimswissprot\d+.pep',
	       'slimtrembl\d+.pep',
	       'yeast\d+.pep',
	      );

  # get the list of files in the BLAST directory
  opendir FH, $blastdir;
  my @list = readdir(FH);
  closedir FH;

  foreach my $regex (@roots) {
    my @files = grep /$regex/, @list;
    if (scalar @files > 1) {
      # sort by creation time
      my @sort = sort {-M "$blastdir/$a" <=> -M "$blastdir/$b"} @files;
      my $youngest = shift @sort; # get the youngest file's release number
      my ($youngest_release) = ($youngest =~ /^[a-zA-Z_]+(\d+)/);
      #print "DONT DELETE release $youngest_release\n";
      #print "dont delete $youngest\n";
      foreach my $file (@sort) {
	if ($file =~ /^[a-zA-Z_]+${youngest_release}/) {
	  #print "dont delete $file\n";
	  next;
	}
	$log->write_to("    Deleting $file*\n");
	$wormbase->run_command("rm -f $blastdir/${file}*", $log);
	#print "rm -f $blastdir/${file}*\n";
      }
    }
  }
}



$log->mail;
exit(0);


##########################################################################################


sub process_human {
    use File::Listing qw(parse_dir);
    use POSIX qw(strftime);

    # determine last update done
    my @files  = glob($blastdir."/ipi_human*");
    my $file = shift(@files);
    my ($m,$d) = $file =~ /ipi_human_(\d+)_(\d+)/;

    my $remote_date;
    my $filename;

    my $ftp;
    my $ipi_file = "/ebi/ftp/pub/databases/IPI/last_release/current/ipi.HUMAN.fasta.gz";
    if (-e $ipi_file) {
      my @a = parse_dir(`ls -l $ipi_file`);
      my @x = localtime($a[0]->[3]);
      $remote_date = strftime("%m_%d",@x);
    } else {
      $ftp = Net::FTP->new("ftp.ebi.ac.uk", Debug => 0) or $log->log_and_die("Cannot connect to ftp.ebi.ac.uk: $@");
      $ftp->login('anonymous','wormbase@sanger.ac.uk');
      $ftp->cwd('pub/databases/IPI/last_release/current/');
      $filename = 'ipi.HUMAN.fasta.gz';
      my $ls = $ftp->dir("$filename");
      my @a = parse_dir($ls);
      $remote_date = strftime("%m_%d",localtime($a[0]->[3]));
    }

    if ($remote_date ne "${m}_${d}") {
	#download the file
	$log->write_to("\tupdating human to $remote_date\n");
	my $target = "/tmp/ipi_human_${remote_date}.gz";

	if (-e $ipi_file) {
	  $wormbase->run_command("cp $ipi_file $target", $log);
	} else {
	  $ftp->binary(); 
	  $ftp->get($filename,$target) or $log->log_and_die("failed getting $filename: ".$ftp->message."\n");
	  $ftp->quit;
	}

	$wormbase->run_command("gunzip $target",$log);
	$target =~ s/\.gz//; #no longer gzipped
	$wormbase->run_script("BLAST_scripts/parse_IPI_human.pl -file $target", $log);
	$wormbase->run_command("rm $target", $log);
	$log->write_to("\tIPI updated\n");
    } else {
	$log->write_to("\tIPI_human is up to date $remote_date\n");
    }
}

sub process_uniprot {
  my ($ver, $do_swiss, $do_trembl) = @_;

  my ($target1, $target2, $target3); 

  if ($do_swiss) {
      $target1 = $ENV{'PIPELINE'}."/blastdb/Supported/uniprot_sprot.fasta.gz";
      my $filename = "pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz";
      if (-e "/ebi/ftp/${filename}") {
	  $wormbase->run_command("cp /ebi/ftp/$filename $target1", $log);
      }
      else {
	  $log->log_and_die("Can't find file $filename\n") if $wormbase->run_command("wget -O $target1 ftp://ftp.ebi.ac.uk/$filename");
      }
  }
  
  if ($do_trembl) {
      $target2 = $ENV{'PIPELINE'}."/blastdb/Supported/uniprot_trembl.fasta.gz";
      my $filename = "pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz";
      if (-e "/ebi/ftp/${filename}") {
	  $wormbase->run_command("cp /ebi/ftp/$filename $target2", $log);
      }
      else {
	  $log->log_and_die("Can't find file $filename\n") if $wormbase->run_command("wget -O $target2 ftp://ftp.ebi.ac.uk/$filename");
      }
      $wormbase->run_command("cp $filename $target2",$log);
  }
  
  $target3 = $ENV{'PIPELINE'}."/blastdb/Supported/uniprot";
  $wormbase->run_command("rm -f $target3",$log);
  $wormbase->run_command("touch $target3", $log);

  if (defined $target1 and -e $target1) {
    $wormbase->run_command("gunzip -c $target1 >> ${target3}.pre",$log);
  }
  if (defined $target2 and -e $target2) {
    $wormbase->run_command("gunzip -c $target2 >> ${target3}.pre",$log);
  }

  # change Pyrrolysine (O) to Lysin (K) as WU BLAST can't deal with 'O' residues
  open (UNI, "<${target3}.pre") || $log->log_and_die("Can't open file ${target3}.pre\n");
  open (UNIOUT, ">${target3}") || $log->log_and_die("Can't open file ${target3}\n");
  while (<UNI>) {
    if (/^>\S+\|(\S+)\|(\S+)(.+?)SV=(\d+)$/) {
      print UNIOUT ">$1.$4 $2 $3\n";
    } else {
      s/O/K/; 
      print UNIOUT $_;
    }
  }
  close (UNIOUT);
  close (UNI);
  unlink "${target3}.pre";
}

sub process_swissprot {

  #swissprot
  my $ver = shift;

  my $target = $swalldir."/uniprot_sprot.dat.gz";

  my $filename = "/ebi/ftp/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz";

  if (-e $filename) {
    $wormbase->run_command("cp $filename $target",$log);
  } else {
    my $login = "anonymous";
    my $passw = 'wormbase@sanger.ac.uk';
    my $ftp = Net::FTP->new("ftp.ebi.ac.uk", Timeout => 18000);
    $ftp->login("anonymous",'wormbase@sanger.ac.uk');
    $ftp->cwd('pub/databases/uniprot/knowledgebase');
  
    my $filename = 'uniprot_sprot.dat.gz';
    $ftp->binary(); 
    $ftp->get($filename,$target) or $log->error("failed getting $filename: ".$ftp->message."\n");
    $ftp->quit;
  }

  $wormbase->run_script("BLAST_scripts/swiss_trembl2dbm.pl -s -file $target", $log);
  $wormbase->run_script("BLAST_scripts/swiss_trembl2slim.pl -s $ver",$log);
  $wormbase->run_script("BLAST_scripts/fasta2gsi.pl -directory $swalldir -f slimswissprot",$log);
  copy ("$swalldir/slimswissprot", "$blastdir/slimswissprot${ver}.pep");
  
}

sub process_trembl {

  my $ver = shift;

  my $tfile = 'uniprot_trembl.dat.gz';

  my $final_target = "$swalldir/$tfile";
  my $okay;
  my $local_target = $wormbase->scratch_area . "/$tfile";

  my $filename = "/ebi/ftp/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz";
  if (-e $filename) {
    $wormbase->run_command("cp $filename $final_target",$log);
    $okay = 1;
  } else {
    my $login = "anonymous";
    my $passw = 'wormbase@sanger.ac.uk';
    
    # fetch the file
    $okay = 1;

    my $ftp = Net::FTP->new("ftp.ebi.ac.uk", Timeout => 300);
    $ftp->login("anonymous",'wormbase@sanger.ac.uk');
    $ftp->cwd('pub/databases/uniprot/knowledgebase');
    $ftp->binary;
    $ftp->get($tfile,$local_target) or do {
      $log->error("process_trembl: Failed getting $tfile: ".$ftp->message."\n");
      $okay = 0;
    };
    $ftp->quit;

    if ($okay) {
      # move it to where it needs to be
      $wormbase->run_command("mv $local_target $final_target", $log) and do {
	$log->error("process_trembl: Could not mv $local_target to $final_target\n");
	$okay = 0;
      };
    }
  }


  if ($okay) {
    $wormbase->run_script("BLAST_scripts/swiss_trembl2dbm.pl -t -file $final_target", $log);
    $wormbase->run_command("rm -f $final_target", $log);
    $wormbase->run_script("BLAST_scripts/swiss_trembl2slim.pl -t $ver",$log);
    
    $wormbase->run_script("BLAST_scripts/blast_kill_list.pl -infile $swalldir/slimtrembl -outfile $blastdir/slimtrembl${ver}.pep -killfile $swalldir/kill_list.txt",$log);
    copy("$blastdir/slimtrembl${ver}.pep","$swalldir/slimtrembl_f");
    $wormbase->run_script("BLAST_scripts/fasta2gsi.pl -directory $swalldir -f slimtrembl_f",$log);
  }
}


sub process_yeast {
  my ($target) = @_;
  my $ver = &determine_last_vers('yeast');
  $ver++;
  my $source_file = "$blastdir/yeast${ver}.pep";
  my $acefile     = "$acedir/yeast.ace";
  # output initally goes to tmp file
  my $pepfile  = "$blastdir/yeast${ver}.pep.tmp"; 


# extract info from main FASTA file and write ace file
    open (SOURCE,"<$target") || $log->log_and_die("Couldn't open $target\n");
    open (PEP,">$pepfile") || $log->log_and_die("Couldn't open $pepfile\n");
    open (ACE,">$acefile") || $log->log_and_die("Couldn't open $acefile\n");

    while (<SOURCE>) {
	if( /\>/ ) { 
	    if (/\>(\S+)\s+(\S+)\s+SGDID:(\w+).+\"(.+)/) {
		my $ID = $1;
		my $GENE = $2;
		my $SGDID = $3;
		my $DESC = $4; 
		$DESC =~ s/\"$//;
		$DESC =~ s/"/'/g; # "

		print ACE "\nProtein : \"SGD:$ID\"\n";
		print ACE "Peptide \"SGD:$ID\"\n";
		print ACE "Species \"Saccharomyces cerevisiae\"\n";
		print ACE "Gene_name  \"$GENE\"\n";
		print ACE "Database \"SGD\" \"SGD_systematic_name\" \"$ID\"\n";
		print ACE "Database \"SGD\" \"SGDID\" \"$SGDID\"\n";
		print ACE "Description \"$DESC\"\n" if ($DESC);

		print ACE "\nPeptide : \"SGD:$ID\"\n"; 	
		
		print PEP ">$ID\n";
	    }
	    else {
		print $_;
	    }
	}
	else { 
	    print PEP $_;
	    print ACE $_;
	}
    }

    close(SOURCE);
    close(PEP);
    close(ACE);

# Now overwrite source file with newly formatted file
    system("mv $pepfile $source_file") && $log->log_and_die("Couldn't write peptide file $source_file\n");


}


sub determine_last_vers {
    my $db = shift;
    my @files  = glob($blastdir."/$db*");
    my $file = shift(@files);
    if (!defined $file || $file eq '') {$log->log_and_die("Can't find file $blastdir/$db*\n");}
    my ($ver) = $file =~ /$db(\d+)\.pep/;
    return $ver ? $ver : '666';
}

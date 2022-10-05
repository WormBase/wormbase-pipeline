#!/usr/bin/env perl
# create a local PFAM database from the FTP site and use it to extract the location of the active sites

use strict;
use Net::FTP;
use Getopt::Long;
use DBI;
use Digest::MD5 qw(md5_hex);

use Bio::SeqIO;

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $acefile, $species, $noload);

my $port= 3478;
my $server='ebiworm-db';
my $dbname = "worm_pfam";
my $tmpDir= $ENV {'WB_SCRATCH'};
my ($user, $pass, $update);
my $pfam_release = "current_release";

GetOptions ('help'       => \$help,
            'debug=s'    => \$debug,
            'test'       => \$test,
	    'verbose'    => \$verbose,
	    'store:s'    => \$store,
	    'species:s'  => \$species,
	    'user:s'     => \$user,
	    'pass:s'     => \$pass,
	    'update'     => \$update,
            'port=s'     => \$port,
            'host=s'     => \$server,
            'dbname=s'   => \$dbname, 
            'tmpdir=s'   => \$tmpDir,
            'acefile=s'  => \$acefile,
            "noload"     => \$noload,
            "pfamrelease=s" => \$pfam_release,
           );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism=> $species
			     );
}

$acefile = $wormbase->acefiles . "/PFAM_active_sites.ace" if not defined $acefile;
my $tax_id = $wormbase->ncbi_tax_id;

my $log = Log_files->make_build_log($wormbase);
$log->write_to("Getting PFAM active sites for $species\n");


$log->write_to("\tconnecting to $server:$dbname as $user\n");
my $DB = DBI->connect("DBI:mysql:$dbname:$server;port=$port;mysql_local_infile=1", 
                      $user,
                      $pass, 
                      {RaiseError => 1}) or
    $log->log_and_die("cannot connect to db, $DBI::errstr");

&update_database if $update;

my $sth_f = $DB->prepare (qq{	
	SELECT pfamseq.pfamseq_id, pfamseq_markup.residue, markup_key.label
	FROM pfamseq,pfamseq_markup,markup_key WHERE
	pfamseq.pfamseq_acc = pfamseq_markup.pfamseq_acc AND
	pfamseq_markup.auto_markup = markup_key.auto_markup AND
 	pfamseq.ncbi_taxid = $tax_id AND
        md5 = ?
} );

open (my $acefh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

my $pepfile = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep".$wormbase->get_wormbase_version;
my $seqs = Bio::SeqIO->new('-file' => $pepfile, '-format' => 'fasta');

my %seen;

while(my $pep = $seqs->next_seq){ 
  my $cds_id = $pep->id;
  my ($pepid) = $pep->desc =~ /^(\S+)/;
  $pepid=~s/wormpep=//; # the C.elegans protein header is slighlty different
  my $protein_id = $pepid;

  next if exists $seen{$protein_id};
  $seen{$protein_id} = 1;

  my ($seq) = $pep->seq;

  my $md5 = md5_hex($seq);

  $sth_f->execute($md5);
  my ($res) = @{$sth_f->fetchall_arrayref};

  next if not defined $res or not @$res;

  my ($seq_id, $residue, $label) = @$res;

  print $acefh "\n// $cds_id $pepid $seq_id $residue $label\n";

  my $method = "Pfam";
  
  if ($label =~ /^(\S+)\s+predicted/i) {
    $method = $1;
  }
  
  my $motif;
  if ($label =~/active site/i) {
    $motif = "Active_site";
  } elsif ($label =~ /metal ion/i) {
    $motif = "Metal_ion_binding_site";
  } else {
    next;
  }

  print $acefh "\nProtein : \"$protein_id\"\n";
  print $acefh "Motif_homol \"$motif\" \"$method\" 0 $residue $residue 1 1\n"; 
}
close($acefh) or $log->log_and_die("Could not close $acefile after writing\n");

if (not $noload) {
  $wormbase->load_to_database($wormbase->autoace, $acefile, 'PFAM_active_sites', $log);
}

$log->mail();
exit(0);


########################################################
sub update_database {
  $log->write_to("\n\nUpdating database from PFAM ftp site\n");
  
  
  my @tables = qw(pfamseq markup_key pfamseq_markup);
  foreach my $table (@tables){
    $log->write_to("\tfetching $table.txt\n");
    
    my $ftp = Net::FTP->new('ftp.ebi.ac.uk',Debug => 0)
        ||$log->log_and_die("Cannot connect to some.host.name: $@\n");
    $ftp->login("anonymous",'-anonymous@')
        ||$log->log_and_die("Cannot login ${\$ftp->message}\n");
    
    my $ftp_folder = "/pub/databases/Pfam/$pfam_release/database_files";

    $ftp->cwd($ftp_folder)
        ||$log->log_and_die("Cannot change working directory ${\$ftp->message}\n");
    $ftp->binary()||$log->log_and_die("cannot change mode to binary ${\$ftp->message}\n");
    
    # pfamseq table is subject to unannounced column re-ordering, so update the schema.
    if ($ftp->get("${table}.sql.gz","$tmpDir/${table}.sql.gz")){
      $log->write_to("\tupdating the $table table schema\n");
      $wormbase->run_command("echo \"SET FOREIGN_KEY_CHECKS=0;\"> $tmpDir/${table}.sql",$log);
      $wormbase->run_command("zcat $tmpDir/$table.sql.gz >> $tmpDir/${table}.sql", $log);
      $wormbase->run_command("echo \"SET FOREIGN_KEY_CHECKS=1;\">> $tmpDir/${table}.sql",$log);
      $wormbase->run_command("mysql -h $server -P$port -u$user -p$pass $dbname < $tmpDir/${table}.sql", $log);
      $wormbase->run_command("rm -f $tmpDir/${table}.sql.gz", $log);
      $wormbase->run_command("rm -f $tmpDir/${table}.sql", $log);
    } else {
      $log->write_to("\tcouldn't update the $table table schema\n");
    }
    
    if ($table eq 'pfamseq') {
      # treat this one specially; it is massive
      if ($ftp->get("${table}.txt.gz","$tmpDir/${table}.txt.gz")){
        $log->write_to("Parsing nematode proteins out of pfamseq...\n");
        $wormbase->run_command("zgrep Nematoda $tmpDir/${table}.txt.gz > $tmpDir/${table}.txt",$log);
      }else{$log->log_and_die("could not download ${table}.txt.gz\n")}
    } else {
      if ($ftp->get("${table}.txt.gz","$tmpDir/${table}.txt.gz")){
        $log->write_to("\tunzippping $tmpDir/$table.txt\n");
        $wormbase->run_command("gunzip -f $tmpDir/$table.txt.gz", $log);
      } elsif ($ftp->get("${table}.txt","$tmpDir/${table}.txt")){
        $log->write_to("\tgzip archive absent....using $table.txt.\n");
      } else {
        $log->log_and_die("Couldn't find $table file to download :(\n");
      }
    }
    
    $ftp->quit;
    
    # flush the table
    $log->write_to("\tclearing table $table\n");
    $DB->do("TRUNCATE TABLE $table") or $log->log_and_die($DB->errstr."\n");
    # load in the new data.
    $log->write_to("\tloading data in to $table\n");
    $DB->do("SET FOREIGN_KEY_CHECKS=0");		
    $DB->do("LOAD DATA LOCAL INFILE \"$tmpDir/$table.txt\" INTO TABLE $table".' FIELDS ENCLOSED BY \'\\\'\'') or $log->log_and_die($DB->errstr."\n");
    $DB->do("SET FOREIGN_KEY_CHECKS=1");
    
    # this will fall to pieces as soon as Rob changes the name of the column again
    if ($table eq 'pfamseq') {
      $log->write_to("\tcleaning quotation marks from $table\n");
      $DB->do("UPDATE pfamseq SET description=REPLACE(description,'\\'','')");
      $DB->do("UPDATE pfamseq SET description=REPLACE(description,'\"','')");
    }
    # clean up files
    
    foreach my $f ("$tmpDir/$table.txt", "$tmpDir/$table.txt.gz") {
      unlink $f if -e $f;
    }
  }
  $log->write_to("Database update complete\n\n");
}

1;

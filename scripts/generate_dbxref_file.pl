#!/usr/bin/env perl
#
# generate_dbxrefs_file.pl
#
# Generates a table of xrefs that should be generally useful for
# Uniprot, ENA and Ensembl. Columns are:
#
# 1  Gene sequence name
# 2  WBGene id
# 3  CGC name
# 4  Transcript name (not CDS name for coding transcripts, but the Coding_transcript name)
# 5  Pep name
# 6  ENA clone accession
# 7  ENA protein_id
# 8  Uniprot accession
#
# 
#  Last updated on: $Date: 2014-11-13 13:48:43 $
#  Last updated by: $Author: klh $

use strict;
use Getopt::Long;
use Storable;
use Net::FTP;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::WormSlurm;

my ($test,
    $debug,
    $store,
    $species,
    $database,
    $wormbase,
    $outfile,
    $no_header,
    $no_coding_transcripts,
    $ebi_upload,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "outfile:s"       => \$outfile,
  "noheader"        => \$no_header,
  "nocodingtrans"   => \$no_coding_transcripts,
  "ebiupload"       => \$ebi_upload,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
                             -autoace  => $database,
      );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
my $tace = $wormbase->tace;
my $full_species_name = $wormbase->full_name;
my $wormbase_version = $wormbase->get_wormbase_version_name;
my $dbdir = ($database) ? $database : $wormbase->autoace;
if (not defined $outfile) {
  if ($ebi_upload) {
    $outfile = $wormbase->reports . "/wormbase_xrefs." .$wormbase->get_wormbase_version_name . ".txt.gz";
  } else {
    $outfile = $wormbase->reports . "/${species}.dbxrefs.txt";
  }
}

my (%wbgene, %gene, %cds, %transcds, %clone2acc, %clone2chr, $out_fh);

$log->write_to("Generating sequence accession table\n");

my $query = &generate_accession_query();
my $command = "Table-maker -p $query\nquit\n";
open(my $tacefh,  "echo '$command' | $tace $dbdir |");
while(<$tacefh>) {
  chomp; s/\"//g;
  my ($seq_name, $accession, $subseq, $parent) = split(/\t/, $_);

  next if not defined $accession;

  $clone2acc{$seq_name} = $accession;
  
  if (defined $parent) {
    $clone2chr{$seq_name} = $parent;
  }
}
close($tacefh) or $log->log_and_die("Could not close tace TM query\n");
unlink($query);

# want to associate every sequence name with the ENA accesssion of its top-level parent
# i.e. for any feature on a clone, we want to give the accession of the chromosome
#  that the clone belongs to
foreach my $k (sort keys %clone2chr) {
  my $seq = $k;
  while(exists $clone2chr{$seq}) {
    if (exists $clone2acc{$clone2chr{$seq}}) {
      $clone2acc{$k} = $clone2acc{$clone2chr{$seq}};
    }
    $seq = $clone2chr{$seq};
  }
}

$log->write_to("Generating protein-coding table\n");

$query = &generate_coding_query($full_species_name);
$command = "Table-maker -p $query\nquit\n";

open ($tacefh, "echo '$command' | $tace $dbdir |");
while (<$tacefh>) {
  chomp; s/\"//g;

  my ($cds, $gene, $trans, $prot, $clone, $pid, $pid_version, $uniprot ) = split(/\t/, $_);

  next if $gene !~ /^WBGene/;

  $wbgene{$gene}->{cds}->{$cds} = 1;
  if ($no_coding_transcripts) {
    $wbgene{$gene}->{transcript}->{$cds} = 1;
    $transcds{$cds} = $cds;
  } else {
    $wbgene{$gene}->{transcript}->{$trans} = 1;
    $transcds{$trans} = $cds;
  }
  
  $prot =~ s/\S+://;
  $cds{$cds}->{protein} = $prot;
  if ($clone and $pid and $pid_version) {
    $cds{$cds}->{pid}->{"$clone:$pid:$pid_version"} = 1;
  }
  $cds{$cds}->{uniprot} = $uniprot if $uniprot;

}
close($tacefh) or $log->log_and_die("Could not cleanly close tace TM query\n");
unlink $query;

foreach my $class ('Transcript', 'Pseudogene') {
  $log->write_to("Generating non-coding $class table\n");

  $query = &generate_noncoding_query($full_species_name, $class);
  $command = "Table-maker -p $query\nquit\n";
  open (my $tacefh, "echo '$command' | $tace $dbdir |");
  while (<$tacefh>) {
    chomp; s/\"//g;
    
    my ($trans, $gene, $parent ) = split(/\t/, $_);
    next if $gene !~ /^WBGene/;
            
    $wbgene{$gene}->{transcript}->{$trans} = 1;
    $wbgene{$gene}->{sequence} = $parent;
  }
  close($tacefh) or $log->log_and_die("Could not cleanly close tace TM query\n");
  unlink $query;
}  


$log->write_to("Generating gene table\n");

$query = &generate_gene_query($full_species_name);
$command = "Table-maker -p $query\nquit\n";

open ($tacefh, "echo '$command' | $tace $dbdir |");
while (<$tacefh>) {
  chomp; s/\"//g;
  
  my ($wbgene, $sequence_name, $cgc_name, $locus_tag) = split(/\t/, $_);
  next if $wbgene !~ /^WBGene/;

  $gene{$sequence_name}->{$wbgene} = 1;
  if ($cgc_name) {
    $wbgene{$wbgene}->{cgc} = $cgc_name;
  }
  if ($locus_tag) {
    $wbgene{$wbgene}->{locus_tag} = $locus_tag;
  }
}
close($tacefh) or$log->log_and_die("Could not cleanly close tace TM query\n");
unlink $query;

if ($outfile =~ /\.gz$/) {
  open($out_fh, "| gzip -c > $outfile") or $log->log_and_die("Could not open $outfile for writing\n");
} else {
  open($out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");
}

&write_header($out_fh) unless $no_header;

foreach my $g (sort keys %gene) {
  foreach my $wbgeneid (keys %{$gene{$g}}) {
    my $cgc_name = (exists $wbgene{$wbgeneid}->{cgc}) ? $wbgene{$wbgeneid}->{cgc} : ".";
    my $locus_tag = (exists $wbgene{$wbgeneid}->{locus_tag}) ? $wbgene{$wbgeneid}->{locus_tag} : ".";
    foreach my $trans (keys %{$wbgene{$wbgeneid}->{transcript}}) {
      my @pid_list;
      my ($cds, $pepid, $uniprot) = (".", ".", ".");
      
      if (exists $transcds{$trans}) {
        # coding
        $cds = $transcds{$trans};        
        $pepid = $cds{$cds}->{protein};
        $uniprot = $cds{$cds}->{uniprot}  if exists $cds{$cds}->{uniprot};
        
        if (exists $cds{$cds}->{pid}) {
          foreach my $str (keys %{$cds{$cds}->{pid}}) {
            my ($clone, $pid) = split(/:/, $str);
            push @pid_list, [$clone2acc{$clone}, $pid];
          }
        } else {
          @pid_list = (['.', '.']);
        }
      } else {
        my ($pid, $clone) = (".", ".");

        if (exists $wbgene{$wbgeneid}->{sequence}) {
          $clone = $wbgene{$wbgeneid}->{sequence};
          if (exists $clone2acc{$clone}) {
            $clone = $clone2acc{$clone};
          } else {
            $clone = "$clone:NOACC";
          }
        }
        @pid_list =  ([$clone, $pid]);
      }

      foreach my $pidpair (@pid_list) {
        my ($clone, $pid) = @$pidpair;
        
        print $out_fh join("\t", 
                           $g, 
                           $wbgeneid, 
                           $cgc_name, 
                           $trans,
                           $pepid,
                           $clone,
                           $locus_tag,
                           $pid, 
                           $uniprot), "\n";
      }
    }
  }
}

close($out_fh) or $log->log_and_die("Could not cleanly close output file\n");

if ($ebi_upload) {
  &upload_to_ebi();
}

$log->mail();
exit(0);


#####################################################
sub upload_to_ebi {

  #
  # First, upload to EBI ENA xref drop-box
  #
  my ($ftp_host, $ftp_user, $ftp_pass, $ftp_dir);

  my $login_details_file = $wormbase->wormpub . "/ebi_resources/ENAXREFFTP.s";
  open(my $infh, $login_details_file)
      or $log->log_and_die("Can't open secure account details file $login_details_file\n");
  while (<$infh>){
    /^HOST:(\S+)$/ and $ftp_host = $1;
    /^USER:(\S+)$/ and $ftp_user = $1;
    /^PASS:(\S+)$/ and $ftp_pass = $1;
    /^DIR:(\S+)$/  and $ftp_dir  = $1;
  }
  close($infh);

  # Establish ftp connection 
  my $ftp = Net::FTP->new($ftp_host, Debug => 0) 
      or $log->log_and_die("Cannot connect to $ftp_host: $@");
  $ftp->login($ftp_user,"$ftp_pass")
      or $log->log_and_die ("Cannot login to $ftp_host using WormBase credentials\n". $ftp->message);
  $ftp->cwd($ftp_dir) 
      or $log->log_and_die ("Cannot change into to_ena dir for upload of files\n". $ftp->message);
  
  $ftp->binary();
  
  # If test dont upload any file to ENA
  unless ($test) {
      $ftp->put($outfile)
      or $log->log_and_die ("FTP-put failed for $outfile: ".$ftp->message."\n");
   }   
  
  $ftp->quit;

  $log->write_to("Successfully uploaded file $outfile to ENA Xref FTP drop-box unless you ran in test mode\n");

  # Secondly, copy the file to a reserved area on the FTP site - for UniProt to pick up xrefs
  # 
  my $uni_file = join(".", "wormbase_xrefs", $wormbase->get_wormbase_version_name, $wormbase->ncbi_tax_id, "txt", "gz");
  my $uni_symlink = join(".", "wormbase_xrefs", "latest", $wormbase->ncbi_tax_id, "txt", "gz");

  my $dest_dir = $ENV{'XREF_FTP_DROPBOX'};
  my $job_id = WormSlurm::submit_job_and_wait("cp $outfile $dest_dir/$uni_file", 'datamover', '200m', '00:15:00', '/dev/null', '/dev/null'); 
  my $exit_code = WormSlurm::get_exit_code($job_id);
  if ($exit_code) {
      $log->error("Copying xref file to $dest_dir failed\n");
  } else {
      $log->write_to("Copied xref file to $dest_dir\n");
  }
  $job_id = WormSlurm::submit_job_and_wait("cd $dest_dir && ln -sf $uni_file $uni_symlink", 'datamover', '200m', '00:01:00', '/dev/null', '/dev/null');
  $exit_code = WormSlurm::get_exit_code($job_id);
  if ($exit_code) {
      $log->error("Creating symlink to xref file in $dest_dir failed\n");
  } else {
      $log->write_to("Created symlink to  xref file in $dest_dir\n");
  }
}


##########################################
sub write_header {
  my ($out_fh) = @_;

  my $header = <<"HERE";
//
// WormBase $full_species_name XREFs for $wormbase_version
//
// Columns (tab separated) are:
//    1. WormBase Gene sequence name
//    2. WormBase Gene accession
//    3. WormBase Gene CGC name
//    4. WormBase Transcript sequence name
//    5. WormPep protein accession
//    6. INSDC parent sequence accession
//    7. INSDC locus_tag id
//    8. INSDC protein_id
//    9. UniProt accession
//
// Missing or not applicable data (e.g. protein identifiers for non-coding RNAs) is denoted by a "."
//
HERE

  print $out_fh $header;

}


##########################################
sub generate_accession_query {

  my $tmdef = "/tmp/clone2acc_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Condition Genomic_canonical OR Source OR Subsequence
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 1 
Tag Database  
Condition EMBL OR NDB
 
Colonne 3 
Width 12 
Mandatory 
Hidden 
Class 
Class Database_field 
Right_of 2 
Tag HERE   
Condition NDB_AC
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 3 
Tag HERE   
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Subsequence  
 
Colonne 6 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Source  
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;

}


##########################################
sub generate_coding_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/cod_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $tablemaker_template = <<"EOF";


Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class CDS 
From 1 
Condition Method = "curated" AND Species = "$full_species"
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Transcript 
From 1 
Tag Corresponding_transcript 
 
Colonne 4
Width 12 
Optional 
Visible 
Class 
Class Protein
From 1 
Tag Corresponding_protein

Colonne 5
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Protein_id 
 
Colonne 6 
Width 12 
Optional 
Visible 
Text 
Right_of 5 
Tag  HERE  
 
Colonne 7 
Width 12 
Optional 
Visible 
Integer 
Right_of 6 
Tag  HERE  
 
Colonne 8 
Width 12 
Optional 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition UniProt
 
Colonne 9 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 8 
Tag  HERE  
Condition UniProtAcc
 
Colonne 10 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 9 
Tag  HERE  

EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}


############################################
sub generate_noncoding_query {
  my ($full_species, $class) = @_;

  my $tmdef = "/tmp/nc_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";
  if ($class eq 'Transcript') {
    $condition = "NOT Method = \"Coding_transcript\" AND NOT Method = \"history_transcript\" AND Species = \"$full_species\"";
  } elsif ($class eq 'Pseudogene') {
    $condition = "(Method = \"Pseudogene\" OR Method = \"rRNA_pseudogene\" OR Method = \"tRNA_pseudogene\") AND Species = \"$full_species\"";
  } else {
    $log->log_and_die("Unrecognised non-coding class: $class\n");
  }

  my $tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class $class
From 1 
Condition $condition 
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Sequence 
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;

}


#######################################
sub generate_gene_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/gene_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";

  my $tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Condition Species = "$full_species"
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Gene_name 
From 1 
Tag Sequence_name 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag CGC_name 

Colonne 4 
Width 12 
Optional 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition NDB
 
Colonne 5 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 4 
Tag  HERE  
Condition locus_tag
 
Colonne 6 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 5 
Tag  HERE  
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}
__END__

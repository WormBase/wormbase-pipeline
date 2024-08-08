#!/usr/bin/env perl

package WBSRS;

use strict;
use LWP::UserAgent;



###########################
# return hash-ref of the following information for each Accession in Uniprot
# key is 'acc' value is hash of 'id', 'des', 'gene', 'key', 'len'
# converted this to use the ENA RESTful interface
# See: http://www.uniprot.org/faq/28
# e.g. https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Ckeyword%2Clength&format=tsv&query=%2A%20AND%20%28model_organism%3AHomo%20sapiens%29

sub uniprot {
  my ($species) = @_;

  if (! defined $species) {die("SRS: Full species name not given\n")};
  $species =~ s/\s+/\%20/g;

  my %acc2ids;

  my $ua       = LWP::UserAgent->new;
  
  my $base     = 'https://rest.uniprot.org/uniprotkb/stream?';
  my $query = "fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Ckeyword%2Clength&format=tsv&query=%2A%20AND%20%28model_organism%3A${species}%29";
  
  my $tmp_file = "/tmp/srs_results.$$.txt";
    
  my $qa2 = $ua->get($base.$query, ':content_file' => $tmp_file);
  die("WBSRS: Could not fetch EMBL entries using EBI Uniprot server") 
      if not $qa2->is_success;
  
  open(my $f, $tmp_file);
  my $count=0;
  while(<$f>) {
    if (++$count == 1) {next;} # skip the title line
    
    my @f = split /\t/;
    next unless @f;
    $acc2ids{$f[0]}{'id'} = $f[1];
    $acc2ids{$f[0]}{'des'} = $f[2];
    my @gene_names = split /\s/, $f[3];
    $acc2ids{$f[0]}{'gene'} = $gene_names[0];
    $acc2ids{$f[0]}{'key'} = $f[5];
    $acc2ids{$f[0]}{'len'} = $f[6];
  }
  unlink $tmp_file;

  return \%acc2ids;
}

sub old_uniprot {
  my ($species) = @_;

  my ($ggenus, $gspecies) = $species =~ /^(\S+)\s+(\S+)/;
  if (! defined $gspecies) {die("WBSRS: Full species name not given\n")};

  my %acc2ids;

  my $ua       = LWP::UserAgent->new;

  my $base     = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession';
  my $query    = "+[uniprot-org:$ggenus]&[uniprot-org:$gspecies]";
  
  my $cResult  = '+-page+cResult+-ascii';
  my $fullview = '+-view+UniprotView+-ascii';
  
  # Doing query: ${base}${query}${cResult}
  # getting the number of entries in the set of results
  
  my $qa1 = $ua->get($base.$query.$cResult);
  die("SRS: Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;

  if($qa1->content =~/^(\d+)/) {
    my $total = $1;
    # print "EBI SRS server Uniprot query returned $total entries; fetching...\n";
    
    my $tmp_file = "/tmp/srs_results.$$.txt";
    
    my $chunksize = 10000;
    my $chunkstart = 1; # start of next chunk
    while ($chunkstart < $total) {
      my $bv = "+-bv+${chunkstart}";
      my $lv = "+-lv+${chunksize}";
      # print "EBI SRS server Uniprot query returned $total entries; fetching no. $chunkstart...\n";
      my $qa2 = $ua->get($base.$query.$fullview.$lv.$bv, ':content_file' => $tmp_file);
      die("WBSRS: Could not fetch EMBL entries using EBI SRS server") 
        if not $qa2->is_success;
      
      open(my $f, $tmp_file);
      while(<$f>) {
	my @f = split /\t/;
	$f[0] =~ s/UNIPROT:(\S+)/$1/;
	$acc2ids{$f[1]}{'id'} = $f[0];
	$acc2ids{$f[1]}{'des'} = $f[3];
	$acc2ids{$f[1]}{'gene'} = $f[4];
	$acc2ids{$f[1]}{'key'} = $f[6];
	$acc2ids{$f[1]}{'len'} = $f[7];
      }
      $chunkstart+=$chunksize;
    }
  } else {
    die("WBSRS: Unexpected content from SRS query\n");
  }
  

  return \%acc2ids;
}


###########################
# return hash-ref of the following information for each Accession in EMBL
# this will probably often fail if done in one go, so the returned result has been chunked
# key is 'acc' value is hash of 'id', 'des', 'len'
sub embl {
  my ($species) = @_;

  my ($ggenus, $gspecies) = $species =~ /^(\S+)\s+(\S+)/;
  if (! defined $gspecies) {die("WBSRS: Full species name not given\n")};

  my %acc2ids;


  my $ua       = LWP::UserAgent->new;

  my $base     = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession';
  my $query    = "+[embl-org:$ggenus]&[embl-org:$gspecies]";

  my $cResult  = '+-page+cResult+-ascii';
  my $fullview = '+-view+EMBLSeqSimpleView+-ascii';

  # Doing query: ${base}${query}${cResult}
  # getting the number of entries in the set of results

  my $qa1 = $ua->get($base.$query.$cResult);
  die("WBSRS: Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;

  if($qa1->content =~/^(\d+)/) {
    my $total = $1;
    # print "EBI SRS server Uniprot query returned $total entries; fetching...\n";
    
    my $tmp_file = "/tmp/srs_results.$$.txt";
    
    my $chunksize = 10000;
    my $chunkstart = 1; # start of next chunk
    while ($chunkstart < $total) {
      my $bv = "+-bv+${chunkstart}";
      my $lv = "+-lv+${chunksize}";
      # print "EBI SRS server Uniprot query returned $total entries; fetching no. $chunkstart...\n";
      my $qa2 = $ua->get($base.$query.$fullview.$lv.$bv, ':content_file' => $tmp_file);
      die("WBSRS: Could not fetch EMBL entries using EBI SRS server") 
        if not $qa2->is_success;
      
      open(my $f, $tmp_file);
      while(<$f>) {
	my @f = split /\t/;
	$f[0] =~ s/EMBL:(\S+)/$1/;
	$acc2ids{$f[1]}{'id'} = $f[0];
	$acc2ids{$f[1]}{'des'} = $f[3];
	$acc2ids{$f[1]}{'len'} = $f[4];
      }
      $chunkstart+=$chunksize;
    }
  } else {
    die("WBSRS: Unexpected content from SRS query\n");
  }
  
  return \%acc2ids;
}

# e.g.
#http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-id+PERMtestSession+[embl-id:AK033621|AK078912]+-f+id%20acc%20div%20mol%20des%20key%20org+-ascii+-view+EMBLSeqSimpleView

###########################
# return text of the specified entry
# http://srs.ebi.ac.uk/cgi-bin/wgetz?-noSession+-e+[EMBL-id:V01512]+-ascii
# example databases: EMBL, SWISSPROT, UNIPROT
# example types of query: 'id', 'acc', (more general types like 'des, 'all' will probably crash your script)
sub entry {

  my ($database, $type, $entry) = @_;

  my $text;
  my $ua       = LWP::UserAgent->new;

  my $base     = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession';
  my $query    = "+-e+[${database}-${type}:${entry}]";
  my $fullview = '+-ascii';

  # Doing query: ${base}${query}${fullview}

  my $qa1 = $ua->get($base.$query.$fullview);
  die("WBSRS: Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;
  $text = $qa1->content;

  return $text;
}
###########################
# return fasta sequences
# http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession+[uniprot-id:pax6_*]+-view+FastaSeqs
# example databases: EMBL, SWALL (same as UNIPROT), UNIPROT, REFSEQ
# example types of query: 'id', 'acc'

sub fasta {

  my ($database, $type, $entry) = @_;

  my $text;
  my $ua       = LWP::UserAgent->new;

  my $base     = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession';
  my $query    = "+[${database}-${type}:${entry}]";
  my $fullview = '+-ascii+-view+FastaSeqs';

  # Doing query: ${base}${query}${fullview}

  my $qa1 = $ua->get($base.$query.$fullview);
  die("WBSRS: Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;
  $text = $qa1->content;

  return $text;

}

1;

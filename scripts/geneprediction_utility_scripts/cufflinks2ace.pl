#!/software/bin/perl -w
#
# Small script to convert GF3 gene predictions to ace
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008/02/14 11:02:17 $      

use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

my ($method, $species);

GetOptions ( "method:s"   => \$method, # method to specify in ace output
             "species:s"  => \$species,
	    )||die(@!);
$species ||= 'Brugia malayi';
$method  ||='cufflinks';

my %clone_child;
my %gene_info;
my %sources;			# keep a note of the sources used

# suck the data in
while (my $line = <>) {
  
  if ($line =~ /^#/) {next;}
  my @fields = split /\s+/, $line;

  # ignore lines that are not from the source we want
  $sources{$fields[1]} = 1;	# keep a note of the sources

  my ($start, $end);

  if($fields[2] eq 'transcript'){
    my $clone = $fields[0];
    $start = $fields[3];
    $end   = $fields[4];
    my $sense = $fields[6]; 
    my ($id) = $line =~/transcript_id \"([^\s;]+)\"/;
    $gene_info{$id}{'sense'}= $sense;
    $gene_info{$id}{'chrom_start'} = $sense eq '+'?$start:$end;
    $gene_info{$id}{'clone'} = $clone;

    $clone_child{$clone}{$id}{'start'}=$start;
    $clone_child{$clone}{$id}{'end'}=$end;
    $clone_child{$clone}{$id}{'sense'}=$sense;
     
  } elsif($fields[2] eq 'exon') {
    my ($id) = $line =~/transcript_id \"([^\s;]+)\"/;

    # get start/end of exon relative to the start of the gene
    my $chrom_start = $gene_info{$id}{'chrom_start'};
    if (! defined $chrom_start) {die "The chrom_start for gene $id is not defined\n";}
    if ($gene_info{$id}{'sense'} eq '+') {
      $start = $fields[3] - $chrom_start + 1;
      $end = $fields[4] - $chrom_start + 1;

      # store the exon
      push @{$gene_info{$id}{'starts'}}, $start;
      push @{$gene_info{$id}{'ends'}}, $end;
    } else { # reverse sense
      $start = $chrom_start - $fields[3] + 1;
      $end = $chrom_start - $fields[4] + 1;

      # store the exon
      push @{$gene_info{$id}{'starts'}}, $end;
      push @{$gene_info{$id}{'ends'}}, $start;
      }
   }
}


# write out the exons of the gene
foreach my $id (keys %gene_info) {
  if ($gene_info{$id}{'sense'} eq '-') {			# reverse sense
    @{$gene_info{$id}{'starts'}} = reverse @{$gene_info{$id}{'starts'}};
    @{$gene_info{$id}{'ends'}} = reverse @{$gene_info{$id}{'ends'}};
  }
  my $clone = $gene_info{$id}{'clone'};
  print "\nTranscript : \"$id\"\n";
  print "Sequence \"$clone\"\n";
  print "Species \"$species\"\n";
  print "Method \"$method\"\n";
  for (my $i = 0; $i < @{$gene_info{$id}{'starts'}}; $i++) {
    my $start = $gene_info{$id}{'starts'}[$i];
    my $end = $gene_info{$id}{'ends'}[$i];
    print "Source_exons $start $end\n";
  }
}

# now write out the CDS_child data
foreach my $clone (keys %clone_child) {
  print "\nSequence : \"$clone\"\n";
  foreach my $id (keys %{$clone_child{$clone}}){
    my $start = $clone_child{$clone}{$id}{'start'};
    my $end = $clone_child{$clone}{$id}{'end'};
    my $sense = $clone_child{$clone}{$id}{'sense'};
    if ($sense eq '+') {
      print "Transcript \"$id\" $start $end\n";
    } else {
      print "Transcript \"$id\" $end $start\n";
    }
  }
}

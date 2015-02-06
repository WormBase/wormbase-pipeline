#!/usr/bin/env perl

# 

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Ace;
use IO::String;
use IO::File;
use Bio::AlignIO;
#use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::Muscle;
use Getopt::Long;
use DBI;

#
# color table based on malign, but changed for the colour blind
#
my %COLOURS = ( 
'*'       =>  '666666',         #mismatch          (dark grey)
'.'       =>  '999999',         #unknown           (light grey)
 A        =>  '33cc00',         #hydrophobic       (bright green)
 B        =>  '666666',         #D or N            (dark grey)
 C        =>  '2c5197',         #cysteine          (st.louis blue)
 D        =>  '0033ff',         #negative charge   (bright blue)
 E        =>  '0033ff',         #negative charge   (bright blue)
 F        =>  '009900',         #large hydrophobic (dark green)
 G        =>  '33cc00',         #hydrophobic       (bright green)
 H        =>  '009900',         #large hydrophobic (dark green)
 I        =>  '33cc00',         #hydrophobic       (bright green)
 K        =>  'cc0000',         #positive charge   (bright red)
 L        =>  '33cc00',         #hydrophobic       (bright green)
 M        =>  '380474',         #hydrophobic       (blue deep)
 N        =>  '6600cc',         #polar             (purple)
 P        =>  '33cc00',         #hydrophobic       (bright green)
 Q        =>  '6600cc',         #polar             (purple)
 R        =>  'cc0000',         #positive charge   (bright red)
 S        =>  '0099ff',         #small alcohol     (dull blue)
 T        =>  '0099ff',         #small alcohol     (dull blue)
 V        =>  '33cc00',         #hydrophobic       (bright green)
 W        =>  '009900',         #large hydrophobic (dark green)
 X        =>  '666666',         #any               (dark grey)
 Y        =>  '009900',         #large hydrophobic (dark green)
 Z        =>  '666666',         #E or Q            (dark grey)
);


my ($debug,
    $store,
    $pepfile,
    $database,
    $test,
    $batch_id,
    $batch_total,
    $wb,
    $aligner_dir,
    $rdb_name,
    $rdb_user,
    $rdb_pass,
    $rdb_host,
    $rdb_port);

GetOptions(
  'debug=s'  => \$debug,
  'store=s'  => \$store,
  'pepfile=s' => \$pepfile,
  'database=s'  => \$database,
  'test'        => \$test,
  'batchid=s'    => \$batch_id,
  'batchtotal=s' => \$batch_total,
  'host=s'      => \$rdb_host,
  'port=s'      => \$rdb_port,
  'user=s'      => \$rdb_user,
  'dbname=s'    => \$rdb_name,
  'password=s'  => \$rdb_pass,
  'alignerdir=s' => \$aligner_dir,
    ) or die "Incorrect invocation\n";


# logfile / wormbase setup
if ($store) {
    $wb = Storable::retrieve($store) or die("cannot restore wormbase from $store");
}
else { 
  $wb = Wormbase->new( -debug => $debug, 
                       -test => $test );
}
my $log = Log_files->make_build_log($wb);
$ENV{MUSCLEDIR} = $aligner_dir;

$database = $wb->autoace if not defined $database;
# database setup
my $db = Ace->connect( -path => $database ) 
    or do { print "cannot connect to $database:", Ace->error; die };

my $pgdb = DBI->connect("DBI:mysql:dbname=${rdb_name};host=${rdb_host};port=${rdb_port}" ,$rdb_user,$rdb_pass);

my $sth = $pgdb->prepare('INSERT INTO clustal(peptide_id,alignment) VALUES (?,?)');
my $qst = $pgdb->prepare('SELECT COUNT(peptide_id) FROM clustal WHERE peptide_id = ?');

# file setup
my $infile= new IO::File $pepfile || $log->write_to("cannot open $pepfile\n");

### here goes the output bit ####

my $entry_no=0;
my %seen;
while (my $line = <$infile>){
    next unless $line=~/^\>/;

    $entry_no++;
    
    my ($cds,$wpid,$gene)=split /\s+/,$line;
    $seen{$wpid}++;

    #
    # For batch_total X:
    #   batch_id 1   => entry 1, X+1, 2X+1 etc 
    #   batch_id 2   => entry 2, X+2, 2X+2 etc
    #   batch_id X-1 => entry X-1, X + (X-1), 2x + (X-1) etc
    #   batch id X   => entry X, X + X, 2x + X etc  

    my $offset = $entry_no % $batch_total;
    $offset = $batch_total if not $offset;
    
    next unless $offset == $batch_id;

    # the same wormpep id may be seen multiple times in the file
    # so we will only process the first one
    if ($seen{$wpid} > 1) {
      $log->write_to("Skipping $cds because $wpid has already been referenced earlier in the file\n");
      next;
    }

    my ($protein)=$db->fetch(Protein => $wb->wormpep_prefix .':'. $wpid);
    $log->write_to("can't find $wpid\n") unless $protein;
    
    next unless $protein;
    next unless length $protein->asPeptide() >=1;
    
    $qst->execute("$protein");
    my ($found)=$qst->fetchrow_array;
    # if in test mode, we are not writing to the db, so go ahead and calculate anyway
    next if $found>0 and not $test;
    
    print "processing $protein\n" if $debug;
    
    my $cl_aln=print_alignment($protein);
    if ($test) {
      print $cl_aln, "\n";
    } else {
      $sth->execute("$protein",$cl_aln);
    }
}
$db->close();

$log->mail();
exit(0);

# main function to generate the alignments
sub print_alignment {
  my @sequences;
  my $protRecord = shift;
  my $peptide = $protRecord->asPeptide();
  $peptide =~ s/^.+?[\s\n]//;
  $peptide =~ s/\s//g;
  
  print "processing $protRecord\n" if $debug;
 
  push( @sequences, Bio::Seq->new(-id => $protRecord, -seq => $peptide));

  my (@candObjs) = Best_BLAST_Objects($protRecord);
  
  my $header; # will get the linkout table
  
  foreach (@candObjs) {
    next if "$_"=~/^MSP/;
    my $pept = $_->asPeptide();
    $pept =~ s/^.+?[\s\n]//;
    $pept =~ s/\s+//g;    # remove any space
    print "  adding $_\n" if $debug;
    next unless length($pept) > 3;
    
    if ($_->Corresponding_CDS){
      my @cds = $_->Corresponding_CDS;
      foreach my $cd (@cds){
        next unless $cd->Gene;
        $header.=sprintf("<tr><td>%s</td><td>(%s)<td>%s <a href=\"/db/gene/gene?name=%s\">%s</a> %s</td></tr>",
                         &protein_url($_),$_->Species,$cd,$cd->Gene,$cd->Gene,($cd->Gene->CGC_name||''));
      }
    }else{
      $header.=sprintf("<tr><td>%s</td><td>(%s)</td><td>%s</td></tr>",&protein_url($_),$_->Species,$_->Description);
    }
    
    push(@sequences, Bio::Seq->new(-id => $_ , -seq => $pept));
  }
  
  if (scalar @sequences <2){
    $log->write_to("$protRecord has no homologs\n");
    return undef;
  }
  # at this point, sequences and alignment should be ready
  $log->write_to("No sequences for $protRecord\n") if (!@sequences);
  
  ######## HaCK from ClustalW bioperltut ################
  use Bio::Tools::Run::Alignment::Muscle;
  my @parameters = (-quiet => 1);
  my $factory = Bio::Tools::Run::Alignment::Muscle->new(@parameters);
  # $factory->save_tempfiles(1);
  my $aln = $factory->align(\@sequences);

 ##############################################################
  my $alignString;
  my $alignIO = IO::String->new($alignString);

  my $aligner = Bio::AlignIO->newFh(
      -interleaved => 0,
      -fh => $alignIO,
      -format => 'clustalw',
      -idlength => 15,
  );
  print $aligner $aln;
  
  #my $calign = _postprocess($alignString);
  $alignString=~s/CLUSTAL.*multiple sequence alignment/MUSCLE (v3.8.31) multiple sequence alignment/;
  return "<table border='0'>$header</table><pre>$alignString</pre>";
}

# from the website
# grab the ids of the best blast hits
sub Best_BLAST_Objects {
  my $protein = shift;  # ace object
  my %matchCtr;
  
  my @pep_homol = $protein->Pep_homol;
  my $length    = $protein->Peptide(2);
  
  # find the best pep_homol in each category
  my %best;
  return () unless @pep_homol;
  
  for my $hit (@pep_homol) {
    next if "$hit" eq "$protein";
    my($method,$score) = $hit->row(1) or next;
    
    # Perl is interpreting integers as strings here (5.8.5, crap need to upgrade)
    my $prev_score     = (!$best{$method}) ? $score : $best{$method}{score};
    $prev_score        = ($prev_score =~ /\d+\.\d+/) ? $prev_score . '0': "$prev_score.0000";
    my $curr_score     = ($score =~ /\d+\.\d+/) ? $score . '0': "$score.0000";
    
    $best{$method}     = {score=>$score,hit=>$hit,adjusted_score=>$curr_score} if !$best{$method} || $prev_score < $curr_score;
  }
  
  my @bestIDs;  # Ace objects; best matches in different species
  
  foreach my $currVal (sort {$b->{adjusted_score}<=>$a->{adjusted_score}}(values %best)) {
    push(@bestIDs, $currVal->{hit});
  }
  return (@bestIDs);
} # end Best_BLAST_Match  

# create protein linkouts
sub protein_url {
  my $p=shift;
  return "<a href=\"http://www.wormbase.org/db/seq/protein?name=$p\">$p</a>";
}

# colour the raw alignment
sub _postprocess{
     my $raw_al=shift;
     my @line=split("\n",$raw_al);
     my $coloured;
     foreach my $l(@line) {
         my @cols=split(//,$l);
         my $flip=0;
         for(my $position=0;$position < scalar(@cols);$position++){
           next if $l=~/CLUSTAL/;
           $flip=1 if $cols[$position]=~/\s/;
           next unless $flip;
           $cols[$position]="<font color=\"#$COLOURS{$cols[$position]}\">$cols[$position]</font>"
            if $COLOURS{$cols[$position]};
         }
         $coloured.=join('',@cols);
         $coloured.="\n";
     }
     return $coloured; 
}

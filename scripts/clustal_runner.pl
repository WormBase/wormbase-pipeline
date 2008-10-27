#!/usr/bin/env perl

# get a inputfile, a database, offset and windowsize and create clustal alignments for them.

$ENV{CLUSTALDIR} = '/software/worm/bin';

use lib $ENV{CVS_DIR};
use lib '/software/worm/lib/site_perl';
use lib '/software/worm/ensembl-51/bioperl-live';
use lib '/software/worm/ensembl-51/bioperl-run';
use Wormbase;
use Ace;
use IO::String;
use IO::File;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Getopt::Long;
use DBI;

my ($debug,$store,$pepfile,$database,$test,$offset,$window,$species,$wb,$user,$pass);
GetOptions(
    'debug=s'  => \$debug,
    'store=s'  => \$store,
    'pepfile=s' => \$pepfile,
    'database=s'  => \$database,
    'test'        => \$test,
    'offset=s'    => \$offset,
    'window=s'    => \$window,
    'species=s'   => \$species,
    'user=s'      => \$user,
    'password=s'  => \$pass,
) or &print_usage();

sub print_usage {
print <<USAGE;
make_clustal.pl options:
            -debug USER_NAME    sets email address and debug mode
            -store FILE_NAME    use a Storable wormbase configuration file
            -pepfile PEP_FILE   input wormpep file
            -database           DATABASE_DIRECTORY use a different AceDB
            -offset NUMBER      chunk number
            -window NUMBER      chunk size
            -user NAME          pg_database user name
	    -password PASSWORD  pg_database password
USAGE

    exit 1;
}

# logfile / wormbase setup
if ($store) {
    $wb = Storable::retrieve($store) or $log->write_to("cannot restore wormbase from $store");
}
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species, -autoace => $database ) }
my $log = Log_files->make_build_log($wb);

# database setup
my $db = Ace->connect( -path => ($database||$wb->autoace) ) 
    || do { print "cannot connect to ${$wb->autoace}:", Ace->error; die };

my $pgdb = DBI->connect('dbi:Pg:dbname=clx;host=deskpro16391.dynamic.sanger.ac.uk',$user,$pass);

my $sth = $pgdb->prepare('INSERT INTO clustal(peptide_id,alignment) VALUES (?,?)');
my $qst = $pgdb->prepare('SELECT COUNT(peptide_id) FROM clustal WHERE peptide_id = ?');

# file setup
my $infile= new IO::File $pepfile || $log->write_to("cannot open $pepfile\n");

### here goes the output bit ####

my $line_no=0;
while (my $line = <$infile>){
    next unless $line=~/^\>/;
    # window magick
    $line_no++;
    next unless $line_no >= $window*($offset-1);
    last if $line_no > $offset*$window;
    
    my ($cds,$wpid,$gene)=split /\s+/,$line;
    my ($protein)=$db->fetch(Protein => $wb->wormpep_prefix .':'. $wpid);
    $log->write_to("can't find $wpid\n") unless $protein;
    
    next unless $protein;
    next unless length $protein->asPeptide() >=1;
    
    $qst->execute("$protein");
    my ($found)=$qst->fetchrow_array;
    next if $found>0;
    
    #print "processing $protein\n";
    
    my $cl_aln=print_alignment($protein);
    $sth->execute("$protein",$cl_aln);
}

$log->mail('All',"LOG: from $0");

# main function to generate the alignments
sub print_alignment{
  my @sequences;
  my $protRecord = shift;
  my $peptide = $protRecord->asPeptide();
  $peptide =~ s/^.+?[\s\n]//;
  $peptide =~ s/\s//g;
  
  #print "processing $protRecord\n";
 
  push( @sequences, Bio::Seq->new(-id => $protRecord, -seq => $peptide));

  my ($candObjs,) = Best_BLAST_Objects($protRecord);
  
  foreach (@$candObjs) {
    my $pept = $_->asPeptide();
    $pept =~ s/^.+?[\s\n]//;
    $pept =~ s/\s+//g;    # remove any space
    #print "  adding $_\n";
    next unless length($pept) > 3;
    push(@sequences, Bio::Seq->new(-id => $_ , -seq => $pept));
   }

  if (scalar @sequences <2){
      $log->write_to("$protRecord has no homologs\n");
      return undef;
  }
  # at this point, sequences and alignment should be ready
  $log->write_to("No sequences for $protRecord\n") if (!@sequences);

  ######## HaCK from ClustalW bioperltut ################
  use Bio::Tools::Run::Alignment::Clustalw;
  my @parameters = (-quiet => 1);
  my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@parameters);
  $factory->save_tempfiles(1);
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

  return $alignString;
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
    return "no hits" unless @pep_homol;

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
    return (\@bestIDs);
} # end Best_BLAST_Match  

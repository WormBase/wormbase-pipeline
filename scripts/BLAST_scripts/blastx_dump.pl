#!/usr/bin/env perl
#
# DESCRIPTION:
#   script to dump single blastx logic_names from an ensembl_database
#   it is possible to set TEST_FEATURE = one_feature_id to just use one
#
# WARNING:
#   the script is a memory hog due to the huge ensembl objects being handled.
#   therefore some shortcuts were used like the while shift(@bla) loops, to
#   reduce the memory footprint. If you change it, DON'T use foreach and a large
#   array of EnsEMBL objects, it invites disaster as it makes a copy of the array.
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2011-01-20 11:09:57 $ 

my $usage = <<USAGE;
blastx_dump.pl options:
      -logicname LOGIC   name a logic_name from the ensembl database
      -database  DB_NAME name of the ensembl database
      -toplevel          use the largest coordinate system for dumping (most times chromosomes)
      -test              use only wormpepX / 1 slice for speeding up testing
      -store STORABLE    use a premade storable
      -selfhits          remove selfhits from the blastx hits
      -outfile FILENAME  print the dumps to a specific file
      -debug USERNAME    send log messages to one user only
      -sequence SEQUENCE only process the sequence specified
USAGE



###################
use lib '/software/worm/ensembl/ensembl/modules/';
use lib '/software/worm/lib/bioperl-live/';
###################

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use IO::File;

use lib $ENV{CVS_DIR};

use Storable;
use Wormbase;
use strict;

# for testing run it only on one clone
# setting the environment variable TEST_FEATURE limits the whole thing to the specified feature

my ($logicname,$test,$dbname,$toplevel,$store,$selfhits,$outfile,$debug);
my @sequences;
GetOptions(
    "logicname=s" => \$logicname,
    'database=s'  => \$dbname,
    'toplevel'    => \$toplevel,
    'test'        => \$test,
    'store=s'     => \$store,
    'selfhits'    => \$selfhits,
    'outfile=s'   => \$outfile,
    'debug=s'     => \$debug,
    'sequence=s'  => \@sequences,
) || die ($usage);


my $wormbase;
if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
        -debug    => $debug,
        -test     => $test,
    );
}


my $outf = new IO::File $outfile,'w';

# evil elegans hack to get the clone id/accession mapping
my %clone2acc;
my %acc2clone;
if ($dbname=~/elegans/){
 $wormbase->FetchData('clone2accession',\%clone2acc);
 %acc2clone=reverse %clone2acc;
 undef %clone2acc;
}

# connect to database

my $database_name = ($dbname ||'worm_ensembl_briggsae_test');
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => 'farmdb1',
        -user   => 'wormro',
        -dbname => $database_name,
        -port   => 3306,
    );

    # actually it should really get a $wormbase object to get db-names and so on ...

my $slice_adaptor = $db->get_SliceAdaptor();
my $feature_adaptor = $db->get_ProteinAlignFeatureAdaptor();

# get all superlinks
my @superlinks;
if ($sequences[0]){
  @sequences = split(/,/,join(',',@sequences));
  foreach my $sequence(@sequences) {
    push @superlinks, $slice_adaptor->fetch_by_region('toplevel',$sequence);
  }
}elsif ($toplevel){
  @superlinks = @{$slice_adaptor->fetch_all('toplevel')} if $toplevel;
}else {
  @superlinks = @{$slice_adaptor->fetch_all('clone')};
  @superlinks = @{$slice_adaptor->fetch_all('superlink')} if scalar(@superlinks)==0;
}

# hardcoded bits
my %logic2type = (
    remaneiX       => 'wublastx_remanei',
    brigpepX       => 'wublastx_briggsae',
    wormpepX       => 'wublastx_worm',
    GadflyX        => 'wublastx_fly',
    ipi_humanX     => 'wublastx_human',
    slimtremblX    => 'wublastx_slimTrEmbl',
    yeastX         => 'wublastx_yeast',
    slimswissprotX => 'wublastx_slimSwissProt',
    ppapepX        => 'wublastx_pristionchus',
    jappepX        => 'wublastx_japonica',
    brepepX        => 'wublastx_brenneri',
);

my %logic2prefix = (
    remaneiX       => 'RP:',
    brigpepX       => 'BP:',
    GadflyX        => 'FLYBASE:',
    ipi_humanX     => 'ENSEMBL:',
    slimtremblX    => 'TR:',
    yeastX         => 'SGD:',
    slimswissprotX => 'SW:',
    wormpepX       => 'WP:',
    ppapepX        => 'PP:',
    jappepX        => 'JA:',
    brepepX        => 'CN:'
);

$logicname='wormpepX' if $test;
die "ERROR: bad logicname\n" unless $logic2type{$logicname};

my %cds2wormpep=%{&read_table()};

# iterate over all superlinks
while (my $link = shift @superlinks){
    print STDERR "iterating over: ",$link->seq_region_name,"\n" if $ENV{TEST};

    # hack around the limitations of the clone slices
    # by projecting them onto the toplevel overlapping genes should be fetched
    if (! $toplevel) {
      my $seq_pairs = $link->project('toplevel');
      my $parent=shift @$seq_pairs;
      $link = $parent->to_Slice();
    }


#######
# split the slice into pieces ?
#######

    # get all ProteinAlignFeatures on it
    my @dafs =  @{ $feature_adaptor->fetch_all_by_Slice($link,$logicname) };
    @dafs    =  ($feature_adaptor->fetch_by_dbID($ENV{TEST_FEATURE}) ) if $ENV{TEST_FEATURE}; # 5970637

    my $type=$logic2type{$logicname} || die "cannot find $logicname\n";
    my $prefix=$logic2prefix{$logicname};

    &clone_ace($link->seq_region_name,$type,$link->length);
    
    my @features=&remove_selfhits(\@dafs,$link);
    undef @dafs;
    @features=&filter_features(\@features,$link->length);
    while (my $f= shift @features){
#       print $outf $link->seq_region_name;
#       printf( "  %s %d-%d (%+d)\t=> %d-%d (%+d) %s %s\n",
#       $_->hseqname, $_->hstart, $_->hend, $_->hstrand,  $_->start,  $_->end, $_->strand,$_->p_value,$_->cigar_string)
        my $hname=$f->hseqname;
        $f->hseqname($prefix.$hname);
        &feature2ace($f,$type); #if $f->p_value < 0.0001;
    }

    print $outf "\n";
#   last if $ENV{TEST};
}



# create the Sequence link
sub clone_ace {
    my ($name,$type,$length)=@_;
    my $tmp_name=$name;
    $tmp_name=~s/\.\d+$//;
    $name = ($acc2clone{$tmp_name}||$name); # because of clones having accession numbers in the db
    print $outf "Sequence : \"$name\"\n";
    print $outf "Homol_data \"$name:$type\" 1 $length\n\n";
    print $outf "Homol_data : \"$name:$type\"\n"
}

# feature2ace
sub feature2ace {
    my ($f,$org)=@_;
    # flip coordinates if - strand
    &_feature_flip($f) if ($f->strand < 0);

    my $blasthit_id=($cds2wormpep{$f->hseqname}||$f->hseqname);

    my $part1 = "Pep_homol\t\"$blasthit_id\" \"$org\"  %.3f %d %d %d %d";
    my $debug_block="// DEBUG: strand(%s) cigar(%s) pvalue(%s) dbID(%s)" if $ENV{TEST};
    my $blocks=&cigar2ace($f);
    my $log_e= &p_value($f->p_value);
    if (scalar @$blocks >1 ){
       foreach my $block (@$blocks) {
            printf $outf ("$part1 %s $debug_block\n",
            $log_e,$f->start,$f->end,$f->hstart,$f->hend,$block,$f->strand,$f->cigar_string,$f->p_value,$f->dbID);
       }
    }
    else {
        printf $outf ("$part1 $debug_block\n",
            $log_e,$f->start,$f->end,$f->hstart,$f->hend,$f->strand,$f->cigar_string,$f->p_value,$f->dbID);
    }
}

# p_value shorthand
sub p_value {
    my ($p)=@_;
    my $log = (eval($p) > 0 ? -(log(eval($p)))/log(10) : 999.9); # if e-value is zero set it to 999.9
    return $log;
}

# convert cigar to ACeDB align pairs
sub cigar2ace {
    my ($f)=@_;

    my @aceblocks;

    my $aligned_pos_pep=0;
    my $aligned_pos_dna=0;
    my $cigar="${\$f->cigar_string}";

    while ($cigar=~m/(\d*)([MID])/g){
        my $size=($1||1);
        my $type=$2;
        if ($type eq 'D'){
            $aligned_pos_pep+=($size/3);  
        }
        elsif ($type eq 'M') {

            my ($d_start);
            if ($f->strand > 0 ){ $d_start= $f->start+$aligned_pos_dna }
            else { $d_start= -(-$f->start+$aligned_pos_dna) }

            my $align= "AlignDNAPep $d_start ".($f->hstart+$aligned_pos_pep)." ".($size/3);

            push @aceblocks, $align;
            $aligned_pos_pep+=($size/3);
            $aligned_pos_dna+=($size);
        }
        elsif ($type eq 'I') {
            $aligned_pos_dna+=($size);
        }
    }

    return \@aceblocks;       
}

# flip start and stop
sub _feature_flip {
    my ($f)=@_;
    my $t=$f->start;
    $f->start($f->end);
    $f->end($t);
}


#################
# only  for C.elegans ?
# else we need some hack to pull the data from different flatfile(formats)
#
# to remove selfhits
sub remove_selfhits {
    my ($features,$link,$wormbase)=@_;
    my (%cds2wormpep,@results);

    #$wormbase->FetchData( 'cds2wormpep', \%cds2wormpep ) if $wormbase;
    %cds2wormpep=%{&read_table()}; #unless $wormbase;
    die('cannot fetch cds2wormpep') unless %cds2wormpep;

    my @index=&build_search_struct($link->get_all_Genes()); # using a dirty search structure

    FEATURE: while (my $feature = shift @$features){
        my $name=($cds2wormpep{$feature->hseqname}||$feature->hseqname );

        if ($selfhits) { # should probably not do this

            print STDERR "processing selfhits for: ${\$feature->hseqname} ($name)\n" if $ENV{TEST};
            foreach my $hit(&search(\@index,$feature)){
                my %ids;
                map {$ids{$cds2wormpep{$_}}=1} @{$hit->{ids}};
                # kick it if it is the same
                map {print STDERR "-- found $_\n"} keys %ids if $ENV{TEST};
                if ($ids{$feature->hseqname} || $ids{$cds2wormpep{$feature->hseqname}}){
                    print STDERR "kicking: ",$feature->dbID," ", $feature->hseqname," ",$feature->start," ",$feature->end," ",$feature->p_value,"\n" if $ENV{TEST};
                    next FEATURE;
                }
            }
        }
        $feature->hseqname($name);
        push @results,$feature;
    }   
    return @results;
}

# build sorted search list
sub build_search_struct {
    my ($genes)=@_;
    
    # lets see what we actually need:
    # start, end, {$_->translation->stable_id} @{$hit->get_all_Transcripts}
    my @slim_genes = map { 
        {
            start => $_->start,
            end => $_->end,
            ids => [map {$_->translation->stable_id} @{$_->get_all_Transcripts}]
        }
    } @$genes;
    my @sorted_genes=sort {$a->{start} <=> $b->{start} || $a->{end} <=>$b->{end}} @slim_genes
    
}

# search in list
sub search {
    my ($list,$feature)=@_;
    my @hits;
    foreach my $g (@{$list}){
        push(@hits,$g) if ($feature->start <= $g->{end} && $feature->end >= $g->{start});
#       last if ($feature->start > $g->end && $feature->end > $g->end); # need to comment that if it returns crap
    }
    return @hits;
}

# remove < 75% of evalue features from 100bp windows
# basic idea:
#   create bins for every 100bp and put all overlapping features into them
#   iterate over the bins and add the best 5 within 25% of the highest bin p-value into the results
sub filter_features {
    my ($_features,$length)=@_;
    my %f_features;
    
    # bin size
    my $size=100;

    # should I bin them instead?
    my @bins;

    # put the feature in all bins between start and stop
    while( my $f= shift @$_features ){
        my ($_start,$_end)=sort {$a <=> $b } ($f->start,$f->end);
#       next if $_start < 0;
#       next if $_end > $length;
        for (my $n=int($_start/$size);$n <=int($_end/$size);$n++){
            push @{$bins[$n]},$f;
        }
    }

    # get the best 5 within 25% of the best hsp and add them to a hash
    foreach my $bin(@bins) {
        next unless $bin; # skip empty bins
        my $best=0;
        my $max_hsp=0;
        my @sorted_bin = sort {$a->p_value <=> $b->p_value} @$bin;
        $best=&p_value($sorted_bin[0]->p_value);
#       map { $f_features{$_->dbID}=$_ if ($max_hsp++ <5)} @$bin; # <= cutoff place, lets try 25%
        map { $f_features{$_->dbID}=$_ if (&p_value($_->p_value) > $best*0.5 && $max_hsp++ <5)} @$bin; # <= cutoff place, lets try 25%
    }
    
    # flatten hash into an array
    my @_filtered=values %f_features;
    return @_filtered;
}

sub get_latest_pep {
   my @species =qw(wormpep remapep brigpep ppapep jappep brepep);
   my @history_files;
   SPECIES: foreach my $s (@species){
       my @files = sort {$b cmp $a} glob("~wormpub/BUILD/WORMPEP/$s*/*history*");
           foreach my $file(@files){
               next if $file=~/666|665/;
                           push (@history_files,$file);
               next SPECIES;
           }
   }
   return @history_files;
}

# create mapping file
sub read_table {
    my ($wormbase)=@_;
    my %mapping;
    # evil hardcoded path. should use wormpep instead.
    my @files=&get_latest_pep;
    foreach my $f(@files){
       my $file= new IO::File $f, 'r' || die(@!);
       while (<$file>){
         my @a=split;
         $mapping{$a[0]}=$a[1];
       }
       $file->close;
    }
    return \%mapping;
}

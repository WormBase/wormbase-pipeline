#!/usr/bin/env perl
# DESCRIPTION:
#   script to dump single blastx logic_names from an ensembl_database
#   it is possible to set TEST_FEATURE = one_feature_id to just use one
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2015-04-27 13:19:48 $ 

my $usage = <<USAGE;
blastx_dump.pl options:
      -logicname LOGIC   name a logic_name from the ensembl database
      -database  DB_NAME name of the ensembl database
      -test              use only wormpepX / 1 slice for speeding up testing
      -store STORABLE    use a premade storable
      -selfhits          remove selfhits from the blastx hits
      -outfile FILENAME  print the dumps to a specific file
      -debug USERNAME    send log messages to one user only
      -sequence SEQUENCE only process the sequence specified
USAGE

###################
use lib $ENV{'WORM_PACKAGES'} . '/ensembl/ensembl/modules/';
use lib $ENV{'WORM_SW_ROOT'} . '/lib/bioperl-live/';
use lib $ENV{CVS_DIR};

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use IO::File;
use Storable;
use Wormbase;
use strict;

my ($logicname,$test,$dbname,$toplevel,$store,$selfhits,$outfile,$debug);
my @sequences;
GetOptions(
    "logicname=s" => \$logicname,
    'database=s'  => \$dbname,
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

# connect to database
my $database_name = ($dbname ||'worm_ensembl_'.$wormbase->species);
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $ENV{'WORM_DBHOST'},
        -user   => 'wormro',
        -dbname => $database_name,
        -port   => $ENV{'WORM_DBPORT'},
    );

my $slice_adaptor = $db->get_SliceAdaptor();
my $feature_adaptor = $db->get_ProteinAlignFeatureAdaptor();

# get all superlinks
my @superlinks;
if ($sequences[0]){
  @sequences = split(/,/,join(',',@sequences));
  foreach my $sequence(@sequences) {
    push @superlinks, $slice_adaptor->fetch_by_region('toplevel',$sequence);
  }
} else {
  @superlinks = @{$slice_adaptor->fetch_all('toplevel')};
}

# hardcoded bits
my %logic2type = (
		  wormpepx       => 'wublastx_worm',
		  ppapepx        => 'wublastx_pristionchus',
		  jappepx        => 'wublastx_japonica',
		  brugpepx       => 'wublastx_brugia',
		  brigpepx       => 'wublastx_briggsae',
		  remapepx       => 'wublastx_remanei',
		  brepepx        => 'wublastx_brenneri',
		  ovolpepx       => 'wublastx_ovolvulus',
                  srapepx        => 'wublastx_sratti',
		  gadflyx        => 'wublastx_fly',
		  ipi_humanx     => 'wublastx_human',
		  yeastx         => 'wublastx_yeast',
		  slimswissprotx => 'wublastx_slimSwissProt',
		 );		 


my %logic2prefix = (
		    wormpepx       => 'WP:',
		    ppapepx        => 'PP:',
		    jappepx        => 'JA:',
		    brugpepx       => 'BM:',
		    ovolpepx       => 'OV:',
                    srapepx        => 'SRP:'
		    brigpepx       => 'BP:',
		    remapepx       => 'RP:',
		    brepepx        => 'CN:',
		    gadflyx        => 'FLYBASE:',
		    ipi_humanx     => 'ENSEMBL:',
		    yeastx         => 'SGD:',
		    slimswissprotx => 'SW:',
		   );


$logicname='wormpepx' if $test;
die "ERROR: bad logicname\n" unless $logic2type{$logicname};

my %cds2wormpep=%{&read_table()};

# iterate over all superlinks
while (my $link = shift @superlinks){
    print STDERR "iterating over: ",$link->seq_region_name,"\n" if $ENV{TEST};

    # get all ProteinAlignFeatures on it
    my @dafs =  @{ $feature_adaptor->fetch_all_by_Slice($link,$logicname) };
    @dafs    =  ($feature_adaptor->fetch_by_dbID($ENV{TEST_FEATURE}) ) if $ENV{TEST_FEATURE}; # 5970637

    my $type=$logic2type{$logicname} || die "cannot find $logicname\n";
    my $prefix=$logic2prefix{$logicname};

    &seq_ace($link->seq_region_name,$type,$link->length);
    
    my @features=&remove_selfhits(\@dafs,$link);
    undef @dafs;
    @features=&filter_features(\@features,$link->length);
    while (my $f= shift @features){
        my $hname=$f->hseqname;
        $f->hseqname($prefix.$hname);
        &feature2ace($f,$type);
    }

    print $outf "\n";
    last if $ENV{TEST};
}

# create the Sequence link
sub seq_ace {
    my ($name,$type,$length)=@_;
    my $tmp_name=$name;
    $tmp_name=~s/\.\d+$//;
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
    my $log = ($p > 0 ? -(log($p))/log(10) : 999.9); # if e-value is zero set it to 999.9
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
            else { $d_start= $f->start-$aligned_pos_dna }

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
# to remove selfhits
sub remove_selfhits {
    my ($features,$link,$wormbase)=@_;
    my (%cds2wormpep,@results);

    %cds2wormpep=%{&read_table()};
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
                    print STDERR "kicking: ",$feature->dbID," ", $feature->hseqname," ",$feature->start," ",$feature->end," ",$feature->p_value,"\n" 
			if $ENV{TEST};
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
    
    my @slim_genes = map { 
        {
            start => $_->start,
            end => $_->end,
            ids => [map {$_->translation->stable_id} @{$_->get_all_Transcripts}]
        }
    } @$genes;
}

# search in list
sub search {
    my ($list,$feature)=@_;
    my @hits;
    foreach my $g (@{$list}){
        push(@hits,$g) if ($feature->start <= $g->{end} && $feature->end >= $g->{start});
    }
    return @hits;
}

# remove < 50% of evalue features from 100bp windows
# basic idea:
#   create bins for every 100bp and put all overlapping features into them
#   iterate over the bins and add the best 5 within 50% of the highest bin p-value into the results
sub filter_features {
    my ($_features,$length)=@_;
    my %f_features;
    
    # bin size
    my $size=100;

    my @bins;

    # put the feature in all bins between start and stop
    while( my $f= shift @$_features ){
        my ($_start,$_end)=sort {$a <=> $b } ($f->start,$f->end);
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
        map { $f_features{$_->dbID}=$_ if (&p_value($_->p_value) > $best*0.5 && $max_hsp++ <5)} @$bin; # <= cutoff place, best 5 and within 50%
    }
    
    # flatten hash into an array
    my @_filtered=values %f_features;
    return @_filtered;
}

sub get_latest_pep {
   my @species =qw(wormpep remapep brigpep ppapep jappep brepep brugpep ovolpep srapep);
   my @history_files;
   SPECIES: foreach my $s (@species){
       my @files = sort {$b cmp $a} glob($ENV{'WORMPUB'}."/BUILD/WORMPEP/$s*/*history*");
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

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Getopt::Long;
use Data::Dumper;
use Coords_converter;
use Wormbase;
use Modules::Transcript;
use Modules::CDS;
use Modules::GFF_feature;
use POSIX;


my @features;     # stores GFF_feature objs of features of interest
my @exclusions;   # stores GFF_feature objs of features to exclude

my @exclude_gff;  # gff files of exclusion data
my @feature_gff;  # gff file of feature of interest

my %feature;     # feature methods
my %exclude;      # exclusion methods

my $config;       #text file to specify parameters

my $database;
my $chromosome;   # so that Coords converter can give clone details

our %score;
my $window = 200;

GetOptions ( 
	    "database:s"    => \$database,
	    "chromosome:s"  => \$chromosome,

	    "feature_gff:s" => \@feature_gff,
	    "exclude_gff:s" => \@exclude_gff,

	    "config:s"      => \$config
	   );

die unless $chromosome;

if ($config) {
  &parse_config($config);
}

#else {
#  @feature_gff  = split(/,/,join(',',@feature_gff));
#  @exclude_gff  = split(/,/,join(',',@exclude_gff));
#  @f_source     = split(/,/,join(',',@f_source));
#  @x_source     = split(/,/,join(',',@x_source));
#  @f_method     = split(/,/,join(',',@f_method));
#  @x_method     = split(/,/,join(',',@x_method));
#}

$database = glob("~wormpub/DATABASES/current_DB") unless $database;
&load_exlusions(\@exclude_gff,\@exclusions,\%exclude);#"curated","Pseudogene"," Transposon","RepeatMasker");

my $method_check;#= join('|',@f_method);
my $source_check;# = join('|',@f_source);

foreach my $feature_file ( @feature_gff ) {
  open( GFF, "<$feature_file") or die "cant open $feature_file GFF file\t$!\n";
  print "Reading feature info from $feature_file\n\n";
  while(<GFF>) {
    my @data = split;
    next unless( defined $feature{$data[1]} and $feature{$data[1]} eq $data[2] );
    # my ($id) = $data[9] =~ /\"$type:(\S+)\"/;
    my $id = $data[9];
    my $feature = GFF_feature->new($id, $data[3], $data[4], $data[5], $data[6] );
    push( @features, $feature);  
    #CHROMOSOME_I    wublastx        protein_match   74626   74943   38.19   +       0       Target "Protein:TR:O77132" 105 206
  }
}

my @ordered_exclusions = sort { $a->start <=> $b->start } @exclusions;
# check overlap
print "checking feature overlaps\n\n";
FEATURE: foreach my $feature (@features) {
  my $start = $feature->start;
  my $end = $feature->end;

  my $match = 0;
  GENE: foreach my $exclude (@ordered_exclusions) {
    if( $start > $exclude->end ) {
      next GENE;
    }
    elsif ( $end < $exclude->start ) {
      last GENE;
    }
    elsif( ( $start > $exclude->start and $start < $exclude->end ) or
	   (   $end > $exclude->start and $end   < $exclude->end ) or 
	   ( $start < $exclude->start and $end   > $exclude->end )
	 ) {
      $match = 1; # matches so report and move on
      last GENE ;
    }
  }
  &addscore($feature) unless ($match == 1);
}

# output hits
print "outputting  . . . \n\n";
my $coords = Coords_converter->invoke($database,1);
my $count = 0;
print "Clone\tstart\tend\tscore\tGFF_coord\n";
foreach my $region ( sort{ $score{$b} <=> $score{$a} } keys %score ) {
  my $scaled_region = ($region * $window) - ($window * 0.5);
  my @clone = $coords->LocateSpan("$chromosome",$scaled_region,$scaled_region + $window);
  print "$clone[0]\t$clone[1]\t$clone[2]\t$score{$region}\t",$scaled_region,"\n";
  last if $count++ > 50;
}


#########################################################
## # ## # #  suboutines


sub load_exlusions 
  {
    my $files      = shift;
    my $exclusions = shift;
    my $exclude    = shift;

    foreach my $file ( @{$files} ) {
      # parse GFF file to get CDS, exon and cDNA info
      open( GFF,"<$file") or die "cant open $file to read gene info\n";
      print "reading exclusion info from $file\n\n";
    LINE: while (<GFF>) {
	my @data = split;
	next unless( defined $$exclude{$data[1]} and $$exclude{$data[1]} eq $data[2] );
	my $gff_exclude = GFF_feature->new($data[9], $data[3], $data[4], $data[5], $data[6]);
	push( @{$exclusions},$gff_exclude);
	next LINE;
      }
    }
  }

sub addscore
  {
    my $feature = shift;
    my $midpoint = ($feature->start + $feature->end) / 2;
    my $mid_region = floor($midpoint / $window);
    unless ( $feature->score =~ /\d/ ) {
      my $span = abs($feature->end - $feature->start);
      $feature->score($span);
    }
    $score{$mid_region} += $feature->score;
  }


sub parse_config 
  {
    my $config = shift;

    my @f_method;
    my @x_method;
    open(CF,"<$config") or die "cant open $config :$!\n";

    my $parameter;
    while( <CF> ) {
      next unless /\w/;
      chomp;
      if( />(\w+)/ ){
	CASE:{
          ($1 eq "feature_gff") && do { $parameter = \@feature_gff; last CASE; };
          ($1 eq "exclude_gff") && do { $parameter = \@exclude_gff; last CASE; };
	  ($1 eq "feature")     && do { $parameter = \@f_method;    last CASE; };
	  ($1 eq "exclude")     && do { $parameter = \@x_method;    last CASE; };
        }
      }
      else {
	push(@{$parameter},"$_");
      }
    }
    %feature = split(/,/, join(",",@f_method) );
    %exclude = split(/,/, join(",",@x_method) );
  }


__END__

=pod

=head2 NAME feature_overlap_checker.pl

=head1 USAGE

=over 4

=item feature_overlap_checker.pl -chromosome X -config config

=item feature_overlap_checker.pl -chromosome X -f_method Genefinder -f_source CDS -x_method gene -x_source gene -exclude_gff CHROMOSOME_X.WBGenes.gff -feature_gff CHROMSOME_X.Genefinder.gff

=back

This script searches for specified features that do not overlap another type of specified features ( exclusion features ).  All features are defined by the method and source as written in the GFF files being used.  Multiple types of features and exclusion features can be processed together 

eg if you want to look for places where waba hits and blastp homologies do not overlap gene predictions you can specify exclusion features such as repeats pseudogenes and features as waba_coding and wublastx specifying all appropriate sources and methods;

This means that command line parameters can get v. long so can be written to a config file like this. .  .

 >feature_gff
 CHROMOSOME_X.Genefinder.gff

 >f_method
 Genefinder

 >f_source
 CDS

 >exclude_gff
 CHROMOSOME_X.WBGene.gff 

 >x_method
 gene

 >x_source
 gene


=cut

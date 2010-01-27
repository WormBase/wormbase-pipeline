#!/software/bin/perl -w

use strict;
use DBI;
use DBD::mysql;
use lib "/software/worm/lib/perl5";
use GD::Graph::lines;

#DATA SOURCE NAME
my $dsn = "dbi:mysql:worm_rgasp:ia64d:3306";

# PERL DBI CONNECT
my $db = DBI->connect($dsn, 'wormro');

my $goi = "AH6.1"; #gene of interest
# PREPARE THE QUERY
my $query = "SELECT transcript,stage,rpkm FROM results where program_id = 1 and transcript = \"$goi\" order by transcript,(field(stage,\"SRX004863\",\"SRX004865\",\"SRX001872\",\"SRX001875\",\"SRX001874\",\"SRX001873\")) limit 200";
my $query_handle = $db->prepare($query);

# EXECUTE THE QUERY
$query_handle->execute();

# BIND TABLE COLUMNS TO VARIABLES
my( $trans,$stage,$rpkm);
$query_handle->bind_columns(undef, \$trans, \$stage, \$rpkm);

my @life_stages = qw (SRX004863 SRX004865 SRX001872 SRX001875 SRX001874 SRX001873);
my %real_stages= ('SRX004863' => 'EEMB',
		  'SRX004865'  => 'LEMB',
		  'SRX001872' => 'L1',
		  'SRX001875' => 'L2',
		  'SRX001874'  => 'L3',
		  'SRX001873' => 'L4',
		  );
my %data;
&reset;
my $gene;

#print "Gene,",join(",",@life_stages),"\n";
# LOOP THROUGH RESULTS
while($query_handle->fetch()) {	
    if(!$gene) {$gene = $trans;}
    if($gene ne $trans) {
	print $gene;
	foreach (@life_stages){ 
	    print ",", $data{$_};
	}
	&draw_graph(\%data, $trans);
	exit;
	print "\n";
	&reset;
	$gene = $trans;
    }
    $data{$stage} = $rpkm;
}
&draw_graph(\%data, $trans);

sub reset {
    foreach (@life_stages){ 
	$data{$_} = 0;
    }
}


sub draw_graph {
    my $data = shift;
    my $trans = shift;
    my @points;
    foreach (@life_stages){ 
	push(@points,$$data{$_});
    }
    my @s = map($real_stages{$_},@life_stages);
    my @ave  = (34.48, 30.20, 35.60, 33.62, 30.62, 27.39 );
    my @graph_data = (\@s,
		      \@ave,
		      \@points
		      );
    
    my $mygraph = GD::Graph::lines->new(600, 300);
    $mygraph->set(
		  x_label     => 'Life stage',
		  y_label     => 'RPKM',
		  title       =>  "RPKM for $trans" ,
		  # Draw datasets in 'solid', 'dashed' and 'dotted-dashed' lines
		  line_types  => [3,1],
		  # Set the thickness of line
		  line_width  => 2,    # Set colors for datasets
		  dclrs       => ['black', 'red'],
		  legend_placement => 'CR'
		  ) or warn $mygraph->error;

    $mygraph->set_legend('stage ave',$trans);
    my $myimage = $mygraph->plot(\@graph_data) or die $mygraph->error;

    my $png;
    open($png,">$trans.png");
    print $png $myimage->png;
    close $png;
}

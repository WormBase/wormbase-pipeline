package rpkm_graph;use strict;
use DBI;
use DBD::mysql;
use GD::Graph::lines;

my @STAGES = qw (SRX004863 SRX004865 SRX001872 SRX001875 SRX001874 SRX001873);
my %STAGE_NAMES = ('SRX004863' => 'EEMB',
		  'SRX004865'  => 'LEMB',
		  'SRX001872' => 'L1',
		  'SRX001875' => 'L2',
		  'SRX001874'  => 'L3',
		  'SRX001873' => 'L4',
		  );


sub new {
    my $class = shift;
    my $self = {};
    bless ($self,$class);

    
    #DATA SOURCE NAME
    my $dsn = "dbi:mysql:worm_rgasp:ia64d:3306";

    # PERL DBI CONNECT
    my $db = DBI->connect($dsn, 'wormro');
    $self->{'db'} = $db;
    return $self;
}

sub submit_genes {
    my $self = shift;
    my $genes = shift; #array ref
	
    # PREPARE THE QUERY
    my $query = "SELECT transcript,stage,rpkm FROM results where program_id = 1 and transcript = ? order by field(stage,\"SRX004863\",\"SRX004865\",\"SRX001872\",\"SRX001875\",\"SRX001874\",\"SRX001873\")";
    my $query_handle = $self->{'db'}->prepare($query);
    my %data;

    foreach my $gene (@$genes) {

    # EXECUTE THE QUERY
	$query_handle->execute($gene);

    # BIND TABLE COLUMNS TO VARIABLES
	my( $trans,$stage,$rpkm);
	$query_handle->bind_columns(undef, \$trans, \$stage, \$rpkm);

   # LOOP THROUGH RESULTS
	while($query_handle->fetch()) {	
	    push(@{$data{$trans}},$rpkm);
	}
    }
    my $image = $self->draw_graph(\%data);
}

sub draw_graph {
    my $self = shift;
    my $data = shift;
    
    #fixed data
    my @s = map($STAGE_NAMES{$_},@STAGES);
    my @ave  = (34.48, 30.20, 35.60, 33.62, 30.62, 27.39 );

    my @legend_list; # keep genes in order for legend

    #array of array refs starting with life stages and average rpkm per stage as point of ref.
    my @graph_data = (\@s,
		      \@ave,
		      );
    #gene data points
    foreach my $gene (keys %$data){ 
	push(@legend_list,$gene);
	push(@graph_data,\@{$$data{$gene}});
    }
    
    my $mygraph = GD::Graph::lines->new(600, 300);
    $mygraph->set(
		  x_label     => 'Life stage',
		  y_label     => 'RPKM',
		  x_label_position => '0.2',
		  title       =>  "RPKM" ,
		  # Draw datasets in 'solid', 'dashed' and 'dotted-dashed' lines
		  line_types  => [3,1,1,1,1,1],
		  # Set the thickness of line
		  line_width  => 2,    # Set colors for datasets
		  dclrs       => ['black', 'red','blue','green','purple','orange'],
		  legend_placement => 'CR',
		  ) or warn $mygraph->error;

    $mygraph->set_legend_font('gdMediumBoldFont',12);
    $mygraph->set_legend('stage ave',@legend_list);
    my $myimage = $mygraph->plot(\@graph_data) or die $mygraph->error;
    
    return \$myimage;
#    my $png;
#    open($png,">/tmp/$trans.png");
#    print $png $myimage->png;
#    close $png;
#    return
}

1;

#!/software/bin/perl
##########################################################
# generates Jeff Holmes allele dumps
# He has his special format wishes (look under print_one)
# 
# AUTHOR: mh6@sanger.ac.uk
#

# MAIN

use Getopt::Long;
use lib $ENV{CVS_DIR};
use lib '/software/worm/lib/site_perl';
use Wormbase;
use IO::File;

use Storable;
# use YAML::Syck;
# $YAML::Syck::ImplicitTyping = 1;
# use CGI qw/:standard/;
use strict;
	
my ($store,$test,$debug,$database,$file,$wormbase);
GetOptions(
		"store=s"    => \$store,
		"test"       => \$test,
		"debug=s"    => \$debug,
		"database=s" => \$database,
		"file=s"     => \$file,
) || die `perldoc $0`;
	

if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database)}

$database = $wormbase->autoace unless $database;
my $log = Log_files->make_build_log($wormbase);

my $jeff = new Jeff $database,$wormbase || die 'cannot create connections to $database\n';
$jeff->get_all_alleles;
$log->write_to("dumping ".(scalar $jeff->get_alleles)." knockout consortium alleles\n");

my $outfile= new IO::File "|/software/worm/bin/tidy -xml|bzip2 -9 -c > $file" if $file;

$jeff->print_alleles($outfile);
$outfile->close;
$jeff->DESTROY;

$log->write_to("finished dumping to $file\n");

my $mailto=$debug?"$debug\@sanger.ac.uk":"All";
$log->mail($mailto,"BUILD REPORT: dump_ko.pl");

# class definitions

package Jeff;

use Ace;
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Coords_converter;

sub new {
	my $class = shift;
	my $db    = shift;
	my $wormbase = shift;
	$wormbase->{'tace'}='/software/worm/bin/acedb/tace';
	my $self = {
		dbh => Ace->connect(-path=>$db),
		conv => Coords_converter->invoke(glob($db),0, $wormbase),
		wormbase => $wormbase,
	};
	bless $self, $class;
	return $self;
}

sub DESTROY {
	my $self= shift;
	$self->{dbh}->close;
	$self->{conv}=undef;
	undef $self;
}

# ACCESSORS 

# accessor for the database handle
sub dbh { 
	my $self = shift;
	return $self->{dbh};
}

# accessor for wormbase
sub wormbase {
	my $self= shift;
	return $self->{wormbase};
}

# accessor for the Coords_converter
sub conv {
	my $self=shift;
	return $self->{conv};
}

# GETTERS and SETTERS

sub set_alleles {
	my $self = shift;
	$self->{allele}=\@_ if (length @_ > 0);
	return @{$self->{allele}};
}

sub get_alleles {
	my $self = shift;
	return @{$self->{allele}}
}

sub get_all_alleles {
	my $self = shift;
	$self->set_alleles($self->dbh->fetch(-query => 'Find Variation WHERE KO_consortium_allele'));
}

# CLASS METHODS

sub print_alleles {
	my ($self,$file)= @_;
	$file = \*STDOUT unless $file;

	my @timelist=localtime();
	my $time=sprintf('%02u/%02u/%u',$timelist[3],$timelist[4],$timelist[5]+1900);
	print $file "<report name=\"KO-Report\" date_generated=\"$time\">\n";
	foreach my $allele($self->get_alleles){
		&print_one($allele,$self,$file);
	}
	print $file "</report>\n";
}

# CLASS FUNCTIONS

sub print_one {
	my ($allele,$self,$file) = @_;
	# name chromosome start stop left_flank right_flank insert genes strain strain_genotype strain_location
	
	# if a value is undef set it to " " instead
	
	my $name = $allele->name;
	my $clone = $allele->Sequence || " ";
	my $left_flank  = $allele->Flanking_sequences||" ";
	my $right_flank = (length $left_flank > 2 ?$left_flank->right:" ");
	my @genes=  $allele->Gene ? map {$_->Sequence_name} $allele->Gene:" ";
	###
	my @strains;
	foreach my $strain($allele->Strain){
		my %strain_;
		$strain_{name}="$strain";
		$strain_{genotype} = length($strain)>2 && $strain->Genotype ? "<![CDATA[${\$strain->Genotype}]]>" : " ";
		$strain_{location} = length($strain)>2 && $strain->Location ? "${\$strain->Location}" : " ";
		push @strains, \%strain_;
        }
	###
	my $start=  $allele->Sequence && $allele->Sequence->at("SMap.S_child.Allele.$name")->right ?
       	$allele->Sequence->at("SMap.S_child.Allele.$name")->right : " ";
	my $stop =  $start=~/\d+/? $allele->Sequence->at("SMap.S_child.Allele.$name")->right->right : " ";
        my $chrom   = " ";
	my $c_start = " ";
	my $c_stop  = " ";
        ($chrom,$c_start,$c_stop) = get_coords($self,$clone,$start,$stop) if $stop=~/\d+/;
	if ($chrom eq " " && $allele->Sequence){
		$chrom=$allele->Sequence->Source->Source;
	}

	my $insert=$allele->Insertion || " ";

	my %output;
	$output{name}="$name";
	if (("$c_start"."$c_stop")=~/[^\d]/) {
		$output{deletion_size}= " ";
	}
	else {
		$output{deletion_size}= abs($c_start-$c_stop);
	}
	$output{flanking_sequences}->{left}="$left_flank";
	$output{flanking_sequences}->{right}="$right_flank";
	$output{flanking_sequences}->{chromosome}="$chrom";
	$output{flanking_sequences}->{left_coord}="$c_start";
	$output{flanking_sequences}->{right_coord}="$c_stop";
	$output{insertion_sequence}="$insert";

	#pcr_product crud
	$output{pcr_product}=&get_pcr($allele);
	
	$output{datasource}=$self->wormbase->get_wormbase_version_name;

	@{$output{affects}{gene}}= map {"$_"} @genes;
	@{$output{strain}} = @strains;
#	print $file YAML::Syck::Dump(\%output); # should be one day replaced with an XML emitter (maybe TT or XML::Writer)

	print $file xml('variation',\%output);
}


sub get_pcr {
	my $allele=shift;
	my %pcr;
	foreach my $product($allele->PCR_product){
		"$product"=~/_(external|internal)/;
		my $type=$1;
		next unless ($type && $product->Oligo);
		foreach my $_primer ($product->Oligo){
			"$_primer"=~/_(f|b)/;
			my $ptype= $1 eq 'f' ? 'left' : 'right';
			next unless $ptype;
			my $psequence=$_primer->Sequence;
			$pcr{"${type}_${ptype}_seq"}="$psequence";
		}
	}
	return \%pcr;
}

sub get_coords {
	my ($self,$clone,$start,$stop)=@_;
	my ($chr,$c_start)=$self->conv->Coords_2chrom_coords($clone,$start);
	my ($chr2,$c_stop)=$self->conv->Coords_2chrom_coords($clone,$stop);
	return ($chr,$c_start,$c_stop);
}


# should be like xml(tag_name,content,attrib)
# attrib is optional
sub xml {
	my ($tag_name,$content,$attrib)=@_;

	my $attrib_long="";
	if ($attrib && ref($attrib) eq 'HASH'){
		while (my($k,$v)=each %$attrib){
			$attrib_long.=" $k=\"$v\"";
		}
	}

	return "" unless $tag_name;

	if (ref($content) eq 'ARRAY'){
		my $line="";
		foreach my $item(@$content){
			$line.=&xml($tag_name,$item);
		}
		return $line;
	}
	elsif (ref($content) eq 'HASH'){
		my $line;
		while (my($k,$v)=each %$content){
			$line.=&xml($k,$v);
			# die(
			#	"$k $v\n".
			#	YAML::Syck::Dump($content)
			# ) unless( defined($v) && defined($k) && defined(&xml($k,$v)));
		}
		return &xml($tag_name,$line);
	}
	else {  
		#chomp $content;
		$content="" unless $content;
		$content = &right_shift($content);

		return "<$tag_name"."$attrib_long>\n$content</$tag_name>\n";
	}
}

sub right_shift {
	my $text=shift;
	chomp $text;
	return (join "\n", (map {"  $_"} (split /\n/,$text)))."\n";
}

=pod

=head1 NAME 

dump_ko.pl

=head1 DESCRIPTION 

script to dump Jeffs knockout consortium data

As Jeff wants a dump of his data to synchonise his data (can we see a potential disaster here?)
The script extracts all the data from ACeDB and converts the clone coordinates to chromsomome.
If the script is called directly, it will use the main block and dump the data.


=head1 USAGE

dump_ko.pl [-store WORMBASE_STORABLE] [-test] [-debug USERNAME] [-database DATABASE_DIRECTORY] [-file OUTFILE]
It will default to BUILD / *STDOUT if no options are set and the outputfile will be bzip2 compressed.

=head1 CONSTRUCTIR

=head2 new($acedb_path)

use like: $j = new Jeff "path_to_acedb"
opens up 2 databases connections for access and conversion.

=head1 CLASS METHODS

=head2 get_all_alleles

populates the alleles attribute with the Knockout Consortium Alleles from the ACeDB. The identification is by tag KO_consortium_allele.
you can access the alleles with $object->get_alleles 

=head2 print_alleles

prints out the allele data in  "name , chromosome , start , stop , left_flank , right_flank , insert , genes , strain , strain_genotype , strain_location" format.
If a value is undefined/not existing it will print it as space.

=cut

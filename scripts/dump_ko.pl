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
use Wormbase;
use IO::File;
use strict;
	
my ($store,$test,$debug,$database,$file,$wormbase);
GetOptions(
		"store=s"    => \$store,
		"test"       => \$test,
		"debug=s"    => \$debug,
		"database=s" => \$database,
		"file=s"     => \$file,
) || die `perldoc $0`;
	

if ($store) { $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database)}

$database = $wormbase->autoace unless $database;
my $log = Log_files->make_build_log($wormbase);

my $jeff = new Jeff $database,$wormbase || die 'cannot create connections to $database\n';
$jeff->get_all_alleles;
$log->write_to("dumping ".(scalar $jeff->get_alleles)." knockout consortium alleles\n");

my $outfile= new IO::File "|bzip2 -9 -c > $file" if $file;

$jeff->print_alleles($outfile);
$outfile->close;
$log->write_to("finished dumping to $file\n");

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
	};
	bless $self, $class;
	return $self;
}

# ACCESSORS 

# accessor for the database handle
sub dbh { 
	my $self = shift;
	return $self->{dbh};
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
	foreach my $allele($self->get_alleles){
		&print_one($allele,$self,$file);
	}
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
	my @genes=$allele->Gene ||(" ");
	my $strain=$allele->Strain || " ";
	my $strain_genotype = length($strain)>2 && $strain->Genotype ? $strain->Genotype : " ";
	my $strain_location = length($strain)>2 && $strain->Location ? $strain->Location : " ";
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

	printf $file ("%s , %s , %s , %s , %s , %s , %s , %s , %s , %s , %s\n",
		$name,$chrom,$c_start,$c_stop,$left_flank,$right_flank,$insert,
		join(" & ",@genes),$strain,$strain_genotype,$strain_location);
}

sub get_coords {
	my ($self,$clone,$start,$stop)=@_;
	my ($chr,$c_start)=$self->conv->Coords_2chrom_coords($clone,$start);
	my ($chr2,$c_stop)=$self->conv->Coords_2chrom_coords($clone,$stop);
	return ($chr,$c_start,$c_stop);
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

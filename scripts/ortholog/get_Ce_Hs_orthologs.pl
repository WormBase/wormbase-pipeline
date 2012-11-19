#! /usr/bin/perl
#
# Perl script to use the TreeFam API to get all orthologs
# between C. elegans, C. briggsae and C. remanei genes:
# Avril Coghlan, 7-Nov-06.

#-----------------------------------------------------------------------------
#
# will take any random pair of human / C.elegans / C.briggsae / C.remanei / C.brenneri
#
use strict;
use lib 'treefam-api';
use Treefam::DBConnection;

# find the query species:
my $query_species = $ARGV[0]; # eg. CAEEL
my $query_species_long;
# find the target species:
my $target_species = $ARGV[1]; # eg. CAEBR
my $target_species_long = "";
if    ($target_species eq 'CAEEL') { $target_species_long = "Caenorhabditis elegans"; }
elsif ($target_species eq 'CAEBR') { $target_species_long = "Caenorhabditis briggsae";}
elsif ($target_species eq 'CAERE') { $target_species_long = "Caenorhabditis remanei"; }
elsif ($target_species eq 'HUMAN') { $target_species_long = 'Homo sapiens'; }
else { print STDERR "ERROR: target_species $target_species.\n"; exit;            }

if    ($query_species eq 'CAEEL') { $query_species_long = "Caenorhabditis elegans"; }
elsif ($query_species eq 'CAEBR') { $query_species_long = "Caenorhabditis briggsae";}
elsif ($query_species eq 'CAERE') { $query_species_long = "Caenorhabditis remanei"; }
elsif ($query_species eq 'HUMAN') { $query_species_long = "Homo sapiens"; }
else { print STDERR "ERROR: target_species $query_species.\n"; exit;            }

my $dbc = Treefam::DBConnection->new();
my $gh = $dbc->get_GeneHandle;

my %cds2wbgene_id=%{get_common_data(glob('~wormpub/DATABASES/current_DB/COMMON_DATA/cds2wbgene_id.dat'))};
#my %cds2wbgene_id=%{get_common_data(glob('~wormpub/BUILD/autoace/COMMON_DATA/worm_gene2geneID_name.dat'))};

# find all query species genes in TreeFam:
my @query_genes;
if    ($query_species eq 'CAEEL') { @query_genes = $gh->get_all_by_species('CAEEL');}
elsif ($query_species eq 'CAEBR') { @query_genes = $gh->get_all_by_species('CAEBR');}
elsif ($query_species eq 'CAERE') { @query_genes = $gh->get_all_by_species('CAERE');}
elsif ($query_species eq 'HUMAN') { @query_genes = $gh->get_all_by_species('HUMAN');}
else { print STDERR "ERROR: query_species $query_species.\n"; exit;                 }
# get all the $target_species orthologs of the $query_species genes:

my %genes;

foreach my $gene (@query_genes) {
	
   my @orthologs = $gene->get_orthologs($target_species_long);
   my $wbgeneid=($cds2wbgene_id{$gene->ID}||$gene->ID);
   next unless scalar(@orthologs)>=1;
   next unless $cds2wbgene_id{$gene->ID};

   foreach my $ortholog (@orthologs) {
	  my $ortho_geneid=($cds2wbgene_id{$ortholog->ID}||$cds2wbgene_id{$ortholog->ID.'.1'}||$ortholog->ID);
	  # next unless ($cds2wbgene_id{$ortholog->ID}||$cds2wbgene_id{$ortholog->ID.'.1'});
	  # next if $ortholog->bootstrap < 99;
	  $genes{$wbgeneid}.="Ortholog_other $ortho_geneid \"$target_species_long\" From_analysis TreeFam // Ortholog of $query_species ${\$gene->ID} is $target_species ${\$ortholog->ID} bootstrap=${\$ortholog->bootstrap}%\n";
   }
}

while (my($k,$v)=each %genes){
	print "Gene : \"$k\"\n";
	print "$v\n";
}

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

sub get_common_data {
	my ($file)=@_;
	use IO::File;
	my $VAR1;
	die "cannot find file: $file" unless -e $file;
	my $infile = new IO::File ("$file","r");
	undef $/;
	my $data = <$infile>;
	eval $data;
	$/="\n";
	
	# some correction for merged genes
	my $merged_genes=get_merged_genes(glob('~mh6/tmp/Merged_genes'));
	
	while(my($k,$v)=each %$VAR1){
		if ($$merged_genes{$v}){
			$$VAR1{$k}=$$merged_genes{$v}
		}
		if ($k=~/[a-z]$/) {
			$$VAR1{"$`"}=$v;
		}
	}	
	
	return $VAR1;
}

sub get_merged_genes {
	my ($f)=@_;
	use IO::File;
	my %var;
	die "cannot find file: $f" unless -e $f;
	my $inf = new IO::File ("$f","r");
	while(<$inf>){
		s/\"//g;
		chomp;
		my @a=split(';');
		$var{$a[0]}=$var{$a[1]};
	}
	return \%var;
}



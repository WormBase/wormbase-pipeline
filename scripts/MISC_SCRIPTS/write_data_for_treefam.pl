#/software/bin/perl -w

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Getopt::Long;

#write out data for treefam
#>WBG	public PFXXXXX
#>seq

#eg 
#>aap-1_CAEEL WBGene00000001 PF00017:20..94; PF00017:349..422

my ($version, $domain);
GetOptions ("version:i" => \$version,
			"domain"	=> \$domain);

my @species = qw(elegans briggsae remanei);
my $wormbase = retrieve("/nfs/disk100/wormpub/DATABASES/current_DB/Elegans.store");
$version = $version or $wormbase->version;
my %gene2cgc = reverse ($wormbase->FetchData('cgc_name2gene'));
my %pep2domain;

my %postfix = (	'elegans' 	=> '_CAEEL', 
				'briggsae' 	=> '_CAEBR',
				'remanei'	=> '_CAERE');

my $tmfile = "/tmp/tm_query.def";
&write_tm_query($tmfile);
&get_domains if $domain;

foreach my $species (@species) {
	my $accessor;
	if($species eq 'elegans') {
		$accessor = $wormbase;
	}else {
		$accessor  = Wormbase->new(
						     -organism => $species  );
	}
	
	my %info;
	#my $seq_file  = $accessor->wormpep."/".$wormbase->pepdir_prefix.'pep'.$version;
	my $seq_file = $accessor->basedir."/WORMPEP/".$accessor->pepdir_prefix.'pep'.$version."/".$accessor->pepdir_prefix.'pep'.$version;
	print STDOUT "reading $seq_file\n";
	open(PEP,"<$seq_file") or die "cant open $seq_file : $!\n";
	my ($gene, $cds, $pep, $pep_seq);
	while (<PEP>){
		chomp;
		if(/>(\S+)\s+(\S+)\s+(\S+)/) {
			if($pep_seq){
				unless($info{$gene}->{'seq'} and length($info{$gene}->{'seq'}) < length($pep_seq)){
					$info{$gene}->{'seq'} = $pep_seq;
					$info{$gene}->{'cds'} = $cds;
					$info{$gene}->{'pep'} = $pep;
				}
				undef $pep_seq;
			}
				
			$cds = $1;
			$pep = $2;
			$gene =$3;
		}else {
			$pep_seq .= $_;
		}
	}
	
	my $tf_file   = $wormbase->wormpub."/treefam_$species";
	open(TF,">$tf_file") or die "cant open $tf_file : $!\n";
	foreach my $wbg (keys %info){
		my $public = $gene2cgc{$wbg} ? $gene2cgc{$wbg} : $info{$wbg}->{'cds'};
		print TF ">$wbg\t".$public.$postfix{$species}."\t".$info{$wbg}->{'cds'};
		print TF "\t".join("; ",@{ $pep2domain{ $info{$wbg}->{'pep'} } }) if (defined $pep2domain{$info{$wbg}->{'pep'}}->[0]);
		print TF "\n".$wormbase->format_sequence($info{$wbg}->{'seq'});
		print TF "\n";
	}
}
			
			

sub get_domains {
	print STDOUT "Getting domain info . . ";
	my $query = $wormbase->table_maker_query($wormbase->database('current'),$tmfile);
	while(<$query>) {
		next unless /PFAM/;
		s/\"//g;#"
		my($pep, $domain,$start, $end) = split;
		$pep =~ s/\w+://;
		$domain =~ s/\w+://;
		
		push(@{$pep2domain{$pep}},"$domain:$start..$end");
	}
	print STDOUT "done\n";	
}




sub write_tm_query {
my $file = shift;
open(TM,">$file");

my $def = <<END;
Sortcolumn 1

Colonne 1
Subtitle protein
Width 12
Optional
Visible
Class
Class Protein
From 1

Colonne 2
Subtitle motif
Width 12
Mandatory
Visible
Class
Class Motif
From 1
Tag Motif_homol

Colonne 3
Subtitle method
Width 12
Mandatory
Hidden
Class
Class Method
Right_of 2
Tag  HERE

Colonne 4
Subtitle score
Width 12
Mandatory
Hidden
Float
Right_of 3
Tag  HERE

Colonne 5
Subtitle start
Width 12
Mandatory
Visible
Integer
Right_of 4
Tag  HERE

Colonne 6
Subtitle end
Width 12
Mandatory
Visible
Integer
Right_of 5
Tag  HERE

END

print TM $def;

close TM;
}






















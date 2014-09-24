#!/usr/bin/perl -w
use strict;

use lib '/home/mh6/ebi_home/project/software/packages/ensembl/ensembl/modules';
use lib '/home/mh6/ebi_home/project/software/packages/bioperl/bioperl-live/';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use YAML;
use IO::File;

my $yfile = '/home/mh6/ebi_home/project/wormbase-pipeline/scripts/ENSEMBL/etc/ensembl_lite.t3.conf';
my $targetdir = '/tmp/';
my $worm;
GetOptions( 'yfile=s'  =>\$yfile,
	    'target=s' => \$targetdir,
            'worm=s'   => \$worm,
)||die(@!);

my $c = YAML::LoadFile($yfile);

while (my($k,$config)=each %$c){

   next unless $config->{taxon_id}; # skip the generics header
   if ($worm){
      next unless $k eq $worm;
   }

   my $outf = IO::File->new("$targetdir/$k.ace",'w')||die(@!);
 
   print STDERR "processing ($k) ${\$config->{species}}\n";

   my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $c->{generics}->{host},
        -user   => $c->{generics}->{user},
        -dbname => $config->{dbname},
        -pass   => $c->{generics}->{password},
        -port   => $c->{generics}->{port},
    );

   $config->{species}=~/^(\w).+\s(\w)/;
   my $prefix = uc("$1$2");

   my @genes = @{$db->get_adaptor('Gene')->fetch_all()};
   foreach my $gene(@genes){

		print $outf <<HERE;
Gene : "${\$config->{taxon_id}}:${\$gene->stable_id}"
Public_name "${\$gene->stable_id}"
Species "${\$config->{species}}"
Method Gene
Remark "automatically created" Inferred_automatically "tier_3stubs.pl"

HERE
		foreach my $trans (@{$gene->get_all_Transcripts()}) {
			unless ($trans->translate()){
		        	print STDERR "cannot translate ${\$trans->stable_id}\n";
				next;
			}
		   	print $outf <<HERE;
CDS : "${\$config->{taxon_id}}:${\$trans->stable_id}"
Species "${\$config->{species}}"
Gene "${\$config->{taxon_id}}:${\$gene->stable_id}"
Method curated
Database WormBase TierIII \"${\$trans->stable_id}"
Remark "automatically created" Inferred_automatically "tier_3stubs.pl"
Brief_identification "Non-core species stub" Inferred_automatically "tier_3stubs.pl"
Coding CDS
From_laboratory RW

Protein : "$prefix:${\$trans->stable_id}"
Species "${\$config->{species}}"
Corresponding_CDS "${\$config->{taxon_id}}:${\$trans->stable_id}"
Live

Peptide : "$prefix:${\$trans->stable_id}"
${\$trans->translate()->seq}

HERE
		}
  }
  $outf->close;
}

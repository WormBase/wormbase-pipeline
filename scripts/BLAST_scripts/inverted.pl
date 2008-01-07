#!/usr/local/ensembl/bin/perl -w

use lib '/nfs/farm/Worms/Ensembl/ensembl-pipeline/modules';
use lib '/nfs/farm/Worms/Ensembl/ensembl/modules';
use lib '/nfs/disk100/humpub/modules/PerlModules';
use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;
use Bio::Root::RootI;
use Bio::Seq;
use Wormbase;

my($dbname, $dbhost, $dbuser, $help);

GetOptions(
    "dbname=s"  => \$dbname,
    "dbhost=s"  => \$dbhost,
    "dbuser=s"  => \$dbuser,
    "help"      => \$help,
);

$dbname = "worm_dna" unless $dbname;
$dbhost = "ia64b" unless $dbhost;
$dbuser = "wormro" unless $dbuser;

unless ($dbname && $dbuser && $dbhost) {
    print STDERR "Must specify all DB parameters\n";
    exit 1;
}


my $dbobj = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
							  '-host'   => $dbhost,
							  '-user'   => $dbuser,
							  '-dbname' => $dbname
							 ) or die "Can't connect to DB $dbname on $dbhost as $dbuser";

my $clone_adaptor = $dbobj->get_CloneAdaptor();
my $clones = $clone_adaptor->fetch_all();

my $root_dir = "/acari/work2a/wormpipe/dumps";
my $out_file = "$root_dir/inverted_repeats.ace";
my %acc2clone;
&FetchData('accession2clone',\%acc2clone);

open ( OUT, ">$out_file") or die "cant open $out_file for writing\t$!\n";

foreach my $clone ( @{$clones} ) {
  my $contigs = $clone->get_all_Contigs();
  foreach my $contig ( @{$contigs} ) {
    my $seq = $contig->seq;
    my $name = $clone->id;
    my $length = $contig->length;
    my $clone_name = $acc2clone{"$name"};
    my $seq_file = "$root_dir/inv_seq.seq";

    open (SEQ,">$seq_file") or die "Cant write $seq_file\t$!\n";
    print SEQ ">$name\n";
    print SEQ "$seq\n";
    close SEQ;

    print OUT "\nSequence : \"$clone_name\"\nFeature_data \"$clone_name:inverted\" 1 $length\n";
    print OUT "\nFeature_data : $clone_name:inverted\n";

    open(INV,"/usr/local/pubseq/bin/inverted $seq_file |" );

    my ($perc_ident,$gaps,$L_start, $R_start, $L_end, $R_end);

    # convert to .ace output
    while ( <INV> ) {

      #  Score 147: 89/116 ( 76%) matches, 1 gaps
      #    4369 ccaaacgtgacgttttgcgattttcgcgctaaaattacagtaagtggggtctcgacacgaca-atttttgtgaaatacaaacgggcgtgtgtctttaagaagtactgtagtttaaaa 4484    
      #            |   |   |   |            ||                                 |      |       | | |      |||      || |       |       
      #    4688 ggtctgagctttaaagcgctataagcccagtttttatgacattaaccccagagctgtgctgtgcaaaaacgctttaagctcgtccgcacgtggaaatttcttatgacattaaagttt 4572    


      if ( /Score\s\d+:\s.*(\d+)%.*(\d+)\sgaps/ ) {
	if ( $perc_ident ) {
	  my $loop = $R_start - $L_start -1;
	  $loop .= ", $gaps gaps" if ( $gaps != 0 );
	  print OUT "Feature inverted $L_start $L_end $perc_ident \"loop $loop\"\n";
	  print OUT "Feature inverted $R_start $R_end $perc_ident \"loop $loop\"\n"; 
	  undef $perc_ident;
	  undef $gaps;
	  undef $L_start;
	  undef $R_start;
	  undef $L_end;
	  undef $R_end;
	}

	$perc_ident = $1;
	$gaps = $2;
      } elsif ( /(\d+)\s[acgt-]+\s(\d+)/ ) {
	if ( $L_start ) {
	  $R_start = $1;
	  $R_end = $2;
	} else {
	  $L_start = $1;
	  $L_end = $2;
	}
      } elsif ( /\|/ ) {
	#print STDERR "alignment\n";
      } elsif ( eof ) {
	  my $loop = $R_start - $L_start -1;
	  $loop .= ", $gaps gaps" if ( $gaps != 0 );
	  print OUT "Feature inverted $L_start $L_end $perc_ident \"loop $loop\"\n";
	  print OUT "Feature inverted $R_start $R_end $perc_ident \"loop $loop\"\n";
      }
    }
    if( $L_start and $R_end ) {
      my $loop = $R_start - $L_start -1;
      $loop .= ", $gaps gaps" if ( defined $gaps and $gaps != 0 );
      print OUT "Feature inverted $L_start $L_end $perc_ident \"loop $loop\"\n";
      print OUT "Feature inverted $R_start $R_end $perc_ident \"loop $loop\"\n";	
      undef $perc_ident;
      undef $gaps;
      undef $L_start;
      undef $R_start;
      undef $L_end;
      undef $R_end;
    }
    close INV;
  }
}



=pod

=head1 inverted.pl

  Extracts clone seq from the mysql database and run '/usr/local/pubseq/bin/inverted' to find inverted repeats.

=head1 OPTIONS ( defaults )

  -dbname  ( worm_dna )
  -dbuser  ( wormro )
  -dbhost  ( ia64b )


Ouput written to /acari/work2a/wormpipe/dumps so needs to be run on a machine that can see here eg ecs4

=cut

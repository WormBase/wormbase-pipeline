#!/usr/bin/perl 

use lib $ENV{'CVS_DIR'};
use strict;
use Ace;
use Wormbase;
use Log_files;
use Storable;
use Getopt::Long;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($affy, $agil, $gsc, $all);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "affy"       => \$affy,
	    "agilent"    => \$agil,
	    "gsc"        => \$gsc,
	    "all"        => \$all
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $database_path = $wormbase->autoace;
my $program       = $wormbase->tace;

print "Connecting to database...";

#my $db = Ace->connect('sace://aceserver.cshl.org:2005') || die "Connection failure: ", Ace->error;                      # uncomment to use aceserver.cshl.org - may be slow
my $db = Ace->connect(-path => $database_path,  -program => $program) or $log->log_and_die( print "Connection failure: ". Ace->error);   # local database

($affy = $agil = $gsc = 1) if $all;
&write_table('affy') if $affy;
&write_table('agil') if $agil;
&write_table('gsc')  if $gsc;

$log->mail;

exit(0);


sub write_table {
  my $type  = shift;
  $log->write_to("writing $type\n");
#  if ($type]) {
#    if ($type=~/affymetrix/i) {
#	$type="Affymetrix";
#    }
#    elsif ($ARGV[1]=~/agilent/i) {
#	$type="Agilent";
#    }
#    elsif ($ARGV[1]=~/gsc/i) {
#	$type="GSC at WashU";
#    }
#    else {
#	print "unrecognized type\n";
#	exit;
#    }
#}
    
my %type_remark = ('affy' => "Affymetrix",
		   'agil' => "Agilent",
		   'gsc'  => "GSC at WashU"
		  );



  my $query;
  if ($type) {
    $query="find oligo_set remark=\"*".$type_remark{$type}."*\"";
  } else {
    $query="find oligo_set";
  }
  my $it=$db->fetch_many(-query=>$query);
  my $count=0;
  my ($gene, $seq_name, %gene_hash);
  my (@cds, @transcript, @pseudogene);
  my ($cds_count, $tr_count, $pseudo_count, $no_count, $coding_count, $non_coding_count) = (0,0,0,0,0,0);

  open(OUT,">".$wormbase->autoace."/${type}_oligo_mapping") or $log->log_and_die("cant open ${type}_oligo_mapping :$!\n");
  print OUT "Oligo_set\tWBGeneID\tGene_sequence_name\tGene_type\tMicroarray_type\tRemark\n";
  while (my $obj=$it->next) {

    my @rem=$obj->Remark;
    if (! $type) {
      $type='';
      foreach (@rem) {
	$type=$type_remark{$type};
      }
    }
    $count++;
    if ($count % 1000 == 0) {
      print "$count objects processed\n";
    }
    @cds=();
    @transcript=();
    @pseudogene=();
    if ($obj->Overlaps_CDS) {
      @cds=$obj->Overlaps_CDS;
      $cds_count++;
    }
    if ($obj->Overlaps_transcript) {
      @transcript=$obj->Overlaps_transcript;
      $tr_count++;
    }
    if ($obj->Overlaps_pseudogene) {
      @pseudogene=$obj->Overlaps_pseudogene;
      $pseudo_count++;
    }
    if (! @cds and ! @transcript and ! @pseudogene) {
      if (scalar @rem == 1) {
	print OUT "$obj\tno overlapping gene\t\t\t$type_remark{$type}\t\n";
      } elsif (scalar @rem > 1) {
	print OUT "$obj\tno overlapping gene\t\t\t$type_remark{$type}\t";

	foreach (@rem) {
	  if (/microarray/) {
	    next;
	  }
	  print OUT "$_ ";
	}
	print OUT "\n";
      } else {
	print OUT "$obj\tno overlapping gene\t\t\t\t\n";
      }

      $no_count++;
    } else {
      %gene_hash=();
      foreach (@cds) {
	$gene=$_->Gene;
	$seq_name=$gene->Sequence_name;
	if (! $seq_name) {
	  $seq_name=$gene->Public_name;
	}
	push (@{$gene_hash{$gene}}, $seq_name, "CDS");
      }
      foreach (@transcript) {
	if ($_->Gene) {
	  $gene=$_->Gene;
	  $seq_name=$gene->Sequence_name;
	  if (! $seq_name) {
	    $seq_name=$gene->Public_name;
	  }
	  push (@{$gene_hash{$gene}}, $seq_name, "Non-coding transcript");
	  $non_coding_count++;
	} else {
	  $gene=$_->Corresponding_CDS->Gene;
	  $seq_name=$gene->Sequence_name;
	  if (! $seq_name) {
	    $seq_name=$gene->Public_name;
	  }
	  push (@{$gene_hash{$gene}}, $seq_name, "CDS");
	  $coding_count++;
	}
      }
      foreach (@pseudogene) {
	$gene=$_->Gene;
	$seq_name=$gene->Sequence_name;
	if (! $seq_name) {
	  $seq_name=$gene->Public_name;
	}
	push (@{$gene_hash{$gene}}, $seq_name, "Pseudogene");
      }
      if (scalar keys %gene_hash == 0) {
	print "$obj does not have a corresponding gene!\n";
	next;
      }
      foreach (sort {$a cmp $b} keys %gene_hash) {
	if (! $gene_hash{$_}[0] or ! $gene_hash{$_}[1]) {
	  print "$obj\t$_\t$gene_hash{$_}[0]\t$gene_hash{$_}[1]\n";
	}
	if (scalar @rem == 1) {
	  print OUT "$obj\t$_\t$gene_hash{$_}[0]\t$gene_hash{$_}[1]\t$type_remark{$type}\t\n";
	} elsif (scalar @rem > 1) {
	  print OUT "$obj\t$_\t$gene_hash{$_}[0]\t$gene_hash{$_}[1]\t$type_remark{$type}\t";

	  foreach (@rem) {
	    if (/microarray/) {
	      next;
	    }
	    print OUT "$_ ";
	  }
	  print OUT "\n";
	} else {
	  print OUT "$obj\t$_\t$gene_hash{$_}[0]\t$gene_hash{$_}[1]\t\t\n";
	}
      }
    }
  }

  close OUT;
  print "$count objects written\n";
  print "$cds_count overlap CDS\n";
  print "$tr_count overlap transcript\n";
  print "$pseudo_count overlap pseudogene\n";
  print "$no_count overlap nothing\n";
}

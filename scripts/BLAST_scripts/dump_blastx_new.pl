#!/usr/local/bin/perl5.6.1 -w

# Marc Sohrmann (ms2@sanger.ac.uk)
# with small amendment by Keith Bradnam (krb@sanger.ac.uk)
# dumps wublastx from ensembl mysql (dna) database to an ace file


use strict;
use DBI;
use Getopt::Long;
use DB_File;

my ($version, $map, $worm, $start_clone, $modifed);
my @analysis;
my @clones;
GetOptions ("version=s"      => \$version,
	    "map"            => \$map,
	    "worm"           => \$worm,
	    "restart=s"      => \$start_clone,
	    "analysis=s"     => \@analysis,
	    "clone=s"        => \@clones,
	    "modified"       => \$modifed
	   );

@analysis = split(/,/,join(',',@analysis));
@clones   = split(/,/,join(',',@clones));

print "dumping for @analysis\n" if @analysis;

my $usage = "dump_blastx.pl\n";
$usage .= "-map  map accessions to names, needed for elegans select if DONT want clone names. Default is to output clone names\n";
$usage .= "-worm  dump only worm matches\n";

die "Please enter version no \n" unless $version;

####################################################################
# set some parameters
####################################################################
# mysql database
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "worm_dna";
my $dbpass = "";

# define the organisms we deal with
my @species = qw(fly human slimswissprot slimtrembl worm yeast briggsae);

# define the AceDB protein prefixes for the organisms used
my %org2acedb;

$org2acedb{worm} = "WP:";
$org2acedb{fly} = "GADFLY:";
$org2acedb{ensembl} = "ENSEMBL:";
$org2acedb{yeast} = "SGD:";
$org2acedb{slimswissprot} = "SW:";
$org2acedb{slimTrEMBL} = "TR:";
$org2acedb{briggsae} = "BP:";
$org2acedb{human} = "never_used"; #this should always be done in sub getPrefix. This is here to check for valid analysis selection


my %worm_dna_processIds = ( 'worm'          => '2',
			    'briggae'       => '3',
			    'human'         => '4',
			    'yeast'         => '5',
			    'fly'           => '6',
			    'slimswissprot' => '7',
			    'slimtrembl_1'  => '8',
			    'slimtrembl_2'  => '9',
			  );

my %anal2db; 
foreach ( keys %worm_dna_processIds ){
  my $anal;
  $anal = $_;
  $anal = "slimtrembl" if ($anal =~ /trembl/);
  $anal2db{ "$worm_dna_processIds{$_}" } = $anal;
}

# set an evalue threshold (as -log base 10)
my $e_threshold = 6;

# size of searchspace used for the calculation of the evalue
my $searchspace = 10000000;


# input files
my $helper_files_dir = glob("~wormpipe/Elegans");
my $agp_file = "$helper_files_dir/WS$version.agp";
my $cds_file = "$helper_files_dir/cds$version.gff";
my $cos_file = "$helper_files_dir/cos$version.gff";
my $wormpep_file = glob("~wormpipe/BlastDB/wormpep$version.pep");

unless ( -e $agp_file and -e $wormpep_file and -e $cds_file and -e $cos_file ) {
  print "missing helper file please check these exist : \n$agp_file\n$cds_file\n$cos_file\n$wormpep_file\n";
  die;
}

# create output files
my $dump_dir = "/acari/work2a/wormpipe/dumps/blastx";
my $log = "$dump_dir/blastx_ensembl.log";
open (LOG, ">$log") || die "cannot create log file $log\n";

# create a hash of output file handles keyed off orgs
# this will allow the outputting by org
my %output_ACE;
my $ace;
$dump_dir .= "/updates" if @clones;
#foreach my $org( keys %org2acedb ) { 
#  $ace = "$dump_dir/${org}_blastx_ensembl.ace";
#  open ($output_ACE{$org}, ">$ace") || die "cannot create ace file $ace\n";
#}

# make ACE and LOG un-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

#$old_fh = select(ACE);
#$| = 1;
#select($old_fh);

print LOG "DUMPing cosmid wublastx data from mysql to ace [".&now."]\n";
print LOG "--------------------------------------------------------------------\n\n";

####################################################################
# read the agp file,
# populating %accession2version
####################################################################
my %accession2version;


print LOG "get the accession 2 version mapping from the agp file [".&now."]\n\n";

open (AGP , "$agp_file") || die "cannot read $agp_file";
while (<AGP>) {
  chomp;
  my @ary = split /\t/;
  $ary[5] =~ /^([^\.]+)\.(\d+)/;
  $accession2version{$1} = $2;
}
close AGP;


####################################################################
# read the appropriate wormpep file (same as BlastDB version),
# populating %name2id (maps gene name to protein id)
####################################################################
my %name2id;

print LOG "get the name 2 id mapping for the wormpep proteins [".&now."]\n\n";

open (WP , "$wormpep_file") || die "cannot read $wormpep_file";
while (<WP>) {
  next unless /^>/;
  /^>(\S+)\s+(\S+)/;
  $name2id{$1} = $2;
}
close WP;


####################################################################
# read the file with cds coordinates
####################################################################
my %cds;

print LOG "get chromosomal cds coordinates [".&now."]\n\n";

open (CDS , "$cds_file") || die "cannot read $cds_file";
while (<CDS>) {
  chomp;
  my @ary = split /\t/;
  my $start = $ary[3];
  my $end = $ary[4];
  $ary[8] =~ /\"(\S+)\"/;
  my $name = $1;
  if ($name) {
    #$name =~ tr/a-z/A-Z/;  -  dont know why it was doing this ! !
    $cds{$name} = [$start, $end];
  } else {
    die "cannot process $cds_file";
  }
}
close CDS;


####################################################################
# read the file with cosmid coordinates
####################################################################
my %cos;
print LOG "get chromosomal cosmid coordinates [".&now."]\n\n";

open (COS , "$cos_file") || die "cannot read $cos_file";
while (<COS>) {
  chomp;
  my @ary = split /\t/;
  my $start = $ary[3];
  my $end = $ary[4];
  $ary[8] =~ /\"(\S+)\"/;
  my $name = $1;
  if ($name) {
    $cos{$name} = [$start, $end];
  } else {
    die "cannot process $cos_file";
  }
}
close COS;


####################################################################
# connect to AceDB using TableMaker,
# populating %accession2name (maps embl accession to contig name)
####################################################################

# this could use Common_data.pm - will investigate -  the file acc2clone.data will need to be copied over 
my %accession2name;

unless ($map) {
  print LOG "loading the accession 2 name mappings for the AceDB contigs\n";
  print LOG "using file ~/dumps/accession2clone.list [".&now."]\n\n";

  open(ACCESSION,"/nfs/acari/wormpipe/dumps/accession2clone.list") || die "Could not open accession2clone file\n";
  while (<ACCESSION>) {
    chomp;
    if (/(\S+)\s+(\S+)/) {
      $accession2name{$1} = $2;
      print "Output: \"$1\"\t\"$2\"\n";
    }
  }
  close ACCESSION;
}


# Open DBM to find out database of HID
dbmopen my %ACC2DB, "/acari/work2a/wormpipe/dumps/acc2db.dbm", 0666 or die "cannot open acc2db \n";
open (IPI_LIST, ">/acari/work2a/wormpipe/dumps/ipi_hits_list_x") or die "cant open hitlist\n";



####################################################################
# connect to the Mysql database
####################################################################
print LOG "connect to the mysql database $dbname on $dbhost as $dbuser [".&now."]\n\n";
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
  || die "cannot connect to db, $DBI::errstr";


####################################################################
# get the mapping of analysisId to organism for each database,
# populating %analysis2org
####################################################################
my %analysis2org;

print LOG "get mapping of organism to analysisId for each database [".&now."]:\n";

my $sth = $dbh->prepare ( q{ SELECT analysis_id, db
			     FROM analysis
			     WHERE program = 'wublastx'
			   } );
$sth->execute;

while (my @row = $sth->fetchrow_array) {
  if ($row[1] =~ /(wormpep|gadfly|ensembl|yeast|slimswissprot|slimtrembl|human|brigpep)/) {
    my $org;
    if ($1 =~ /wormpep/) {
      $org = "worm";
    } elsif ($1 =~ /gadfly/) {
      $org = "fly";
    } elsif ($1 =~ /ensembl/) {
      $org = "human";
    } elsif ($1 =~ /yeast/) {
      $org = "yeast";
    } elsif ($1 =~ /slimswissprot/) {
      $org = "slimswissprot";
    } elsif ($1 =~ /slimtrembl/) {
      $org = "slimtrembl";
    } elsif ($1 =~ /human/) {
      $org = "human";
    } elsif ($1 =~ /brigpep/) {
      $org = "briggsae";
    }
    $analysis2org{$row[0]} = $org;
  }
}
foreach my $anal (sort {$a cmp $b} keys %analysis2org) {
  print LOG "\t$anal\t=> $analysis2org{$anal}\n";
}
print LOG "\n";

####################################################################
# prepare sql queries
####################################################################
# contig table
my $sql = "SELECT contig.contig_id, contig.name, contig.length, clone.embl_version FROM contig, clone WHERE contig.contig_id = clone.clone_id";
$sql .= " and contig.contig_id >= ?" if ($start_clone);
my $x = join(" OR contig.contig_id = ",@clones);
$sql .= " AND (contig.contig_id = $x)" if ( @clones );
$sql .= " ORDER BY contig_id";

my $sth_c;
$sth_c = $dbh->prepare ( $sql );


#if ($start_clone) {
#  $sth_c = $dbh->prepare ( q{ SELECT contig.contig_id, contig.name, contig.length, clone.embl_version
#			      FROM contig, clone
#			      WHERE contig.contig_id = clone.clone_id and contig.contig_id >= ? 
#			      ORDER BY contig_id
#			    } );
#} else {
#  $sth_c = $dbh->prepare ( q{ SELECT contig.contig_id, contig.name, contig.length, clone.embl_version
#			      FROM contig, clone
#			      WHERE contig.contig_id = clone.clone_id
#			      ORDER BY contig_id
#			    } );
#}
# this was an attempt to spee it up ! 
#$sql = " SELECT dna_align_feature_id, analysis_id, contig_start, contig_end, hit_name, hit_start, hit_end, score, evalue,contig_strand, cigar_line FROM dna_align_feature WHERE contig_id = ?";
#$x = join("OR analysis_id = ",@analysis);
#$sql .= " AND ( analysis_id = $x )", if ( @analysis );
#$sql .= " ORDER BY contig_start and contig_end";

#my $sth_f = $dbh->prepare ( $sql );


# feature table
my $sth_f = $dbh->prepare ( q{ SELECT dna_align_feature_id,
			       analysis_id,
			       contig_start, contig_end,
			       hit_name, hit_start, hit_end,
			       score, evalue,contig_strand, cigar_line
			       FROM dna_align_feature
			       WHERE contig_id = ?
			       ORDER BY contig_start and contig_end
			     } );


####################################################################
# loop over all the contigs, selecting only the following matches:
#   - evalues treated as -log(10)
#   - do the evaluation for both + and - strand, and every 25bp interval
#   - keep best of fly, human, yeast and slimSwall + within 25% of evalue
#   - if above => keep all better worm, or if none the best within 25% of evalue,
#     else => keep best worm + within 25% of evalue
####################################################################
# loop over all the contigs

print LOG "get all the features (per contig) [".&now."]:\n";

if ($start_clone) {
  $sth_c->execute($start_clone);
} else {
  $sth_c->execute;
}
my $ref = $sth_c->fetchall_arrayref;
foreach my $aref (@$ref) {
  # get sequence info
  my ($internal_id, $id, $length, $version) = @$aref;
  my $fragment_number = (int($length/50))+1;
  # process id's from ensembl names
  if ($id =~ /^(.+)\.[^\.]+\.[^\.]+\.[^\.]+/) {
    $id = $1;
  }
  my $name;
  # map accession 2 name
  unless ($map) {
    if (exists $accession2name{$id}) {
      $name = $accession2name{$id};
      print LOG "\n\tprocessing $name $id $version ($internal_id) of length $length , $fragment_number 50bp fragments [".&now."]\n";
    } else {
      print LOG "\nno mapping of accession $id to clone name -> skip it\n\n";
      next;
    }
  } else {
    $name = $id;
    print LOG "\n\tprocessing $name $version ($internal_id) of length $length , $fragment_number 50bp fragments [".&now."]\n";
  }
  # check that we have the current version
  #   if ($opt_a) {
  if ($version < $accession2version{$id}) {
    print LOG "\nold version ($version) of $id, current one is $accession2version{$id} -> skip it\n\n";
    next;
    #       }
  }
  # %hsp:      stores all mysql column ref's keyed by org, hid and feature_id
  # %accepted: stores all accepted hsp's, keyed by org, hid and feature_id
  # %sum:      stores all hsp's that are part of the best sum statistics, keyed by org, hid and feature_id
  my %accepted;
  my %hsp;
  my %strand;
  my %sum;
  my %min_sum;    
  # feature table query
  $sth_f->execute ($internal_id);
  my $ref = $sth_f->fetchall_arrayref;
  my %count;
  # loop over all returned hsp columns
  # (feature_id , analysis , seq_start, seq_end, hid, hstart, hend, score, evalue, strand, cigar)
  foreach my $aref (@$ref) {
    # check that analysis2org mapping exists
    my $org;
    unless ($org = $analysis2org{$aref->[1]}) {
      next;
    }
    # worm: change gene name to protein id (e.g. AH6.2 to CE01456), and reject self-matches
    if ($org eq "worm") {
      my $gene_name = $aref->[4];
      $gene_name =~ /^(\S+)\.\S+$/;
      my $parent = $1;
      my $match_start = $cos{$name}->[0]+$aref->[2]+1;
      my $match_end = $cos{$name}->[0]+$aref->[3]+1;

      if (($match_start > $cds{$gene_name}->[0] && $match_start < $cds{$gene_name}->[1]) || ($match_end > $cds{$gene_name}->[0] && $match_end < $cds{$gene_name}->[1])) {
	next;
      }

      unless ($aref->[4] = $name2id{$aref->[4]}) {
	die "no mapping of protein name ".$aref->[4]." to id";
      }
    }
    # store the hsp's that are part of the best sum statistics (strand independent)
    unless (exists $min_sum{$org}->{$aref->[4]}) {
      $min_sum{$org}->{$aref->[4]} = 1;
    }
    if ($aref->[8] < $min_sum{$org}->{$aref->[4]}) {
      $sum{$org}->{$aref->[4]} = ();
      $sum{$org}->{$aref->[4]}->{$aref->[0]} = $aref->[2];
      $min_sum{$org}->{$aref->[4]} = $aref->[8];
    } 
    elsif ($aref->[8] == $min_sum{$org}->{$aref->[4]}) {
      $sum{$org}->{$aref->[4]}->{$aref->[0]} = $aref->[2];
    }

    # recalculate the evalue, not using sum statistics
    $aref->[8] = ($searchspace*$length)*2**(-($aref->[7]));

    $count{$aref->[1]}++;

    # set evalue to -log(10)
    if ($aref->[8] != 0) {
      $aref->[8] = sprintf ("%.3f", -(log($aref->[8]))/log(10));
    }
    else {
      $aref->[8] = 999.999;
    }

    # populate %hsp
    $hsp{$org}->{$aref->[4]}->{$aref->[0]} = $aref;

    # populate %plus or %minus, splitting hsp' into intervals
    my $first_fragment = (int($aref->[2]/50));
    my $last_fragment = (int($aref->[3]/50));

    for (my $n = $first_fragment ; $n <= $last_fragment ; $n++) {
      if ($aref->[9] > 0) {
	push (@{$strand{plus}->[$n]}, $aref);
      } elsif ($aref->[9] < 0) {
	push (@{$strand{minus}->[$n]}, $aref);
      } else {
	die "wrong strand assignment ".$aref->[9]." for feature ".$aref->[0];
      }
    }
  }


  foreach (sort {$a <=> $b} keys %count) {
    print LOG "\t\t\t$analysis2org{$_} analysis $_ ($count{$_} HSP's)\n";
  }
  # make the feature table querie again, this time per 10kb fragments
  # (this speeds up the analysis of every 50 bp interval quite a lot)
  my @strands = qw(plus minus);
  foreach my $strand (@strands) {
    print LOG "\t\t$strand strand\n\t\t\t";
    for (my $n = 0 ; $n < $fragment_number ; $n++) {
      print LOG "." if $n % 100 == 0;

      next unless $strand{$strand}->[$n];

      my $worm_switch = 0;
      my @worm_tmp;
      my %accept_per_org;

      # initialise hash to zero
      $accept_per_org{"worm"}  = 0;
      $accept_per_org{"fly"}   = 0;
      $accept_per_org{"human"} = 0;
      $accept_per_org{"yeast"} = 0;
      $accept_per_org{"briggsae"} = 0;
      $accept_per_org{"slimswissprot"} = 0;
      $accept_per_org{"slimtrembl"} = 0;

      my $max_ev = 0;
      # loop over all hsp columns, ordered by increasing evalue,
      # parsing the results
    LOOP:foreach my $aref (sort {$b->[8] <=> $a->[8]} @{$strand{$strand}->[$n]}) {
	my ($fid , $analysis, $start, $end, $hid, $hstart, $hend, $score, $ev, $strand, $cigar) = @$aref;

	next LOOP unless $ev > $e_threshold;
	my $org = $analysis2org{$analysis};
	# worm
	if ($org eq "worm") {
	  # if there are no better matches from others
	  if ($max_ev == 0) {
	    if ($accept_per_org{$org} >= 5) {

	      next LOOP;
	    }
	    else {
	      push (@worm_tmp , [$fid, $hid, $ev]);
	      $worm_switch = 1;
	      $accept_per_org{$org}++;

	    }
	  }
	
	  # if the best match is from others , and worm is within 75%
	  elsif (($max_ev != 0) && ($ev > 0.75*$max_ev) && ($worm_switch == 0)) {
	    if ($accept_per_org{$org} >= 2) {

	      next LOOP;
	    }
	    else {
	      $accepted{$org}->{$hid}->{$fid} = 1;
	      $accept_per_org{$org}++;

	    }
	  }
	}
	# others
	else {
	  if ($ev > $max_ev) {
	    $max_ev = $ev;
	    if ($accept_per_org{$org} >= 2) {

	      next LOOP;
	    } else {
	      $accepted{$org}->{$hid}->{$fid} = 1;
	      $accept_per_org{$org}++;

	    }
	  } 
	  elsif ($ev > 0.75*$max_ev) {
	    if ($accept_per_org{$org} >= 2) {
	      next LOOP;
	    }
	    else {
	      $accepted{$org}->{$hid}->{$fid} = 1;
	      $accept_per_org{$org}++;

	    }
	  }
	  else {
	    last LOOP;
	  }
	}
      }

      # go through the worm again:
      if ($max_ev == 0) {
	# there were no matches from others at all, keep within 75% of best worm match
	foreach my $aref (@worm_tmp) {

	  if ( (defined $worm_tmp[1]->[2]) and $aref->[2] > 0.75*$worm_tmp[1]->[2]) {
	    $accepted{worm}->{$aref->[1]}->{$aref->[0]} = 1 if (defined $aref->[1]);

	  }
	}
      }
      else {
	# keep all
	foreach my $aref (@worm_tmp) {
	  $accepted{worm}->{$aref->[1]}->{$aref->[0]} = 1;
	}
      }
    }
    print LOG "\n";
  }

  # loop over all the matches, printing the accepted ones to an ACE file
  print LOG "\t\toutput\n";

  # if specific analysis selected just have them in the species @ to dump
  # put in hash then @ coz we dont want slimtrembl there twice.
  if (@analysis) {
    my %species;
    foreach ( @analysis ) {
      $species{$_} = $anal2db{$_};
    }
    @species = values %species;
  }
  foreach my $org (@species) {
    my $output;

    # if filehandle no already active create & store it
    unless( $output_ACE{$org} ) {
      my $ace = "$dump_dir/${org}_blastx.ace";
      open($output_ACE{$org}, ">$ace") ;
    }
    $output = $output_ACE{$org};

    # write the sequence info to the text file
    # if this is just updating specific clones then we need to -D the exisiting data
    if( @clones ) {
      print $output "\n\n";
      print $output "Sequence : \"$name\"\n";
      print $output "-D Homol_data";
      print $output "\n";
    }

    # output the data
    print $output "\n\n";
    print $output "Sequence : \"$name\"\n";
    print $output "Homol_data \"$name:wublastx_$org\" 1 $length\n";
    print $output "\n";
    print $output "Homol_data : \"$name:wublastx_$org\"\n";

    my $count_target = 0;
    my $count_accepted = 0;
    my $count_hsp = 0;
  TARGET:foreach my $target (keys %{$hsp{$org}}) {
      $count_target++;
      my %written;
      unless (exists $accepted{$org}->{$target}) {
	next TARGET;
      }
      $count_accepted++;
    FID:foreach my $fid (sort { $hsp{$org}->{$target}->{$a}->[2] <=> $hsp{$org}->{$target}->{$b}->[2] } keys %{$hsp{$org}->{$target}}) {
	unless (exists $accepted{$org}->{$target}->{$fid}) {
	  next FID;
	}
	if (exists $written{$fid}) {
	  next FID;
	}
	my @write;
	my ($feature_id , $analysis, $start, $end, $hid, $hstart, $hend, $score, $e, $hsp_strand, $cigar) = @{$hsp{$org}->{$target}->{$fid}};
	if (exists $sum{$org}->{$target}->{$fid}) {
	  foreach my $fid (keys %{$sum{$org}->{$target}}) {
	    push (@write, $hsp{$org}->{$target}->{$fid});
	    $count_hsp++;
	    $written{$fid} = 1;
	  }
	} else {
	  push (@write, $hsp{$org}->{$target}->{$fid});
	  $count_hsp++;
	  $written{$fid} = 1;
	}
	foreach my $aref (sort { $a->[2] <=> $b->[2] } @write) {
	  my ($fid , $analysis, $start, $end, $hid, $hstart, $hend, $score, $e, $hsp_strand, $cigar) = @$aref;
	  # check strand
	  if ($hsp_strand < 0) {
	    my $tmp = $end;
	    $end = $start;
	    $start = $tmp;
	  }
	  # write the HSP info to the ACE file
	  my @cigars = split (/:/, $cigar);
	  my $prefix;
	  if ("$org" eq "human") {
	    $prefix = &getPrefix($hid);
	  } else {
	    $prefix = ${org2acedb{$org};
	  }
	  ;
	}
	if (@cigars == 1) {
	  print $output "Pep_homol\t\"$prefix$hid\" \"wublastx_$org\" ";
	  print $output "$e $start $end $hstart $hend\n";
	} else {
	  print $output "Pep_homol\t\"$prefix$hid\" \"wublastx_$org\" ";
	  print $output "$e $start $end $hstart $hend ";
	  print $output "align $start $hstart\n";
	  shift @cigars;
	  foreach my $string (@cigars) {
	    my ($coor, $hcoor) = split (/,/, $string);
	    print $output "Pep_homol\t\"$prefix$hid\" \"wublastx_$org\" ";
	    print $output "$e $start $end $hstart $hend ";
	    print $output "align $coor $hcoor\n";
	  }
	}
      } 
    }
  }
  print LOG "\t\t\t$org -> $count_target targets, $count_accepted accepted ($count_hsp HSP's)\n";
}
}

  # close the mysql database and the filehandles
  $sth->finish;
$sth_c->finish;
$sth_f->finish;
$dbh->disconnect;

close IPI_LIST;

foreach ( keys %output_ACE ) {
  close $output_ACE{$_};
}

print LOG "cat'ing individual files in to ${dbname}_blastx.ace\n";
system( "cat $dump_dir/*_blastx.ace > /acari/work2a/wormpipe/dumps/${dbname}_blastx.ace");

close LOG;
print "\nEnd of dump \n";



exit 0;
sub getPrefix 
  {
    my $name = shift;
    if ( $ACC2DB{$name} ) {
      print IPI_LIST "$name\n";
      my $n = $ACC2DB{$name}.":";
      return $n; 
    }
    # NOTE this is only the prefix - not the method (it will look like wublastp_ipi_human ENSEMBL:ENS00342342 etc)
    if ( $name =~ /ENS\w+/ ) {
      return $org2acedb{'ensembl'};
    } else {
      if (length $name > 6 ) {
	return $org2acedb{'slimswissprot'};
      } else {
	return $org2acedb{'slimtrembl'};
      }
    }
  }

# subroutine to get the current time...
sub now {
  return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
		  sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

#!/usr/local/bin/perl5.6.0 -w

# Marc Sohrmann (ms2@sanger.ac.uk)
# with small amendment by Keith Bradnam (krb@sanger.ac.uk)
# dumps wublastx from ensembl mysql (dna) database to an ace file
# must currently be run on wormsrv2, since it needs access to autoace

use strict;
use DBI;
use Getopt::Std;
use vars qw($opt_w $opt_a $opt_g $opt_c $opt_m $opt_o);

getopts ("w:a:g:c:mo");

my $usage = "dump_blastx.pl\n";
$usage .= "-w [wormpep file, same version as database used in blast searches]\n";
$usage .= "-a [agp file], checks for up-to-date sequence version number before dumping\n";
$usage .= "-g [gff file of cds coordinates relative to chromosomes], needed to exclude self matches\n";
$usage .= "-c [gff file of cosmid coordinates relative to chromosomes], needed to exclude self matches\n";
$usage .= "-m map accessions to names, needed for elegans\n";
$usage .= "-o dump only worm matches\n";

unless ($opt_w) {
    die "$usage";
}

####################################################################
# set some parameters
####################################################################
# mysql database
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "worm01";
my $dbpass = "";

# AceDB database
my $ace_dir = "/wormsrv2/current_DB";
my $wquery_dir = "/wormsrv2/autoace/wquery";
my $tace = "/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4/tace";

# define the organisms we deal with
my @species = qw(fly human slimSwissProt slimTrEMBL worm yeast);

# define the AceDB protein prefixes for the organisms used
my %org2acedb;

$org2acedb{worm} = "WP:";
$org2acedb{fly} = "GADFLY:";
$org2acedb{human} = "ENSEMBL:";
$org2acedb{yeast} = "SGD:";
$org2acedb{slimSwissProt} = "SW:";
$org2acedb{slimTrEMBL} = "TR:";

# set an evalue threshold (as -log base 10)
my $e_threshold = 6;

# size of searchspace used for the calculation of the evalue
my $searchspace = 10000000;

# subroutine to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
open (LOG, ">blastx_dump_resume.log") || die "cannot create log file";
open (ACE, ">blastx_dump_resume.ace") || die "cannot create ace file";
#*LOG = *STDOUT;

# make ACE and LOG un-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

$old_fh = select(ACE);
$| = 1;
select($old_fh);

print LOG "DUMPing cosmid wublastx data from mysql to ace [".&now."]\n";
print LOG "--------------------------------------------------------------------\n\n";

####################################################################
# read the agp file,
# populating %accession2version
####################################################################
my %accession2version;

if ($opt_a) {
    print LOG "get the accession 2 version mapping from the agp file [".&now."]\n\n";

    open (AGP , "$opt_a") || die "cannot read $opt_a";
    while (<AGP>) {
        chomp;
        my @ary = split /\t/;
        $ary[5] =~ /^([^\.]+)\.(\d+)/;
        $accession2version{$1} = $2;
    }
    close AGP;
}

####################################################################
# read the appropriate wormpep file (same as BlastDB version),
# populating %name2id (maps gene name to protein id)
####################################################################
my %name2id;

print LOG "get the name 2 id mapping for the wormpep proteins [".&now."]\n\n";

open (WP , "$opt_w") || die "cannot read $opt_w";
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

if ($opt_g) {
    print LOG "get chromosomal cds coordinates [".&now."]\n\n";

    open (CDS , "$opt_g") || die "cannot read $opt_g";
    while (<CDS>) {
        chomp;
        my @ary = split /\t/;
        my $start = $ary[3];
        my $end = $ary[4];
        $ary[8] =~ /\"(\S+)\"/;
        my $name = $1;
        if ($name) {
            $name =~ tr/a-z/A-Z/;
            $cds{$name} = [$start, $end];
        }
        else {
            die "cannot process $opt_g";
        }
    }
    close CDS;
}

####################################################################
# read the file with cosmid coordinates
####################################################################
my %cos;

if ($opt_c) {
    print LOG "get chromosomal cosmid coordinates [".&now."]\n\n";

    open (COS , "$opt_c") || die "cannot read $opt_c";
    while (<COS>) {
        chomp;
        my @ary = split /\t/;
        my $start = $ary[3];
        my $end = $ary[4];
        $ary[8] =~ /\"(\S+)\"/;
        my $name = $1;
        if ($name) {
            $cos{$name} = [$start, $end];
        }
        else {
            die "cannot process $opt_c";
        }
    }
    close COS;
}

####################################################################
# connect to AceDB using TableMaker,
# populating %accession2name (maps embl accession to contig name)
####################################################################
my %accession2name;

if ($opt_m) {
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

my $sth = $dbh->prepare ( q{ SELECT analysisId, db
                               FROM analysisprocess
                              WHERE program = 'wublastx'
			   } );
$sth->execute;

while (my @row = $sth->fetchrow_array) {
    if ($row[1] =~ /(wormpep|gadfly|ensembl|yeast|slimSwissProt|slimTrEMBL)/) {
        my $org;
        if ($1 eq "wormpep"){$org = "worm";}
        elsif ($1 eq "gadfly"){$org = "fly";}
        elsif ($1 eq "ensembl"){$org = "human";}
        elsif ($1 eq "yeast"){$org = "yeast";}
        elsif ($1 eq "slimSwissProt"){$org = "slimSwissProt";}
        elsif ($1 eq "slimTrEMBL"){$org = "slimTrEMBL";}
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
my $sth_c = $dbh->prepare ( q{ SELECT contig.internal_id, contig.id, contig.length, clone.embl_version
                                 FROM contig, clone
                                WHERE contig.internal_id = clone.internal_id and 
contig.internal_id = 2575 
                             ORDER BY internal_id
                             } );

# feature table
my $sth_f = $dbh->prepare ( q{ SELECT id,
                                      analysis,
                                      seq_start, seq_end,
                                      hid, hstart, hend,  
                                      score, evalue, strand, cigar
                                 FROM feature
                                WHERE contig = ?
                             ORDER BY seq_start and seq_end
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

$sth_c->execute;
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
    if ($opt_m) {
        if (exists $accession2name{$id}) {
            $name = $accession2name{$id};
            print LOG "\n\tprocessing $name $id $version ($internal_id) of length $length , $fragment_number 50bp fragments [".&now."]\n";
        }
        else {
	    print LOG "\nno mapping of accession $id to clone name -> skip it\n\n";
            next;
        }
    }
    else {
        $name = $id;
        print LOG "\n\tprocessing $name $version ($internal_id) of length $length , $fragment_number 50bp fragments [".&now."]\n";
    }
    # check that we have the current version
    if ($opt_a) {
        if ($version < $accession2version{$id}) {
            print LOG "\nold version ($version) of $id, current one is $accession2version{$id} -> skip it\n\n";
            next;
        }
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
            if ($opt_g && $opt_c) {
                my $gene_name = $aref->[4];
                $gene_name =~ /^(\S+)\.\S+$/;
                my $parent = $1;
                my $match_start = $cos{$name}->[0]+$aref->[2]+1;
                my $match_end = $cos{$name}->[0]+$aref->[3]+1;

#                print $aref->[0]." $gene_name $parent  ".$cos{$parent}->[0]." $match_start  $match_end \n";

                if (($match_start > $cds{$gene_name}->[0] && $match_start < $cds{$gene_name}->[1]) || ($match_end > $cds{$gene_name}->[0] && $match_end < $cds{$gene_name}->[1])) {

#                    print LOG "reject self match: ".$aref->[0]." $gene_name $parent\n";

                    next;
	        }
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

#        foreach my $org (keys %sum) {
#            foreach my $hid (keys %{$sum{$org}}) {
#                foreach my $fid (keys %{$sum{$org}->{$hid}}) {
#                    print "$org  $hid  $fid\n";
#		}
#	    }
#	}

        # recalculate the evalue, not using sum statistics
        $aref->[8] = ($searchspace*$length)*2**(-($aref->[7]));
        #
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
   	    }
            elsif ($aref->[9] < 0) {
                push (@{$strand{minus}->[$n]}, $aref);
	    }
	    else {
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
#            print LOG "$n of $fragment_number\n";
            next unless $strand{$strand}->[$n];
#            print LOG "\t$n exists\n";
            my $worm_switch = 0;
            my @worm_tmp;
            my %accept_per_org;
            my $max_ev = 0;
            # loop over all hsp columns, ordered by increasing evalue,
            # parsing the results
            LOOP:foreach my $aref (sort {$b->[8] <=> $a->[8]} @{$strand{$strand}->[$n]}) {
                my ($fid , $analysis, $start, $end, $hid, $hstart, $hend, $score, $ev, $strand, $cigar) = @$aref;

#                print LOG "N: $n I: $i  $hid  $analysis  $analysis2org{$analysis}  $ev\n";

                next LOOP unless $ev > $e_threshold;
                my $org = $analysis2org{$analysis};
                # worm
                if ($org eq "worm") {
                    # if there are no better matches from others
                    if ($max_ev == 0) {
                        if ($accept_per_org{$org} >= 5) {

#                            print LOG "\t $i failed $hid $analysis $org} $accept_per_org{$org}\n";

                            next LOOP;
 			}
                        else {
                            push (@worm_tmp , [$fid, $hid, $ev]);
                            $worm_switch = 1;
                            $accept_per_org{$org}++;

#                            print LOG "\t $i accept as best $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

			}
		    }
                    # if the best match is from others , and worm is within 75%
                    elsif (($max_ev != 0) && ($ev > 0.75*$max_ev) && ($worm_switch == 0)) {
                        if ($accept_per_org{$org} >= 2) {

#                            print LOG "\t $i failed $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

                            next LOOP;
			}
                        else {
                            $accepted{$org}->{$hid}->{$fid} = 1;
                            $accept_per_org{$org}++;

#                            print LOG "\t $i accept $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

			}
                    }
		}
                # others
                else {
                    if ($ev > $max_ev) {
                        $max_ev = $ev;
                        if ($accept_per_org{$org} >= 2) {

#                            print LOG "\t $i failed $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

                            next LOOP;
			}
                        else {
                            $accepted{$org}->{$hid}->{$fid} = 1;
                            $accept_per_org{$org}++;

#                            print LOG "\t $i accept $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

			}
	            }
          	    elsif ($ev > 0.75*$max_ev) {
                        if ($accept_per_org{$org} >= 2) {

#                            print LOG "\t $i failed $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

                            next LOOP;
			}
                        else {
                            $accepted{$org}->{$hid}->{$fid} = 1;
                            $accept_per_org{$org}++;

#                            print LOG "\t $i accept $hid count $ev  $analysis $org} $accept_per_org{$org}\n";

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
                    if ($aref->[2] > 0.75*$worm_tmp[1]->[2]) {
                        $accepted{worm}->{$aref->[1]}->{$aref->[0]} = 1;

#                        print LOG "\t ACCEPT finally ".$aref->[0]." ".$aref->[1]."\n"; 

		    }
		}
	    }
            else {
                # keep all
                foreach my $aref (@worm_tmp) {
                    $accepted{worm}->{$aref->[1]}->{$aref->[0]} = 1;

#                    print LOG "\t ACCEPT finally all ".$aref->[0]." ".$aref->[1]."\n"; 

 		}
	    }
	}
        print LOG "\n";
    }
    # loop over all the matches, printing the accepted ones to an ACE file
    print LOG "\t\toutput\n";
    if ($opt_o) {
        @species = qw(worm);
    }
    foreach my $org (@species) {
        # write the sequence info to the text file
        print ACE "\n\n";
        print ACE "Sequence : \"$name\"\n";
        print ACE "Homol_data \"$name:wublastx_$org\" 1 $length\n";
        print ACE "\n";
        print ACE "Homol_data : \"$name:wublastx_$org\"\n";
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
		}
                else {
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
                    if (@cigars == 1) {
                        print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastx_$org\" ";
                        print ACE "$e $start $end $hstart $hend\n";
                    }
                    else {
                        print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastx_$org\" ";
                        print ACE "$e $start $end $hstart $hend ";
                        print ACE "align $start $hstart\n";
                        shift @cigars;
                        foreach my $string (@cigars) {
                            my ($coor, $hcoor) = split (/,/, $string);
                            print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastx_$org\" ";
                            print ACE "$e $start $end $hstart $hend ";
                            print ACE "align $coor $hcoor\n";
			}
		    }
		}

   	    }
	}
        print LOG "\t\t\t$org -> $count_target targets, $count_accepted accepted ($count_hsp HSP's)\n";
    }
}
print LOG "End of BlastX dump\n";
print "\n"End of BlastX dump\n";

# close the mysql database and the filehandles
$sth->finish;
$sth_c->finish;
$sth_f->finish;
$dbh->disconnect;

close LOG;
close ACE;



exit 0;


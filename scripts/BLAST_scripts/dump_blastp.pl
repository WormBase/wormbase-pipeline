#!/usr/local/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)
# dumps wublastp from ensembl mysql (protein) database to an ace file

use strict;
use DBI;
use Getopt::Std;
use vars qw($opt_w $opt_s);

getopts ("w:s");

my $usage = "dump_blastp.pl\n";
$usage .= "-w [wormpep file, same version as used for blast searches]\n";
$usage .= "-s only print HSP's that are part of the sum statistics\n";

unless ($opt_w) {
    die "$usage";
}

# mysql database parameters
my $dbhost = "ecs1f";
#my $dbuser = "wormadmin";
my $dbuser = "wormro";
my $dbname = "wormprot";
my $dbpass = "";

# define the organisms we deal with (except slimSwissProt and slimTrEMBL)
my @species = qw(worm fly human yeast);

# set a E-value threshold (used in the sql query)
my $e_threshold = 1.e-6;

# set max. of target id's to be dumped per query sequence
my $target_threshold = 10;
my $others_target_threshold = 1;

# to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
open (LOG, ">blastp_dump.log") || die "cannot create log file";
open (ACE, ">blastp_dump.ace") || die "cannot create ace file";
#*ACE = *STDOUT;

# make the LOG filehandle line-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

my $old_fh = select(ACE);
$| = 1;
select($old_fh);


print LOG "DUMPing protein wublastp data from mysql to ace [".&now."]\n";
print LOG "---------------------------------------------------------------\n\n";
print LOG "parameters:\n";
print LOG "\tE-value threshold = $e_threshold\n";
print LOG "\tmax. number of targets per query = $target_threshold\n";
print LOG "\tmax. number of targets from SwissProt and TrEmbl = $others_target_threshold\n";
if ($opt_s) {
    print LOG "\tonly HSP's that are part of the sum statistics were accepted\n\n";
}

# read the appropriate wormpep file (same as BlastDB version),
# and keep the id2name mapping in a hash
my %id2name;
my %name2id;
open (WP , "$opt_w") || die "cannot read $opt_w";
while (<WP>) {
    next unless /^>/;
    /^>(\S+)\s+(\S+)/;
    $name2id{$1} = $2;
    push (@{$id2name{$2}} , $1);
}

# connect to the mysql database
print LOG "connect to the mysql database $dbname on $dbhost as $dbuser [".&now."]\n\n";
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || die "cannot connect to db, $DBI::errstr";

# get the mapping of org 2 analysis
print LOG "get mapping of database organism to analysis id [".&now."]:\n";
my $sth = $dbh->prepare ( q{ SELECT analysisId, db
                               FROM analysisprocess
                              WHERE program = 'wublastp'
			   } );
$sth->execute;
my %org2analysis;
while (my @row = $sth->fetchrow_array) {
   if ($row[1] =~ /(wormpep|gadfly|ensembl|yeast|slimswissprot|slimtrembl)/) {
        my $org;
        if ($1 =~ /wormpep/){$org = "worm";}
        elsif ($1 =~ /gadfly/){$org = "fly";}
        elsif ($1 =~ /ensembl/){$org = "human";}
        elsif ($1 =~ /yeast/){$org = "yeast";}
        elsif ($1 =~ /slimswissprot/){$org = "slimSwissProt";}
        elsif ($1 =~ /slimtrembl/){$org = "slimTrEMBL";}
        # assume that there is only one analysis per organism
        $org2analysis{$org} = $row[0];
    }
}
foreach my $org (sort {$a cmp $b} keys %org2analysis) {
    print LOG "\t$org\t=> $org2analysis{$org}\n";
}
print LOG "\n";

# set the mapping of org 2 acedb abbreviations
my %org2acedb;
$org2acedb{worm}          = "WP:";
$org2acedb{fly}           = "GADFLY:";
$org2acedb{human}         = "ENSEMBL:";
$org2acedb{yeast}         = "SGD:";
$org2acedb{slimSwissProt} = "SW:";
$org2acedb{slimTrEMBL}    = "TR:";

# prepare sql queries:
# protein table
my $sth_p = $dbh->prepare ( q{ SELECT proteinId, length
                                 FROM protein
                             ORDER BY proteinId
                             } );

# protein_feature table
my $sth_f = $dbh->prepare ( q{ SELECT protein_featureId,
                                      start, end,
                                      hId, hstart, hend,  
                                      strand, evalue, cigar
                                 FROM protein_feature
                                WHERE proteinId = ? and analysis = ? and evalue < ?
                             ORDER BY hId
	  	  	     } );

# loop over all the proteins
print LOG "get all the features (per protein) [".&now."]:\n";
my %written;
# get all the proteins plus their length from the protein table
$sth_p->execute;
my $ref = $sth_p->fetchall_arrayref;
foreach my $aref (@$ref) {
    # get id and length of sequence (id = CE.....)
    my ($id, $length) = @$aref;
    # dump only the currently active entries
   unless (exists $id2name{$id}) {
        print LOG "\n$id is inactive, skip it...\n";
        next;
    }
    # write the sequence info to the ACE file
    print ACE "\n//\n";
    print ACE "Protein : \"WP:$id\"\n";
    print LOG "\n\tprocessing $id [".&now."]\n";
    my %hit;
    my @reverse_lines;
    my $min_cutoff_ev;
    my %processed;
    # loop over all the org's
    foreach my $org (@species) {
        unless (exists $org2analysis{$org}) {
            print LOG "\t\tno mapping of org to analysis for $org\n";
            next;
	}
        else {
            print LOG "\t\t$org\n";
	}
        # make the sql query on the protein_feature table
        $sth_f->execute ($id, $org2analysis{$org}, $e_threshold);
        my $ref = $sth_f->fetchall_arrayref;

#this is what is returned
#   0  131770976
#   1  1
#   2  407
#   3  'CE00003'
#   4  1
#   5  407
#   6  1
#   7  '5.600e-220'
#   8  '1,1,0'
#1  ARRAY(0x140ffb5d8)
#   0  131771615
#   1  1
#   2  407
#   3  'CE00075'
#   4  4496
#   5  4902
#   6  1
#   7  '4.300e-218'
#   8  '1,4496,0'

        # we would like to sort the columns:
        #   order the target sequences by their minimal E-value,
        #   and the HSP's per target sequence by their E-value (increasing)
        # (it would be better to do this in the sql query, but since the
        #  evalue is stored as string, it cannot be done that way)
        my %tmp;
        my %ev;
        my $min_ev = 0;
        LOOP:foreach my $aref (@$ref) {
            # get rid of self-matches (including splice variants)
            if ($org eq "worm") {
	      if ("$id" eq "$aref->[3]" ){ next;}
                my $gene_name = $aref->[3];  #CE00003
                my $trunc_target;
                if ($aref->[3] =~ /([^\.]+\.\d+)[a-zA-Z]$/) {
                    $trunc_target = $1;
                }
                else {
                    $trunc_target = $aref->[3];  #CE00003
	        }
                if (exists $id2name{$id}) {
                    foreach my $name (@{$id2name{$id}}) {
                        my $trunc_name;
                        if ($name =~ /([^\.]+\.\d+)[a-zA-Z]$/) {
                            $trunc_name = $1;
	                }
                        else {
                            $trunc_name = $name;
    	                }
                        if ($trunc_name eq $trunc_target) {
                            next LOOP;
                        }
	            }
		}
#i commented this bit out to dump WS88 as it was causing a crash
# i think this was because I'd used the wrong wormpep file that had headers
# CE00001 AH6.1a rather than the other way round
               # convert target name 2 id
		unless ($aref->[3] = $name2id{$aref->[3]}) {
		  die "problem renaming  name to id (fid ".$aref->[0].")";
		}
                # get rid of duplicated matches due to several gene names
                # mapping to the same protein id
                if (scalar(@{$id2name{$aref->[3]}}) > 1) {
		  unless ($gene_name eq $id2name{$aref->[3]}->[0]) {
		      next LOOP;
		    }
		}
######  end of bit I commented #######################

	      }
            # take -log(10) of evalues
            if ($aref->[7] != 0) {
	      $aref->[7] = sprintf ("%.2f", -(log($aref->[7]))/log(10));
	    }
            else {
	      $aref->[7] = 999.99;
	    }
            # store the lowest evalue
            if ($aref->[7] > $min_ev) {
	      $min_ev = $aref->[7];
	    }
            # store the lowest evalue per target sequence as hash value
            if (exists $ev{$aref->[3]}) {
	      if ($aref->[7] > $ev{$aref->[3]}) {
                    $ev{$aref->[3]} = $aref->[7];
  	        }
            }
            else {
                $ev{$aref->[3]} = $aref->[7];
      	    }
            # push the rows onto an array, stored as hash value
            # keyed by the target name
            push (@{$tmp{$aref->[3]}}, $aref);
	}
        # calculate evalue cut-off
        # (currently set to 66% of -log(10) of evalue)
        my $cutoff_ev = 0.66*$min_ev;
        if ($min_cutoff_ev && $min_ev != 0) {
            if ($cutoff_ev < $min_cutoff_ev) {
                $min_cutoff_ev = $cutoff_ev;
	    }
	}
        elsif ($min_ev != 0) {
            $min_cutoff_ev = $cutoff_ev;
	}
        # loop over all the targets
        my $target_count = 0;
        my $xref_line;
        foreach my $target (sort {$ev{$b} <=> $ev{$a}} keys %ev) {
            # check if this match target passes the cutoff_evalue
            unless ($ev{$target} > $cutoff_ev) {
                last;
            }            
            # check if we have dumped the max. number of targets
            if (++$target_count > $target_threshold) {
                last;
	    }
            # keep track of the worm query-target pairs already processed,
            # in order to prevent writing duplicates
            if ($org eq "worm") {
                my $key = "$id"."_"."$target";
                my $reverse_key = "$target"."_"."$id";
                if (exists $written{$reverse_key}) {
                    next;
                }                
                else {
                    $written{$key} = 1;
		}
	    }
            #
            $xref_line = "Protein : \"${org2acedb{$org}}$target\"\n";
            # loop over all the HSP's (accepting only the hsp that are part of the sum stats if opt_s)
            foreach my $aref (sort {$b->[7] <=> $a->[7] or $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @{$tmp{$target}}) {
                # $target and $hid are the same (or we are in trouble...)
                my ($fid , $start, $end, $hid, $hstart, $hend, $strand, $e, $cigar) = @$aref;
                # check if the hsp is part of the sum stats
                if ($opt_s) {
                    unless ($e == $ev{$hid}) {
                        last;
	   	    }
		}
                # mask the regions of the query sequence that have been matched
	        for (my $i = $start ; $i <= $end ; $i++) {
                    $hit{$i} = 1;
		}
                # check the strand, and set start/end accordingly
                # (not really required for blastp, since all matches are "forward")
	        if ($strand > 0) {}
                elsif ($strand < 0) {
                    my $tmp = $end;
                    $end = $start;
                    $start = $tmp;
  	        }
                else {
                    print LOG "incorrect strand assignment $strand for $id, feature $fid -> skip it\n";
                    next;
   	        }
                # write the HSP info to the ACE file
                my @cigars = split (/:/, $cigar);
                if (@cigars == 1) {
                    print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                    print ACE "$e $start $end $hstart $hend\n";
                    # write only the Xref's for active wormpep id's (CE.....)
                    $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                    $xref_line .= "$e $hstart $hend $start $end\n";
   		}
                else {
                    print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                    print ACE "$e $start $end $hstart $hend ";
                    print ACE "align $start $hstart\n";
                    $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                    $xref_line .= "$e $hstart $hend $start $end ";
                    $xref_line .= "align $hstart $start\n";
                    shift @cigars;
                    foreach my $string (@cigars) {
                        my ($coor, $hcoor) = split (/,/, $string);
                        print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                        print ACE "$e $start $end $hstart $hend ";
                        print ACE "align $coor $hcoor\n";
                        $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                        $xref_line .= "$e $hstart $hend $start $end ";
                        $xref_line .= "align $hcoor $coor\n";
		    }
		}
	    }
            push (@reverse_lines , $xref_line);
	}
    }
    # look at the "other" proteins (slimSwissProt and slimTrEMBL)
    my @orgs = qw(slimSwissProt slimTrEMBL);
    my $target_count = 0;
    foreach my $org (@orgs) {
        if ($target_count > $others_target_threshold) {
	    last;
	}
        unless (exists $org2analysis{$org}) {
            print LOG "\t\tno mapping of org to analysis for $org\n";
            next;
        }
        else {
            print LOG "\t\t$org\n";
        }	
        # make the sql query on the protein_feature table
        $sth_f->execute ($id, $org2analysis{$org}, $e_threshold);
        my $ref = $sth_f->fetchall_arrayref;
        # we would like to sort the columns:
        #   order the target sequences by their minimal E-value,
        #   and the HSP's per target sequence by their E-value (increasing)
        # (it would be better to do this in the sql query, but since the
        #  evalue is stored as string, '>' works but 'order by' doesn't)
        my %tmp;
        my %ev;
        foreach my $aref (@$ref) {
            # take -log(10) of evalues
            if ($aref->[7] != 0) {
                $aref->[7] = sprintf ("%.2f", -(log($aref->[7]))/log(10));
            }
            else {
                $aref->[7] = 999.99;
	    }
            # store the lowest evalue per target sequence as hash value
            if (exists $ev{$aref->[3]}) {
                if ($aref->[7] > $ev{$aref->[3]}) {
                    $ev{$aref->[3]} = $aref->[7];
   	        }
            }
            else {
                $ev{$aref->[3]} = $aref->[7];
       	    }
            push (@{$tmp{$aref->[3]}}, $aref);
        }
        # loop over all the targets
        my $xref_line;
        foreach my $target (sort {$ev{$b} <=> $ev{$a}} keys %ev) {
            # check if this match target passes the cutoff_evalue
            unless ($ev{$target} > $min_cutoff_ev) {
                last;
            }                    
            # check if any of the HSP's of this target hit a previously unmatched region
            my $switch = 0;
            foreach my $aref (sort {$b->[7] <=> $a->[7]} @{$tmp{$target}}) {
                my ($fid , $start, $end, $hid, $hstart, $hend, $strand, $e) = @$aref;
                for (my $i = $start ; $i <= $end ; $i++) {
                    if ($hit{$i} != 1) {
                        $switch = 1;
                        last;
		    }
   	        }
	    }
            unless ($switch == 1) {
                next;
            }       
            # check if we have dumped the max. number of targets
            if (++$target_count > $others_target_threshold) {
                last;
            }
            $xref_line = "Protein : \"${org2acedb{$org}}$target\"\n";
            # loop over all the HSP's
            foreach my $aref (sort {$b->[7] <=> $a->[7]} @{$tmp{$target}}) {
                # $target and $hid are the same (or we are in trouble...)
                my ($fid , $start, $end, $hid, $hstart, $hend, $strand, $e, $cigar) = @$aref;
                # check if the hsp is part of the sum stats
                if ($opt_s) {
                    unless ($e == $ev{$hid}) {
                        last;
	   	    }
		}
                # check the strand, and set start/end accordingly
                # (not really required for blastp, since all matches are "forward")
                if ($strand > 0) {}
                elsif ($strand < 0) {
                    my $tmp = $end;
                    $end = $start;
                    $start = $tmp;
     	        }
                else {
                    print LOG "incorrect strand assignment $strand for $id, feature $fid -> skip it\n";
                    next;
      	        }
                # write the HSP info to the ACE file
                my @cigars = split (/:/, $cigar);
                if (@cigars == 1) {
                    print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                    print ACE "$e $start $end $hstart $hend\n";
                    $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                    $xref_line .= "$e $hstart $hend $start $end\n"
   	        }
                else {
                    print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                    print ACE "$e $start $end $hstart $hend ";
                    print ACE "align $start $hstart\n";
                    $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                    $xref_line .= "$e $hstart $hend $start $end ";
                    $xref_line .= "align $hstart $start\n";
                    shift @cigars;
                    foreach my $string (@cigars) {
                        my ($coor, $hcoor) = split (/,/, $string);
                        print ACE "Pep_homol\t\"${org2acedb{$org}}$hid\" \"wublastp_$org\" ";
                        print ACE "$e $start $end $hstart $hend ";
                        print ACE "align $coor $hcoor\n";
                        $xref_line .= "Pep_homol\t\"WP:$id\" \"wublastp_$org\" ";
                        $xref_line .= "$e $hstart $hend $start $end ";
                        $xref_line .= "align $hcoor $coor\n";
		    }
   	        }
	    }
            push (@reverse_lines , $xref_line);
        }
    }
    foreach my $entry (@reverse_lines) {
        print ACE "\n";
        print ACE "$entry";
    }
}


print LOG "End of blastp dump [".&now."]:\n";
print "End of blastp dump [".&now."]:\n"; 

$sth->finish;
$sth_p->finish;
$sth_f->finish;
$dbh->disconnect;

close LOG;
close ACE;

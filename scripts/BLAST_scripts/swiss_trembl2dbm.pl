#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

# parses either the swissprot or trembl .dat flat files,
# and writes three DBM files: org, description, key words

use strict;
use Getopt::Std;
use GDBM_File;
use vars qw($opt_s $opt_t $opt_v);
#swissprot trembl verbose

getopts ("stv");
my $verbose = $opt_v;

my $output_dir = "/acari/work2a/wormpipe/swall_data";

my %ORG;
my %DES;
my %KEY;
my $id;
my $org;
my $total_org;
my $des;
my $key;
my $switch = 0;

my $usage = "cat swissprot/trembl .dat file | swiss_trembl2dmb.pl\n";
$usage .= "-s for swissprot\n";
$usage .= "-t for trembl\n";

if ($opt_s && $opt_t) {
    die "$usage";
}
elsif ($opt_s) {
  
  `rm $output_dir/swissprot2org` if (-e "$output_dir/swissprot2org" );
  `rm $output_dir/swissprot2des` if (-e "$output_dir/swissprot2des" );
  `rm $output_dir/swissprot2key` if (-e "$output_dir/swissprot2key" );
  
  tie %ORG,'GDBM_File',"$output_dir/swissprot2org",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
  tie %DES,'GDBM_File',"$output_dir/swissprot2des",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
  tie %KEY,'GDBM_File',"$output_dir/swissprot2key",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
}
elsif ($opt_t) {
  
  `rm $output_dir/trembl2org` if (-e "$output_dir/trembl2org" );
  `rm $output_dir/trembl2des` if (-e "$output_dir/trembl2des" );
  `rm $output_dir/trembl2key` if (-e "$output_dir/trembl2key" );
  
  tie %ORG,'GDBM_File', "$output_dir/trembl2org",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
  tie %DES,'GDBM_File', "$output_dir/trembl2des",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
  tie %KEY,'GDBM_File', "$output_dir/trembl2key",&GDBM_WRCREAT, 0666 or die "cannot open DBM file";
}
else {
    die "$usage";
}

while (my $line = <>) {
    if ($line =~ /^AC\s+(\S+)\;/) {
        $id = $1;
        $switch = 1;
    }
    elsif ($line =~ /^OS\s+/) {
        $line =~ s/^OS\s+//;
        $total_org .= $line;
    }
    elsif ($line =~ /^DE\s+/) {
        $line =~ s/^DE\s+//;
        $des .= $line;
    }
    elsif ($line =~ /^KW\s+/) {
        $line =~ s/^KW\s+//;
        $key .= $line;
    }

    elsif ($line =~ /^SQ\s+/ && $switch == 1) {
        # this gets rid of the tremblnew entries that are
        # in the same file with the swissprot entries
        # (currently, swall-1 = swissprot + tremblnew, swall-2 = trembl)
        if (($opt_s) || ($opt_t)) {
            #
            $total_org =~ s/\n/ /g;
            $total_org =~ s/\.//g;
            $total_org =~ s/\,/()/g;
            $total_org =~ s/\s+and\s+/()/g;
            chomp $total_org;
            my @split;
            my @tmp;
            @split = split (/\([^\(\)]*\)/, $total_org);
            my %seen;
            foreach (@split) {
                s/^\s+//;
                s/\s+$//;
                if ($_ eq "") {
                    next;
	        }
                if (exists $seen{$_}) {
                    next;
	        }
                else {
                    $seen{$_} = 1;
                    push (@tmp, $_);
                }                
	    }
            $org .= join (";", @tmp);
            if (exists $ORG{$id}) {
                print "ORG PRESENT\t$id\t($org)\n" if $verbose;
            }
            else {
                $ORG{$id} = $org;
                print "ORG ADDED\t$id\t($org)\n" if $verbose;
            }
            #
            chomp $des;
            $des =~ s/\n/\s/g;
	    $des =~ s/\"//g;
            if (exists $DES{$id}) {
                print "DES PRESENT\t$id\t($des)\n" if $verbose;
            }
            else {
                $DES{$id} = $des;
                print "DES ADDED\t$id\t($des)\n" if $verbose;
            }
            #
            chomp $key;
            $key =~ s/\n/\s/g;
            if (exists $KEY{$id}) {
                print "KEY PRESENT\t$id\t($key)\n" if $verbose;
            }
            else {
                $KEY{$id} = $key;
                print "KEY ADDED\t$id\t($key)\n" if $verbose;
            }
	}
        else {
            print "\tINCORRECT ID: $id\n";
	}
        $id = "";
        $org = "";
        $total_org = "";
        $des = "";
        $key = "";
        $switch = 0;
    }
}

untie %ORG;
untie %DES;
untie %KEY;






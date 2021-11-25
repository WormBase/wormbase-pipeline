#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

unless (@ARGV >0) {
        &USAGE;
}

sub USAGE {

die '
perl ~/bin/perl/FTP_checker.pl  ftp_species_dir1 ftp_species_dir2
Put in the full paths to the ftp locations you would like to compare between. The directory should be on the
level containing a species dir per species i.e. 
';

}

my @infs = @ARGV;

foreach my $in (@infs) {
    my $ino=$in;
    $ino=~s/\.//;
    open (OUT, "> $ino.fol.sync") || die "\nCannot write to file $ino.fol.sync\n";
    open (OUT2, "> $ino.fil.sync") || die "\nCannot write to file $ino.fil.sync\n";


    # Find all folders and files in the system
    # get the file sizes

    #my $inf = '/nfs/ftp/pub/databases/wormbase/releases/' .$in . '/species/c_elegans/PRJNA13758';
    #my $inf= '/nfs/production/panda/ensemblgenomes/wormbase/BUILD/TestFTP/releases/' . $in;
    #my $inf= '/nfs/production/panda/ensemblgenomes/wormbase/BUILD/FAKE_FTP/releases/' . $in;
    my $inf = '/nfs/ftp/pub/databases/wormbase/releases/' .$in;

    my @folders = `find $inf -type d -exec du -sk \'{}\' \\;`;
    my @files = `find $inf -type f -exec ls -l \'{}\' \\;`;

    foreach my $fol (@folders) {
        my @ar = split(/\s+/, $fol);
        $fol = "$ar[1]\t$ar[0]";
    }

    foreach my $fil (@files) {
        my @ar = split(/\s+/, $fil);
        $ar[4]=~s/ //g;
    #    $ar[8]=~s/^ //;
        $fil = "$ar[8]\t$ar[4]";
    }


    foreach my $fol (@folders) {
        print OUT "$fol\n";
    }

    foreach my $fil (@files) {
        print OUT2 "$fil\n";
    }

}


# Now compare file sizes and existance between all

my $i = 1;
my @a = (1..(scalar(@infs)-1));
#print "@a\n";

foreach(@a){

my $ino=$infs[0];
$ino=~s/\.//;

open (OUT, ">$ino.$infs[$i].fol.sync") || die "\nCannot write to file $ino.$infs[$i].fol.sync\n";
open (OUT2, ">$ino.$infs[$i].fil.sync") || die "\nCannot write to file $ino.$infs[$i].fil.sync\n";
#open (OUT2, ">$infs[0].$infs[$i].fol.sync.newcoms") || die "\nCannot write to file $infs[0].$infs[$i].fol.sync.newcoms\n";


    # Read in both lists of folders
    open (IN, "<$ino.fol.sync") || die "\nCannot write to file $ino.fol.sync\n";
    open (IN2, "<$infs[$i].fol.sync") || die "\nCannot write to file $infs[$i].fol.sync\n";

    my %in;
    my %in2;
    my %in2fols;

    while (<IN>) {
        chomp;
	my @ar = split(/\t/, $_);
	$in{$ar[0]}=$ar[1];
	#$in2fols{$ar[0]}=$ar[1];
    }
    while (<IN2>) {
        chomp;
	$_=~s/$infs[$i]/$infs[0]/g;
	my @ar = split(/\t/, $_);
	$in2{$ar[0]}=$ar[1];
    }

    # Get shared folders
    foreach my $key (keys %in2) {
        #print "#$key#\n";
        if (exists $in{$key}) {
		# Check if is empty
		if ($in{$key}<30 or $in2{$key}<30) {
			print OUT "WARNING EMPTY\t$key\t$in{$key}\t$in2{$key}\n";
		}
		# Check if it is simlar value or not
		# Skip folders we know might be different
		if ($key=~/MULTI_SPECIES\/hub/ or $key=~/PATCHES/ or $key=~/acedb\/database/) {
		}
		elsif ($in{$key} == $in2{$key} ) {
            		print OUT "IDENTICAL\t$key\t$in{$key}\t$in2{$key}\n";
    		}
		elsif ( abs(($in{$key} - $in2{$key})/$in{$key}) < 0.05 ) {
            		print OUT "SIMILAR\t$key\t$in{$key}\t$in2{$key}\n";
    		}
		else {
			print OUT "WARNING\t$key\t$in{$key}\t$in2{$key}\n";
		}

        }
        else {
		print OUT "ERROR Missing $key\n";

            }
        }



    close (IN);
    close (IN2);

    #__END__

    # Now lets deal with the files

    #my $i = 1;

    # Read in both lists of files
    open (IN, "<$ino.fil.sync") || die "\nCannot write to file $ino.fil.sync\n";
    open (IN2, "<$infs[$i].fil.sync") || die "\nCannot write to file $infs[$i].fil.sync\n";

    #print "\n#Reading in file-lists\n\n";
    my %fin;
    my %fin2;

    while (<IN>) {
        chomp;
	$_=~s/\.\./\./g;
	my @ar = split(/\t/, $_);
	$fin{$ar[0]}=$ar[1];
	#$in2fols{$ar[0]}=$ar[1];
    }
    while (<IN2>) {
        chomp;
	$_=~s/$infs[$i]/$infs[0]/g;
	$_=~s/\.\./\./g;
	my @ar = split(/\t/, $_);
	$fin2{$ar[0]}=$ar[1];
    }

    # Get shared folders
    foreach my $key (keys %fin2) {
        #print "#$key#\n";
        if (exists $fin{$key}) {
		# Check if is empty
		if ($fin{$key}<30 or $fin2{$key}<30) {
			print OUT2 "WARNING EMPTY\t$key\t$fin{$key}\t$fin2{$key}\n";
		}
		# Check if it is simlar value or not
		# Skip folders we know might be different
		if ($key=~/MULTI_SPECIES\/hub/ or $key=~/PATCHES/ or $key=~/acedb\/database/) {
		}
		elsif ($fin{$key} == $fin2{$key} ) {
            		print OUT2 "IDENTICAL\t$key\t$fin{$key}\t$fin2{$key}\n";
    		}
		elsif ( abs(($fin{$key} - $fin2{$key})/$fin{$key}) < 0.05 ) {
            		print OUT2 "SIMILAR\t$key\t$fin{$key}\t$fin2{$key}\n";
    		}
		else {
			print OUT2 "WARNING\t$key\t$fin{$key}\t$fin2{$key}\n";
		}

        }
        else {
		print OUT2 "ERROR Missing $key\n";

            }
        }

	$i++;
}

__END__
system "cat $in2.fol.sync.newcoms | sort > $in2.fol.sync.newcoms.sh";
system "rm -f $in2.fol.sync.newcoms ";
print "\n\n#Good news, Syncing is ready now.\n";
print "#Please go to the target folder and execute this command:\n\n";
print "bash $in2.fol.sync.newcoms.sh\n\n";
print "Then take this list of files and copy over:\n\n";
print "wc -l $in2.fol.sync.files_to_sopy\n\n";
print "These files exist in both locations but have different file sizes:\n\n";
print "wc -l $in2.fol.sync.partial_files\n\n";
}
else {
    &USAGE;
}
exit;
__END__
Or rsync:
rsync -av -e "ssh -p12345" clusteruser@localhost:/location/on/the/cluster .
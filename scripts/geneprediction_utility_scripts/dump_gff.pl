#!/usr/bin/env perl
# dump an acedb in GFF2/GFF3

use strict;
use Getopt::Long;
use Ace;

my ($targetDB,$sace);

GetOptions(
    'sourceDB=s'   => \$targetDB,
    'acedb=s'      => \$sace,
)||die($!);

die("need a target\n") unless ($targetDB);

$sace ||= '/scratch/acedb/giface';

my @seqs = dump_get_Sequences($targetDB);

dump_gff($targetDB,@seqs);

#################################################### Functions ###########################


# dumps the curated featuers on sequences as GFF2
# arguments:
#   (database, list of sequences to dump)
sub dump_gff{
        my ($database,@seqs)=@_;
       
        foreach my $s(@seqs){
           my $cmd = "gif seqget $s; seqfeatures -version 2 -file $database/$s.gff2"; 
           open (WRITEDB,"echo '$cmd' | $sace $database |")||die "$!\n";
           while (my $line = <WRITEDB>) {
	       print "ERROR detected while GFF dumping $s:\n\t$line\ngiface: $cmd\n" if $line=~/ERROR/;
	   }
	   close WRITEDB;
           system("cat $database/$s.gff2 >> $database/genome.gff2") && die($!);
           unlink "$database/$s.gff2"||die($!);

           my $cmd = "gif seqget $s; seqfeatures -version 3 -file $database/$s.gff3"; 
           open (WRITEDB,"echo '$cmd' | $sace $database |")||die "$!\n";
           while (my $line = <WRITEDB>) {
	       print "ERROR detected while GFF dumping $s:\n\t$line\ngiface: $cmd\n" if $line=~/ERROR/;
	   }
	   close WRITEDB;
           system("cat $database/$s.gff3 >> $database/genome.gff3") && die($!);
           unlink "$database/$s.gff3"||die($!);

        }
}

# find a  list of genomic_canonical sequences
# arguments:
#  (target database path)
sub dump_get_Sequences {
	my $dbpath = shift;

        # sequence bit
        my @sequenceNames;
        my $acedb = Ace->connect(-path => $dbpath) ||die(Ace->Error);
        my @sequences = $acedb->fetch(-query => 'Find Sequence; Genomic_canonical');
        foreach my $s(@sequences){
            push @sequenceNames,"$s";
        }

        return @sequenceNames;
}

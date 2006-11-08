package GSI;

use strict;

################################################################
# GSI - generic sequence index file support
# SRE, Mon Sep  7 07:56:57 1998
# Adapted from squid's gsi.c
#
# package variables available by $GSI::
#
# NOTE: Only one GSI file can be open at a time.
#
# Modified from S. Eddy source by Lorenzo Cerutti 2000/05/10
#                              by Marc Sohrmann (ms2@sanger.ac.uk) 2001/04/24
################################################################

# Size of the records to store the indexed information
my $RECSIZE = 70; # A64 + n + N

# fmt numbers (file format)
my $fmt_genbank = 2;
my $fmt_embl    = 4;
my $fmt_fasta   = 7;
my $fmt_pir     = 12;

# some setup variables
my $nfiles;
my $nkeys;

###########################################################
# Functions to write an indexed file                      #
###########################################################

###########################################################
# writeHeader: create the header
sub  writeHeaderRecord
{
    local *FH = shift;
    my ($nfiles, $nkeys) = @_;
    print FH pack "a64nN", "GSI", $nfiles, $nkeys;
}

###########################################################
# writeFile: write a file name
sub  writeFileRecord
{
    local *FH = shift;
    my ($file, $nfile, $fmt) = @_;
    print FH pack "a64nN", $file, $nfile, $fmt;
}

###########################################################
# writeKey: write a key => the same as above
*writeKeyRecord=\&writeFileRecord;
#sub  writeKeyRecord
#{
#    local *FH = shift;
#    my ($key, $nfile, $offset) = @_;
#    print FH pack "a64nN", $key, $nfile, $offset;
#}

###########################################################
# Functions to search a gsi indexed file                  #
###########################################################

###########################################################
# Open and setup variables
sub openGSI
{
	my $file = shift;
        my $record;
        my $check;
	open (GSI, "$file") || die "Failed to open GSI file $file";
	read(GSI, $record, $RECSIZE);
	($check, $nfiles, $nkeys) = unpack "A64nN", $record;
	die "File $file is not in the GSI format" if ($check ne "GSI");
}

###########################################################
# Close file
sub closeGSI
{
	close GSI;
}

###########################################################
# Get offset
sub getOffset
{
	my $key = shift;
	my $record;
        my ($name, $file_number, $offset);
        my ($seqfile, $fmt);
	my $left  = $nfiles + 1;
	my $right = $nfiles + $nkeys;
	my $mid   = int(($left + $right)/2);
	seek(GSI, $mid * $RECSIZE, 0);
	
	while (read(GSI, $record, $RECSIZE))
	{
		($name, $file_number, $offset) = unpack "A64nN", $record;
#                print "\t:$key  :$name  :$left  :$right\n";
		if    ($key  eq $name ) { last }
		elsif ($left >= $right) { return (-1, -1, -1) }
		elsif ($key  gt $name ) { $left  = int(($left+$right)/2) + 1 }
		else                    { $right = int(($left+$right)/2) - 1 }
	
		$mid = int(($left+$right)/2);
		seek(GSI, $mid * $RECSIZE, 0);
	}

        seek(GSI, $file_number * $RECSIZE, 0);
	read(GSI, $record, $RECSIZE);
	($seqfile, $file_number, $fmt) = unpack "A64nN", $record;
	return ($seqfile, $fmt, $offset);
}

1;

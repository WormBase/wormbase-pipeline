#!/usr/local/bin/perl
#
# check_protein_ID.pl
# v0.1
# dl
# 2000-07-12
#
# Consistency script to check camace for 'Protein_id' tags.
#  (a) warns if Predicted_gene does not have a Protein_id tag
#  (b) warns if Protein_id tag is corrupted
#
# Usage : check_protein_ID.pl 
#
# History:
#
# v0.1
# 000712 : dan : PP version
# 020116 : dan : moved to /wormsrv1/camace as acedb database

#use strict;
use Socket;
use vars qw ($debug $seq_len $sv_acc $sv_ver);
use lib '/wormsrv2/scripts';
use Wormbase;

$|=1;

#####################################################################################################
# Declare variables                                                                                 # 
#####################################################################################################

my $CDS = "";                                           # CDS name
my $genome_sequence = "";                                           # CDS name
my %clone2acc = "";
my %clone2id = "";
my %acc2ver = "";

my $tace = "/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace"; # tace executable path
my $db = "/wormsrv1/camace";                            # Database path
my $exec="$tace $db";     
my $command=<<EOF;
find Predicted_gene
show -a Protein_id
quit
EOF
my $command_db=<<EOF;
find Genome_sequence
show -a Database
quit
EOF






#####################################################################################################
# MAIN LOOP                                                                                         # 
#####################################################################################################
open(textace_db, "echo '$command_db' | $exec | ");
while (<textace_db>) {
#    print "$_";
    chomp;
    next if ($_ eq "");
    next if (/acedb/);
    next if (/\/\//);
    s/\"//g;
    if (/^Sequence : (\S+)/) {
	$genome_sequence = $1;
#	print "processing $genome_sequence\n";
	next;
    }
    if (/Database\s+EMBL\s+(\S+)\s+(\S+)\s+(\d+)/) {
	$clone2acc{$genome_sequence} = $2;
	$clone2id{$genome_sequence} = $1;
	$acc2ver{$2} = $3;
    }
    
}   
close textace_db;  


#####################################################################################################
# Count the No. of CDS objects with Protein_ID tags                                                 # 
#####################################################################################################

open(textace, "echo '$command' | $exec | ");
while (<textace>) {
#    print "$_";
    chomp;
    next if ($_ eq "");
    next if (/acedb/);
    next if (/\/\//);
    s/\"//g;
    if (/^Sequence : (\S+)/) {
	$CDS = $1;
	print "processing $CDS\n";
	next;
    }
    # correct syntax for the protein_id fields
    if (/Protein_id\s+(\S+)\s+(\S+)\s+(\d+)/) {
	my ($cam_clone,$cam_pid,$cam_ver) = ($1,$2,$3);
	$PID_count++;
	print "[acedb]  Protein_id $CDS in $cam_clone $cam_pid.$cam_ver\n";
	@embl_data = &getseq($CDS,$cam_pid);

	my ($EMBL_clone,$EMBL_pid,$EMBL_ver);
	
	foreach (@embl_data) {
	    ($EMBL_clone,$EMBL_pid,$EMBL_ver) = split (/ /,$_);
#	    print "[embl]   Protein_id $CDS in $EMBL_clone $EMBL_pid.$EMBL_ver\n";
	    
	    # if accessions match then check protein_id
	    if ($EMBL_clone eq $clone2acc{$cam_clone}) {
		print "[embl]   Protein_id $CDS in $EMBL_clone $EMBL_pid.$EMBL_ver\n";
		print "---> checking protein_id for clone $cam_clone [$clone2acc{$cam_clone}]\n";
		if ($EMBL_pid ne $cam_pid) {
		    print "Protein_ID for $CDS identifier is incorrect [$cam_clone $cam_pid.$cam_ver  EMBL: $EMBL_pid.$EMBL_ver]\n\n";
		} 
		elsif ($EMBL_ver ne $cam_ver) {
		    print "Protein_ID for $CDS version is incorrect [$cam_clone $cam_pid.$cam_ver  EMBL: $EMBL_pid.$EMBL_ver]\n\n";
		}
		else {
		    print "Protein_ID for $CDS is synchronised [$cam_clone $cam_pid.$cam_ver]\n\n";
		}
		print "End loop: Checked for genome_sequence $cam_clone\n\n";
		last;
	    }
	    else {
		print "Discard data: Incorrect clone\n";
	    }
	}
    }
    
#    print "\n";
    # eeek
    
}    
close textace;     

exit(0);

#####################################################################################################
#################

sub getseq {

    my ($cds,$protein_id) = @_;
    my ($EMBL_clone,$EMBL_pid,$EMBL_ver,$EMBL_misc);
    my @data;

    my $querycontent = "-e+[{SWALL_SP_REMTREMBL}-prd:'$protein_id']";
    my $request      = "/srs6bin/cgi-bin/wgetz?$querycontent";
    my $server       = "srs.ebi.ac.uk";

    if (!defined(open_TCP(*G,$server,80))) {
        print "Error connecting to server \n";
        exit(-1);
    }
    
    print G "GET $request HTTP/1.0\n\n";
    print G "Accept: */*\n";
    print G "User-Agent: socketsrs/1.0\n\n";
    
# Parsing annotation
    while (my $return_line=<G>) {
	
	if ($return_line =~ /embl-ProteinID/) {
#	    print "$return_line";
	    ($EMBL_clone) = $return_line =~ (/EMBL-acc\:(\S+)\]/);
	    ($EMBL_pid,$EMBL_ver) = $return_line =~ (/AMP_gt;parent\)\+\-e\"\>(\S+)\<\/A\>\.(\d+)\;/);
	    ($EMBL_misc) = $return_line =~ (/\d\;(.+)$/);
	    $EMBL_misc  =~ s/^\s+//g;
	    $EMBL_misc  =~ s/\.$//g;
#	    ($EMBL_pid,$EMBL_ver) = $return_line =~ (/AMP_gt;parent\)\+\-e\"\>(\S+)\<\/A\>\.(\d+)\;.[0,2]\-/);
	    print "[getseq] Protein_id $cds in $EMBL_clone $EMBL_pid.$EMBL_ver '$EMBL_misc'";
	    if ($EMBL_misc eq "JOINED") {
		print "\t !! Discard !!\n";
		next;
	    }
	    else {
		print "\n";
	    }
	    push (@data,"$EMBL_clone $EMBL_pid $EMBL_ver");
#	    return ($EMBL_clone,$EMBL_pid,$EMBL_ver);
#	    last;
	}
    }
    close G;
    
    return (@data);

}


 ########################################################
 # Output: successful network connection in file handle #
 ########################################################

sub open_TCP {
    my ($FS,$dest,$port) = @_;
    my $proto = getprotobyname ('tcp');
    socket ($FS,PF_INET,SOCK_STREAM,$proto);
    my $sin = sockaddr_in($port,inet_aton($dest));
    connect($FS,$sin) || return undef;
    my $old_fh = select($FS);
    $| = 1;
    select($old_fh);
}


#####################################################################################################
# Tidy up                                                                                           # 
#####################################################################################################

unlink '/tmp/gene_list';
unlink '/tmp/protein_ID_list';

exit(0);

__END__

=pod

=head1 NAME - check_protein_ID.pl

=head2 AUTHOR

Dan Lawson dl1@sanger.ac.uk

=head2 USAGE

check_protein_ID.pl 

check_protein_ID.pl uses tace to query an ACeDB database for predicted 
gene objects (default database is camace). It will then check the 
'Protein_id' field for the presence/absence of the tag and perform a
simple consistency check (3 tag format - '?Sequence Text Int').

The script will report objects with corrupted tags and in the 
case of a discrepency between the number of prediced CDS and the
number with Protein_id tags (i.e. missing data), those CDS which
do not yet have a Protein_ID in the ACeDB database.

=head2 MANDATORY ARGUMENTS

=over 2

=item *
none

=back
=cut

=head2 REVISIONS

=over 2

=item *
v0.1
2000-07-12 : PP version

=back
=cut




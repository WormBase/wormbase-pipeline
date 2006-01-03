#!/usr/local/bin/perl5.8.0 -w
#
# getProteinID.pl
#
# Parses the weekly mail from Nadeem to get the Protein_ID and
# SWALL accession for C.elegans entries in EMBL 
#
# written by Dan Lawson
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2006-01-03 15:24:38 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Data::Dumper;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $file;    # specify protein ID file to use (defaults to ~wormpub/protein_ID.mail)
my $load;    # option specifies whether resulting acefile will be loaded to autoace


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
            "file=s"     => \$file,
	    "load"       => \$load,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES


##########################
# MAIN BODY OF SCRIPT
##########################


# fetch hashes made by other scripts
my %acc2clone = &FetchData('accession2clone');
my %gene2CE   = &FetchData('cds2wormpep');
my %Ip2Go     = &FetchData('interpro2go');
my %cds2gene  = &FetchData('cds2wbgene_id');
our %swall;

# get the InterPro => GOterm mapping ( space separated ie 'IPR004794' => 'GO:0008703 GO:0008835 GO:0009231 ')

our %databases = (
		  'SW' => 'SWISSPROT',
		  'TR' => 'TREMBL',
		  'TN' => 'TREMBLNEW'
		  );

my %db_ids_acc = (
		  'SW_id' => 'SwissProt_ID',
		  'SW_ac' => 'SwissProt_AC',
		  'TR_ac' => 'TrEMBL_AC',
		  'TN_ac' => 'TrEMBLNEW_AC'
		 );


# get swall data
&getswalldata;

# set default file if -file not specified on command line
($file = "/nfs/disk100/wormpub/protein_ID.mail") if (!defined($file));

my $ace_file = "$ace_dir/acefiles/WormpepACandIDs.ace";
$ace_file = "/tmp/WormpepACandIDs.ace" if ($debug);

$log->write_to("Using $file as protein_id file - writing to $ace_file\n");

open (OUT, ">$ace_file");

my @f = "";
my ($dbxref_ac,$dbxref_id,$dbxref_db);

# File format for protein_ID file
# <acc>          <ver>    <proteinID>   <ver>     <checksum>     <gene_name>  <swall> <standard_name>
# non-genome entry (e.g. from a mRNA) 
# AB000913        2       BAA21715        1       2898126014      unc-14      O15940   
# genome entry
# AC006633        10      AAK68374        1       3313183671      F35B3.3     Q966K1   F35B3.3

open (FILE, "<$file") || die "Can't open the protein_ID file\n";
while (<FILE>) {
    chomp;
    next if ($_ eq "");
    next if ( (/^From/) || (/^Date/) || (/To/) || (/^Subject/) );

    @f = split /\t/;

    # discard non genome_sequence entries
    next unless (defined $acc2clone{$f[0]});
    
    my @interpro = "";

    # Standard name

    print "$f[7]\t" if ($verbose);
    my $protein = $gene2CE{$f[7]};                      # incremented array slice to handle new SWALL column
    
    if (!defined( $swall{$f[2]} ) or ( $swall{$f[2]}{Accession} ne $f[6]) ){ 
	my $ebi = $f[6] ? $f[6] : "-";
	my $getz = $swall{$f[2]}{Accession} ? $swall{$f[2]}{Accession} : "-";
	$log->write_to("ERROR:  mismatch between getz and EBI for $protein [EBI:$ebi|GETZ:$getz)\n");
	next;
    }

 
    # Protein entries

    if ($protein and $swall{$f[2]}{Database}) {
	print "$protein" if ($verbose);
	print OUT "\nProtein : \"WP:$protein\"\n";
	
	if ( ($swall{$f[2]}{Identifier}) and  ($swall{$f[2]}{Identifier} ne $swall{$f[2]}{Accession}) ) {
	    print OUT "Database $databases{ $swall{$f[2]}{Database} } ",$db_ids_acc{ $swall{$f[2]}{Database}."_id" }," $swall{$f[2]}{Identifier}\n"
	} #creates tons of uninitialized values

	print OUT "Database $databases{ $swall{$f[2]}{Database} } ",$db_ids_acc{ $swall{$f[2]}{Database}."_ac" }," $swall{$f[2]}{Accession}\n";
	
        if ( $swall{$f[2]} ) {
	    foreach (@{$swall{$f[2]}{Interpro}}) {
		next if ($_ eq "");
		print OUT "Motif_homol\t\"INTERPRO:$_\"\n"
		}	
	}
	else {
	    print "ERROR: gene $f[7] has no protein (has common_data been updated ?)\n";
	    next;
	    }
    }
    
    print "\n" if ($verbose);

    # CDS entry

    print OUT "\nCDS : \"$f[7]\"\n";
    print OUT "Protein_id \"$acc2clone{$f[0]}\" $f[2] $f[3]\n";
    if ( $swall{$f[2]}{Database} ) {
	print OUT "Database $databases{ $swall{$f[2]}{Database} } ",$db_ids_acc{ $swall{$f[2]}{Database}."_id" }," $swall{$f[2]}{Identifier}\n"	
                     if ( ($swall{$f[2]}{Identifier}) and  ($swall{$f[2]}{Identifier} ne $swall{$f[2]}{Accession}) ); #creates tons of uninitialized values
	print OUT "Database $databases{ $swall{$f[2]}{Database} } ",$db_ids_acc{ $swall{$f[2]}{Database}."_ac" }," $swall{$f[2]}{Accession}\n";
    }

    if (@{$swall{$f[2]}{Interpro}}) {
	# assign GO terms based on InterPro Motifs
	foreach my $ip (@{$swall{$f[2]}{Interpro}}) {
	    if ( exists $Ip2Go{$ip} ) {
		my @GOterms = split(/\s/,$Ip2Go{$ip});
		foreach my $go (@GOterms) {
		    print OUT "GO_term\t\"$go\" \"IEA\" Inferred_automatically\n";
		}
	    }
	}
	# connect Gene to GO_term as well.
	my $gene = $cds2gene{$f[7]};
	unless ( $gene ) {
	    print STDERR "no gene_id for $f[7]\n";
	    next;
	}
	
	# Gene entry
	print OUT "\nGene : $gene\n";
	foreach my $ip (@{$swall{$f[2]}{Interpro}}) {
	    if ( exists $Ip2Go{$ip} ) {
		my @GOterms = split(/\s/,$Ip2Go{$ip});
		foreach my $go (@GOterms) {
		    print OUT "GO_term\t\"$go\" \"IEA\" Inferred_automatically\n";
		}
	    }
	}
    }
}


# now load to autoace if -load specified
if ($load) {
    my $command = "autoace_minder.pl -load $ace_dir/acefiles/WormpepACandIDs.ace -tsuser wormpep_IDs";
    my $status = system($command);
    if (($status >>8) != 0) {
	die "ERROR: Loading WormpepACandIDs.ace file failed \$\? = $status\n";
    }
}

close FILE;
close OUT;

$log->mail;
exit(0);

##############################################################
# Subroutines
##############################################################


sub getswalldata {
    
# ID   1433_CAEEL     STANDARD;      PRT;   248 AA. 
# AC   P41932; Q21537; 
# DR   EMBL; U05038; AAA61872.1; -.
# DR   EMBL; Z73910; CAA98138.1; -.
# DR   PIR; JC2581; JC2581.
# DR   PIR; T23759; T23759.
# DR   HSSP; P29312; 1A38.
# DR   WormPep; M117.2; CE06200.
# DR   InterPro; IPR000308; 14-3-3.
# DR   Pfam; PF00244; 14-3-3; 1.
# DR   PRINTS; PR00305; 1433ZETA.
# DR   SMART; SM00101; 14_3_3; 1.
# DR   PROSITE; PS00796; 1433_1; 1.
# DR   PROSITE; PS00797; 1433_2; 1. 

    my $acc; 
    my $id; 
    my $db;
    my $text;
    my @interpro;
    my @proteinID;

    $/ = "ID";

    open (LOOK, "/usr/local/pubseq/bin/getzc -f 'id acc dbxref'  \"[SWALL-organism:Caenorhabditis elegans]\" |");
    while (<LOOK>) {
	$text = "ID" . $_;         # add the leading 'ID'
	chop $text;                # remove the trailing 'ID'
	chop $text;

	next if ($text eq "ID");     # shortcircuit loop for empty set at beginning

#	print "\n$text\n";

	if ($text =~ /ID\s+(\S+)/)  {                                        # ID  Q95YB2  PRELIMINARY;  PRT;  409 AA.
	    $id = $1;
#	    print "assign ID\n";
	}
	if ($text =~ /AC\s+(\S+);/) {                                        # AC  Q95YB2;
	    $acc = $1;
	    # assign SwissProt : TrEMBL based on diff between ID and ACC eg SW = YPC1_CAEEL Q11178 vs TR = Q52GY9_CAEEL Q52GY9
	    $db = $id =~ /$acc/ ? 'TR' : 'SW';
	}
	while ($text =~ /DR\s+EMBL;\s\S+\s(\S+)\.\d+; \-\;/g) {            # DR   EMBL; U05038; AAA61872.1; -.
	    push (@proteinID, $1);
#	    print "assign proteinID\n";
	}
	while ($text =~ /DR\s+InterPro;\s(\w+)/g ) {                      # DR  InterPro; IPR006089; Acyl-CoA_dh.
	    push (@interpro, $1);
#	    print "assign InterPro\n";
	}
	
	foreach my $i (@proteinID) {
	    next if ($i eq "");
	    $swall{$i}{Identifier} = $id;
	    $swall{$i}{Accession}  = $acc;
	    $swall{$i}{Database}   = $db;
	    foreach my $j (@interpro) {
		push (@{$swall{$i}{Interpro}},$j);
	    }
	    
	    print "// Assign to proteinID $i [ID:$id\tAC:$acc\tDB:$db\tInterPro: " . (join ' ',@interpro) . "]\n" if $verbose;
	}

	@proteinID = "";
	@interpro = "";
	$acc = "";
	$id  = "";
	$text = "";

    }

    $/ = "\n";

}




##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


__END__

=pod

=head2 NAME - getProteinID.pl

=head1 USAGE

=over 4

=item getProteinID [-options]

=back

=head4 MANDATORY arguments:

=over 4

=back

-head4 OPTIONAL arguments:

=over 4

=item -file <file>

Extracted e-mail from EBI containing the Protein_ID data

=item -help

This help

=item -verbose

Extra command line output, shows status of each WormPep protein (runs slower!)

=item -load

if specified will load resulting acefile to autoace

=back

=head1 AUTHOR

=over 4

=item Dan Lawson (dl1@sanger.ac.uk)

=back

=cut

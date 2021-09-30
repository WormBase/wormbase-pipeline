#!/usr/local/bin/perl5.8.0 -w
use strict;
use Getopt::Long;
use lib -e $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

my ($filein, $fileout, $wormbase, $USER, $species, $debug, $test, $old);

GetOptions (
    'filein=s'  => \$filein,
    'fileout=s' => \$fileout,
    'user:s'    => \$USER,
    'debug:s'   => \$debug,
    'old'       => \$old,
    );

$species = "elegans";
$wormbase = Wormbase->new("-organism" =>$species, -debug => $debug, -test => $test);
my $db;
my $tace            = $wormbase->tace;
my $database = "/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/geneace";
my $ace = Ace->connect (-path => $database,
			-program => $tace) || die "cannot connect to database at $database\n";

unless ((defined $filein) && ($fileout)) {
    print "You need to specif both -filein and -fileout\n";
}
print "Using $filein\n" if defined($filein);
print "Creating $fileout\n" if defined($fileout);


my $curator;
if (defined $USER){ 
    if ($USER eq 'pad') {
        $curator = 'WBPerson1983';
    } elsif ($USER eq 'skd') {
        $curator = 'WBPerson51134';
    } elsif ($USER eq 'mz3') {
        $curator = 'WBPerson21950';
    }
}
open (IN, "<$filein") or die("Failed to open input file\n");
open (OUT, ">$fileout") or die("Failed to open output file\n");

#0 - Strain - RG3235
#1 - Allele - ve735
#2 - Gene - Y39B6A.3
#3 - Del size  - 852 bp
#4 - allele description - Y39B6A.3(ve735[LoxP + myo-2p::GFP::unc-54 3' UTR + rps-27p::neoR::unc-54 3' UTR + LoxP]) V.
#5 - Phenotype - Homozygous early larval arrest
#6 - Upstream - tttattagcattttttctagaatgtacacg
#7 - Downstream - tttttttctgtaaattttttacgaaaatat
#8 - Strain - RG3235
#9 - Genotype - Y39B6A.3(ve735[LoxP + myo-2p::GFP::unc-54 3' UTR + rps-27p::neoR::unc-54 3' UTR + LoxP])/sC4(s2172) [dpy-21(e428)] V.
#10 - Description - Homozygous early larval arrest. Deletion of 852 bp with Calarco/Colaiacovo selection cassette conferring myo-2 GFP and G418 resistance inserted at break in parental strain N2. Heterozygotes are wild-type GFP+, and segregate wild-type GFP+, GFP+ arrested larvae (ve7#35 homozygotes) and arrested non-GFP (stage unknown) (sC4 homozygotes). Maintain by picking wild-type GFP+. Left flanking Sequence: tttattagcattttttctagaatgtacacg ; Right flanking sequence: tttttttctgtaaattttttacgaaaatat. sgRNA #3: CGTCACCGATAAGCTATCGT; sgRNA #4: ggtaaactacacgcgtggcc. Please reference Au et al., G3 9(1): 135-144 2019 in any work resulting from use of this mutation.
#11 - Chromosome - V
#12 - sgRNA #1 CGTCACCGATAAGCTATCGT
#13 - sgRNA #2 ggtaaactacacgcgtggcc



my ($wbgene_id, $wbstrain_id, $maptg, $type);
my $count = "0"; 
while (<IN>) {
    my @obj;
    my $lab;
    chomp;
    if (/Deletion Size/) {
	next;
    }
    my @f;
    if ($old) {
	@f=split",";
    }
    else {
	@f=split"\t";
    }
    if ($count =~ /0$/){
	print "Processed $count\n";
    }


    print "Processing : $f[0]\n" if ($debug);
    $count++;
    unless ($old){
	if ($f[0] ne $f[8]) {
	    print "STRAIN Missmatch $f[0] $f[8]\n";
	}
    }
    #code currently only takes a single target gene, if there are large deletions these will currently not get connected to a Gene object.
    # Only a single example in the 240 data points we have.
    
    if ($old) {
	my $ngene;
	my $var;
	my $wbvar_id;
#	$f[3] = s/Please reference Au et al.; G3 9(1): 135-144 2019 in any work resulting from use of this mutation.//;
	my $strain_name = $f[0];
	my @strobj = $ace->fetch(-query => "FIND Strain where Public_name = $strain_name ");
	if (defined $strobj[0]) {
	    $wbstrain_id = $strobj[0]->name;
	    print "Found Strain : \"$wbstrain_id\"\n" if ($debug);
	}
	else {
	    print "\nCan't find Strain : \"$strain_name\"\n" if ($debug);
	    $wbstrain_id = $strain_name;
	}
	    print OUT "\nStrain $wbstrain_id\nPublic_name $f[0]\nLive\nSpecies \"Caenorhabditis elegans\"\nLocation RG\nEvidence Curator_confirmed WBPerson1983\nGenotype \"$f[1]\"\nRemark $f[3]\n";
#	if ($f[3] =~ /Deletion/){print OUT "Type_of_mutation Deletion\n";}
	if ($f[3] =~ /135-144 2019/){print OUT "Reference WBPaper00055671\n";}
	if ($f[1] =~ /(ve\d+)/){
	    $var = $1;
	    my @varobj = $ace->fetch(-query => "FIND Variation where Public_name = $var");
	    if (defined $varobj[0]) {
		$wbvar_id = $varobj[0]->name;
		print "Found Variation : \"$wbvar_id\"\n" if ($debug);
	    }
	    else {
		print "\nCan't find Variation : \"$var\"\n" if ($debug);
		$wbvar_id = $var;
	    }
	    print OUT "\nVariation : $wbvar_id\nEvidence Curator_confirmed WBPerson1983\nPublic_name $var\n";}
	if ($f[3] =~ /Left flanking Sequence: ([a-zA-Z]+)/){print OUT "Flanking_sequences $1"}
	if ($f[3] =~ /Right flanking sequence: ([a-zA-Z]+)/){print OUT " $1\n";}
	my $mapping_target;
	if ($f[1] =~ /(\S+)\(/){$ngene = $1;
				if ($ngene =~ /\S+\-\d+/) {
				    $type = "Public_name";
				}
				elsif ($ngene =~ /(\S+)\.\d+/) {
				    $type = "Sequence_name";
				    $mapping_target = $1;
				    print OUT "Mapping_target $mapping_target\n";
				}
				@obj = $ace->fetch(-query => "FIND Gene where $type = $ngene");
				if (defined $obj[0]) {
				    $wbgene_id = $obj[0]->name;
				    print "Found Gene : \"$wbgene_id\"\n" if ($debug);
				    print OUT "Gene $wbgene_id\n";
				    unless (defined $mapping_target){
					my $seq_name = $obj[0]->Sequence_name->name;
					if ($seq_name =~ /(\S+)\.\d+/) {
					    $mapping_target = $1;
					    print OUT "Mapping_target $mapping_target\n";
					}
				    }
				}
				else {
				    print OUT "Gene $ngene\n";
				    print OUT "Mapping_target $ngene\n"; 
				}
	}
	print OUT "Deletion\nSequenced\nEngineered_allele\nSpecies \"Caenorhabditis elegans\"\nStrain $f[0]\nCRISPR_Cas9\nLive\nRemark $f[3]\nMethod Engineered_allele\n\n";
    }
    
    else {
    my $gene = $f[2];
    if ($gene =~ /,/){print "Multi gene deletion this will need manual curation\n";}
    if ($gene =~ /\S+\-\d+/) {
	$type = "Public_name";
    }
    elsif ($gene =~ /\S+\.\d+/) {
	$type = "Sequence_name";
    }   
    @obj = $ace->fetch(-query => "FIND Gene where $type = $gene");
    #Dead genes with a Sequence_based Public_name exception catch
    unless (defined $obj[0]) {
	@obj = $ace->fetch(-query => "FIND Gene where Public_name = $gene");
    }
    if (defined $obj[0]) {
        $wbgene_id = $obj[0]->name;
        print "Found Gene : \"$wbgene_id\"\n" if ($debug);
	#check status of the gene as dead genes cause issues.
	if ($obj[0]->Status->name =~ 'Dead') {
	    if ($obj[0]->Merged_into) {
		my $mergeg = $obj[0]->Merged_into->name;
		#switch to the gene the dead has been merged into.
		print "Switching Gene $wbgene_id for $mergeg as there has been a gene merge even.\n";
		@obj = $ace->fetch(-query => "FIND Gene $mergeg");
	    }
	    else {
		print "\nCan't retrieve a live gene for this one $f[0] $gene.\n";
		if (defined $f[11]) {
		    $maptg = "CHROMOSOME_$f[11]";
		}
		else {
		    print "\nCan't even infer mapping target from supplied chromosome\n";
		}
	    }
	}
	unless ($obj[0]->Status->name =~ 'Dead') {
	    $maptg = $obj[0]->Sequence_name->name;
	    if ($maptg =~ /(\S+)\.\d+/) {
		$maptg = $1;
	    }
	}
    }
    else {
	print "Can't find Gene : \"$gene\"\n";
	$wbgene_id = $gene;
    }
    
    my $strain_name = $f[0];
    my @strobj = $ace->fetch(-query => "FIND Strain where Public_name = $strain_name ");
    if (defined $strobj[0]) {
	$wbstrain_id = $strobj[0]->name;
	print "Found Strain : \"$wbstrain_id\"\n" if ($debug);
    }
    else {
	print "\nCan't find Strain : \"$strain_name\"\n" if ($debug);
	$wbstrain_id = $strain_name;
    }
    
    if ($f[0] =~ /([A-Z]+)/) {
	$lab = $1;
    } 
    my $wbvar_id;
    my $var_name = $f[1];
    my @varobj = $ace->fetch(-query => "FIND Variation where Public_name = $var_name ");
    if (defined $varobj[0]) {
        $wbvar_id = $varobj[0]->name;
        print "Found Variation : \"$wbvar_id\"\n" if ($debug);
    }
    else {
        print "\nCan't find Variation : \"$var_name\"\n" if ($debug);
        $wbvar_id = $var_name;
    }

    
    #Print the Variation Information    
    print OUT "\nVariation \: \"$wbvar_id\"\nPublic_name $var_name\nStrain $wbstrain_id\n";
    print OUT "Evidence Curator_confirmed $curator\n" if (defined $curator);
    my ($lflank, $rflank);
    if (defined $f[6]) {
	$lflank = lc $f[6];
    }
    if (defined $f[7]){
	$rflank= lc $f[7];
    }
    print OUT "Flanking_sequences $lflank $rflank\n";
    if (defined $maptg) {print OUT "Mapping_target $maptg\n";}
    print OUT "Deletion\nSequenced\nEngineered_allele\nSpecies \"Caenorhabditis elegans\"\nCRISPR_Cas9\nLive\n";
    print OUT "Gene $wbgene_id\nRemark \"Whole gene deletions made by the Rougvie, Moerman, Hutter and Sternberg labs that replace the coding sequence with a [LoxP + myo-2p::GFP::unc-54 3Prime UTR + rps-27p::neoR::unc-54 3Prime UTR + LoxP] cassette.\"\nMethod Engineered_allele\n";
    print OUT "Remark \"Allele Description\:$f[4]\"\n";
    if (defined $f[12]){print OUT "Remark \"sgRNA1:$f[12] sgRNA2: $f[13]\"\n\n";}

    #Print out the Strain information
    print OUT "\nStrain $wbstrain_id\nPublic_name $strain_name\n";
    print OUT "Evidence Curator_confirmed $curator\n" if (defined $curator);
    print OUT "Genotype \"$f[9]\"\nLocation $lab\nLive\nSpecies \"Caenorhabditis elegans\"\nVariation $wbvar_id\nRemark \"$f[10]\"\n";   
    }
}
print "Final Number $count\n";
print "Diaskeda same Poli\n"; #we had alot of fun#
close (IN);
close (OUT);
exit(0);
__END__

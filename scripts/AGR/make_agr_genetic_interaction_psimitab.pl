#!/usr/bin/perl
# instructions
# https://docs.google.com/spreadsheets/d/1DE1Ba3XPd0L2yXYwSq2rY1cmE79VrvDNqCV1Fr-_tS4/edit#gid=0

use Getopt::Long;
use IO::File;
use Ace;
use strict;
#breaks Ace.pm use diagnostics;

my ($dbpath,$test,$outfile);
GetOptions(
	'database=s' => \$dbpath,
	'test=s'     => \$test,
	'outfile=s'  => \$outfile,
)||die(@!);


my $of = IO::File->new($outfile,'w');
die("cannot open $outfile\n") unless $of;

my %interactorType;
$interactorType{TABTAB}{any}                                          = '-';
$interactorType{TABTABDiverging}{any}                                 = '-';
$interactorType{TABTABNeutral}{any}                                   = '-';
$interactorType{TABAll_suppressingTAB}{Effector}                      = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{TABAll_suppressingTAB}{Affected}                      = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{TABAll_suppressingTAB}{Non_directional}               = '-';
$interactorType{TABEnhancingTAB}{Effector}                            = 'psi-mi:"MI:2351"(enhancer gene)';
$interactorType{TABEnhancingTAB}{Affected}                            = 'psi-mi:"MI:2352"(enhanced gene)';
$interactorType{TABEnhancingTAB}{Non_directional}                     = '-';
$interactorType{TABMaskingTAB}{Effector}                              = 'psi-mi:"MI:2353"(epistatic gene)';
$interactorType{TABMaskingTAB}{Affected}                              = 'psi-mi:"MI:2354"(hypostatic gene)';
$interactorType{TABMaskingTAB}{Non_directional}                       = '-';
$interactorType{TABSub_suppressingTAB}{Effector}                      = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{TABSub_suppressingTAB}{Affected}                      = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{TABSub_suppressingTAB}{Non_directional}               = '-';
$interactorType{TABSuper_suppressingTAB}{Effector}                    = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{TABSuper_suppressingTAB}{Affected}                    = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{TABSuppressingTAB}{Effector}                          = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{TABSuppressingTAB}{Affected}                          = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{TABSuppressingTAB}{Non_directional}                   = '-';
$interactorType{A_phenotypicTABTABDiverging}{any}                     = '-';
$interactorType{A_phenotypicTABTABNeutral}{any}                       = '-';
$interactorType{Cis_phenotypicTABTABNeutral}{any}                     = '-';
$interactorType{Cis_phenotypicTABAll_suppressingTAB}{Non_directional} = '-';
$interactorType{Cis_phenotypicTABCo_suppressingTAB}{Non_directional}  = '-';
$interactorType{Cis_phenotypicTABEnhancingTAB}{any}                   = '-';
$interactorType{Cis_phenotypicTABEnhancingTABDiverging}{any}          = '-';
$interactorType{Cis_phenotypicTABEnhancingTABNeutral}{any}            = '-';
$interactorType{Cis_phenotypicTABInter_suppressingTAB}{any}           = '-';
$interactorType{Cis_phenotypicTABMaskingTAB}{Effector}                = 'psi-mi:"MI:2353"(epistatic gene)';
$interactorType{Cis_phenotypicTABMaskingTAB}{Affected}                = 'psi-mi:"MI:2354"(hypostatic gene)';
$interactorType{Cis_phenotypicTABMaskingTAB}{Non_directional}         = '-';
$interactorType{Cis_phenotypicTABSemi_suppressingTAB}{Effector}       = 'psi-mi:"MI:2353"(epistatic gene)';
$interactorType{Cis_phenotypicTABSemi_suppressingTAB}{Affected}       = 'psi-mi:"MI:2354"(hypostatic gene)';
$interactorType{Cis_phenotypicTABSuper_suppressingTAB}{any}           = '-';
$interactorType{Cis_phenotypicTABSuppressingTAB}{any}                 = '-';
$interactorType{Iso_phenotypicTABMaskingTAB}{any}                     = '-';
$interactorType{Mono_phenotypicTABTABNeutral}{any}                    = '-';
$interactorType{Mono_phenotypicTABAll_suppressingTAB}{Effector}       = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{Mono_phenotypicTABAll_suppressingTAB}{Affected}       = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{Mono_phenotypicTABEnhancingTAB}{Effector}             = 'psi-mi:"MI:2351"(enhancer gene)';
$interactorType{Mono_phenotypicTABEnhancingTAB}{Affected}             = 'psi-mi:"MI:2352"(enhanced gene)';
$interactorType{Mono_phenotypicTABEnhancingTABDiverging}{Effector}    = 'psi-mi:"MI:2351"(enhancer gene)';
$interactorType{Mono_phenotypicTABEnhancingTABDiverging}{Affected}    = 'psi-mi:"MI:2352"(enhanced gene)';
$interactorType{Mono_phenotypicTABSub_suppressingTAB}{Effector}       = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{Mono_phenotypicTABSub_suppressingTAB}{Affected}       = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{Mono_phenotypicTABSuppressingTAB}{Effector}           = 'psi-mi:"MI:0581"(suppressor gene)';
$interactorType{Mono_phenotypicTABSuppressingTAB}{Affected}           = 'psi-mi:"MI:0582"(suppressed gene)';
$interactorType{Trans_phenotypicTABAll_suppressingTAB}{any}           = '-';
$interactorType{Trans_phenotypicTABAll_suppressingTABNeutral}{any}    = '-';
$interactorType{Trans_phenotypicTABEnhancingTAB}{Effector}            = 'psi-mi:"MI:2351"(enhancer gene)';
$interactorType{Trans_phenotypicTABEnhancingTAB}{Affected}            = 'psi-mi:"MI:2352"(enhanced gene)';
$interactorType{Trans_phenotypicTABMaskingTAB}{Effector}              = 'psi-mi:"MI:2353"(epistatic gene)';
$interactorType{Trans_phenotypicTABMaskingTAB}{Affected}              = 'psi-mi:"MI:2354"(hypostatic gene)';
$interactorType{Trans_phenotypicTABSuppressingTAB}{any}               = '-';
$interactorType{Trans_phenotypicTABSub_suppressingTAB}{any}           = '-';

my %modToPsi = (
	TABTAB 					=> 'psi-mi:"MI:2402"(genetic interaction)',
	TABTABDiverging				=> 'psi-mi:"MI:0794"(synthetic)',
	TABTABNeutral				=> 'psi-mi:"MI:0932"(neutral multigenic phenotype result)',
	TABAll_suppressingTAB	 		=> 'psi-mi:"MI:1290"(genetic suppression (complete))',
	TABEnhancingTAB				=> 'psi-mi:"MI:2380"(genetic enhancement)',
	TABMaskingTAB				=> 'psi-mi:"MI:0797"(genetic epistasis (sensu Bateson))',
	TABSub_suppressingTAB 			=> 'psi-mi:"MI:1291"(genetic suppression (partial))',
	TABSuper_suppressingTAB			=> 'psi-mi:"MI:2388"(genetic over-suppression)',
	TABSuppressingTAB 			=> 'psi-mi:"MI:2379"(genetic suppression)',
	A_phenotypicTABTABDiverging 		=> 'psi-mi:"MI:0794"(synthetic genetic interaction)',
	A_phenotypicTABTABNeutral 		=> 'psi-mi:"MI:2398"(aphenotypic neutral multigenic phenotype result)',
	Cis_phenotypicTABTABNeutral 		=> 'psi-mi:"MI:2392"(mutual genetic enhancement (expected))',
	Cis_phenotypicTABAll_suppressingTAB 	=> 'psi-mi:"MI:2396"(cisphenotypic genetic suppression (complete))',
	Cis_phenotypicTABCo_suppressingTAB 	=> 'psi-mi:"MI:2397"(cisphenotypic co-suppressing genetic interaction)',
	Cis_phenotypicTABEnhancingTAB 		=> 'psi-mi:"MI:2401"(mutual genetic enhancement)',
	Cis_phenotypicTABEnhancingTABDiverging 	=> 'psi-mi:"MI:1278"(mutual genetic enhancement (sensu unexpected))',
	Cis_phenotypicTABEnhancingTABNeutral 	=> 'psi-mi:"MI:2392"(mutual genetic enhancement (expected))',
	Cis_phenotypicTABInter_suppressingTAB 	=> 'psi-mi:"MI:1283"(cisphenotypic inter-suppressing genetic interaction)',
	Cis_phenotypicTABMaskingTAB 		=> 'psi-mi:"MI:1273"(maximal genetic epistasis)',
	Cis_phenotypicTABSemi_suppressingTAB 	=> 'psi-mi:"MI:1274"(minimal genetic epistasis)',
	Cis_phenotypicTABSuper_suppressingTAB 	=> 'psi-mi:"MI:1287"(mutual genetic over-suppression)',
	Cis_phenotypicTABSuppressingTAB		=> 'psi-mi:"MI:1280"(cisphenotypic genetic suppression)',
	Iso_phenotypicTABMaskingTAB 		=> 'psi-mi:"MI:0795"(asynthetic genetic interaction)',
	Mono_phenotypicTABTABNeutral 		=> 'psi-mi:"MI:2400"(monophenotypic neutral multigenic phenotype result)',
	Mono_phenotypicTABAll_suppressingTAB 	=> 'psi-mi:"MI:1293"(monophenotypic genetic suppression (complete))',
	Mono_phenotypicTABEnhancingTAB 		=> 'psi-mi:"MI:1279"(monophenotypic genetic enhancement)',
	Mono_phenotypicTABEnhancingTABDiverging => 'psi-mi:"MI:1279"(monophenotypic genetic enhancement)',
	Mono_phenotypicTABSub_suppressingTAB 	=> 'psi-mi:"MI:1294"(monophenotypic genetic suppression (partial))',
	Mono_phenotypicTABSuppressingTAB 	=> 'psi-mi:"MI:1292"(monophenotypic genetic suppression)',
	Trans_phenotypicTABAll_suppressingTAB 	=> 'psi-mi:"MI:2391"(transphenotypic genetic suppression (complete))',
	Trans_phenotypicTABAll_suppressingTABNeutral => 'psi-mi:"MI:2393"(transphenotypic genetic suppression (expected))',
	Trans_phenotypicTABEnhancingTAB		=> 'psi-mi:"MI:1288"(transphenotypic enhancing genetic interaction)',
	Trans_phenotypicTABMaskingTAB 		=> 'psi-mi:"MI:1276"(opposing genetic epistasis)',
	Trans_phenotypicTABSuppressingTAB 	=> 'psi-mi:"MI:2390"(transphenotypic genetic suppression)',
	Trans_phenotypicTABSub_suppressingTAB 	=> 'psi-mi:"MI:2390"(transphenotypic genetic suppression)',
);

my $db = Ace->connect(-path => $dbpath)||die(Ace->Error);

my %intValidGiNoOther;


my $it = $test ? $db->fetch_many(Interaction => $test) : $db->fetch_many(-query => 'Find Interaction; Interaction_type="Genetic" AND COUNT(Interactor_overlapping_gene)=2 AND !Other_interactor AND !Molecule_interactor AND !Rearrangement AND Paper');

while (my $int = $it->next) {

  my %intGene;
  my %intVar;
  my %varGene;
  my %intValidTransgene;
  my %intPhen;
  my %ginName;
  my %ginSeq;
  my %papPmid;
  my %intGeneType;
  my %intSuppress;
  my %intValidGiRnai;
    
  $intValidGiNoOther{"$int"}++;

  my ($rnai)=$int->Interaction_RNAi;
  my ($gene,$effectorGene)=$int->Interactor_overlapping_gene;
  my ($variation) = $int->Variation_interactor;

  next unless $gene->Species eq 'Caenorhabditis elegans';
  map {$papPmid{"$_"}=get_pmid($_)} $int->Paper; # get the pmID => MEDLINE PMID XYZ
  next unless values %papPmid;

  if ($rnai){
    my %interaction_genes = map {$_->name => 1} $int->Interactor_overlapping_gene;
    foreach my $r ($int->Interaction_RNAi){
          map {
		  if ($interaction_genes{"$_"}){
			  $intValidGiRnai{"$int"}{"$_"}{"$r"}++
		  }else{
			  $intSuppress{$int}++
		  }
	  } $r->Gene;
    }
  }
  if ($int->Variation_interactor){
	  my @interactors = $int->Variation_interactor;

	  foreach my $v (@interactors){
		  my @g = $v->Gene;
                  $intSuppress{"$int"}++ if scalar(@g) > 1;
		  map {$varGene{"$v"} = "$_" if ("$_" eq "$gene" ||"$_" eq "$effectorGene")} @g;
	  }
  }
  if (grep {/Transgene/} $int->Detection_method){
	  map {
		  my ($t) = $_->asAce =~/WBTransgene\d+/g;
		  $intValidTransgene{"$int"}{"$_"}{"$t"}++ if $t;
	  } $int->Interactor_overlapping_gene;
  }
  if ($variation){
          map{$intVar{"$int"}{"$_"}++}$int->Variation_interactor;
  }
  map {$intPhen{"$int"}{"$_"}++} $int->Interaction_phenotype;
  map { 
	  if($intGeneType{"$int"}{"$_"}){ 
		  $intSuppress{"$int"}++ 
	  }else{
		  $intGeneType{"$int"}{"$_"} = "${\$_->right(2)}";
	  }
          $ginSeq{"$_"}="${\$_->Sequence_name}" if $_->Sequence_name;
	  $ginName{"$_"}="${\$_->Public_name}" if $_->Public_name;$intGene{"$int"}{"$_"}++;
  } $int->Interactor_overlapping_gene;

  ##########################
  my $type  = $int->Genetic;
  my $mod1  = $int->GI_module_one||'';
  my $mod2  = $int->GI_module_two||'';
  my $mod3  = $int->GI_module_three||'';
  my $pap   = $int->Paper;
  my $brief = $pap->Brief_citation;

  next if $intSuppress{$int};
  next if ($mod3 eq 'Neutral');

  my $key = $mod1 . 'TAB' . $mod2 . 'TAB' . $mod3;
  unless ($modToPsi{$key}) { 
	  print STDERR "ERR: $key does not map to psi-mi\n";
	  next;
  }

  my @tab;
  for(0 .. 41){ push @tab, '-'} # initialize tabs

  my $psi = $modToPsi{$key};
  $tab[11] = $psi;
  my $ginCount = scalar keys %{ $intGene{$int} };
  next if ($ginCount > 2); 
  my ($gene1, $gene2) = sort keys %{ $intGeneType{$int} };
  next unless $gene2; 
  my $type1 = $intGeneType{$int}{$gene1};
  my $type2 = $intGeneType{$int}{$gene2};

  if (($type2 eq 'Affected') && ($type1 eq 'Effector')) { # if the order is inverted 
	    ($gene2,$gene1) = ($gene1,$gene2); 
	    ($type2,$type1) = ($type1,$type2); 
  }

  $tab[0] = 'wormbase:' . $gene1; 
  $tab[1] = 'wormbase:' . $gene2;

  my @gin1Name; 
  my @gin2Name;
  if ($ginName{$gene1}){push @gin1Name, 'wormbase:' . $ginName{$gene1} . '(public_name)'  }
  if ($ginSeq{$gene1}) {push @gin1Name, 'wormbase:' . $ginSeq{$gene1}  . '(sequence_name)'}
  if ($ginName{$gene2}){push @gin2Name, 'wormbase:' . $ginName{$gene2} . '(public_name)'  }
  if ($ginSeq{$gene2}) {push @gin2Name, 'wormbase:' . $ginSeq{$gene2}  . '(sequence_name)'}

  if (scalar @gin1Name > 0) { $tab[4] = join"|", @gin1Name }
  if (scalar @gin2Name > 0) { $tab[5] = join"|", @gin2Name }

  $tab[6] = 'psi-mi:"MI:0254"(genetic interference)';

  if ($brief =~ m/^\s+/){
	  $brief =~ s/^\s+// 
  }

  if($brief =~ m/^(.*?\))/){
	  $tab[7] = $1 
  }elsif($brief =~ m/^(.*?et al)/){
   	  $tab[7] = $1 
  }elsif($brief =~ m/^([\S+]\s[\S+])\s/){
	  $tab[7] = $1 
  }elsif($brief =~ m/^([\S+])\s/){        
	  $tab[7] = $1
  }elsif($brief =~ m/^([\S+])/){
	  $tab[7] = $1
  }else{
	  $tab[7] = '-'; 
	  print ERR "$int brief citation does not match : $brief\n" 
  }

  if ($papPmid{$pap}) { 
      $tab[8] = $papPmid{$pap}
  } else { 
      next
  }

  $tab[9]  = 'taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)';
  $tab[10] = 'taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)';
  $tab[12] = 'psi-mi:"MI:0487"(wormbase)';
  $tab[13] = 'wormbase:' . $int;

  unless ($interactorType{$key}) { 
	  print ERR "$key does not map to interactorType\n";
	  next
  }

  if($interactorType{$key}{'any'}) { 
      $tab[18] = $interactorType{$key}{'any'}; 
      $tab[19] = $interactorType{$key}{'any'}
  }else {
      if($interactorType{$key}{$type1}){ 
	      $tab[18] = $interactorType{$key}{$type1}
      }else{
	      print ERR "gene1 $gene1 has type $type1 does not map to interactorType\n"
      }

      if($interactorType{$key}{$type2}){
	      $tab[19] = $interactorType{$key}{$type2}
      }else{
	      print ERR "gene2 $gene2 has type $type2 does not map to interactorType\n"
      }
  }

  $tab[20] = 'psi-mi:"MI:0250"(gene)';
  $tab[21] = 'psi-mi:"MI:0250"(gene)';

  my $skipBecauseVar = 0;

  my @tab25; 
  my @tab26;
  if ($intVar{$int}){

   my $varCount = scalar keys %{ $intVar{$int} };
   $skipBecauseVar++  if ($varCount > 2);

   foreach my $var (sort keys %{ $intVar{$int} }){
      my $varGene = $varGene{$var};
      print STDERR "$var doesn't have an attached gene\n" unless $varGene;
      $skipBecauseVar++ if $varGene eq 'skip';
      next unless $varGene;
      push @tab25, "wormbase:$var" if $varGene eq $gene1;
      push @tab26, "wormbase:$var" if $varGene eq $gene2;
   }
  }
  next if $skipBecauseVar;

  unless ($intValidGiNoOther{$int}){
    print ERR "$int no intValidGiNoOther\n";
    next;
  }

  if ($intValidGiRnai{$int}{$gene1})    { push @tab25, 'wormbase:rnai'}
  if ($intValidGiRnai{$int}{$gene2})    { push @tab26, 'wormbase:rnai'}
  if ($intValidTransgene{$int}{$gene1}) { push @tab25, 'wormbase:transgene'}
  if ($intValidTransgene{$int}{$gene2}) { push @tab26, 'wormbase:transgene'}

  if (scalar keys @tab25 > 1){
	  print STDERR qq(too many types in column 26 $int $gene1 @tab25\n);
  }
  $tab[25] = $tab25[0];
  $tab[25] ||= '-';

  if (scalar keys @tab26 > 1){
	  print STDERR qq(too many types in column 27 $int $gene2 @tab26\n);
  }
  $tab[26] = $tab26[0];
  $tab[26] ||= '-';


  my (@phens) = map {"wormbase:$_"} sort keys %{ $intPhen{$int} };

  if(scalar @phens > 0){
	  $tab[27] = join"|", @phens
  }

  $tab[35]='false'; # because Chris

  my $line = join "\t", @tab;
  print $of "$line\n";
}
$of->close;
$db->close;
exit(0);

sub get_pmid {
  my ($paper) = @_;
  my $pmid;
  foreach my $db ($paper->Database) {
	    if ($db->name eq 'MEDLINE') {
	      $pmid = $db->right->right->name;
	      last;
	    }
	}
  $pmid = $pmid ? "pubmed:$pmid" : undef;
  return $pmid;
}

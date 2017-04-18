
#!/software/bin/perl -w
#
# generate_RRID_data.pl
#
# Usage : generate_RRID_data.pl [-options]
#
# Can be run against the build or geneace via the -database option.
#
# Target file unless specified - autoace/acefiles/WS<wormbase_version>_RRIDs.dat
# Soirce database unless specified - autoace
#
use strict;
use Getopt::Long;
use Storable;
use lib $ENV{CVS_DIR};
use Wormbase;

my ($test,
    $debug,
    $store,
    $species,
    $database,
    $wormbase,
    $outfile,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "outfile:s"       => \$outfile,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
      );
}

#################################
# Set up some useful stuff      #
#################################
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
my $tace = $wormbase->tace;
my $full_species_name = $wormbase->full_name;
my $wormbase_version = $wormbase->get_wormbase_version_name;
my $dbdir;

if (defined $database) {
  $dbdir = $database;
}
unless (defined $database) {
  $dbdir = $wormbase->autoace;
}

if (not defined $outfile) {
  $outfile = $wormbase->acefiles . "/${wormbase_version}_RRIDs.dat";
}

##########################
# MAIN BODY OF SCRIPT
##########################

my (%genes,
    %geneid2gene,
    %WBvar2name,
    %WBPaper2pmid,
    %straingenotypes,
    %strainspecies,
    %strainoutcross,
    %strainmutagen,
    %strainmadeby,
    %Straingenes,
    %Strainremarks,
    %Strainvars,
    %wbpersonhuman,
    %Strainrefs);
my @Unique_strains;
my $out_fh;
$log->write_to("Generating RRID Query\n");
my $query = &generate_RRDID_query();
$log->write_to("Retrieving RRID data from $dbdir\n");
my $command = "Table-maker -p $query\nquit\n";
my $count2;
open(my $tacefh,  "echo '$command' | $tace $dbdir |");
## Strain Genotype Species VarID Variation GeneID Gene outcross Mutagen Remark Reference PersonID Name ##
print "Table_maker_finished\n";
while(<$tacefh>) {
  $count2++;

  unless (/\"/) {next;}
  chomp; 
  s/\"//g;
  s/\\//g;
  #print "$_\n";
  my $F;
  my @F = split"\t",$_,-1;
  if ($test){print "$F[0] Main $count2\n";}

  #Store WBID :: Human readables in associated arrays.
  if ((exists $F[4])&&(!grep (/$F[4]/, @{$geneid2gene{$F[3]}}))) {push @{$geneid2gene{$F[3]}}, $F[4];}
  #if (!grep (/$F[4]/, @{$WBvar2name{$F[3]}}))  {push @{$WBvar2name{$F[3]}}, $F[4];}
  if (!grep (/$F[0]/, @Unique_strains))        {push @Unique_strains, $F[0];}

  #store Strain: data
  if (exists $F[1]) {push @{$straingenotypes{$F[0]}}, $F[1];}
  if ((exists $F[2])&&($F[2] =~ /\S+/))   {push @{$strainspecies{$F[0]}}, $F[2];}
  if ((exists $F[3])&&(!grep (/$F[3]/, @{$Straingenes{$F[0]}})))     {push @{$Straingenes{$F[0]}}, $F[3];}
  if ((exists $F[5])&&($F[5] =~ /\S+/))  {push @{$strainoutcross{$F[0]}}, $F[5];}
  if ((exists $F[6])&&($F[6] =~ /\S+/))   {push @{$strainmutagen{$F[0]}}, $F[6];}
  #\Q\E makes it ignore special characters in the data.
  if ((exists $F[7])&&(!grep (/\Q$F[7]\E/, @{$Strainremarks{$F[0]}})&&($F[7] =~ /\S+/)))   {push @{$Strainremarks{$F[0]}}, $F[7];}
  if ((exists $F[8])&&($F[8] =~ /\S+/)&&(!grep (/$F[8]/, @{$Strainrefs{$F[0]}})))     {push @{$Strainrefs{$F[0]}}, $F[8];}
  if ((exists $F[9])&&($F[9] =~ /\S+/))   {push @{$strainmadeby{$F[0]}}, $F[9];}
  if ((exists $F[10])&&($F[10] =~ /\S+/)) {push @{$wbpersonhuman{$F[9]}}, $F[10];}
  if ($F[11] =~ /\S+/) {
    if ((exists $F[11])&&(!grep (/$F[11]/, @{$WBPaper2pmid{$F[8]}})))  {
      push @{$WBPaper2pmid{$F[8]}}, $F[11];
    }
  }
}
$log->write_to("Finished Retrieving RRID data\n");


my $query2 = &generate_RRDIDvar_query();
$log->write_to("Retrieving RRID var  data from $dbdir\n");
my $command2 = "Table-maker -p $query2\nquit\n";
my $count3;
open(my $tacefh2,  "echo '$command2' | $tace $dbdir |");

## Strain Genotype Species VarID Variation GeneID Gene outcross Mutagen Remark Reference PersonID Name ##
print "Var Table_maker_finished\n";
while(<$tacefh2>) {
  $count3++;
  unless (/\"/) {next;}
  chomp; s/\"//g;
  #print "$_\n";
  my $G;
  my @G = split"\t";
  if ($test) {print "$count3 variation $G[0]\n"}
  if ((exists $G[2])&&(!grep (/$G[2]/, @{$WBvar2name{$G[1]}})))  {push @{$WBvar2name{$G[1]}}, $G[2];}
  if ((exists $G[1])&&!grep (/$G[1]/, @{$Strainvars{$G[0]}}))      {push @{$Strainvars{$G[0]}}, $G[1];}

}
$log->write_to("Finished Retrieving RRID Var data\n");


open($out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");
print "Output file opened $outfile\n";

my $count;
print $out_fh "#Version $wormbase_version\n#\n#\n# This file contains information on all non-wild isolate strains in the WS$wormbase_version release of WormBase.\n";
print $out_fh "# ID\tSpecies\tGenotype\tVariants\tGenes modified\tStrain Info\tMutagen\tComments\tReferences\tMade By\tURL\n";
foreach my $strain (@Unique_strains) {
  $count++;
  if ($test) {print "$count Processing $strain\n";}
  # spit out array of genes - print "$strain: @{ $Straingenes{$strain} }\n";

  #(Identifying number/id perorganism,)
  my $strainws;
  if ($strain =~ /\s+/) {
    ($strainws = $strain) =~ s/\s/_/g;
  print $out_fh "WB-STRAIN\:$strainws\t";
  }
  else {
  print $out_fh "WB-STRAIN\:$strain\t";
  }
  #(name of organism,)
  if (defined$strainspecies{$strain}) {
    print $out_fh "\"$strainspecies{$strain}[0]\"\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }

  #Genotype
  if (exists $straingenotypes{$strain}) {
    if ($straingenotypes{$strain}[0] =~ /\S+/) {
      my $straingeno = $straingenotypes{$strain}[0];
      print $out_fh "\"$straingeno\"\t";
    }
    else {
      print $out_fh "EMPTY\t";
    }
  }
  else {
    print $out_fh "EMPTY\t";
  }

  #Variations genetic modification of organism if modified, 
  if (exists $Strainvars{$strain}) {
    if (defined  $Strainvars{$strain}) {
      print $out_fh "$Strainvars{$strain}[0]\($WBvar2name{$Strainvars{$strain}[0]}[0]\)";
      foreach my $i ( 1 .. $#{ $Strainvars{$strain} } ) {
	print $out_fh "|$Strainvars{$strain}[$i]\($WBvar2name{$Strainvars{$strain}[$i]}[0]\)";
      }
      print $out_fh "\t";
    }
  }
  else {
    print $out_fh "EMPTY\t";
  }
  #Genes (gene(s) modified if any (list is better than array))
  if ((exists $Straingenes{$strain}[0]) && (($Straingenes{$strain}[0]) =~ /\w+/)){
      print $out_fh "$Straingenes{$strain}[0]\($geneid2gene{$Straingenes{$strain}[0]}[0]\)";
      foreach my $i ( 1 .. $#{ $Straingenes{$strain} } ) {
	print $out_fh "|$Straingenes{$strain}[$i]\($geneid2gene{$Straingenes{$strain}[$i]}[0]\)";
      }
      print $out_fh "\t";
    }
  else {
    print $out_fh "EMPTY\t";
  }
  #inbred strain information (outcrossed etc), 
  if (exists $strainoutcross{$strain}) {
    print $out_fh "$strainoutcross{$strain}[0]\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }
  if (exists $strainmutagen{$strain}[0]) {
    print $out_fh "$strainmutagen{$strain}[0]\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }
  # possible multiple remarks
  if (exists $Strainremarks{$strain}[0]) {
    print $out_fh "\"$Strainremarks{$strain}[0]\"";
    foreach my $i ( 1 .. $#{ $Strainremarks{$strain} } ) {
      print $out_fh "|\"$Strainremarks{$strain}[$i]\"";
    }
    print $out_fh "\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }

  if (defined $Strainrefs{$strain}) {
    my @refarray = @{$Strainrefs{$strain}}; 
    foreach my $tmppaper (@refarray) {
      if ((exists $WBPaper2pmid{$tmppaper}) && ($WBPaper2pmid{$tmppaper}[0] =~ /\S+/)) {
	if($tmppaper == $refarray[-1]  ) {
	  print $out_fh "$tmppaper(PMID:$WBPaper2pmid{$tmppaper}[0])";
	}
	else {
	  print $out_fh "$tmppaper(PMID:$WBPaper2pmid{$tmppaper}[0])|";
	}
      }
      else {
	if($tmppaper == $refarray[-1]  ) {
	  print $out_fh "$tmppaper(PMID:EMPTY)";
	}
	else {
	  print $out_fh "$tmppaper(PMID:EMPTY)|";
	}
      }
    }
    print $out_fh "\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }


  if (defined $strainmadeby{$strain}) {
    print $out_fh "$strainmadeby{$strain}[0]\($wbpersonhuman{$strainmadeby{$strain}[0]}[0]\)\t";
  }
  else {
    print $out_fh "EMPTY\t";
  }

  #URL http://www.wormbase.org/db/get?name=$species;class=STRAIN
  print $out_fh "http://www.wormbase.org/db/get?name=${strain};class=STRAIN";


  print $out_fh "\n";
}
close($out_fh);
$log->write_to("\nRRID data generated for $count strains........ ;) \n\n");
$log->mail();
exit(0);

sub generate_RRDIDvar_query {
  my ($full_species) = @_;

  my $tmdef2 = "/tmp/gene_tmqueryvar.$$.def";
  open my $qfhv, ">$tmdef2" or 
      $log->log_and_die("Could not open $tmdef2 for writing\n");  

  my $condition = "";

  my $tablemaker_template2 = <<"EOF";

Sortcolumn 1

Colonne 1 
Subtitle Strain 
Width 12 
Optional 
Visible 
Class 
Class Strain
Condition !Wild_isolate
From 1 

//wildisolate

Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Variation 
From 1 
Tag Variation  
 
Colonne 3 
Subtitle Variation 
Width 12 
Optional 
Visible 
Class 
Class Variation_name 
From 2 
Tag Public_name  

EOF

  print $qfhv $tablemaker_template2;
  return $tmdef2;
}



sub generate_RRDID_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/gene_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";


#Colonne 1 
#Subtitle Strain 
#Width 12 
#Optional 
#Visible 
#Class 
#Class Strain
#Condition "N2*"
#From 1

  my $tablemaker_template = <<"EOF";
 
Sortcolumn 1

Colonne 1 
Subtitle Strain 
Width 12 
Optional 
Visible 
Class 
Class Strain
Condition Species
From 1 
 
Colonne 2 
Subtitle Genotype 
Width 25 
Optional 
Visible 
Text 
From 1 
Tag Genotype  
 
Colonne 3 
Subtitle Species 
Width 12 
Optional 
Visible 
Class 
Class Species 
From 1 
Tag Species  
  
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene  
 
Colonne 5 
Subtitle Gene 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 4 
Tag Public_name  
 
Colonne 6 
Subtitle Outcrossed 
Width 12 
Optional 
Visible 
Text 
From 1 
Tag Outcrossed  
 
Colonne 7 
Subtitle Mutagen 
Width 12 
Optional 
Visible 
Text 
From 1 
Tag Mutagen  
 
Colonne 8 
Subtitle Remark 
Width 12 
Optional 
Visible 
Class 
Class Text 
From 1 
Tag Remark  
 
Colonne 9 
Subtitle Reference 
Width 12 
Optional 
Visible 
Class 
Class Paper 
From 1 
Tag Reference  
 
Colonne 10 
Width 12 
Optional 
Visible
Class 
Class Person 
From 1 
Tag Made_by  
 
Colonne 11 
Subtitle Made_by 
Width 12 
Optional 
Visible 
Class 
Class Person_name 
From 10 
Tag Standard_name 

Colonne 12 
Width 12 
Optional 
Hidden 
Class 
Class Database 
From 9 
Tag Database 
 
Colonne 13 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 12 
Tag  HERE  
Condition PMID

Colonne 14 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 13 
Tag  HERE
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}
__END__

=pod

=head2 NAME - EMBL_Sequencefetch_species.pl

=head1 USAGE

=over 4

=item EMBL_Sequencefetch_species.pl  [-options]

=back

This script gets sequence data for all species (currently hacked to store the taxid within the script)

EMBL_Sequencefetch_species.pl

MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4

=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.

=back

=over 4

=item -test, Test mode.

=back

=over 4

=item -species, only performs syncronisation for given species

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut

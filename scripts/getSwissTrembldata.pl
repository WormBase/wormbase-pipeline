#!/usr/local/bin/perl5.6.1 -w
#
# getSwissTrembldata.pl
#
# dl
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2002-09-05 12:29:58 $


use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;
use Ace;

##################
# Paths and stuff
##################

my $tace      = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $dbdir     = "/wormsrv2/autoace";
$ENV{'ACEDB'} = $dbdir;

our %databases = (
	      'SW' => 'SWISSPROT',
	      'TR' => 'TREMBL',
	      'TN' => 'TREMBLNEW'
	      );


open (ACE_GSC,  ">/nfs/disk100/wormpub/analysis/SWALL/output_stl.ace");
open (ACE_WTSI, ">/nfs/disk100/wormpub/analysis/SWALL/output_cam.ace");
open (OUTPUT,   ">/nfs/disk100/wormpub/analysis/SWALL/output_autoace");


# Grab accessions and sequence versions from autoace
my $command=<<EOF;
Table-maker -p "$dbdir/wquery/accession2clone.def"
quit
EOF

open (TACE, "echo '$command' | $tace | ") || die "Could not open pipe to tace\n";

while (<TACE>) {
  print;
  chomp;
  s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
  next if ($_ eq "");
  next if (/\/\//);
  s/\"//g;

  my ($acename,$acc) = (/^(\S+)\s+(\S+)/);
  #    print "\n// Parsing genome sequence $acename [$acc]\n";
  my $GSC = "";
  my $CDS_xref_count = 0;
  my $CDS_found_count = 0;
  
  my ($CDS_dbxref,$CDS_gene,$CDS_protein,$CDS_prod,$CDS_name,$pid,$ver,$db,$acc,$EMBL_acc,$CDS_xref_count);
  my ($CDS_on,$CDS_dbxref_ac,$CDS_dbxref_id,$CDS_dbxref_db);
  my $carryover = 0;
  
  open (PFETCH, "/usr/local/pubseq/bin/pfetch -F $acc |");
  while (<PFETCH>) {
    chomp;
    if (/^SV\s+(\S+)\.\d+/) {$EMBL_acc = $1; next;}
    if (/^DR/)              {$CDS_xref_count++; next;}
    #	if (/^DR/)              {print "DR lines $_\n"; $CDS_xref_count++; next;}
    if (/^FH\s+Key/)        {next;}
    #	if (/^FH\s+Key/)        {print "\nExpecting $CDS_xref_count CDS in this entry\n\n";next;}
    
    # begin CDS loop
    if (/^FT\s+CDS/) {
      $CDS_on = 1;
      ($CDS_dbxref,$CDS_gene,$CDS_protein,$CDS_prod,$CDS_name) = "";
      ($pid,$ver,$db,$acc) = "";
      next;
    }
    
    if (/^FT\s+\/gene=\"(\S+)\"/) {
      $CDS_gene = $1;
      next;
    }
    
    if (/^FT\s+\/product=\"(\S+.+)/) {
      $CDS_prod = $1;
      
      # Set file on
      $GSC = 1;
      
      # line ends in a speech mark (i.e. full entry)
      if ($CDS_prod =~ /\"/) {
	chop $CDS_prod;
	next;
      }
      else {
	$carryover = 1;
	next;
      }
    }
    if ($carryover == 1) {
      (/^FT\s+(\S+.+)/);
      $CDS_prod .= " $1";
      if ($CDS_prod =~ /\"/) {
	chop $CDS_prod;
	$carryover = 0;
	next;
      }
      else {
	$carryover = 1;
      }
    }
    
    if (/^FT\s+\/protein_id=\"(\S+)\"/) {
      $CDS_protein = $1;
      next;
    }
    
    if (/^FT\s+\/translation/) {
      if (substr($CDS_prod,-1) eq ")") {chop $CDS_prod};
      
      ($CDS_name) = $CDS_prod =~ (/(\S+)$/);
      
      # check the entry from the protein_id directly
      
      #	    print "Checking with protein_id [$CDS_protein]\n";
      ($CDS_dbxref_ac,$CDS_dbxref_id,$CDS_dbxref_db) = &get_from_protein_id($CDS_protein);
      
      ($pid,$ver) = split(/\./,$CDS_protein);
      
      if ($GSC) {
	print OUTPUT "EMBL [$acename|$EMBL_acc] CDS [gene:$CDS_name protein_id:$CDS_protein DB_xref:$CDS_dbxref_ac]\n";
	print ACE_GSC "\nSequence : \"$CDS_name\"\nProtein_id $acename $pid $ver\n";
	print ACE_GSC "Database $databases{$CDS_dbxref_db} $CDS_dbxref_id $CDS_dbxref_ac\n";
      }
      else {
	print OUTPUT "EMBL [$acename|$EMBL_acc] CDS [gene:$CDS_gene protein_id:$CDS_protein DB_xref:$CDS_dbxref_ac]\n";
	print ACE_WTSI "\nSequence : \"$CDS_gene\"\nProtein_id $acename $pid $ver\n";
	print ACE_WTSI "Database $databases{$CDS_dbxref_db} $CDS_dbxref_id $CDS_dbxref_ac\n";
      }
      $CDS_found_count++;
      next;
    }
    
    if (/^SQ/) {
      #	    print "\nFound $CDS_found_count CDS in this entry\n\n";
      
      next;
    }
    
  }
  close PFETCH;
}

close(TACE);
close ACE_GSC;
close ACE_WTSI;
close OUTPUT;

exit(0);



sub get_from_protein_id {
    
  my $protein_id = shift;
  my $acc; 
  my $id; 
  my $db;
  
  open (LOOK, "/usr/local/pubseq/bin/pfetch -F $protein_id |");
  while (<LOOK>) {
    
    if (/^AC\s+(\S+);/) {
      $acc = $1;
    }
    if (/^ID\s+(\S+)/) {
      $id  = $1;
    }
  }
  close LOOK;
  
  
  if ((length ($acc)) == 8) {
    $db = "TN";
  }
  elsif ($id =~ /_CAEEL/) {
    $db = "SW";
  }
  else {
    $db = "TR";
  }
  
  return ($acc,$id,$db);
}



#!/usr/local/bin/perl5.6.0
#
# getSwissTrembldata.pl
#
# dl
#
# 

%databases = (
	      'SW' => 'SWISSPROT',
	      'TR' => 'TREMBL',
	      'TN' => 'TREMBLNEW'
	      );

$file = shift;

open (ACE_GSC,  ">/nfs/disk100/wormpub/analysis/SWALL/output_stl.ace");
open (ACE_WTSI, ">/nfs/disk100/wormpub/analysis/SWALL/output_cam.ace");
open (OUTPUT,   ">/nfs/disk100/wormpub/analysis/SWALL/output_autoace");

open (FILE, "<$file");
while (<FILE>) {
    s/\"//g;
    ($acename,$acc) = (/^(\S+)\s+(\S+)/);

#    print "\n// Parsing genome sequence $acename [$acc]\n";
    $CDS_xref_count  = 0;
    $CDS_found_count = 0;
    
    undef ($GSC);
    
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
close FILE;

close ACE_GSC;
close ACE_WTSI;
close OUTPUT;

exit(0);



sub get_from_protein_id {
    
    my $protein_id = shift;
    my $acc; my $id; 
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



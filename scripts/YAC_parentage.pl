#!/usr/local/bin/perl
#
# YAC_parentage.pl v1.0
# dl
# 1999-11-29
#
# Aceperl script to check YAC parent-child relationship. 
#
###########################################################################################
#
# 991129 dl : PP version


########################################
# iniatialise                          #
########################################

use Ace;

########################################
# command-line parsing                 #
########################################

while ($ARGV[0] =~ /^-/) {
    $_=shift;
    if (/^-c(.*)/) {
       $ace=1;
    }
    elsif (/^-s(.*)/) {
        $ace=2;
    }
    elsif (/^-d(.*)/) {
        $debug=1;
    }
    elsif (/^-w(.*)/) {
        $html=1;
    }
    else {
        &usage;
    }
}

########################################
# usage subroutine                     #
########################################

sub usage {
    print "Usage: YAC_parentage.pl [-options]/n/n";
    print "Options:\n";
    print "-c   Use Cambridge data  \n";
    print "-s   Use St Louis data   \n";
    print "-d   Debug/Verbose mode\n\n";
    print "-w   Output HTML table\n\n";
    exit;
}

########################################
# Connect with acedb database          #
########################################

my $stlacepath="/nfs/disk100/wormpub/acedb/ace4/stl";
my $camacepath="/nfs/disk100/wormpub/acedb/ace4/cam";

$|=1;

if ($ace == 1) {
    if ($debug == 1) {print  "Opening camace ....\n";}
    $db = Ace->connect(-path=>$camacepath) || do { print "Connection failure: ",Ace->error; die();};
}
elsif ($ace == 2) {
    if ($debug == 1) {print  "Opening stlace ....\n";}
    $db = Ace->connect(-path=>$stlacepath) || do { print "Connection failure: ",Ace->error; die();};
}
if ($debug == 1) {print "Connection OK.\n\n";}

########################################
# Main Loop                            #
########################################

if ($html == 1) {
    print "<TABLE WIDTH=100% BORDER=1 CELLPADDING=0>\n<TR>\n<TH>Project</TH><TH>Status</TH>";
    print "<TH>Canonical</TH><TH>EMBL_ID</TH><TH>EMBL_AC</TH><TH>2nd AC</TH><TH>Remarks</TH></TR>\n";
}

if ($debug == 1) {
    $count = $db->fetch(-query=> 'find Sequence "Y*" &! "y*.*"');
    print "checking $count YAC sequences\n\n";
}

$i = $db->fetch_many(-query=> 'find Sequence "Y*" &! "y*.*"');  
while ($obj = $i->next) {
    $project = $obj;
	
########################################
# !! NIGHTMARE !! errors in camace     #
########################################
    
    if ($project eq "Y47E12")  {next;}  # empty object -> add DB_info etc
    if ($project eq "Y47H10B") {next;}  # empty object -> add DB_info etc
    if ($project eq "Y53B1L")  {next;}  # empty object -> delete ?
    if ($project eq "Y53G8")   {next;}  # empty object -> St Louis
    if ($project eq "Y53G8A")  {next;}  # empty object -> St Louis
    if ($project =~ /^yk/)     {last;}  # erroneous empty EST objects
    
########################################
# Recognise parent YAC by '\d$'        #
########################################

    if ($project =~ /\d$/) {
	if ($html == 1) {
	    print "<TR><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;";
	    print "</TD><TD>&nbsp;</TD><TD>&nbsp;</TD></TR>\n";
	}
	else {
	    print "\n";
	}
    }

    if (($project eq "Y9C2U") && ($odd_dump == 0)) {
	if ($html == 1) {
	    print "<TR><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;</TD><TD>&nbsp;";
	    print "</TD><TD>&nbsp;</TD><TD>&nbsp;</TD></TR>\n";
	}
	else {
	    print "\n";
	}
	$odd_dump++;
    }
        
########################################
# Database information                 #
########################################
    
    ($tag1,$tag2,$db,$db_ID,$db_AC) = $obj->DB_info->row();
    
    ($db_AC_2) = $obj->Secondary_accession;
    if ($db_AC_2 eq "") {$db_AC_2 = "n/a";}
    
########################################
# Genomic_canonical                    #
########################################
    
    $canonical=$obj->Genomic_canonical;
    
########################################
# Finished                             #
########################################
    
    $finished = $obj->Finished;

########################################
# Remarks                              #
########################################
    
    @remarks = $obj->Remark(1);

########################################
# output                               #
########################################
    
    if (($finished eq "") && ($canonical ne "")) {
	if ($html == 1) {
	    print "<TR ALIGN=\"middle\"> <TD>$project</TD> <TD>Unfinished</TD> <TD>&nbsp;</TD> ";
	    print "<TD>$db_ID</TD> <TD>$db_AC</TD> <TD>$db_AC_2</TD>";
	    &filter_remarks;
	}
	else {
	    print "$project\tUnfin\t \t\t\t$db_ID  \t$db_AC  \t$db_AC_2\t";
	    &filter_remarks;
	}
    }
    else {
	if ($canonical eq "") {
	    if ($html == 1) {
		print "<TR ALIGN=\"middle\"> <TD>$project</TD> <TD>&nbsp;</TD> <TD>&nbsp;</TD> ";
		print "<TD>$db_ID</TD> <TD>$db_AC</TD> <TD>$db_AC_2</TD>";
		&filter_remarks;
	    }
	    else {
		print "$project\t\t\t\t\t$db_ID  \t$db_AC  \t$db_AC_2\t";
		&filter_remarks;
	    }
	}
	else {
	    if ($html == 1) {
		print "<TR ALIGN=\"middle\"> <TD>$project</TD> <TD>Finished</TD> <TD>Yes</TD> ";
		print "<TD>$db_ID</TD> <TD>$db_AC</TD> <TD>$db_AC_2</TD>";
		&filter_remarks;
	    }
	    else {
		print "$project\tFin\t$canonical\t$db_ID  \t$db_AC  \t$db_AC_2\t";
		&filter_remarks;
	    }
	}
    }
} # end object loop    


if ($html == 1) {
    print "</TABLE>\n";
}


sub filter_remarks {

    foreach (@remarks) {
	if (/was genomic_canonical/) {
	    $gen_can =  $_;
	}
	elsif (/replaced by/) {
	    $replaced = $_;
	}
    }

    $remark = $gen_can . " " . $replaced;
    if (($remark eq " ") && ($html == 1)) {
	$remark = "&nbsp;";
    }
     

    if ($html == 1) {
	print "<TD ALIGN=\"left\"><FONT SIZE=-1>$remark</FONT></TD> </TR>\n";
    }
    else {
	print "$remark\n";
    }
    $gen_can = ""; $replaced = "";
}


exit;

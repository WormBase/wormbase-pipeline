#!/usr/local/bin/perl

$file = shift;

# variables

$operon_num = "";                     # Operon accession number
$operon_status = "";                  # Operon status {predicted/confirmed}
$num_genes = "";                      # Number of genes in the operon
$genes = "";                          # Constituent genes in the operon (list delimited by ','


&header;

open (FILE, $file);
while (<FILE>) {

    # header line : "| OPERON CEOP001 N=3 predicted [B0205.7, B0205.6, B0205.5]                                                  |"
    if (/OP (\S+);/) {$operon_num = $1;    next;}
    if (/NG (\d+);/) {$num_genes = $1;     next;}
    if (/ST (\S+);/) {$operon_status = $1; next;}
    if (/CH (\S+);/) {$chromosome = $1;    next;}
    if (/GL (\S+.+);/) {
	$genes = $1; 
	print "<IMG SRC=\"/Projects/C_elegans/IMAGES/key.gif\" HEIGHT=90 WIDTH=90 ALIGN=right>";
	print "<H2><FONT COLOR=\"purple\">Operon $operon_num : $operon_status on Chromosome $chromosome</FONT></H2>\n";
	print "</P>";
 
	push (@out, "<TABLE WIDTH=\"100%\" CELLPADDING=\"2\" CELLSPACING=\"0\" BORDER=\"0\">");
	push (@out, "<TR BGCOLOR=\"darkblue\"><TH><FONT COLOR=\"white\">Gene</FONT></TH> <TH><FONT COLOR=\"white\">Locus</FONT></TH> <TH><FONT COLOR=\"white\">Operon position</FONT></TH> <TH><FONT COLOR=\"white\">SL2 Microarray</FONT></TH> <TH><FONT COLOR=\"white\">SL1 clones</FONT></TH> <TH><FONT COLOR=\"white\">SL2 clones</FONT></TH><TH><FONT COLOR=\"white\">Brief_ID</FONT></TH> </TR>");
	next;}
    
    if (/PO (\S+);/) {$operon_pos = $1; next;}
    if (/GN (\S+);/) {$current_gene = $1; next;}
    if (/LC (\S+);/) {$current_locus = $1; next;}
    if (/MY (\S+);/) {$current_array = $1; next;}
    if (/S1 (\d+);/) {$current_SL1 = $1; next;}
    if (/S2 (\d+);/) {$current_SL2 = $1; next;}
    if (/DE \"(\S+.+)\";/) {$current_id = $1; next;}
    if (/SL (\S+);\s+(\S+);\s+\[(\S+)\];\s+\"(\S+.+)\";/) {
	($TSL_seq,$CDS_seq) = split (/ /,$4);
	$SL_type = $1; 
	$pushline = "$current_gene $1 $2 $3 $TSL_seq $CDS_seq";
	push (@TSL, $pushline);
	next;
    }

    if (/GS (\d+)/) {$genic_distance = $1; next;}

    if (/IS (\d+)/) {
	$intergenic_distance = $1; 
	if ($SL_type eq "") {$SL_type = "SLx";}
	$span =  "$current_gene $SL_type $genic_distance $intergenic_distance";
	push (@map,$span);
	$SL_type = "";

	
	# gene details
	if ($current_locus eq "") {$current_locus = "&nbsp;";}
	if ($current_id eq "")    {$current_id = "&nbsp;";}
	
	if ($line_count % 2 == 1) {
	    push (@out, "<TR BGCOLOR=\"lightblue\" ALIGN=\"center\">");
	}
	else {
	    push (@out, "<TR ALIGN=\"center\">");
	}
	
	push (@out, "<TD><A href=\"http://wormbase.sanger.ac.uk/perl/ace/elegans/seq/sequence?name=$current_gene\">$current_gene</A></TD> <TD>$current_locus</TD> <TD>$operon_pos</TD> <TD>$current_array</TD> <TD>$current_SL1</TD> <TD>$current_SL2</TD> <TD ALIGN=\"left\">$current_id</TD></TR>\n");
	$line_count++;
	$current_locus = "";
	next;
    }

}
close (FILE);

&map;
&CDS_table;
&TSL_table;

&footer;

exit (0);

sub CDS_table {

    foreach (@out) {
	print $_;
    }
    print "</TABLE></P>\n";
}

sub TSL_table {

    $line_count = 0;
    print "<H2><FONT COLOR=\"purple\">TSL sequences</FONT></H2>\n";
    print "<TABLE WIDTH=\"100%\" CELLPADDING=\"0\" CELLSPACING=\"0\" BORDER=\"0\">\n";
    print "<TR BGCOLOR=\"darkblue\" ALIGN=\"center\"><TH><FONT COLOR=\"white\">Gene</FONT></TH> <TH><FONT COLOR=\"white\">TSL</FONT></TH> <TH><FONT COLOR=\"white\">Clone</FONT></TH> <TH><FONT COLOR=\"white\">Reference</FONT></TH> <TH><FONT COLOR=\"white\">Sequence</FONT></TH> </TR>\n";

    foreach (@TSL) {
	($current_gene,$tsl,$est,$ref,$TSL_seq,$CDS_seq) = split (/\s+/,$_);

	if ($line_count % 2 == 1) {
	    print "<TR BGCOLOR=\"lightblue\" ALIGN=\"center\">";
	}
	else {
	    print "<TR ALIGN=\"center\">";
	}	
	print "<TD ALIGN=\"center\"><FONT SIZE=\"2\">$current_gene</FONT></TD> <TD ALIGN=\"center\"><FONT SIZE=\"2\">$tsl</FONT></TD> <TD ALIGN=\"center\"><A href=\"http://www.sanger.ac.uk/srs6bin/cgi-bin/wgetz?-e+[EMBL-ACC:'$est']\"><FONT SIZE=\"2\">$est</FONT></A></TD>";


	unless ($ref =~ /EST/) {
	  print "<TD><A href=\"http://wormbase.sanger.ac.uk/perl/ace/elegans/misc/paper?name=\[$ref\]\">[<FONT SIZE=\"2\">$ref</FONT>]</A></TD>";
	}
	else {
	  print "<TD>[EST]</TD>";
	}
	
	$display_seq = substr($TSL_seq,-30);
	print " <TD ALIGN=\"left\"> <FONT COLOR=\"red\" FACE=\"courier\" SIZE=\"-1\">$display_seq</FONT> <FONT COLOR=\"black\" FACE=\"courier\">$CDS_seq</FONT></TD></TR>\n";
	$line_count++;
    }
    print "</TABLE>\n";

}

sub map {

    $sizeline = " intergenic length (bp) : [";
    
    # 5' region
    print "<IMG SRC=\"/Projects/C_elegans/IMAGES/blank.gif\" HEIGHT=25 WIDTH=50>";
    
    foreach (@map) {
	
	#print $_;
	($gene,$SL_type,$genlen,$intergene)  = split (/ /, $_);
	$genewidth = $genlen / 40;
	$intergenewidth = $intergene / 40 ;

       

	if ($intergene == 0) {

	    # CDS object
	    print "<A href=\"http://wormbase.sanger.ac.uk/perl/ace/elegans/seq/sequence?name=$gene\" onMouseOver=\"status='$gene';return true;\">";
	    if ($SL_type eq "SLx") {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	    elsif ($SL_type eq "SL1") {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene_SL1.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	    elsif ($SL_type =~ /SL2/) {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene_SL2.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	}
	else {
	    # intergenic region
	    $sizeline .=  $intergene . ", ";
	    print "<IMG SRC=\"/Projects/C_elegans/IMAGES/blank.gif\" HEIGHT=25 WIDTH=$intergenewidth>";
	    
	    # CDS object
	    print "<A href=\"http://wormbase.sanger.ac.uk/perl/ace/elegans/seq/sequence?name=$gene\" onMouseOver=\"status='$gene';return true;\">";
	    if ($SL_type eq "SLx") {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	    elsif ($SL_type eq "SL1") {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene_SL1.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	    elsif ($SL_type =~ /SL2/) {
		print "<IMG SRC=\"/Projects/C_elegans/IMAGES/gene_SL2.gif\" HEIGHT=25 WIDTH=$genewidth></A>";
	    }
	    
	}
    }
    
    # 3' region
    print "<IMG SRC=\"/Projects/C_elegans/IMAGES/blank.gif\" HEIGHT=25 WIDTH=50>";

    chop ($sizeline);
    chop ($sizeline);
    $sizeline .= "]";
    print "<PRE>$sizeline</PRE>\n\n";

}

sub header {

    print "<!--#include virtual=\"/Projects/C_elegans/menu_other.shtml\" -->\n\n";

#    print "<style>\n";
#    print "P.indent {text-indent: -1em; margin-left: 1em}\n";
#    print "</style>\n\n";


    print "<TABLE CELLSPACING=0 CELLPADDING=0 WIDTH=\"100%\" >\n";
    print "<TR VALIGN=TOP CLASS=\"h2bg\">\n";
    print "<TD WIDTH=\"100%\">\n";
    print "<CENTER>\n";
    print "<BR>\n";
    print "<H2>Operons in the <I>C.elegans</I> genome</H2></CENTER>\n";
    print "</TD>\n";
    print "</TR>\n";
    print "</TABLE>\n";
    print "</P>\n\n";

}

sub footer {
    
    print "</TABLE>\n"; 
    print "<!--#exec cgi=\"/cgi-bin/footer\" -->\n"; 
}






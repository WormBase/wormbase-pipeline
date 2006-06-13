#!/usr/local/bin/perl -wT
#author ar2
use lib "../lib";
use strict;

use vars qw($USER $PASS $DB $VALID_USERS $VALID_API_USERS $VALID_CGCNAME_USERS $SSO_USER $MAIL_NOTIFY_LIST $LIVE);

use SangerPaths qw(core);
use SangerWeb;
use NameDB_handler;
use Data::Dumper;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use Website::Utilities::Mail;
$| = 1;
$DB 	= 'wbgene_id;mcs2a';
$PASS = "wormpub";
$USER	= "wormpub";

## a list of valid  IDs to use this resources
$VALID_USERS = {
					'avc' 		=> 1,
					'ar2' 		=> 1,
					'pad' 		=> 1,
					'mt3' 		=> 1,
					'gw3' 		=> 1,
					'mh6' 		=> 1,
		
					
					'tbieri' 	=> 1,
					'jspieth' 	=> 1,
					'dblasiar' 	=> 1,
					'pozersky' 	=> 1,
					
					'stlouis' 	=> 1,
					'caltech' 	=> 1,
					'cshl' 		=> 1,
					'sanger' 	=> 1,
			   };

## a list of valid SSO login names for each DB operation
$VALID_API_USERS = {
		'query'		=> [qw(avc ar2 pad gw3 mh6 mt3 tbieri jspieth dblasiar pozersky stlouis caltech cshl sanger)],

		'merge_genes'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'split_gene'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'new_gene'		=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'kill_gene'		=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'add_name'		=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'remove_name'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'change_class'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'load_file'		=> [qw(ar2 tbieri dblasiar)],
		'dump_all'     => [qw(ar2 tbieri dblasiar)],

};

## a list of valid SSO login names able to add GCG name
$VALID_CGCNAME_USERS = {
		'mt3'			=> 1,
		'ar2'			=> 1,
};

## a list of people to mail when a DB operation occurs
$MAIL_NOTIFY_LIST = [qw(ar2 mt3)];

&main();
1;


#################################################################
sub main {
	my $web = 1;
 	my $path = SangerWeb->document_root();
	my $sw = SangerWeb->new( {
				'title'  => "WormBase NameDB",
				'banner' => "Wormbase Gene Name Management<br><font color=red> TEST ONLY - changes will not result in any changes to data in WormBase. You must make requests for new geneIDs, mergers, splits etc in the usual way</font>",
				'inifile'=> "$path/Projects/C_elegans/header.ini",
				'author' => 'avc',
				'onload' => 'init()',
			   });


	$SSO_USER = $sw->username(); ## for SSO
		
	if(!$SSO_USER) {
		my $url  = "http://$ENV{'HTTP_X_FORWARDED_HOST'}$ENV{'SCRIPT_NAME'}?$ENV{'QUERY_STRING'}";
		$sw->cookie($sw->cgi->cookie(	
									'-name'    => 'sssodestination',
                            		'-value'   => $url,
                	   				'-expires' => '+10m',
                					'-domain'  => '.sanger.ac.uk'));

		$sw->redirect("https://enigma.sanger.ac.uk/sso/login");
		$sw->redirect_delay(5);
		print $sw->header();
		print qq(<b>You need to log in to use this resource. You will shortly be redirected to a log in page...</b>);
		print $sw->footer();
		return;
	}
	
	#print @INC;
		
	print $sw->header();
	
  	if (!defined $VALID_USERS->{$SSO_USER}) {
		print qq(<h3>Sorry, you are not authorised to access this resource. Please contact the Wormbase team</h3>);
		print $sw->footer();
		return;
	} else {
		print qq(Authenticated database user: "$SSO_USER"<BR>

		);
	}

	## get CGI parameters
	my $action   = $sw->cgi->param('action') || "query";
	my $type     = $sw->cgi->param('type');
	my $name     = $sw->cgi->param('new_name') || $sw->cgi->param('delete_name');
	my $gene_id  = $sw->cgi->param('id');
	my $merge_id = $sw->cgi->param('id_2');
	my $lookup   = $sw->cgi->param('gene');
	my $user     = $sw->cgi->param('user') || "wormbase";
	my $ispublic = $sw->cgi->param('ispublic');
	my $remark   = $sw->cgi->param('remark');
	my $upfile   = $sw->cgi->param('upfile');
	my $class    = $sw->cgi->param('class');

	print_javascript();
	print_selector($action);

	if($action eq "query") {
		&query($lookup);

	} elsif($action =~ /new_gene/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&new_gene($name,$type);
		}

	} elsif($action =~ /add_name/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&add_name($name,$type,$gene_id,$ispublic);
		}

	} elsif($action =~ /kill_gene/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&kill_gene($gene_id, $remark);
		}

	} elsif($action =~ /split_gene/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&split_gene($name,$type,$gene_id);
		}

	} elsif($action =~ /merge_genes/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&merge_genes($merge_id,$gene_id);
		}
	} elsif($action =~ /remove_name/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&remove_name($name,$type,$gene_id,$ispublic);
		}

	} elsif($action =~ /load_file/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&load_file($upfile);
		}

	}  elsif($action =~ /change_class/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&change_class($name,$class);
		}
	}  elsif($action =~ /dump_all/) {
		if( is_authorised($SSO_USER,$action) == 1) {
			&dump_all('dump');
		}
	} else {
		&query();
	}

	print $sw->footer();

}
#################################################################
sub send_mail {
	my ($from, $to, $subject, $message) = @_;

	#send_mail("webserver",
	#$MAIL_NOTIFY_LIST,
	#"Merged gene $gene_id : $merge_id",
	#"GENE MERGE\nUSER : $SSO_USER\nLIVE:retained geneID $gene_id\nDEAD: killed geneID $merge_id");
	
	Website::Utilities::Mail->new({
    'to'      => $MAIL_NOTIFY_LIST,
    'from'    => $from,
    'subject' => "NAMEDB: $subject",
    'message' => $message,
	})->send(); 
}
#################################################################
sub print_javascript {

	print qq(		
		<link rel="stylesheet" type="text/css" href="/Projects/C_elegans/javascript/autosuggest.css" />
		<script type="text/javascript" src="/Projects/C_elegans/javascript/autosuggest2.js"></script>
		<script type="text/javascript" src="/Projects/C_elegans/javascript/suggestions2.js"></script>
	) if (0);

	##  var oTextbox = new AutoSuggestControl(document.getElementById("autoq"), new StateSuggestions());

	print qq(		
		<script language="javascript">
			function init()
			{

			}
			function formSubmit()
			{
				document.forms.doAction.submit()
			}
			function validate_merge()
			{
				var agree = confirm("Are you sure you want to merge these genes?")
				if (agree) {
					return(true);
				} else {
					return(false)
				}

			}
			function validate_delete()
			{
				var agree = confirm("Are you sure you want to delete this gene?")
				if (agree) {
					return(true);
				} else {
					return(false)
				}

			}

		</script>
	);
}
#################################################################
sub is_authorised {
	my ($user, $action) = @_;
	
	if (!grep {/$user/} @{$VALID_API_USERS->{$action}} ){
		print qq(
		<b>Sorry, you are not authorised to perform the following action "$action" on the gene names database</b>
		<BR><BR>
		);
		return(0);
	} else {
		return(1)
	}
}
#################################################################
sub merge_genes 
  {
    my ($merge_gene,$gene_gene) = @_;

    unless ($merge_gene && $gene_gene){
      print qq(
		<h3>Merge two existing WBGenes</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Gene to stay alive  
		<INPUT TYPE="text" NAME="id" SIZE="20" MAXLENGTH="14"<br>
		Gene to remove after merge   
		<INPUT TYPE="text" NAME="id_2" SIZE="20" MAXLENGTH="14">
		<INPUT TYPE="hidden" NAME="action" VALUE="merge_genes">
		<br><br>
		<INPUT TYPE="submit" VALUE="Merge" onClick="return validate_merge()">
		<input type="reset" value="Clear Form" />		
    	</form>
		);
      print "<BR><font color=red>Please enter genes in both fields</font>" if($merge_gene or $gene_gene);
 
    } else {
      print  "attempting to merge -$merge_gene | $gene_gene-<br>";
      my $db = get_db_connection();
		if( $db->merge_genes($gene_gene, $merge_gene) ) {
			print "Merge complete, $merge_gene is DEAD and has been merged into gene $gene_gene <br>";
			#notify
			send_mail("webserver",$MAIL_NOTIFY_LIST,"Merged gene $gene_gene : $merge_gene", "GENE MERGE\nUSER : $SSO_USER\nLIVE:retained geneID for $gene_gene\nDEAD: killed geneID $merge_gene");
      }
     	else {
			print "Sorry, the gene merge failed<br>";
      }
    }
  }


#################################################################
sub split_gene {
	my ($name,$type,$gene_id) = @_;

	unless ($type && $name && $gene_id){

		print qq(
		<h3>Split an existing WBGene</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Enter the WBGene_id of the gene to split 
		<INPUT TYPE="text" NAME="id" SIZE="13" MAXLENGTH="14">
		<BR>and the new CDS to create
		<INPUT TYPE="text" NAME="new_name" SIZE="13" MAXLENGTH="10" VALUE="new CDS">
		<INPUT TYPE="hidden" NAME="action" VALUE="split_gene">
		<INPUT TYPE="hidden" NAME="type" VALUE="CDS">
		<br><br>
		<INPUT TYPE="submit" VALUE="Split">
		</form>
		);
 
	} else {
	
		my $db = get_db_connection();
		my $new_id =$db->split_gene($name, $type, $gene_id);
		if( $new_id ){
			send_mail("webserver",$MAIL_NOTIFY_LIST,"Split gene $gene_id", "USER : $SSO_USER\nACTION : Split $gene_id\nNEW geneId : $new_id\nNEW CDS : $type $name");

			#report to screen
			$new_id = 'secret' unless ( $LIVE or defined $$VALID_CGCNAME_USERS{$SSO_USER});
			print qq(Split $gene_id creating $new_id with CDS name "$name"<br>);
		}
		else {
			print qq(ERROR : Cant Split $gene_id);		
		}
	}
}

#################################################################
sub kill_gene {
	my ($gene_id,$remark) = @_;

	unless ($gene_id and $remark){

		print qq(
		<h3>Delete an existing WBGene</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Remove the WBGene <INPUT TYPE="text" NAME="id" SIZE="20" MAXLENGTH="14">
		from the database<BR>
                Please give reason for removal (<font color=red>required</font>)<br> <INPUT TYPE="text" NAME="remark" SIZE="100" MAXLENGTH="200"><BR>
		<INPUT TYPE="hidden" NAME="action" VALUE="kill_gene">
		<br><br>
		<INPUT TYPE="submit" VALUE="Kill" onClick="return validate_delete()">
		<input type="reset" value="Clear Form" />		
		</form>
		);
 
	} else {
	
		my $db = get_db_connection();
		if ($db->kill_gene($gene_id)){
			print qq(Gene "$gene_id" has been killed because \"$remark\"<br>);
			
			send_mail("webserver",$MAIL_NOTIFY_LIST,"Kill request $gene_id", "USER :$SSO_USER\nWBGeneID : $gene_id\nAction : KILL\nRemark : $remark");
			#send to caltech to
			#send_mail("webserver",'ar2',"Kill request $gene_id", "CALTECH MAIL\nGene killed please update your annotation\n:\nUSER :$SSO_USER\nWBGeneID : $gene_id\nAction : KILL\nRemark : $remark");
		}
		else {
			print qq(FAILED to kill "$gene_id" <br>);
		}
	}
 }
#################################################################
sub remove_name {
  
    my ($name,$type,$gene_id,$ispublic) = @_;
    unless ($type && $name && $gene_id){
	
      ## print the query form
      print qq(
		<h3>Remove name from existing gene</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Remove the name 
		<SELECT NAME="type">
		  <OPTION SELECTED>CDS
		  <OPTION>CGC
		  <OPTION>other
		</SELECT>
		<INPUT TYPE="text" NAME="delete_name" SIZE="15" MAXLENGTH="10">
		from gene: 
		<INPUT TYPE="text" NAME="id" SIZE="13" MAXLENGTH="14">
		<INPUT TYPE="hidden" NAME="action" VALUE="remove_name">
		<br><br>
		<INPUT TYPE="submit" VALUE="Remove Name">
		</form>
   		);
    } 
    else {
	 	my $db = get_db_connection();
	 	$gene_id = &convert_to_geneID($db,$gene_id);
	 	print "$gene_id $name $type";
      if($db->remove_name($gene_id, $type, $name)) {
      	print "removed $name from $gene_id";
      	#print remaining names
      	$db->printAllNames($gene_id);
		   send_mail("webserver",$MAIL_NOTIFY_LIST,"Remove gene name", "$SSO_USER removed a gene name (gene $gene_id of type $type had name $name removed)");
      }
      else {
	  		print "name removal failed<br>Names for $gene_id";
		}
	 }
  }


#################################################################
sub add_name {
	my ($name,$type,$gene_id,$ispublic) = @_;

	unless ($type && $name && $gene_id){
	
		## print the query form
		print qq(
			<h3>Add name to existing gene</h3>
			<form action="$ENV{'SCRIPT_NAME'}" method="GET">
			<BR><BR>
			Enter the details of the additional gene name<br><br>  
			Name type to be added: <SELECT NAME="type">
			  <OPTION SELECTED>CDS
		);
		## only specified users can add a CGC name
		if (exists $VALID_CGCNAME_USERS->{$SSO_USER}){
		  print qq(
			  <OPTION>CGC
			);
		}
		print qq(
			</SELECT><br>
			Name to add: <INPUT TYPE="text" NAME="new_name" SIZE="15" MAXLENGTH="10" VALUE=""><br>
			for gene (WBGeneID or Sequence_name): <INPUT TYPE="text" NAME="id" SIZE="13" MAXLENGTH="14" VALUE=""><br>
			<INPUT TYPE="hidden" NAME="action" VALUE="add_name">
			<br><br>
    		<INPUT TYPE="submit" VALUE="Add Name">
		);

	} else {
	
		## if people try to hack the GET URL parameters
		if ( ($type eq "CGC") and (!exists $VALID_CGCNAME_USERS->{$SSO_USER}) ){
			print qq(
			<b>Sorry, you are not authorised to add a GCG gene name to the database</b>
			<BR><BR>
			);
		} else {	
			my $db = get_db_connection();
		 	$gene_id = &convert_to_geneID($db,$gene_id);
			if( $db->add_name($gene_id,$name, $type) ) {
				print qq(Added "$name" as $type name for gene $gene_id<br>);
				send_mail("webserver",$MAIL_NOTIFY_LIST,"$type added to $gene_id", "USER : $SSO_USER\nWBGeneID : $gene_id\nName added : $name $type\n");
			}
			else {
				print "<p>FAILED: cant add name $name to $gene_id</p>";
			}
		}
	}
}
#################################################################
sub new_gene {
	my ($name,$type) = @_;

	#print "type: $type, name: $name";
	
	unless ($type && $name){
	
		## print the query form
		print qq(
    	<h3>Make a new gene</h3>
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Enter the details of the new gene<br>  
		<SELECT NAME="type">
		  <OPTION SELECTED>CDS
                );	
		if (exists $VALID_CGCNAME_USERS->{$SSO_USER}){
		  print qq(
			  <OPTION>CGC
			);
		}
		print qq(
		</SELECT>
		<INPUT TYPE="text" NAME="new_name" SIZE="15" MAXLENGTH="10" VALUE="">
		<INPUT TYPE="hidden" NAME="action" VALUE="new_gene">
		<br><br>
    	<INPUT TYPE="submit" VALUE="Create">
    	</form>
		);

	} else {
	
		my $db = get_db_connection();
		
		if(my $id = $db->new_gene($name, $type) ) {
			send_mail("webserver",$MAIL_NOTIFY_LIST, " WBGeneID request $name $SSO_USER","User : $SSO_USER\nType : $type\nName : $name\nID : $id");
			#report to screen
			$id = 'secret' unless ( $LIVE or defined $$VALID_CGCNAME_USERS{$SSO_USER});
			print "The new id for $type: $name is $id\n";
			print qq(
			<hr align="left">
			<BR><BR><BR>
			);
		}
		else {
			print "<p>FAILED : cant create GeneID for $type $name</p>";
		}
	}
}

#################################################################
sub query {
	my $lookup = shift;
	unless ($lookup){
		## print the query form
		print qq(
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
    	  <h3>Retreive gene info</h3>
    	  Gene to retreive<br>  
    	  <INPUT TYPE="text" NAME="gene" SIZE="15" MAXLENGTH="14" VALUE="">
    	  <INPUT TYPE="hidden" NAME="action" VALUE="query">
		<br><br>
    	  <INPUT TYPE="submit" VALUE="Search">
    	</form>
		);
	
	} else {
		
		##run the query
		my $db = get_db_connection();
		
		my $gene = $lookup;
		#$db->validate_name($gene, $type);
		($gene) = $db->idGetByAnyName($lookup);
		unless ($gene and $gene =~ /\w+/) {
		  print($lookup." does not exist in database<br>");
		  return;
		}
		
		print qq(
		The Gene Name database currently holds the following information about $lookup:
		<BR><BR><BR>
		);
		$db->print_history($gene);
		$db->printAllNames($gene);
		print qq(
		<hr align="left">
		<BR><BR><BR>
		);
	}

}


sub load_file {
	my $file = shift;
	unless ( $file ) {
	## print the query form
		print qq(
    	<form action="$ENV{'SCRIPT_NAME'}" method="POST" enctype="multipart/form-data" >
		<BR><BR>
    	  <h3>Upload a file</h3>
    	  File to upload: <INPUT TYPE=FILE NAME="upfile"><BR>
			<INPUT TYPE=SUBMIT VALUE="Submit">
			<INPUT TYPE="hidden" NAME="action" VALUE="load_file">
    	</form>
		);		
	}
	else {
		my $db = get_db_connection();
		#untaint file
#		unless( $file =~ m/^(.+)$/ ) {
#//			$db->dienice("tainted file\n");
#//		}
#//		my $tfile = $1;
#		open(UL, "<$file") or $db->dienice("cant open file",$!);
		print "<TABLE border=1><THEAD>
			<TR>
				<TH>Gene_id
				<TH>CGC
				<TH>Sequence
				<TH>CDS(s)
			</THEAD>
			<TBODY>";
		open FILE,"<$file" or die "cant open $file : $!\n";
		while ( <FILE> ) {
			chomp;
			my $id = $_;
			next unless /WBGene\d{8}/;
			my( $cgc, $seq, @isoforms);
			my %names = $db->idAllNames($_);
			$cgc = ($names{'CGC'} or '-');
			$seq = ($names{'Sequence'} or '-');
			
			#CDS may be single scalar or array of isoforms
			if( ref($names{'CDS'}) eq 'ARRAY' ) {
				@isoforms = @{$names{'CDS'}};
			}
			else {
				push(@isoforms, $names{'CDS'});
			}
			print "<TR><TD>$id<TD>$cgc<TD>$seq<TD>",join(", ",@isoforms);
		}
		print "</TBODY></TABLE>";
		close FILE;
	}
 }
 
sub change_class {
	my ($cds, $class) = @_;
	unless( $cds and $class) {
		## print the query form
		print qq(
    	<h3>Change CDS class</h3>
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Enter the details of the CDS to change<br>Convert CDS 
		<INPUT TYPE="text" NAME="new_name" SIZE="15" MAXLENGTH="10" VALUE="">
		<p>Change to . . <p>
		<input type="radio" name="class" value="CDS"> CDS<br>
		<input type="radio" name="class" value="pseudogene" checked> Pseudogene<br>
		<input type="radio" name="class" value="transcript"> non-coding transcript<br>
		<input type="radio" name="class" value="transposon"> transposon<br>
		
		
		<INPUT TYPE="hidden" NAME="action" VALUE="change_class">
		<br><br>
    	<INPUT TYPE="submit" VALUE="Change Class">
    	</form>
		);

	}
	else {
		my $db = get_db_connection();
		my $gene = $db->idGetByTypedName('CDS', $cds);
		if ( $gene ) {
			# just send a mail - no change to database
			#$self->send_mail('Nameserver',
			send_mail("webserver",$MAIL_NOTIFY_LIST,"$cds to $class", "$cds converted to $class - $SSO_USER\nThis does not affect Nameserver DB");
		}
		print "changing $cds to $class<br>";
	}
}


sub dump_all {
	my $dump = shift;
	if ( $dump ) {
		# iterate over all genes and print details
		my $db = get_db_connection();
#		my $all = $db->allLiveIds;
#		my $tmpfile ="/tmp/nameDB.$$";
#		open (TMP,">$tmpfile") or croak ("cant open file $tmpfile\t:$!");
#		foreach my $gene (@$all) {
#			print TMP "$gene\n";
#		}
#		close TMP;
#		&load_file($tmpfile);
		my $query =<<END;
SELECT primary_identifier.object_public_id, name_type_name, secondary_identifier.object_name 
FROM primary_identifier, secondary_identifier, name_type 
WHERE secondary_identifier.object_id = primary_identifier.object_id 
		AND secondary_identifier.name_type_id = name_type.name_type_id 
		AND primary_identifier.object_live = 1
ORDER BY object_public_id;
END
		my $sth = $db->dbh->prepare($query);
		$sth->execute or die "Unable to execute query: $db->errstr\n";
		my $row;
		while($row = $sth->fetchrow_arrayref) {
    		print "$row->[0]\t$row->[1]\t$row->[2]<br>";
    	}

		$sth->finish;
	}
	else {
		print qq( 
			<h3>Print page of all live genes</h3>
			<form action="$ENV{'SCRIPT_NAME'}" method="GET">
			<INPUT TYPE="hidden" NAME="action" VALUE="dump_all">
			<br><br>
			<INPUT TYPE="submit" VALUE="Dump"
		</form>
		);
	}
}
#################################################################
sub print_selector {
	my ($action) = @_;

	my $qs = "";
	my $ng = "";
	my $an = "";
	my $kg = "";
	my $sg = "";
	my $mg = "";
	my $rn = "";
	my $lf = ""; my $cc = ""; my $da = "";
	if ($action eq "query") 		{ $qs = " selected" };
	if ($action eq "new_gene") 	{ $ng = " selected" };
	if ($action eq "add_name") 	{ $an = " selected" };
	if ($action eq "kill_gene") 	{ $kg = " selected" };
	if ($action eq "split_gene") 	{ $sg = " selected" };
	if ($action eq "merge_genes")	{ $mg = " selected" };
	if ($action eq "remove_name")	{ $rn = " selected" };
	if ($action eq "load_file") 	{ $lf = " selected" };
	if ($action eq "change_class"){ $cc = " selected" };
	if ($action eq "dump_all") 	{ $da = " selected" };


	print qq(
    What do you want to do ?  <br><br>
    <form name="doAction" action="$ENV{'SCRIPT_NAME'}" method="GET">
    Select action to perform:
  	<SELECT NAME="action" onChange="formSubmit()">
    	<OPTGROUP LABEL="Gene information">
    	  <OPTION $qs LABEL="Find gene information" value="query">Search for gene information</OPTION>
    	</OPTGROUP>
	);

	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	<OPTGROUP LABEL="Request a new gene ID">
    	  <OPTION $ng LABEL="new_gene" value="new_gene">Request a new WBGeneID</OPTION>
    	</OPTGROUP>
		);
	}

	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	<OPTGROUP LABEL="Delete an existing gene">
    	  <OPTION $kg LABEL="kill_gene" value="kill_gene">Kill a WBGeneID</OPTION>
    	</OPTGROUP>
		);
	}

	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	<OPTGROUP LABEL="Update gene names">
    	  <OPTION $an LABEL="add_name" value="add_name">Add a gene name</OPTION>
    	  <OPTION $rn LABEL="remove_name" value="remove_name">Remove a gene name</OPTION>
    	</OPTGROUP>
		);
	}

	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	<OPTGROUP LABEL="Merge/split genes">
    	  <OPTION $mg LABEL="merge_genes" value="merge_genes">Merge two genes</OPTION>
    	  <OPTION $sg LABEL="split_gene" value="split_gene">Split a gene</OPTION>
    	</OPTGROUP>
		);
	}
	
	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	<OPTGROUP LABEL="Load file">
    	  <OPTION $lf LABEL="load file" value="load_file">Load a file</OPTION>
    	</OPTGROUP>
		);
	}
		
	if (grep {/$SSO_USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
   	<OPTGROUP LABEL="Miscellaneous">
  		   <OPTION $cc LABEL="change class" value="change_class">Change class</OPTION>
  		   <OPTION $da LABEL="dump all" value="dump_all">Dump all</OPTION>
   	</OPTGROUP>
		);
	}		
	print qq(
    </SELECT>
    <INPUT TYPE="submit" VALUE="Submit">
    </form>
	<BR><BR>
	<HR align="left">
	);
}

#make sure gene names are correct case
sub fix_case
  {
    my ($type, $name) = @_;
    if( $type eq 'CDS' ) {
      my @data = split(/\./,$name);
      return (uc $data[0]).".$data[1]";
    }
    elsif( $type eq 'CGC' ) {
      return lc $name;
    }
    elsif( $type eq 'Gene' ){
      /(\d+$)/;
      return "WBGene$1";
    }
    else {
      #$db->dienice($type is)
    }
  }

#################################################################
sub get_db_connection {
	
	my $DOMAIN = 'Gene';
	my $db = NameDB_handler->new($DB,$USER,$PASS,1); #1 is for web output
	
	$db || return(undef);
	$db->setDomain($DOMAIN);
	
	set_web_output($db); # turn on web reporting
	return($db);
}

#################################################################
sub set_web_output {
	my ($conn) = @_;
	$conn->web(1); # turn on web reporting
	
}
#################################################################

#This is to allow the web page to accept any name for a gene but pass the geneID to NameDB_handler.
sub convert_to_geneID {
	my $db = shift;
	my $id = shift;
	my $id_search = $db->idGetByAnyName($id)->[0];
	my $gene_id = $id_search ? $id_search : $id; #allow any name type to be added to
	return $gene_id;
}








#!/usr/local/bin/perl -wT
#author ar2

use lib "../lib";
use strict;

use vars qw($USER $PASS $DB $VALID_USERS $VALID_API_USERS $VALID_CGCNAME_USERS $SSO_USER $MAIL_NOTIFY_LIST $LIVE);

use SangerPaths qw(core);
use SangerWeb;
use NameDB_handler;
use Data::Dumper;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Mail::Mailer;

$| = 1;
$DB 	= 'wbgene_id;mcs2a';
$PASS 	= "wormpub";
$USER 	= "wormpub";

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
		'new_gene'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'kill_gene'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'add_name'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'remove_name'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],

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

	} else {
		&query();
	}


	
	print $sw->footer();

}
#################################################################
sub send_mail {
	my ($from, $to, $subject, $message) = @_;

	my $server = 'mail.sanger.ac.uk';
	my $mailer = new Mail::Mailer 'smtp', Server => $server;
    #$mailer = new Mail::Mailer;
    #$mailer = new Mail::Mailer $type, @args;

	my $headers = {
		'To'		=>	$to,
		'From'		=>	$from,
		'Subject'	=>	"NameDB: " . $subject,
	};

    $mailer->open($headers);
    print $mailer $message;
    $mailer->close;
	
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
      print  "-$merge_gene | $gene_gene-";
      my $db = get_db_connection();
      my ($gene_id, $merge_id);
      $gene_id  = ($db->idGetByAnyName($gene_gene)->[0]  or $gene_gene);  #allow use of any name
      $merge_id = ($db->idGetByAnyName($merge_gene)->[0] or $merge_gene); #allow use of any name
      $db->validate_id($gene_id);
      $db->validate_id($merge_id);
		
      #remove all gene names for eaten gene
      $db->remove_all_names($merge_id);
	
      #do the merge
      if ($db->idMerge($merge_id,$gene_id)) {
	print "Merge complete, $merge_gene ($merge_id) is DEAD and has been merged into gene $gene_gene ($gene_id)<br>";
	#notify
	send_mail("webserver",$MAIL_NOTIFY_LIST,"Merged gene $gene_id : $merge_id", "GENE MERGE\nUSER : $SSO_USER\nLIVE:retained geneID $gene_id\nDEAD: killed geneID $merge_id");
      } else {
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
		$db->validate_id($gene_id);
		$db->validate_name($name, $type);
		$db->check_pre_exists($name, $type);
		my $id = $db->idSplit($gene_id);		## hmm, does this work?
		#$type = "CDS";
		$db->add_name($id, $name, $type);

		#mail
		send_mail("webserver",$MAIL_NOTIFY_LIST,"Split gene $gene_id", "USER : $SSO_USER\nACTION : Split $gene_id\nNEW geneId : $id\nNEW CDS : $type $name");

		#report to screen
		$id = 'secret' unless ( $LIVE or defined $$VALID_CGCNAME_USERS{$SSO_USER});
		print qq(Split $gene_id creating $id with CDS name "$name"<br>);
		
		
	}
}
#################################################################
sub kill_gene {
	my ($gene_id,$remark) = @_;

	unless ($gene_id){

		print qq(
		<h3>Delete an existing WBGene</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Remove the WBGene <INPUT TYPE="text" NAME="id" SIZE="20" MAXLENGTH="14">
		from the database<BR>
                Please give reason for removal<br> <INPUT TYPE="text" NAME="remark" SIZE="100" MAXLENGTH="200"><BR>
		<INPUT TYPE="hidden" NAME="action" VALUE="kill_gene">
		<br><br>
		<INPUT TYPE="submit" VALUE="Kill" onClick="return validate_delete()">
		<input type="reset" value="Clear Form" />		
		</form>
		);
 
	} else {
	
		my $db = get_db_connection();
		$db->validate_id($gene_id);
		if ($db->idKill($gene_id)){
			print qq(Gene "$gene_id" has been killed because \"$remark\"<br>)
		}
		
		send_mail("webserver",$MAIL_NOTIFY_LIST,"Kill request $gene_id", "USER :$SSO_USER\nWBGeneID : $gene_id\nAction : KILL\nRemark : $remark");
		
	}
 }
#################################################################
sub remove_name
  {
  
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

    } else {
	
      my $db = get_db_connection();
      $db->validate_id($gene_id);
      $db->validate_name($name, $type);
      my ($exist_id) = $db->idGetByTypedName($type,$name);
      if ( "$exist_id" ne "$gene_id" ) {
	$db->dienice("$name is not a name of $gene_id\n");
      }
      if ( $db->delName($gene_id,$type,$name) ) {
	print "Removed $name as name for $gene_id<br><br>Remaining names are <br>";
	#need to update the public name
	if ($type = 'CGC') {
	  my $seq_name = $db->idTypedNames($gene_id,'Sequence');
	  $db->add_name($gene_id,$seq_name->[0],'Public_name');
	} else {
	  print "name removal failed<br>Names for $gene_id";
	}
      }
      #print remaining names
      $db->printAllNames($gene_id);
		
      send_mail("webserver",$MAIL_NOTIFY_LIST,"Remove gene name", "$SSO_USER removed a gene name (gene $gene_id of type $type had name $name removed)");
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
			$db->validate_name($name, $type);
			unless( $gene_id =~ /^WBGene/) {
			  my $ids = $db->idGetByTypedName('Sequence',$name);
			  $db->dienice("$gene_id is not a name of a valid gene\n") unless $ids->[0];
			  $gene_id = $ids->[0];
			}			  
			$db->validate_id($gene_id);
			$db->check_pre_exists($name, $type);
			$db->add_name($gene_id,$name, $type);
			$db->add_name($gene_id,$name, 'Public_name') if ($type eq 'CGC');
			print qq(Added "$name" as $type name for gene $gene_id<br>);

			send_mail("webserver",$MAIL_NOTIFY_LIST,"$type added to $gene_id", "USER : $SSO_USER\nWBGeneID : $gene_id\nName added : $name $type\n");
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
		
		$db->validate_name($name, $type);
		$db->check_pre_exists($name, $type);
		if( my $id = $db->isoform_exists($name, $type) ) {
		  $db->dienice("$name is an isoform of $id<BR>Please confirm and add $name as another name for $id<BR>")
		};
		my $id = $db->make_new_obj($name, $type);

		send_mail("webserver",$MAIL_NOTIFY_LIST, " WBGeneID request $name $SSO_USER","User : $SSO_USER\nType : $type\nName : $name\nID : $id");

		#report to screen
		$id = 'secret' unless ( $LIVE or defined $$VALID_CGCNAME_USERS{$SSO_USER});
		print "The new id for $type: $name is $id\n";
		print qq(
		<hr align="left">
		<BR><BR><BR>
		);
		
		
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
		unless ($gene =~ /\w+/) {
		  $db->dienice($lookup." does not exist in database<br>");
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
	if ($action eq "query") 		{ $qs = " selected" };
	if ($action eq "new_gene") 		{ $ng = " selected" };
	if ($action eq "add_name") 		{ $an = " selected" };
	if ($action eq "kill_gene") 	{ $kg = " selected" };
	if ($action eq "split_gene") 	{ $sg = " selected" };
	if ($action eq "merge_genes") 	{ $mg = " selected" };
	if ($action eq "remove_name") 	{ $rn = " selected" };


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
    elsif( $type = 'Gene' ){
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
	my $db = NameDB_handler->new($DB,$USER,$PASS);
	
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










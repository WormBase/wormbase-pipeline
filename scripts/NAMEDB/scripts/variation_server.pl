#!/usr/local/bin/perl -wT
#author ar2
use lib "../lib";
use strict;

use vars qw($SSO_USER $PASS $DB $VALID_USERS $VALID_API_USERS $VALID_CGCNAME_USERS $USER $MAIL_NOTIFY_LIST $MAILS $LIVE);

use SangerPaths qw(core);
use SangerWeb;
use NameDB_handler;
use Data::Dumper;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use Website::Utilities::Mail;
use File::Path;
$| = 1;

# determine if live or WWWdev and select database accordingly
if( -e '/nfs/WWWdev/SANGER_docs/htdocs'){
	$DB = 'test_wbgene_id;mcs2a';
}
else {
	$DB = 'wbgene_id;mcs2a';
}

## a list of valid  IDs to use this resources
$VALID_USERS = {
					# these are the users WBPerson id
					'ar2' 		=> 1847,
					'pad' 		=> 1983,
					'mt3' 		=> 2970,
					'gw3' 		=> 4025,
					'mh6' 		=> 4036,
					'tbieri' 	=> 1849,
					'jspieth' 	=> 615,
					'dblasiar' 	=> 1843,
					'pozersky' 	=> 1867,
					'stlouis' 	=> 1,
					'caltech' 	=> 1,
					'cshl' 		=> 1,
					'sanger' 	=> 1,
			   };

## a list of valid SSO login names for each DB operation
$VALID_API_USERS = {
		'query'		=> [qw(avc ar2 pad gw3 mh6 mt3 tbieri jspieth dblasiar pozersky stlouis caltech cshl sanger)],
   	'merge_var'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'new_var'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
		'kill_var'	=> [qw(avc ar2 pad gw3 mt3 tbieri jspieth dblasiar pozersky)],
};

## a list of valid SSO login names able to add GCG name
$VALID_CGCNAME_USERS = {
		'mt3'			=> 1,
		'ar2'			=> 1,
};

$MAILS = {
	'ar2'			=>	'ar2@sanger.ac.uk',
	'pad'			=>	'pad@sanger.ac.uk',
	'gw3'			=>	'gw3@sanger.ac.uk',
	'mh6'			=>	'mh6@sanger.ac.uk',
	'mt3'			=>	'mt3@sanger.ac.uk',
	'tbieri'		=>	'tbieri@watson.wustl.edu',
	'jspieth'	=>	'jspieth@watson.wustl.edu',
	'dblasiar'	=>	'dblasiar@watson.wustl.edu',
	'pozersky'	=>	'pozersky@watson.wustl.edu',
	'stlouis'	=>	'stlouis@wormbase.org',
	'caltech'	=>	'caltech@wormbase.org',
	'cshl'		=>	'cshl@wormbase.org',
	'sanger'		=>	'sanger@wormbase.org',
	'cgc'       => 'mt3@sanger.ac.uk'
};

## a list of people to mail when a DB operation occurs
$MAIL_NOTIFY_LIST = [qw(ar2)];

&main();
1;


#################################################################
sub main {
	my $web = 1;
 	my $path = SangerWeb->document_root();
	my $sw = SangerWeb->new( {
				'title'  => "WormBase Variation ID Server $DB",
				'banner' => "Wormbase Variation Name Management<br><font color=red>This system is TEST only</font>",
				'inifile'=> "$path/Projects/C_elegans/header.ini",
				'author' => 'avc',
				'onload' => 'init()',
			   });


	$SSO_USER = $sw->username(); ## for SSO
	if( $SSO_USER =~ /^(\w+)@/) {
		$PASS = $1;
	}
	else {
		$PASS = $SSO_USER;
	}
	$USER = $PASS;
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
			
	print $sw->header();
	
  	if (!defined $VALID_USERS->{$USER}) {
		print qq(<h3>Sorry, you $USER are not authorised to access this resource. Please contact the Wormbase team</h3>);
		print $sw->footer();
		return;
	} else {
		print qq(Authenticated database user: "$USER"<BR>		);
		#send_mail('webserver', $MAILS->{'ar2'}, "user $USER", "is using nameserver");
	}

	## get CGI parameters
	my $action   = $sw->cgi->param('action')|| "query";
	my $name     = $sw->cgi->param('name');
	my $varid    = $sw->cgi->param('varid');
	my $lookup   = $sw->cgi->param('lookup');
	my $keep_id   = $sw->cgi->param('keep_id');
	my $kill_id   = $sw->cgi->param('kill_id');


	print_javascript();
	print_selector($action);

	if($action eq "query") {
		&query($lookup);

	} elsif($action =~ /new_var/) {
		if( is_authorised($USER,$action) == 1) {
			&new_var($name);
		}

	} elsif($action =~ /merge/) {
		if( is_authorised($USER,$action) == 1) {
			&merge($keep_id, $kill_id);
		}
	} elsif( $action =~ /kill_var/) {
		if( is_authorised($USER,$action) == 1) {
			&kill_var($kill_id);
		}
	} else {
		&query();
	}

	print $sw->footer();

}
#################################################################
sub send_mail {
	my ($from, $to, $subject, $message,) = @_;

	Website::Utilities::Mail->new({
    'to'      => $to,
    'from'    => $from,
    'subject' => ($DB =~ 'test' ? 'TEST : ' : " ")."NAMEDB: $subject",
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
				var agree = confirm("Are you sure you want to merge these ?")
				if (agree) {
					return(true);
				} else {
					return(false)
				}

			}
			function validate_delete()
			{
				var agree = confirm("Are you sure you want to delete this ?")
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
		<b>Sorry, you are not authorised to perform the following action "$action" on the variation names database</b>
		<BR><BR>
		);
		return(0);
	} else {
		return(1)
	}
}
#################################################################
sub merge
  {
    my ($keep,$kill) = @_;

    unless ($keep && $kill){
      print qq(
		<h3>Merge two existing Variations</h3>
		<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
		Variation to stay alive  
		<INPUT TYPE="text" NAME="keep_id" SIZE="20" MAXLENGTH="14"<br>
		Variation to remove after merge   
		<INPUT TYPE="text" NAME="kill_id" SIZE="20" MAXLENGTH="14">
		<INPUT TYPE="hidden" NAME="action" VALUE="merge_var">
		<br><br>
		<INPUT TYPE="submit" VALUE="Merge" onClick="return validate_merge()">
		<input type="reset" value="Clear Form" />		
    	</form>
		);
      print "<BR><font color=red>Please enter variation ids in both fields</font>" if($keep or $kill);
 
    } else {
      print  "attempting to merge $keep | $kill<br>";
      my $db = get_db_connection();
      #get ids if names passed
      
      my ($keep_id, $kill_id);
      $keep_id = ($keep =~ /WBVariation\d+/) ? $keep : $db->idGetByAnyName("$keep")->[0];
      $kill_id = ($kill =~ /WBVariation\d+/) ? $kill : $db->idGetByAnyName("$kill")->[0];
      
      #check ids are valid
      foreach ($keep_id, $kill_id) {
      	unless( defined $db->idExists($_) ) {
      		print qq(Variation does not exist<br>);
      		exit(0);
      	}
      }
      
      #do the merge
		if( my $id = $db->idMerge($kill_id,$keep_id) ) {
			print "Merge complete, $kill_id is DEAD and has been merged into variation $keep_id <br>";
			#notify
			send_mail("webserver",[$MAILS->{$USER},$MAILS->{'cgc'}],"Merged Variations $keep ($keep_id) and $kill ($kill_id)", "VARIATION MERGE\nUSER : $USER\nLIVE:retained WBVarID for $keep_id\nDEAD: killed VarID $kill_id \n");
      }
     	else {
			print "Sorry, the variation merge failed<br>";
      }
    }
  }
  
sub new_var {
	my $public = shift;
	unless ($public) {
		## print the query form
		print qq(
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
    	  <h3>Request a new Variation Id</h3>
    	  Enter Public name of Variation<br>  
    	  <INPUT TYPE="text" NAME="name" SIZE="10" MAXLENGTH="10" VALUE="">
    	  <INPUT TYPE="hidden" NAME="action" VALUE="new_var">
		<br><br>
    	  <INPUT TYPE="submit" VALUE="Search">
    	</form>
		);
	}
	else {
		my $db = get_db_connection();
  		my $var = $db->idGetByTypedName('CGCvar'=>$public)->[0];
		if ( $var ) {	
			print "$public already exists as $var"
		}
		else {
			my $id = $db->idCreate;
			$db->addName($id,'CGCvar'=>$public);
			print "$id created with name $public<br>";
		}
	}
}
  
sub kill_var {
	my $kill_id = shift;
	unless ($kill_id) {
		## print the query form
		print qq(
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
    	  <h3>Kill a Variation Id</h3>
    	  Enter Public name or Id of Variation<br>  
    	  <INPUT TYPE="text" NAME="kill_id" SIZE="10" MAXLENGTH="10" VALUE="">
    	  <INPUT TYPE="hidden" NAME="action" VALUE="kill_var">
		<br><br>
    	  <INPUT TYPE="submit" VALUE="Search">
    	</form>
		);
	}
	else {
		my $db = get_db_connection();
		my $death = ($kill_id =~ /WBVariation\d+/) ? $kill_id : $db->idGetByAnyName("$kill_id")->[0];
		if ( $death ) {	
			if( $db->idKill($death) ) {
				print qq( $death killed<br>);
			}
		}
		else {
			print "$kill_id doesn't exists"
		}
	}
}			
  		
  		
sub query {
	my $lookup = shift;
	unless ($lookup){
		## print the query form
		print qq(
    	<form action="$ENV{'SCRIPT_NAME'}" method="GET">
		<BR><BR>
    	  <h3>Retreive Variation info</h3>
    	  Variation to retreive<br>  
    	  <INPUT TYPE="text" NAME="lookup" SIZE="21" MAXLENGTH="20" VALUE="">
    	  <INPUT TYPE="hidden" NAME="action" VALUE="query">
		<br><br>
    	  <INPUT TYPE="submit" VALUE="Search">
    	</form>
		);
	
	} else {
		
		##run the query
		my $db = get_db_connection();
		my $var = $db->idGetByAnyName($lookup)->[0];
		if( $var ) {		
			print qq(
			The Variation Name database currently holds the following information about $lookup:
			<BR><BR><BR>
			);
			&printAllNames($db,$var);
			print qq(
			<hr align="left">
			<BR><BR><BR>
			);
			send_mail("webserver",[$MAILS->{$USER}],"query", "$USER queried for $var");
		}
		else {
			print qq( $lookup is not known to the Variation nameserver<br>);
		}
	}

}




#################################################################
sub print_selector {
	my ($action) = @_;

	my $nv = "";
	my $qv = "";
	my $mv = "";
	my $kv = "";
	if ($action eq "query") 		{ $qv = " selected" };
	if ($action eq "new_var") 		{ $nv = " selected" };
	if ($action eq "merge_var")   { $mv = " selected" };
	if ($action eq "kill_var") 	{ $kv = " selected" };

	print qq(
    What do you want to do ?  <br><br>
    <form name="doAction" action="$ENV{'SCRIPT_NAME'}" method="GET">
    Select action to perform:
  	<SELECT NAME="action" onChange="formSubmit()">
    	<OPTGROUP LABEL="Variation information">
    	  <OPTION $qv LABEL="Find variation information" value="query">Search for variation information</OPTION>
	);

	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
		  <OPTION $nv LABEL="new_var" value="new_var">Request a new WBVarID</OPTION>
		);
	}
	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	 <OPTION $mv LABEL="merge_var" value="merge_var">Merge two variations Ids</OPTION>
		);
	}
	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	 <OPTION $kv LABEL="kill_var" value="kill_var">Kill a Variation Id</OPTION>
		);
	}
	print qq(    
		</OPTGROUP>
    </SELECT>
    <INPUT TYPE="submit" VALUE="Submit">
    </form>
	<BR><BR>
	<HR align="left">
	);
}

#################################################################
sub get_db_connection {
	
	my $DOMAIN = 'Variation';
	my $db = NameDB->connect($DB,$USER,$PASS,1); #1 is for web output
	
	$db || return(undef);
	$db->setDomain($DOMAIN);
	
	#set_web_output($db); # turn on web reporting
	return($db);
}

#################################################################
sub set_web_output {
	my ($conn) = @_;
	$conn->web(1); # turn on web reporting
	
}
#################################################################



sub printAllNames
  {
    my $db = shift;
    my $obj = shift;
    my $names = $db->idAllNames($obj);
    print "<br>Current names for <b>$obj</b> <br>";
    foreach (keys %{$names} ) {
      my $name_str;
      if (ref $names->{"$_"} eq 'ARRAY') {
			$name_str =  join(" ",@{$names->{"$_"}})
      }else {
			$name_str = $names->{$_};
      }
      print "$_ : $name_str<br>";
    }
  }





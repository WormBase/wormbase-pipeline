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
		'bastiani'  => 1,
	       };

## a list of valid SSO login names for each DB operation
$VALID_API_USERS = {
		    'query'		=> [qw(avc ar2 pad gw3 mh6 mt3 tbieri bastiani jspieth dblasiar pozersky stlouis caltech cshl sanger)],
		    'merge_features'	=> [qw(ar2 pad gw3 mt3 mh6 tbieri jspieth dblasiar pozersky)],
		    'new_feature'	=> [qw(ar2 pad gw3 mt3 mh6 tbieri bastiani jspieth dblasiar pozersky)],
		    'kill_feature'	=> [qw(ar2 pad gw3 mt3 mh6 tbieri jspieth dblasiar pozersky)],
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
	  'cgc'       => 'mt3@sanger.ac.uk',
	  'bastiani'  => 'bastiani@its.caltech.edu'
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
				
			    'banner' => "Wormbase Feature Name Management<br><font color=red>This system is TEST only</font>",
			    'inifile'=> "$path/Projects/C_elegans/header.ini",
			    'author' => 'avc',
			    'onload' => 'init()',
			   });
  if ($sw->is_dev()) {
    $DB = 'test_wbgene_id;mcs2a';
  } else {
    $DB = 'test_wbgene_id;mcs2a';
  }

  $SSO_USER = $sw->username();	## for SSO
  if ( $SSO_USER =~ /^(\w+)@/) {
    $PASS = $1;
  } else {
    $PASS = $SSO_USER;
  }
  $USER = $PASS;
  if (!$SSO_USER) {
    my $url  = "http://$ENV{'HTTP_X_FORWARDED_HOST'}$ENV{'SCRIPT_NAME'}?$ENV{'QUERY_STRING'}";
    $sw->cookie($sw->cgi->cookie(	
				 '-name'    => 'sssodestination',
				 '-value'   => $url,
				 '-expires' => '+10m',
				 '-domain'  => '.sanger.ac.uk'));

    $sw->redirect("https://enigma.sanger.ac.uk/sso/login");
    $sw->redirect_delay(5);
    print $sw->header({'title'  => "WormBase Feature ID Server $DB"});
    print qq(<b>You need to log in to use this resource. You will shortly be redirected to a log in page...</b>);
    print $sw->footer();
    return;
  }
			
  print $sw->header({'title'  => "WormBase Feature ID Server $DB"});
	
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
  my $new      = $sw->cgi->param('new');
  my $keep_id  = $sw->cgi->param('keep_id');
  my $kill_id  = $sw->cgi->param('kill_id');
  my $remark   = $sw->cgi->param('remark');

  print_javascript();
  print_selector($action);

  if ($action eq "query") {
    &query();

  } elsif ($action =~ /new_feature/) {
    if ( is_authorised($USER,$action) == 1) {
      &new_feature($new);
    }

  } elsif ($action =~ /merge_features/) {
    if ( is_authorised($USER,$action) == 1) {
      &merge_features($keep_id, $kill_id);
    }
  } elsif ( $action =~ /kill_feature/) {
    if ( is_authorised($USER,$action) == 1) {
      &kill_feature($kill_id, $remark);
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
	
    if (!grep {/$user/} @{$VALID_API_USERS->{$action}} ) {
      print qq(
	       <b>Sorry, you are not authorised to perform the following action "$action" on the feature database</b>
	       <BR><BR>
	      );
      return(0);
    } else {
      return(1)
    }
  }
  #################################################################

   sub merge_features
    {
      my ($keep,$kill) = @_;

      unless ($keep && $kill){
	print qq(
		 <h3>Merge two existing Features</h3>
		 <form action="$ENV{'SCRIPT_NAME'}" method="GET">
		 <BR><BR>
		 Feature to stay alive  
		 <INPUT TYPE="text" NAME="keep_id" SIZE="10" MAXLENGTH="10"<br>
		 Feature to remove after merge   
		 <INPUT TYPE="text" NAME="kill_id" SIZE="10" MAXLENGTH="10">
		 <INPUT TYPE="hidden" NAME="action" VALUE="merge_features">
		 <br><br>
		 <INPUT TYPE="submit" VALUE="Merge" onClick="return validate_merge()">
		 <input type="reset" value="Clear Form" />		
		 </form>
		);
	print "<BR><font color=red>Please enter feature ids in both fields</font>" if($keep or $kill);
 
      } else {
	print  "attempting to merge $keep | $kill<br>";
	my $db = get_db_connection();
	#get ids if names passed
      
	my ($keep_id, $kill_id);
	$keep_id = ($keep =~ /WBsf\d+/) ? $keep : $db->idGetByAnyName("$keep")->[0];
	$kill_id = ($kill =~ /WBsf\d+/) ? $kill : $db->idGetByAnyName("$kill")->[0];
      
	#check ids are valid
	foreach ($keep_id, $kill_id) {
	  unless( defined $db->idExists($_) ) {
	    print qq(Feature does not exist<br>);
	    exit(0);
	  }
	}
      
	#do the merge
	if ( my $id = $db->idMerge($kill_id,$keep_id) ) {
	  print "Merge complete, $kill_id is DEAD and has been merged into feature $keep_id <br>";
	  #notify
	  send_mail("webserver",[$MAILS->{$USER},$MAILS->{'cgc'}],"Merged features $keep ($keep_id) and $kill ($kill_id)", "FEATURE MERGE\nUSER : $USER\nLIVE:retained Feature id for $keep_id\nDEAD: killed Feature id $kill_id \n");
	} else {
	  print "Sorry, the merge failed<br>";
	}
      }
    }
  
  sub new_feature {
    my $new = shift;
    unless ($new) {
      ## print the query form
      print qq(
	       <form action="$ENV{'SCRIPT_NAME'}" method="GET">
	       <BR><BR>
	       <h3>Request a new Feature Id</h3> 
	       <INPUT TYPE="hidden" NAME="action" VALUE="new_feature">
	       <INPUT TYPE="hidden" NAME="new" VALUE="1">
	       <INPUT TYPE="submit" VALUE="Get ID">
	       </form>
	      );
    } else {
      my $db = get_db_connection();
	my $id = $db->idCreate;
	print "$id created<br>";
			
	my $mail_msg = "Feature ID $id created\n";
	send_mail("webserver",
		  [$MAILS->{$USER},$MAILS->{'cgc'}],
		  "New Feature : $id", 
		  "$mail_msg");
    }
  }

  
  sub kill_feature {
    my $kill_id = shift;
    my $remark = shift;
    unless ($kill_id) {
      ## print the query form
      print qq(
	       <form action="$ENV{'SCRIPT_NAME'}" method="GET">
	       <BR><BR>
	       <h3>Kill a Feature Id</h3>
	       Enter ID to kill<br>  
	       <INPUT TYPE="text" NAME="kill_id" SIZE="10" MAXLENGTH="10" VALUE=""><BR>
	       Please give reason for removal (<font color=red>required</font>)<br> 
	       <INPUT TYPE="text" NAME="remark" SIZE="100" MAXLENGTH="200"><BR>
	       <INPUT TYPE="hidden" NAME="action" VALUE="kill_feature">
	       <br><br>
	       <INPUT TYPE="submit" VALUE="Kill" onClick="return validate_kill()">
	       </form>
	      );
    } else {
      if ( $remark ) {
	my $db = get_db_connection();

	if ( $db->idExists($kill_id) ) {	
	  if ( $db->idKill($kill_id) ) {
	    print qq( $kill_id killed<br>);
	    send_mail("webserver",
		      [$MAILS->{$USER},$MAILS->{'cgc'}],
		      "Feature killed - $kill_id", 
		      "\nUSER : $USER\nID : $kill_id\nREASON : $remark"
		     );
	  }
	} else {
	  print "$kill_id doesn't exists"
	}
      } else {
	print "The remark is a <font color=red>required field</font><br>";
      }
    }
  }			

  

  		
  sub query {
      ##run the query
      my $db = get_db_connection();
      my $query = " select object_public_id from primary_identifier where domain_id=3 order by object_public_id desc limit 1";
      my $last_object = $db->dbh->selectcol_arrayref($query,undef)->[0];

      print "the last Feature ID used was $last_object";

  }




  #################################################################
  sub print_selector {
    my ($action) = @_;

    my $nf = "";
    my $qf = "";
    my $mf = "";
    my $kf = "";

    if ($action eq "query") {
      $qf = " selected";
    }
    ;
    if ($action eq "new_feature") {
      $nf = " selected";
    }
    ;
    if ($action eq "merge_features") {
      $mf = " selected";
    }
    ;
    if ($action eq "kill_feature") {
      $kf = " selected";
    }
    ;


    print qq(
	     What do you want to do ?  <br><br>
	     <form name="doAction" action="$ENV{'SCRIPT_NAME'}" method="GET">
	     Select action to perform:
	     <SELECT NAME="action" onChange="formSubmit()">
	     <OPTGROUP LABEL="Feature information">
	     <OPTION $qf LABEL="Find Feature information" value="query">Show last used Feature id</OPTION>
	    );

    if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ) {
      print qq(
	       <OPTION $nf LABEL="new_feature" value="new_feature">Request a new Feature ID</OPTION>
	      );
    }
    if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ) {
      print qq(
	       <OPTION $mf LABEL="merge_features" value="merge_features">Merge two Feature IDs</OPTION>
	      );
    }

    if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ) {
      print qq(
	       <OPTION $kf LABEL="kill_feature" value="kill_feature">Kill a Feature ID</OPTION>
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
	
    my $DOMAIN = 'Feature';
    my $db = NameDB->connect($DB,$USER,$PASS,1); #1 is for web output
	
    $db || return(undef);
    $db->setDomain($DOMAIN);
	
    #set_web_output($db); # turn on web reporting
    return($db);
  }

  #################################################################
  sub set_web_output {
    my ($conn) = @_;
    $conn->web(1);		# turn on web reporting
	
  }
  #################################################################




  sub print_javascript {

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
	     function validate_kill(id)
	     {
	       var agree = confirm("Are you sure you want to kill this feature?")
		 if (agree) {
		   return(true);
		 } else {
		   return(false)
		 }

	     }
	  
	     </script>
	    );
  }

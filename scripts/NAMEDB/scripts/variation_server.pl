#!/usr/local/bin/perl -wT
#author ar2
#use lib "../lib";
use strict;
use vars qw($SSO_USER $PASS $DB $VALID_USERS $VALID_API_USERS $VALID_CGCNAME_USERS $USER $MAIL_NOTIFY_LIST $MAILS $LIVE);
use SangerPaths qw(core celegans);
use SangerWeb;
use NameDB_handler;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use Website::Utilities::Mail;
use File::Path;
$| = 1;

## a list of valid  IDs to use this resources
$VALID_USERS = {
  # these are the users WBPerson id
  'klh'       => 3111,
  'ar2'       => 1847,
  'mt3'       => 2970,
  'mh6'       => 4036,
  'xdwang'    => 1760,
  'kyook'     => 712,
  'jolenef'   => 2021,
  'wen'       => 101,
  'ranjana'   => 324,
  'vanauken'  => 1843,
  'pad'       => 1983,
  'stlouis'   => 1,
  'caltech'   => 1,
  'cshl'      => 1,
  'sanger'    => 1,
};

## a list of valid SSO login names for each DB operation
$VALID_API_USERS = {
		'query'		=> [qw( pad gw3 mh6 mt3 klh stlouis caltech cshl sanger kyook jolenef xdwang wen ranjana vanauken)],
		'dump'		=> [qw( pad gw3 mh6 mt3 klh stlouis caltech cshl sanger kyook jolenef xdwang wen ranjana vanauken)],
		'merge_var'	=> [qw( pad gw3 mt3 mh6 klh kyook jolenef xdwang wen ranjana vanauken)],
		'new_var'	=> [qw( pad gw3 mt3 mh6 klh kyook jolenef xdwang wen ranjana vanauken)], 
		'kill_var'	=> [qw( pad gw3 mt3 mh6 klh kyook jolenef xdwang wen ranjana vanauken)],
		'change_name'	=> [qw( pad gw3 mt3 mh6 klh kyook jolenef xdwang wen ranjana vanauken)],
};

## a list of valid SSO login names able to add GCG name
$VALID_CGCNAME_USERS = {
		'mt3'			=> 1,
};

$MAILS = {
        'klh'           =>      'klh@sanger.ac.uk',
	'mh6'		=>	'mh6@sanger.ac.uk',
	'mt3'		=>	'mt3@sanger.ac.uk',
	'kyook'         =>      'karen@wormbase.org',
	'pad'           =>      'paul.davis@wormbase.org', 
        'xdwang'        =>      'xdwang@caltech.edu',
        'wen'           =>      'wchen@its.caltech.edu',
        'ranjana'       =>      'ranjana@its.caltech.edu',
        'vanauken'      =>      'vanauken@caltech.edu',	
        'stlouis'	=>	'stlouis@wormbase.org',
	'caltech'	=>	'caltech@wormbase.org',
	'cshl'		=>	'cshl@wormbase.org',
	'sanger'	=>	'sanger@wormbase.org',
	'cgc'           =>      'mt3@sanger.ac.uk',
	'gw3'           =>      'gw3@sanger.ac.uk'
};

## a list of people to mail when a DB operation occurs
$MAIL_NOTIFY_LIST = [qw(gw3 klh)];

&main();
1;


#################################################################
sub main {
    my $web = 1;
    my $path = SangerWeb->document_root();
    my $sw = SangerWeb->new( {
	'banner' => "Wormbase Variation Name Management",
	'inifile'=> "$path/Projects/C_elegans/header.ini",
	'author' => 'mt3',
	'onload' => 'init()',
    });
    if ($sw->is_dev()) {
	$DB = 'test_wbgene_id;utlt-db;3307';
	$sw->banner("This is the test server");
    } else {
	$DB = 'nameserver_live;web-wwwdb-core-02;3449';
	$sw->banner("This is the LIVE server");
    }


    $SSO_USER = $sw->username(); ## for SSO
    if( $SSO_USER =~ /^(\S+)\@wormbase.org/) {
      $PASS = $1;
    }
    elsif( $SSO_USER =~ /^(\S+)\@caltech.edu/){
     $PASS = $1; 
    }
    else {
      $PASS = $SSO_USER;
    }

    # This is a hard coded SSO - mysql conversion.
    if ($PASS eq "karen") {
      $PASS = "kyook";
    }
    if ($PASS eq "paul.davis") {
      $PASS = "pad";
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
    print $sw->header({'title'  => "WormBase Variation ID Server $DB"});
    #print $sw->header();
    
    if (!defined $VALID_USERS->{$USER}) {
	print qq(<h3>Sorry, you $USER are not authorised to access this resource. Please contact the Wormbase team</h3>);
	print $sw->footer();
	return;
    } else {
	print qq(Authenticated database user: "$USER"<BR>		);
    }

    ## get CGI parameters
    my $action   = $sw->cgi->param('action')|| "query";
    my $name     = $sw->cgi->param('name');
    my $varid    = $sw->cgi->param('varid');
    my $lookup   = $sw->cgi->param('lookup');
    my $keep_id  = $sw->cgi->param('keep_id');
    my $kill_id  = $sw->cgi->param('kill_id');
    my $change   = $sw->cgi->param('change_name');


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
    
    }elsif( $action =~ /change_name/) {
	if( is_authorised($USER,$action) == 1) {
	    &change_name($varid,$change);
	}
    
    } elsif( $action eq "dump") {
	&dump_all("dump");
    } 
    elsif( $action eq "last") {
	&last("last");
    }
    else {
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
    'subject' => ($DB =~ 'test' ? 'TEST : ' : " ")."Var_NAMEDB: $subject",
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
			send_mail("mt3",
                                  [$MAILS->{$USER},
                                   $MAILS->{'cgc'}],
                                  "Merged Variations $keep ($keep_id) and $kill ($kill_id)", 
                                  "VARIATION MERGE\nUSER : $USER\nLIVE:retained WBVarID for $keep_id\nDEAD: killed VarID $kill_id \n");
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
	if(&_check_name($public) == 0) {
	    my $db = get_db_connection();
	    my $id = $db->idCreate;
	    $db->addName($id,'Public_name'=>$public);
	    print "$id created with name $public<br>";
	    send_mail('mt3', 
                      $MAILS->{'cgc'}, 
                      "WBVarID request $public $USER", 
                      "WBVarID request $public : $id $USER");
	}
	else {
	    print "not a good Var name\n";
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
    	  <INPUT TYPE="text" NAME="kill_id" SIZE="13" MAXLENGTH="13" VALUE="">
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
				send_mail('mt3', 
                                          $MAILS->{'cgc'}, 
                                          "WBVarID $death killed by $USER", 
                                          "WBVarID $death (user entererd:$kill_id) killed by  $USER");
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
			&print_history($db,$var);
			print qq(
			<hr align="left">
			<BR><BR><BR>
			);
		}
		else {
			print qq( $lookup is not known to the Variation nameserver<br>);
		}
	}

}

sub _check_name {
    my $name = shift;
    my $db = get_db_connection();
    my $var = $db->idGetByTypedName('Public_name'=>$name)->[0];

    if($var) {	
	print "$name already exists as $var";
	return 1;
    }
    elsif ($var =~ /(\w+)[a-z]+$/) {
	my $short_var = $db->idGetByTypedName('Public_name'=>$1)->[0];
	if($short_var) {	
	    print "$var looks like a poorly name version of $short_var";
	    return 1;
	}
    }
    return 0;
}



#################################################################
sub print_selector {
	my ($action) = @_;

	my $nv = "";
	my $qv = "";
	my $mv = "";
	my $kv = "";
	my $da = "";
	my $cn = "";
	if ($action eq "query") 		{ $qv = " selected" };
	if ($action eq "new_var") 		{ $nv = " selected" };
	if ($action eq "merge_var")   { $mv = " selected" };
	if ($action eq "kill_var") 	{ $kv = " selected" };
	if ($action eq "dump") 		{ $da = " selected" };
	if ($action eq "changename") 		{ $cn = " selected" };

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
	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	 <OPTION $kv LABEL="change_name" value="change_name">Change name of existing Variation Id</OPTION>
		);
	}
	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	 <OPTION $da LABEL="dump" value="dump">Dump all Variation Ids</OPTION>
		);
	}
	if (grep {/$USER/} @{$VALID_API_USERS->{$action}} ){
		print qq(
    	 <OPTION $da LABEL="last" value="last">Show last used Variation Id</OPTION>
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

sub last {
    my $db = get_db_connection();
    my $domain    = $db->getDomain;
    my $domain_id = $db->getDomainId($domain);
    my $query = " select primary_identifier.object_public_id from primary_identifier where domain_id=? order by object_id desc limit 1";
    my $last_object = $db->dbh->selectcol_arrayref($query,undef, $domain_id);
    my $id = $last_object->[0];

    print "Last used $domain ID is <b>$id<b><br>";
}



sub dump_all {
  my $dump = shift;
  if ( $dump ) {
    # iterate over all genes and print details
      print qq( <meta http-equiv="REFRESH" content="0; URL=make_vars.txt.pl?user=$USER">);
      exit(0);
      my $db = get_db_connection();
      my $query =<<END;
    SELECT primary_identifier.object_public_id, secondary_identifier.object_name
	FROM primary_identifier,secondary_identifier 
	WHERE secondary_identifier.object_id = primary_identifier.object_id and primary_identifier.domain_id = 3
	ORDER by object_public_id
	LIMIT 10;
END
    my $sth = $db->dbh->prepare($query);
    $sth->execute or die "Unable to execute query: $db->errstr\n";
    my $row;
    print "query executed<br>";
    #write a new version of webpage for stl
    my $path = SangerWeb->document_root()."/tmp/Projects/C_elegans/LOCI";
    mkpath $path unless -e $path;
    my $IDpage = "$path/variation_ids.txt";
    print "$IDpage<br>";
    open (ID,">$IDpage") or die "cant open $IDpage : $!\n";

    while ($row = $sth->fetchrow_arrayref) {
      print ID "$row->[0]\t$row->[1]<br>";
      print "$row->[0]\t$row->[1]<br>";
    }
    $sth->finish;

    close ID;	
    print qq(<a href="http://www.sanger.ac.uk/tmp/Projects/C_elegans/LOCI/variation_ids.txt">Updated here</a>);
   # print qq( <meta http-equiv="REFRESH" content="0; URL=http://www.sanger.ac.uk/tmp/Projects/C_elegans/variation_ids.txt">);
  } else {
    print qq( 
			<h3>Print page of all variation ids</h3>
			<form action="$ENV{'SCRIPT_NAME'}" method="GET">
			<INPUT TYPE="hidden" NAME="action" VALUE="dump">
			<br><br>
			<INPUT TYPE="submit" VALUE="Dump">
		</form>
		);
  }
}

sub change_name {
    my $id = shift;
    my $name = shift;

     unless ($id) {
	## print the query form
	print qq(
		 <form action="$ENV{'SCRIPT_NAME'}" method="GET">
		 <BR><BR>
		 <h3>Change the name of an existing Variation</h3>
		 Enter Id of Variation<br>  
		 <INPUT TYPE="text" NAME="varid" SIZE="14" MAXLENGTH="14" VALUE="">
		 <INPUT TYPE="text" NAME="change_name" SIZE="14" MAXLENGTH="14" VALUE="">
		 <INPUT TYPE="hidden" NAME="action" VALUE="change_name">
		 <br><br>
		 <INPUT TYPE="submit" VALUE="Search">
		 </form>
		 );
    }else {
	print  "attempting to rename $id to $name<br>";
	my $db = get_db_connection();
	if($db->addName($id,'Public_name'=>$name)) {
	    print "successfully changed $id to $name<br>";
	}
	else {
	    print "<font color=red>FAILED to rename</font>";
	}
    }
}
    
#These are the bits that are stopping the use of NameDB_handler.pm - because of NameDB_handler::validate_name.
sub print_history {
  my $db = shift;
  my $id = shift;

  &validate_id($db,$id);

  my $history = $db->idGetHistory($id);
  print "$id<br>";
  foreach my $event (@{$history}) {
    print $event->{version}," ";
    print $event->{date}," ";
    print $event->{event}," ";
    print $event->{name_type}," " if (defined $event->{name_type});
    print $event->{name_value}," " if (defined $event->{name_value});
    print $event->{user},"<br>";
  }
  return;
}
sub validate_id
  {
    my $db = shift;
    my $id = shift;
    unless ( $db->idExists($id) ) {
      $db->dienice("$id does not exist");
      return undef ;
    }
    unless( $db->idLive($id) == 1) {
      $db->dienice("$id is not live");
      return undef;
    }
    return $id;
  }

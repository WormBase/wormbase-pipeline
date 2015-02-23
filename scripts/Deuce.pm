#!/software/bin/perl -w
#
# Deuce.pm
# 
# by Gary Williams
#
# Some methods for querying the WormBAse datomic database.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-02-23 09:29:53 $      

=pod

=head1 NAME

 Deuce

=head1 SYNOPSIS

  use Deuce;

  my $db = Deuce->new('db.wormbase.org:8120/colonnade', $user, $password);
  my $results = $db->$query($datomic_querystr);
  my $schema = $db->schema;
  my $object = $db->fetch_object($class, $id);
  my @ID_names = $db->fetch($class, $class_id_pattern);
  my @ID_names = $db->fetch($class, $class_id_pattern, $tace_command_filter);

=head1 DESCRIPTION
Methods for doing various searches in the datomic WormBase databases

=head1 CONTACT

Gary gw3@ebi.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut




package Deuce;

use strict;             
use lib $ENV{'CVS_DIR'};
use Carp;
use warnings;
use Data::Dumper;
use HTTP::Tiny;
use HTTP::CookieJar;
use edn;


=head2 

    Title   :   new
    Usage   :   my $db = Deuce->new('db.wormbase.org:8120/colonnade', $user, $password);
    Function:   initialises the connection to the datomic REST server
    Returns :   Deuce object;
    Args    :   url - string - URL of datomic REST server
                user - string - user name to authenticate as
                password - string - password to authenticate the user

=cut
sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  $self->{'url'} = shift; 
  my $user = shift;
  my $password = shift;

  # default
  if (!defined $self->{'url'}) {
    $self->{'url'} = 'db.wormbase.org:8120/colonnade';

    my $login_details_file = glob("~/.deuce/COLONNADE.s");
    open(my $infh, $login_details_file)
      or die("Can't open secure account details file $login_details_file\n");
    while (<$infh>){
      /^USER_ID:(\S+)$/ and $user = $1;
      /^PASSWD:(\S+)$/ and $password = $1;
    }
    close($infh);

    die("Could not find both user name and password in login details file\n")
      if not defined $user or not defined $password;
  }
  $self->{'url'} =~ s/http\:\/\///;

  my $cookie_jar = HTTP::CookieJar->new;
  $self->{'http'} = HTTP::Tiny->new('cookie_jar' => $cookie_jar);
  my $url = $self->{'url'};
  my $http = $self->{'http'};
  my $response = $self->{'http'}->get("http://${user}:${password}\@${url}");
  #print "token content: ".$response->{'content'}."\n";

  # parse out token
  my ($token) = ($response->{'content'} =~ /trace_token = \'(\S+)';/);

  my %headers;
  $headers{'X-XSRF-Token'} = $token;
  $headers{'content-type'} = 'application/edn';
  my %options;
  $options{'headers'} = \%headers;
  $self->{'options'} =  \%options;

  return $self;
}

#############################################################################

=head2 

    Title   :   query
    Usage   :   my $results = $db->$query($querystr, $max_rows)
    Function:   get the results of a query from the REST server
    Returns :   int count of total number of results, not constrained by the $max_rows
                array of arrays holding the result of the first $max_rows results of the query
    Args    :   querystr - string - datomic query
                max_rows - int - maximum number of rows to return - if undef, then set to 10 million

=cut


sub query {
 my ($self, $querystr, $max_rows) = @_;

 my $max_rows_str = '';
 if (!defined $max_rows) {
   $max_rows = 10000000; # Set to more than the max number of objects in any class.
 }
 $max_rows_str = ", :max-rows $max_rows";
 $querystr = '{:query ' . $querystr . ", :rules [] ${max_rows_str}}"; # encode query as EDN

 $self->{'options'}->{'content'} = $querystr;

 my $response = $self->{'http'}->post(
				  'http://' . $self->{'url'} . '/query',
				  $self->{'options'}
				 );
 
  if (!$response->{success}) {
    warn "Query failure\nContent: $response->{'content'}\nStatus: $response->{status}\nReason: $response->{reason}\n";
    return (0, undef);
  }

 my $content = edn::read($response->{'content'}); 

#    print "Count: ". $content->{count}."\n"; # number of results in total, not limited by max-rows
#    foreach my $i (@{$content->{results}}) { # 'results' hash is array of arrays
#      print "@{$i}\n";
#    }

 return ( $content->{count}, $content->{results});
}

#############################################################################

=head2 

    Title   :   get_schema
    Usage   :   my $schema = $db->get_schema;
    Function:   returns a reference to the complete schema of the database (hash-ref)
    Returns :   the complete schema (hash-ref)
    Args    :   

=cut

sub get_schema {
 my ($self) = @_;

 if (!defined $self->{schema}) {

   print "Getting schema ...\n";
   
   my $schema_url = 'http://db.wormbase.org:8120/schema';
   my $response = HTTP::Tiny->new->get($schema_url);
   
   if (!$response->{success}) {
     die "In schema()\nStatus: $response->{status} Reason: $response->{reason}\nContent: $response->{'content'}\n";
   }
   
   $self->{schema} = edn::read($response->{'content'});

   if (!exists $self->{schema}{classes} || !exists $self->{schema}{attributes}) {
     die "In schema()\nThe expected data is not in the schema\n";
   }
 }

 return $self->{schema};
}

#############################################################################

=head2 

    Title   :   get_classes
    Usage   :   my $classes = $db->get_classes();
    Function:   returns the complete hash-ref of classes keyed on ACE name of the class
            :   useful data includes:
            :     'db/ident' => name of the class in the datomic database (e.g. 'transposon/id')
            :     'pace/is-hash' => boolean, true is a hash (e.g. ?Evidence)
    Returns :   the complete classes - hash keyed by name of class
    Args    :   

=cut

sub get_classes {
 my ($self) = @_;

 if (!defined $self->{classes}) {

   my $schema = $self->get_schema();
   my %classes;
   my %ace_classes;
   foreach my $class (@{$schema->{classes}}) {
     my $ace_class_name = $class->{'pace/identifies-class'};
     $classes{$ace_class_name} = $class;
     my $datomic_class_name = ${$class->{'db/ident'}};
     # simple hash to aid finding the ACE class name from the datomic class name
     # NB we key by both with and without the '/id' on the end of the datomic class name
     $ace_classes{$datomic_class_name} = $ace_class_name; 
     $datomic_class_name =~ s/\/id//;
     $ace_classes{$datomic_class_name} = $ace_class_name; 
   }
   $self->{classes} = \%classes;
   $self->{ace_classes} = \%ace_classes;
 }

 return $self->{classes};
}


#############################################################################

=head2 

    Title   :   get_attributes
    Usage   :   my $attributes = $db->get_attributes();
    Function:   returns the complete attributes of the classes
    Returns :   the complete attributes of the classes - array-ref of hashes
    Args    :   

=cut

sub get_attributes {
 my ($self) = @_;

 if (!defined $self->{attributes}) {

   my $schema = $self->get_schema();
  
   my $attributes = $schema->{attributes};
   $self->{attributes} = $attributes
 }

 return $self->{attributes};
}

#############################################################################

=head2 

    Title   :   get_class_attributes
    Usage   :   my $class_desc = $db->get_class_attributes($class);
    Function:   returns an array holding the attributes for this class
    Returns :   array of hash holding the attributes of the class
    Args    :   class - string - ACE name of the class

=cut

sub get_class_attributes {
 my ($self, $class) = @_;

  if (!defined $self->{class_desc}{$class}) {
    my @desc;
    my $attributes = $self->get_attributes();

    my $class_name_in_datomic = $self->get_class_name_in_datomic($class);

    # look for this classes' attributes
    foreach my $attr (@{$attributes}) {
      my $ident = ${$attr->{'db/ident'}};
      if ($ident =~ /^${class_name_in_datomic}[\/\.]/i) {
	push @desc, $attr;
      }
    }

    # sort them by 'db/id' event number
    my @desc_sorted = sort {$a->{'db/id'} <=> $b->{'db/id'}} @desc;

    $self->{class_desc}{$class} = \@desc_sorted;
  }
  return $self->{class_desc}{$class};
}

#############################################################################

=head2 

    Title   :   get_class_name_in_datomic
    Usage   :   my $class_desc = $db->get_class_name_in_datomic($class);
    Function:   input an ACE class name and it returns the datomic equivalent of the ACE class name
    Returns :   returns the datomic equivalent of the ACE class name
    Args    :   class - string - ACE class name

=cut

sub get_class_name_in_datomic {
  my ($self, $class) = @_;
  
  my $classes = $self->get_classes();
  if (!exists  $classes->{$class}) {
    warn "Can't find datomic name of $class, assuming we already have a datomic name\n";
    return $class;
  }
  my $ident = ${$classes->{$class}{'db/ident'}};

  $ident =~ s/\/id//; # remove the '/id' from the end
  return $ident;
}

#############################################################################

=head2 

    Title   :   get_class_name_in_ace
    Usage   :   my $class_desc = $db->get_class_name_in_ace($class);
    Function:   input an datomic class name and it returns the ACE equivalent of the datomic class name
    Returns :   returns the ACE equivalent of the datomic class name
    Args    :   class - string - class datomic name

=cut

sub get_class_name_in_ace {
  my ($self, $class) = @_;
  
  my $classes = $self->get_classes(); # ensure the ace_classes hash is populated
  my %ace_classes = %{$self->{ace_classes}};
  if (!exists  $self->ace_classes->{$class}) {
    warn "Can't find Ace name of $class, assuming we already have an Ace name\n";
    return $class;
  }
  return $ace_classes{$class};
}

#############################################################################

=head2 

    Title   :   is_class_a_hash
    Usage   :   my $boolean = $db->is_class_a_hash($class);
    Function:   returns whether the class is a hash class, like ?Evidence or not
    Returns :   boolean
    Args    :   class - string - class datomic name

=cut

sub is_class_a_hash {
  my ($self, $class) = @_;
  
  my $classes = $self->get_classes();
  return $classes->{$class}{'pace/is-hash'};
}


#############################################################################

=head2 

    Title   :   output_class_model_as_ace
    Usage   :   my $model = $db->output_class_model_as_ace($class);
    Function:   returns the class model as an text ACE file
    Returns :   text ACE file class format
    Args    :   class - string - ACE name of the class

=cut

sub output_class_model_as_ace {
  my ($self, $class) = @_;

  my $text;

  # walk along the tag structure
  my $objstr = $self->model($class);

#  my $text = $self->display_class_text($objstr);

  # +++ need to add the class name at the start with an indent
  my $element = 0;
  my $level=0;
  my @level;
  my $model = "?$class ";
  my $indent = length $model;
  my $right;
  my $down;

  $level[$level] = $element;

  do {
    do {
      do {
	$level[$level] = $element;
	my $ace_tags = $objstr->[$element]->{ace_tags};
	my $type = $objstr->[$element]->{type};
	print "$ace_tags ($type) ";
	
	$right = $objstr->[$element]->{right};
	if (defined $right) {
	  $element = $right;
	  #print " > $element ";
	  $level++;
	  $level[$level] = $right;
	}
      } while (defined $right);
      print "\n"; # end of line
      
      $down = $objstr->[$element]->{down};
      if (defined $down) {
	$element = $down;
	$level[$level] = $down;
	#print " V $element \n";
      }
      
    } while (defined $down);

    # no more nodes in this tree, go back some levels and down
    my $backtracking = 1;
    while ($backtracking) {
      $level--;
      $element = $level[$level];
      #print " < $element level=$level\n";
      $down = $objstr->[$element]->{down};
      if (defined $down) {
	$element = $down;
	#print " VV $element ";
	$level[$level] = $element;
	#print "Want to exit loop level=$level elem=$element down=$down\n";
      }
      #print "down = $down\n";
      #print "level = $level\n";
      if (defined $down || $level < 0) {$backtracking = 0}
    }
    #print "out of loop level=$level\n";
  } while ($level >= 0);

  return $text;
}


#############################################################################

=head2 

    Title   :   model
    Usage   :   my $model = $db->model($class);
    Function:   returns the class model as an text ACE file
    Returns :   array-ref of class tags and values
    Args    :   string - datomic class name

=cut

sub model {
  my ($self, $class) = @_;

  if (!defined $self->{class_object}{$class}) {

    my $class_attributes = $self->get_class_attributes($class); # attributes of class
    my $classes = $self->get_classes(); # hash of all classes' data
    my @class_desc_array = $classes->{$class};
    my %class_desc = %{$class_desc_array[0]};
    
    my @objstr = ();
    $self->parse_class_object($class_attributes, \%class_desc, 0, \@objstr);
    $self->{class_object}{$class} = \@objstr;
  }
  return $self->{class_object}{$class};
}

#############################################################################

=head2 

    Title   :   parse_class_object
    Usage   :   my $model = $db->parse_class_object($class_attributes, $class_data, $level, $objstr);
    Function:   returns the class model as a parsed object
    Returns :   array of class elements
    Args    :   class_attributes = array of attributes that have an ident that matches this class
                class_desc = hash-ref of class data for this class
                level = int level of tags
                objstr = ref to array of hashes - populated in this method - one hash for each tag-and-value

=cut

sub parse_class_object {
  my($self, $class_attributes, $class_desc, $level, $objstr) = @_;

  my $prev_ident_base = 'dummy ident';

  my @level_names;       # tag names of the levels
  my @level_elements;    # element numbers currently at level[n] for 'down' pointers

  foreach my $element (@{$class_attributes}) {

    my $ace_tags = $element->{'pace/tags'};
    my @ace_tags = split /\s+/, $ace_tags;
    my $ident = ${$element->{'db/ident'}}; # datomic name
    # 'db/cardinality' is not defined for enum tags
    my $unique = undef;
    if (exists $element->{'db/cardinality'}) {(${$element->{'db/cardinality'}} eq 'db.cardinality/one') ? $unique = 1 : $unique = 0};
    my $type = (exists $element->{'db/valueType'}) ? (${$element->{'db/valueType'}}) : undef;
    if (!defined $type || $type eq '') {$type = 'enum'}
    $type =~ s/db.type\///;
    my $hash = '';
    if (exists $element->{'pace/use-ns'}) {$hash = $element->{'pace/use-ns'}[0]} # e.g. ('evidence') - assume we only have one hash object!
    my $xref = (exists $element->{'pace/obj-ref'}) ? (${$element->{'pace/obj-ref'}}) : undef;
    my $order =  $element->{'pace/order'};
    my $name = $ace_tags[$#ace_tags];
    if (!defined $name) {($name) = ($ident =~ /\/(\S+)/); $name = '^'.$name}

    # is this element part of the same line?
    if ($ident =~ /^$prev_ident_base/) { #  ('cds.cds/start' =~ /^'cds.cds/'/)
      print "$ident $type / ";
      my $prev = scalar @{$objstr} - 1;
      push @{$objstr}, {
			attribute => $element,
			ident     => $ident,
			level     => $level, # display level for indentation
			ace_tags  => $ace_tags, # tag or object name
			name      => $name, # last tag
			xref      => $xref, # xref class name
			type      => $type, # type of value
			isUnique  => $unique,
			isHash    => $hash, # want a hash at the end of the line
			down      => undef, # navigation 
			right     => undef, # navigation
		       };
      
      if (scalar @{$objstr} - 1 >= 0) {
	$objstr->[$prev]->{right} = scalar @{$objstr} - 1;
      }

    } else {

      # see if this is the first tag of a line
      if ($ident !~ /\./) { 
	print "\n$ident $type / ";

	push @{$objstr}, {
			  attribute => $element,
			  ident     => $ident,
			  level     => $level, # display level for indentation
			  ace_tags  => $ace_tags, # tags or object name
			  name      => $name, # last tag
			  xref      => $xref, # class name
			  type      => $type, # type of value
			  isUnique  => $unique,
			  isHash    => $hash,
			  down      => undef, # navigation 
			  right     => undef, # navigation
			 };

	# put the new line at the end of the other tags which start lines
	if (scalar @{$objstr} > 1) {
	  my $elem = 0;
	  while (defined $objstr->[$elem]->{down}) {
	    $elem = $objstr->[$elem]->{down};
	  }
	  $objstr->[$elem]->{down} = scalar @{$objstr} - 1;
	}


	
      } else {
	# is this element part of some previous line? like the enum tags at the end of: 
	# ?CDS Properties Coding Prediction_status Predicted
	#                                          Partially_confirmed
	#                                          Confirmed
	# which are at the end of the CDS attributes
	# go through the objstr first column of tags matching the idents
	my ($ident_start) = ($ident =~ /(\S+?)\//);
	$ident_start =~ s/\./\//; # change . to / so it can match a previous start-of-line ident
	my $found = 0;
	my $elem = 0;
	while (defined $elem) {
	  if ($ident_start eq $objstr->[$elem]->{ident}) {
	    print "\n$ident ($type) ";

	    $found = 1;
	    
	    push @{$objstr}, {
			      attribute => $element,
			      ident     => $ident,
			      level     => $level, # display level for indentation
			      name      => $name, # last tag
			      ace_tags  => $ace_tags, # tags or object name
			      xref      => $xref, # class name
			      type      => $type, # type of value
			      isUnique  => $unique,
			      isHash    => $hash,
			      down      => undef, # navigation 
			      right     => undef, # navigation
			     };
	    
	    if (defined  $objstr->[$elem] && !defined $objstr->[$elem]->{right}) {
	      if (scalar @{$objstr} - 1 >= 0) {
		print " right / ";
		$objstr->[$elem]->{right} = scalar @{$objstr} - 1;
	      } else {print "Shouldn't get here"}
	    } else {
	      if ($type eq 'enum') { # multiple enum fields - go down
		$elem = $objstr->[$elem]->{right};
		while(defined $objstr->[$elem] && defined $objstr->[$elem]->{down}) {
		  $elem = $objstr->[$elem]->{down};
		}
		if (scalar @{$objstr} - 1 >= 0) {
		  print " down / ";
		  $objstr->[$elem]->{down} = scalar @{$objstr} - 1;
		} else {print "Shouldn't get here"}

	      } else { # multiple fields on the same line - go right
		$elem = $objstr->[$elem]->{right};
		while(defined $objstr->[$elem] && defined $objstr->[$elem]->{right}) {
		  $elem = $objstr->[$elem]->{right};
		}
		if (scalar @{$objstr} - 1 >= 0) {
		  print " right / ";
		  $objstr->[$elem]->{right} = scalar @{$objstr} - 1;
		} else {print "Shouldn't get here"}
		
	      }
	    }
	    
	    last;
	  }
	  $elem = $objstr->[$elem]->{down};
	}
	if (!$found) {
	  print "Not found a match of the class element $ident with any lines in the attributes\n";
	}
      }


    }    

    # change / to . and add '/' so 'cds.cds/start' will match with regex /^'cds/cds'/
    $prev_ident_base = $ident;
    $prev_ident_base =~ s/\//\./; 
    $prev_ident_base .= '/';
  }

  # now add the tags and xrefs from other classes that point in to this class
  my @incoming_xref = @{$class_desc->{'pace/xref'}};
  foreach my $pace_xref (@incoming_xref) {

    my $ident = ${$pace_xref->{'pace.xref/attribute'}}; # e.g. 'feature.associated-with-cds/cds'
    my $xref = ${$pace_xref->{'pace.xref/obj-ref'}}; # e.g. 'feature/id'
    my $tags = $pace_xref->{'pace.xref/tags'}; # e.g. 'Associated_feature'
    my $type = 'ref';
    if (!defined $type || $type eq '') {$type = 'enum'}
    $type =~ s/db.type\///;

    print "$ident  incoming xref / ";

    # put the new line at the end of the other tags which start lines
    my $elem = 0;
    while (defined $objstr->[$elem] && defined $objstr->[$elem]->{down}) {
      $elem = $objstr->[$elem]->{down};
    }

    # how do we know whether to add an evidence hash or not?

    push @{$objstr}, {
		      class     => $pace_xref,
		      ident     => $ident,
		      level     => $level, # display level for indentation
		      name      => $tags, # last tag
		      ace_tags  => $tags, # tags or object name
		      xref      => $xref, # class name
		      type      => $type, # type of value
		      isUnique  => undef, # how do we find if this is unique or not ? Thomas says that it is always assumed the this is NOT unique but he will think of a way to specify it later. Maybe hard-code the Sequence->Sequence xrefs to be UNIQUE?
		      isHash    => 0,
		      down      => undef, # navigation 
		      right     => undef, # navigation
		     };
	    
    
  }

  return;
}

#############################################################################

=head2 

    Title   :   fetch_object
    Usage   :   my $object = $db->fetch_object($class, $id);
    Function:   gets the complete object data
    Returns :   the object data
    Args    :   string - datomic class name
                string - ID name to fetch

=cut

sub fetch_object {	       
 my ($self, $class, $id) = @_;
 
 my $object_url = 'http://db.wormbase.org:8120/raw2/';
 $object_url .= "/$class;
 $object_url .= "/$id;

 my $response = HTTP::Tiny->new->get($object_url);

 #print "content: ".$response->{'content'}."\n";
 
 if (!$response->{success}) {
   die "In schema()\nStatus: $response->{status} Reason: $response->{reason}\nContent: $response->{'content'}\n";
 }

 return edn::read($response->{'content'});

}

#############################################################################

=head2 

    Title   :   ace_to_querystr
    Usage   :   my $querystr = $db->$ace_to_querystr($class, $class_pattern, $querystr)
    Function:   convert an ACE-style query string into a datomic-like query string
    Returns :   datomic query string or undef if error
    Args    :   $class = datomic name of class
                $class_pattern = optionally wildcarded pattern of ID in class to search for first
                $criteria = any subsequent criteria to search for 
                $classref - ref to string - returns the class of object returned if defined

=cut

sub ace_to_querystr {
  my ($self, $class, $pattern, $criteria, $classref) = @_;

  # [:find ?col1-id = ?col-id${class_var} to get a unique entity col-id that changes with the output class
  # $class_var = name of the class to output
  # col2-id, col3-id = col${criterion_count}-id  local criteria variable
  # ?col1 = ?${class_var} if we consistently use the class name as also being the variable for output then this keeps things easy
  # ?col2, ?col3 etc = ?col${criterion_count} local criteria variable

  my $result=undef; # the result of the query
  my $criterion_count = 0; # the next numeric suffix to use on variable identifiers

  if (!defined $class || $class eq '') {
    return $result;
  }

  if (!defined $pattern || $pattern eq '') {
    $pattern = '*';
  }

  if (!defined $criteria) {
    $criteria = '';
  }

  my $class_var = lc "${class}"; # initialise the class_variable to be the current lowerclassed class name (use wherever colonnade has '?col1, ?col6 etc.)
  if (defined $classref) {$$classref = $class_var}
  my $in;# = ':in $ ';
  my $where = ":where ";
  my $params = '';
  my $tag_name_element;
  my $op;

  print "In ace_to_querystr, Criteria: $criteria\n";

  # start off by constructing a query to get the class and pattern
  # is the pattern a regular expression?
  my $find = ":find [?col-id${class_var} ...] "; # initialise the col-id - return a collection of names
  if ($pattern =~ /\*/) {
# [:find ?col1-id :where [?col1 :cds/id ?col1-id] [(re-pattern "AC3.*") ?col1-regex] [(re-matches ?col1-regex ?col1-id)]] 
    $pattern = $self->protect($pattern);

    $where .= "[?${class_var} :${class_var}/id ?col-id${class_var}]";
    $where .= "[(re-pattern \"(?i)$pattern\") ?var${criterion_count}-regex]"; # case-independent regular expression
    $where .= "[(re-matches ?var${criterion_count}-regex ?col-id${class_var})]";
  } else {
# [:find ?col1-id :where [?col1 :cds/id ?col1-id] [(ground "AC3.3") ?col1-id]]
    # a single non-regex value
    $where .= "[?${class_var} :${class_var}/id ?col-id${class_var}]";
    $where .= "[(ground \"$pattern\") ?col-id${class_var}]";
  }

  my @criteria = split /\;/, $criteria;

  foreach my $criterion (@criteria) {
    if ($criterion =~ /^\s*$/) {next}

    my $tokens = $self->parse_cmd($criterion);
    if (!defined $tokens) {return undef} # syntax error found
    $criterion_count++; # for making variables that are local to this criterion
    
    # get description of first word of command 
    #       Name      Word Expect Type      Start End
    # e.g. ('NEGATE', '!', 'TAG', 'NEGATE', 1,    0)
    my $next_token = 0;
    my ($name, $word, $expect_next, $type, $can_start, $can_end) = @{$tokens->[$next_token]};
    my ($next_name, $next_word, $next_expect_next, $next_type, $next_can_start, $next_can_end);
    if ($next_token+1 < scalar @{$tokens}) {
      ($next_name, $next_word, $next_expect_next, $next_type, $next_can_start, $next_can_end) = @{$tokens->[$next_token+1]};
    }
    if ($can_start) {

      if ($type eq 'TAG') { # fnd objects where the TAG matches various conditions

	$tag_name_element = $self->get_tag_name($class, $tokens, \$next_token);
	my $objstr;
	my $ident;
	if (defined $tag_name_element) {
	  $objstr = $self->model($class);
	  $ident = $objstr->[$tag_name_element]->{ident};
	} else {
	  warn "Unknown type of tag '$word' in '$criterion'\n";
	  return undef;
	}
	
	if ($next_token+1 >= scalar @{$tokens}) { 	# at end of command, just the TAG so check it exists
	  
	  # TAG  (NEXT) (tag exists)
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/source-exons ?col2]]
	  # [:find ?name :where [?e :cds/id ?name][?e :cds/from-laboratory _]]
	  # FIND CDS; Prediction-status (enum )
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/prediction-status ?col2] [?col2 :db/ident ?col2-ident]]
	  
	  $where .= "[?${class_var} :${ident} _]"; # check the ident exists and we don't care what the contents are, so use '_'
	  
	} elsif (defined $next_type && $next_type eq 'OP') { # is there an OP after the tag?
	  
	  # TAG (NEXT) OP STRING|NUM
	  # FIND CDS; Species = "Caenorhabditis elegans" (ID of object)
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/species ?col2] [?col2 :species/id ?col2-id] [(ground "Caenorhabditis elegans") ?col2-id]]
	  # FIND CDS; Remark = "gw3" (simple string)
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/remark ?col2] [(ground "gw3") ?col2]]
	  # FIND CDS; Prediction-status = "confirmed" (enum value)
	  # *Doesn't work* [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/prediction-status ?col2] [?col2 :db/ident ?col2-ident] [(ground "confirmed") ?col2]]
	  
	  # FIND CDS; Sequence = "AC3"
	  # TAG value = ID name of the XREF
	  
	  $next_token++;
	  if ($next_token >= scalar @{$tokens}) {
	    warn "Expected something after '$next_word' in '$criterion'\n";
	    return undef;
	  }
	  my ($last_name, $last_word, $last_expect_next, $last_type, $last_can_start, $last_can_end) = @{$tokens->[$next_token]};
	  print "parsed command is: $word $next_word $last_word\n";
	  
	  # if the 'db/valueType' is 'db.type/ref', then look for the next tag/value on the line by going to the ->{right} in the objstr
	  # unless it has 'db/valueType' = 'db.type/ref' and 'pace/obj-ref' pointing to an outgoing xref
	  my $obj_type = $objstr->[$tag_name_element]->{'type'};
	  my $obj_xref = $objstr->[$tag_name_element]->{'xref'};
	  
	  if (!defined $obj_xref || $obj_type ne 'ref') {
	    # go right to the next TAG
	    $tag_name_element = $objstr->[$tag_name_element]->{'right'};
	    if (!defined $tag_name_element) {
	      warn "What? We seem to have fallen off the end of the model when trying to find the tag with the value to use in '$criterion'\n";
	    }
	    $obj_type = $objstr->[$tag_name_element]->{'type'};
	    $obj_xref = $objstr->[$tag_name_element]->{'xref'};
	    $ident = $objstr->[$tag_name_element]->{ident};
	  }
	  
	  
	  if (defined $obj_xref) { # tag is out-going xref
	    # +++
	    
	    
	    
	    
	  } elsif ($obj_type eq 'string') { # tag is string value - GE, LE should really do lexicographic order comparison
	    if ($last_name eq 'Value_') {$last_word = '"'.$last_word.'"';} # change the raw value into a string
	    
	    if ($next_name eq 'EQ') {
	      # normal string or regex? - Only NOT or EQ are sensible if regex
	      if ($last_word =~ /\*/) {
		#  [?col1 :cds/brief-identification ?col6] [(re-pattern ".*e") ?col6-regex] [(re-matches ?col6-regex ?col6)]
		# [:find ?col1-id :where [?col1 :cds/id ?col1-id] [(re-pattern "AC3.*") ?col1-regex] [(re-matches ?col1-regex ?col1-id)]] 
		$last_word = $self->protect($last_word);

		$where .= "[?${class_var} :${ident} ?var${criterion_count}]";
		$where .= "[(re-pattern $last_word) ?var${criterion_count}-regex]";
		$where .= "[(re-matches ?var${criterion_count}-regex ?var${criterion_count})]";
	      } else {
		# a single non-regex value
		# [:find ?col1-id :where [?col1 :cds/id ?col1-id] [(ground "AC3.3") ?col1-id]]
		$where .= "[?${class_var} :${ident} ?var${criterion_count}]";
		$where .= "[(ground $last_word) ?var${criterion_count}]";
		
	      }
	    } elsif ($next_name eq 'NOT') { 
	      
	      #[:find (count ?artist) .	      :where [?artist :artist/name]       (not-join [?artist]         [?release :release/artists ?artist]         [?release :release/year 1970])]
	      if ($last_word =~ /\*/) {
		$last_word = $self->protect($last_word);

		$where .= "(not-join [?${class_var}]";
		$where .= "[?${class_var} :${ident} ?var${criterion_count}]";
		$where .= "[(re-pattern $last_word) ?var${criterion_count}-regex]";
		$where .= "[(re-matches ?var${criterion_count}-regex ?var${criterion_count})]";
		$where .= ")";
		
	      } else {
		# a single non-regex value
		# [:find ?col1-id :where [?col1 :cds/id ?col1-id] [(ground "AC3.3") ?col1-id]]
		$where .= "(not-join [?${class_var}]";
		$where .= "[?${class_var} :${ident} ?var${criterion_count}]";
		$where .= "[(ground $last_word) ?var${criterion_count}]";
		$where .= ")";
	      }
	      
	    } elsif ($next_name eq 'GE') {
	      warn "Not yet implemented '$next_word', in '$criterion'\n";
	      return undef;
	    } elsif ($next_name eq 'LE') {
	      warn "Not yet implemented '$next_word', in '$criterion'\n";
	      return undef;
	    } elsif ($next_name eq 'GT') {
	      warn "Not yet implemented '$next_word', in '$criterion'\n";
	      return undef;
	    } elsif ($next_name eq 'LT') {
	      warn "Not yet implemented '$next_word', in '$criterion'\n";
	      return undef;
	    } else {
	      warn "Found unexpected token '$next_word', type $next_type after TAG in '$criterion'\n";
	      return undef;
	    }
	    
	  } elsif ($obj_type eq 'long' || $obj_type eq 'double') { # tag is numeric value
	    # [(/ ?f-32 1.8) ?celsius]
	    if ($next_name eq 'EQ') {
	      $where .= "[?${class_var} :${ident} ?var${criterion_count}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      
	    } elsif ($next_name eq 'NOT') { # how do we test for NE string?
	      $where .= "(not-join [?${class_var}]";
	      $where .= "[?${class_var} :${ident} ?var${criterion_count}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      $where .= ")";
	      
	    } elsif ($next_name eq 'GE') {
	      $where .= "[(>= ?var{criterion_count} ${ident}) ?${class_var}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      
	    } elsif ($next_name eq 'LE') {
	      $where .= "[(<= ?var{criterion_count} ${ident}) ?${class_var}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      
	    } elsif ($next_name eq 'GT') {
	      $where .= "[(> ?var{criterion_count} ${ident}) ?${class_var}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      
	    } elsif ($next_name eq 'LT') {
	      $where .= "[(< ?var{criterion_count} ${ident}) ?${class_var}]";
	      $where .= "[(ground $last_word) ?var${criterion_count}]";
	      
	    } else {
	      warn "Found unexpected token '$next_word', type $next_type after TAG in '$criterion'\n";
	      return undef;
	    }
	    
	    
	    
	  } elsif (0) { # tag is enum +++
	    
	  } elsif (0) { # tag is in-coming xref +++
	    
	  } elsif (0) { # tag is hash object tag +++
	    
	  }
	  
	  
	} else { # invalid command
	  warn "Found token '$next_word', type $next_type after TAG in '$criterion'\n";
	  return undef;
	}
	
      } elsif ($type eq 'COUNT') {
	# COUNT TAG (NEXT) OP NUM END - the normal case
	# COUNT TAG (NEXT) OP String END - weird - throw an error for this
	# COUNT TAG (NEXT) END - defaults to finding where the TAG exists?

# Thomas's example of finding gene IDs with exactly 2 GO terms
#	  [:find ?gene-id
#          :where [(datomic.api/q
#              '[:find ?gene-id (count ?go-term)
#                :where [?gene :gene/id ?gene-id]
#                       [?gene :gene/go-term ?go-term]]
#               $)
#            [[?gene-id ?go-count]]]
#           [(= ?go-count 2)]] 

	$tag_name_element = $self->get_tag_name($class, $tokens, \$next_token);
	my $objstr;
	my $ident;
	if (defined $tag_name_element) {
	  $objstr = $self->model($class);
	  $ident = $objstr->[$tag_name_element]->{ident};
	} else {
	  warn "Unknown type of tag '$word' in '$criterion'\n";
	  return undef;
	}
	
	if ($next_token+1 >= scalar @{$tokens}) { 	# at end of command, just the TAG so check it exists
	  
	  # COUNT TAG  (NEXT) (tag exists)
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/source-exons ?col2]]
	  # [:find ?name :where [?e :cds/id ?name][?e :cds/from-laboratory _]]
	  # FIND CDS; COUNT Prediction-status (enum )
	  # [:find ?col1-id :where [?col1 :cds/id ?col1-id] [?col1 :cds/prediction-status ?col2] [?col2 :db/ident ?col2-ident]]
	  
	  $where .= "[?${class_var} :${ident} _]"; # check the ident exists and we don't care what the contents are, so use '_'
	  
	} elsif (defined $next_type && $next_type eq 'OP') { # is there an OP after the tag?
	  
	  # COUNT TAG (NEXT) OP NUM
	  
	  $next_token++;
	  if ($next_token >= scalar @{$tokens}) {
	    warn "Expected something after '$next_word' in '$criterion'\n";
	    return undef;
	  }
	  my ($last_name, $last_word, $last_expect_next, $last_type, $last_can_start, $last_can_end) = @{$tokens->[$next_token]};
	  print "parsed command is: $word $next_word $last_word\n";
	  if ($last_type ne 'INT') {
	    warn "Expected integer after '$next_word' in '$criterion'\n";
	    return undef;	    
	  }
	  
	  # if the 'db/valueType' is 'db.type/ref', then look for the next tag/value on the line by going to the ->{right} in the objstr
	  # unless it has 'db/valueType' = 'db.type/ref' and 'pace/obj-ref' pointing to an outgoing xref
	  my $obj_type = $objstr->[$tag_name_element]->{'type'};
	  my $obj_xref = $objstr->[$tag_name_element]->{'xref'};
	  
	  if (!defined $obj_xref || $obj_type ne 'ref') {
	    # go right to the next TAG
	    $tag_name_element = $objstr->[$tag_name_element]->{'right'};
	    if (!defined $tag_name_element) {
	      warn "What? We seem to have fallen off the end of the model when trying to find the tag with the value to use in '$criterion'\n";
	    }
	    $obj_type = $objstr->[$tag_name_element]->{'type'};
	    $obj_xref = $objstr->[$tag_name_element]->{'xref'};
	    $ident = $objstr->[$tag_name_element]->{ident};
	  }
	  
	  if ($next_name eq 'EQ') {
	    $where .= "[(datomic.api/q ";
	    $where .= " '[:find ?col-id${class_var} (count ?countvar${criterion_count})";
	    $where .= "  :where [?local${criterion_count} :${class_var}/id ?col-id${class_var}]";
	    $where .= "         [?local${criterion_count} :${ident} ?countvar${criterion_count}]]";
	    $where .= "  $)";
	    $where .= " [[?col-id${class_var} ?var${criterion_count}-count]]]"; 
            $where .= " [(= ?var${criterion_count}-count $last_word)]]";

	    
	  } elsif ($next_name eq 'NOT') { 
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    
	  } elsif ($next_name eq 'GE') {
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    
	  } elsif ($next_name eq 'LE') {
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    
	  } elsif ($next_name eq 'GT') {
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    
	  } elsif ($next_name eq 'LT') {
	    $where .= "";
	    $where .= "";
	    $where .= "";
	    
	  } else {
	    warn "Found unexpected token '$next_word', type $next_type after TAG in '$criterion'\n";
	    return undef;
	  }
	  
	  
	} elsif ($type eq 'MISSING') { # '!' or 'NOT' TAG - find objects where is doesn't exist
	  # NOT TAG (NEXT) END
	  #[:find ?name :where [?artist :artist/name ?name] [(missing? $ ?artist :artist/startYear)]]
	  #[:find ?name :where [?cds :cds/id ?name][(missing? $ ?cds :cds/species)]
	  $tag_name_element = $self->get_tag_name($class, $tokens, \$next_token);
	  if (defined $tag_name_element) {
	    my $objstr = $self->model($class);
	    my $ident = $objstr->[$tag_name_element]->{ident};
	    $where .= "[?var${criterion_count} :${class_var}/id ?${class_var}]";
	    $where .= "[(missing? \$ ?var${criterion_count} :${ident})]";
	    
	  } else {
	    warn "Unknown type of tag '$word' in '$criterion'\n";
	    return undef;
	  }
	}

      } elsif ($type eq 'FOLLOW') { # expect a incoming or outgoing xref and change the class to match it and find instances of the new class that xref the existing class
	# FOLLOW should cause the $class_var to change
	# FOLLOW TAG (NEXT)
	$tag_name_element = $self->get_tag_name($class, $tokens, \$next_token);
	# +++ want to update the ${class_var} entity variable links all criteria for the class we now expect +++
	# get the new class      
	# $class = whatever 
	# initialise the class_variable to be the new lowerclassed class name
	# $class_var = lc "${class}"; 
	if (defined $classref) {$$classref = $class_var} # return the updated type of class that is now being returned
	$find = ":find [?col-id${class_var} ...] "; # update the col-id
	# +++ to be done
	
      } else {
	warn "Command '$word' not recognised\n";
	return undef;
      }
    }
  }




# +++ removed $in from the query strng
  my $query = '[' . $find . $where . ']' ;

  return $query;
}
#################################################################
# parse a simple constraint on the current search
sub parse_cmd {
  my ($self, $input) = @_;

# the syntax is very simple, so we can parse and lex in one go
# key => regex, ignore, expect after, expect next, type, can start, can end
# this allows the novel command "FOLLOW TAG OP (NUM|STRING) 
# the word 'NEXT' (repeated as necessary) modifies TAG to get subsequent tags on a line
# COUNT TAG
# COUNT TAG NEXT...
# COUNT TAG OP STRING|NUM
# COUNT TAG NEXT... OP NUM
# TAG
# TAG NEXT...
# NOT TAG
# NOT TAG NEXT...
# TAG OP STRING|NUM
# TAG NEXT... OP STRING|NUM
# FOLLOW TAG
# FOLLOW TAG NEXT...

  my @token_def = (
#                  Name           regex              ignore  expect after                 expect next       type      can start can end
		   [Whitespace => qr{\s+},                1, 'ALL',                       'ALL',            '',       0,        1],
		   [COUNT      => qr{\bCOUNT\b}i,         0, 'START',                     'TAG',            'COUNT',  1,        0],
		   [MISSING    => qr{(\!|\bNOT\b)}i,      0, 'START',                     'TAG',            'MISSING',1,        0],
		   [FOLLOW     => qr{\bFOLLOW\b}i,        0, 'START',                     'TAG',            'FOLLOW', 1,        0],
		   [NEXT       => qr{\bNEXT\b}i,          0, 'TAG',                       'NEXT OP END',    'NEXT',   0,        1],
		   [GE         => qr{\>\=},               0, 'TAG NEXT',                  'NUM',            'OP',     0,        0],
		   [LE         => qr{\<\=},               0, 'TAG NEXT',                  'NUM',            'OP',     0,        0],
		   [GT         => qr{\>},                 0, 'TAG NEXT',                  'NUM',            'OP',     0,        0],
		   [LT         => qr{\<},                 0, 'TAG NEXT',                  'NUM',            'OP',     0,        0],
		   [EQ         => qr{\=},                 0, 'TAG NEXT',                  'NUM STRING',     'OP',     0,        0],
		   [NOT        => qr{(\!=|\!|\bNOT\b|\bNE\b)}i, 0, 'TAG NEXT',            'NUM STRING',     'OP',     0,        0],
		   [INT        => qr{\b[0-9]+\b},         0, 'OP',                        'END',            'NUM',    0,        1],
		   [FLOAT      => qr{\b[0-9]+[\.0-9]*\b}, 0, 'OP',                        'END',            'NUM',    0,        1], 
		   [Value_dq   => qr{\".+?\"},            0, 'OP',                        'END',            'STRING', 0,        1],
		   [Value_sq   => qr{\'.+?\'},            0, 'OP',                        'END',            'STRING', 0,        1],
		   [Value_     => qr{\b[\w\*]+\b},        0, 'OP',                        'END',            'STRING', 0,        1],
		   [TAG        => qr{\b[\w]+\b},          0, 'START FOLLOW COUNT MISSING','OP END',         'TAG',    1,        1],
		  );
  
  
  my @tokens;

  print "In parse_cmd: $input\n";
  
  pos($input) = 0;
  my $last = 'START';
  
  while (pos($input) < length $input) {
    my $matched = 0;
    for my $t (@token_def) {
      my ($name, $re, $ignore_flag, $expect_after, $expect_next, $type, $can_start, $can_end) = @$t;
      if (index($expect_after, $last) != -1 || $expect_after eq 'ALL') {
	if ($input =~ m/\G($re)/gc) {
	  $matched = 1;
	  print "in criterion: $input, match $1 as a $name type $type\n";
	  next if $ignore_flag;
	  if ($last eq 'START' && !$can_start) {
	    warn "parse_cmd(): There is an invalid token at the start of the command: '$input'\n'\n";
	    return undef;
	  }
	  $last = $type;
	  my $word = $1;
	  push @tokens, [$name, $word, $expect_next, $type, $can_start, $can_end];
	}
      }
    }
    if (!$matched) {
      warn "parse_cmd(): Syntax error at position " . pos($input) . " of '$input'\n";
      return undef;
    }
  }

  if ($#tokens == -1) {
    warn "parse_cmd(): There are no recognised tokens in the command: '$input'\n'\n";
    return undef;
  }
  if ($tokens[-1][-1] != 1) {
    warn "parse_cmd(): There is no valid end to the command: '$input'\n'\n";
    return undef;
  }

  return \@tokens;
   
}

#############################################################################

=head2 

    Title   :   get_tag_name
    Usage   :   $tag_name_element = $db->get_tag_name($class, $tokens, \$next_token);
    Function:   get the class schema @{$objstr} element number of the required tag in a class
    Returns :   int array element of the class schema @{$objstr} which has the required tag
                or 'undef' if tag is not found
    Args    :   $class = name of class
                $tokens = array-ref describing tokens in the search criteria
                $next_token = ref of int of next array element in $tokens to look at

+++ Doesn't yet navigate to Evidence hash tags 


=cut


sub get_tag_name {

  my ($self, $class, $tokens, $next_token) = @_;

  my $element = undef;

  my $objstr = $self->model($class);

  # check there are more things in the $tokens array-ref starting at $next_token
  if (scalar @{$tokens} <= $$next_token) {
    warn "Expected a TAG name - it appears to be missing\n";
    return undef;
  }

  # get first tag name
  my ($tag, $tag_name) = @{$tokens->[$$next_token++]};

  # search through @{$objstr} for a match of the tag name to either the Ace or the Datomic tag names
  # valid tag names are: 
  #  the last word of a 'pace/tags' list
  #  the 'ident' name
  #  either or the above with a virtual tag name appended (including the '^'), e.g. 'Database^field'

  my $tag_name_virtual_half = '';
  if ($tag_name =~ /(\S+)\^(\S+)/) {
    $tag_name = $1;
    $tag_name_virtual_half = "/".$2;
  }
  my $datomic_tag_name = $class.'.'.lc($tag_name).$tag_name_virtual_half;

  my $found = 0;
  for ($element = 0; $element < scalar @{$objstr}; $element++) {
    if ($tag_name eq $objstr->[$element]->{name}) {$found = 1; last} 
    if ($tag_name eq $objstr->[$element]->{ident}) {$found = 1; last} 
  }
  if (!$found) {
    warn "Can't find the tag $tag_name in $class\n";
    return undef;
  }

  # repeat:   if next token is 'NEXT' then move 'right' in objstr to the next virtual tag
  while (scalar @{$tokens} > $$next_token && $tokens->[$$next_token]->[0] eq 'NEXT') {
    $$next_token++;
    if (defined $objstr->[$element]->{right}) {
      $element = $objstr->[$element]->{right};
    } else {
      warn "Off end of class model looking NEXT after tag $tag_name in $class\n";
      return undef;
    }
    # +++ will probably have to adjust the element to go to the next {right} if we are on a xref with the real xref following it
  }


  print "\n\nFound tag '$tag_name' in class '$class' model. It has ident: '$objstr->[$element]->{ident}'\n";

  return $element;

}

#############################################################################

=head2 

    Title   :   protect
    Usage   :   $pattern = $db->protect($pattern);
    Function:   change the wildcards in the pattern into a perl or java.util.regex.Pattern regex
    Returns :   perl regex
    Args    :   $pattern - string - optionally containing '.', '?', '*' to be protected

=cut

sub protect {

  my ($self, $pattern) = @_;

  $pattern =~ s/\./\\\\./g; # quote dots
  $pattern =~ s/\*/\.\*/g; # change * to .*
  $pattern =~ s/\?/\.\?/g; # change ? to .?

  return $pattern;
}

########################################################################################################################
# ACE - like routines - these try to stick fairly closely to the old Ace methods - See: /software/worm/lib/site_perl/5.16.2/Ace.pm
########################################################################################################################

=head2 

    Title   :   fetch
    Usage   :   @results = $db->fetch('CDS', 'AC3.3*', 'Method = curated', \$count, \$returned_class)
    Function:   gets a list of the names of the members of the specified class that match the patetrn and the optional query criteria
    Returns :   number of results or sorted ID names of matching objects in array.
    Args    :  $class - string - the class of object to get (datomic or Ace names can be used)
               $pattern - string - pattern to match to ID name of the required object
               $query - string - criteria to restrict the class/pattern results
               $countref - ref to int - returns the count of objects (optional)
               $classref - ref to string - returns the datomic name of the class of object returned (optional)
                
=cut

sub fetch {

  my ($self, $class, $pattern, $query, $countref, $classref) = @_;

  $pattern ||= '*';

  my ($count, @results);
  if (defined $countref) {$$countref = 0};
  if (defined $classref) {$$classref = $class};

  if ($class eq 'Class') { # should this be a case-independent match?
    $self->get_classes();
    my @classes = keys %{$self->{classes}};
    $pattern = $self->protect($pattern);
    @results = sort {$a cmp $b} grep /$pattern/i, @classes;
    $count = scalar @results;
  } else {

    print "in fetch() class name = $class";
    my $datomic_class = $self->get_class_name_in_datomic($class);
    if (!defined $datomic_class) {die "Can't find datomic name of $class\n"}

    my $querystr = $self->ace_to_querystr($datomic_class, $pattern, $query, $classref);
    if (defined $querystr) {
      print "In fetch() query:\n$querystr\n";
      ($count, @results) = $self->query($querystr);
      @results =  sort {$a cmp $b} @{$results[0]};
      if (defined $countref) {$$countref = $count}
    }
  }

  return $count if !wantarray;
  if ($count) {
    return @results;
  } else {
    return ();
  }
}
















1;

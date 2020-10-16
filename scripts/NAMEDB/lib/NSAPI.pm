#!/software/bin/perl -w

# NSAPI - nameserver API
# methods for talking to the datomic Nameserver REST API
# Gary Williams 2019-02-12

# see https://names.wormbase.org/api-docs/index.html
# for the REST documentation

package NSAPI;
use lib '/software/worm/lib/perl';
use Carp;
use warnings;
use strict;
use LWP::UserAgent;
use JSON qw(encode_json decode_json);
#use URI::Encode qw(uri_encode uri_decode);
use CGI; # for making ISO-format dates
use Data::Dumper;


# useful command for Matt to debug the output error messages:
# curl -X "PUT" -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": {"cgc-name": "new-1"}, "prov": { "why": "Testing" }}' https://names.wormbase.org/api/gene/abc-20 | ~/.local/bin/jq ".problems" | xargs -0 echo -e


# my equivalent in perl
#curl -X GET -H "Authorization: Token $TOKEN" -H 'Accept: application/json' -v 'https://names.wormbase.org/api/gene/WBGene00111112' | perl -ne '{$_ =~ /biotype\"\:\"(\S+?)\"/; print $1,"\n"}'

# get info on a gene - WBGeneID, Sequence-name or CGC-name
#curl -X GET -H "Authorization: Token $TOKEN" -H 'Accept: application/json' -v 'https://names.wormbase.org/api/gene/WBGene00027400'


# create an uncloned gene with a CGC-name - CGC-name and Species are required
#curl -X POST -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": [{"cgc-name": "abc-31", "species": "Caenorhabditis elegans"}], "prov": {"why": "Testing"}}' -v https://names.wormbase.org/api/batch/gene

# create a cloned gene with a sequence-name - Sequence-name, Biotype and Species are required
#curl -X POST -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": [{"sequence-name": "AC3.18", "biotype": "cds","species": "Caenorhabditis elegans"}], "prov": {"why": "Testing"}}' -v https://names.wormbase.org/api/batch/gene

# delete a gene by sequence-name
#curl -X DELETE -H "Authorization: Token $TOKEN" -H "content-type: application/json" -v --data '{"prov": {}}' "https://names.wormbase.org/api/gene/CCCC.1"

# delete a cgc-name, leaving the gene untouched.
#curl -X DELETE -H "Authorization: Token $TOKEN" -H "content-type: application/json" -v --data '{"data": ["CCCC.1"], "prov": {}}' "https://names.wormbase.org/api/batch/gene/cgc-name"


# delete a batch of genes
#curl -X DELETE -H "Authorization: Token $TOKEN" -H "content-type: application/json" -v --data '{"data": [{"sequence-name":"CCCC.1"},{"sequence-name":"AC3.3"}], "prov": {}}' "https://names.wormbase.org/api/batch/gene"

# merge into a gene with no CGC-name from a gene with a CGC-name ++++ Get HTTP/1.1 500 Server Error
#curl -X POST -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": {"biotype": "cds"}, "prov": {"why": "Testing"}}' -v https://names.wormbase.org/api/gene/AC3.1/merge/AC3.5

# update gene to change biotype to 'transcript'  - should work but doesn't  - GitHub Issue #177
#curl -X PUT -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": {"biotype": "transcript", "species": "Caenorhabditis elegans"}, "prov": {"why": "Testing"}}' -v https://names.wormbase.org/api/gene/ZK1320.15

# create some Features
#curl -X POST -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": {"n": 10}}' -v https://names.wormbase.org/api/batch/sequence-feature

# recent variation changes in the console data
#curl -X GET -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"prov": {}}' "https://names.wormbase.org/api/recent/variation/console?from=2018-07-24&until=2019-07-26" 

# delete a variation
#curl -X DELETE -H "Authorization: Token $TOKEN" -H "content-type: application/json" -v --data '{"prov": {}}' "https://names.wormbase.org/api/entity/variation/WBVar02149574"

# update the public name of a variation
#curl -X PUT -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"data": {"name": "lazy1"},"prov": {}}'  https://names.wormbase.org/api/entity/variation/WBVar02149581

#resurrect a dead variation
#curl -X POST -H "Authorization: Token $TOKEN" -H "content-type: application/json" --data '{"prov": {}}'  https://names.wormbase.org/api/entity/variation/WBVar02149574/resurrect

# typical gene values:
#$VAR1 = {
#         'history' => [
#                         {
#                           'who' => {
#                                                 'id' => 'WBPerson1971',
#                                                 'name' => 'Keith Bradnam'
#                                               },
#                           'when' => '2004-04-07T11:29:19Z',
#                           'changes' => [],
#                           'what' => 'new-gene',
#                           'how' => 'importer'
#                         }
#                       ],
#          'status' => 'live',
#          'id' => 'WBGene00000001',
#          'sequence-name' => 'Y110A7A.10',
#          'biotype' => 'cds',
#          'cgc-name' => 'aap-1',
#          'species' => 'Caenorhabditis elegans'
#        };

#
# The valid {entity-type} values for the batch API entry points are 'gene' 'variation' and 'sequence-feature'
# 
#

=over 4

=item $db = NameDB->ping

Show the system is working.

=cut

sub ping  { return "Bonjour!" }

=over 4

=item $db = NameDB->connect($test)

Connect to the database.
Returns a NSAPI object if successful.

=cut


#--------------------- connect, disconnect -----------------
# get the Google authentication token ready to make REST calls

# there is an optional 'test' argument - if this is true, then the
# 'test' NameServer database will be used instead of the live server.

sub connect {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  my $test = shift;
  $self->test($test);

  $self->{'TOKEN'} = $self->get_token();

  $self->noise(0); # turn off the debugging output of curl

  return $self;
}

=item $db->disconnect

A null method, just in case someone feels the need to balance the connect() with a disconnect() call.

=cut

sub disconnect {
  my $self = shift;
}

=item $db->noise(0)

Set or reset and return the level of debugging output

=cut

sub noise {
  my $self = shift;
  my $level = shift;
  $self->{'noise'} = $level if (defined $level);
  return $self->{'noise'};
}

=item $db->test(0)

Set or reset and return whether we wish to use the test server or not.

=cut

sub test {
  my $self = shift;
  my $level = shift;
  if (defined $level) {
    $self->{'test'} = $level;
  }
  return $self->{'test'};
}

#======================================================================
# SINGLE GENE operations
#======================================================================
# merge_genes($from, $into, $biotype, $why)
# merges one pair of genes 
# /api/gene/{identifier}/merge/{from-identifier}
# Args:
# $from - string, gene identifier of gene to die
# $into - string, gene identifier of gene to live
# $biotype - string, biotype of resulting merged gene, defaults to 'CDS' if undefined, one of 	'cds'	'transcript'	'pseudogene'	'transposon'
# $why - string, optional reason for the merger, ignored if undefined

# This gets an error if you try to specify a gene to die that has a CGC-name

sub merge_genes {
  my ($self, $from, $into, $biotype, $why) = @_;

  if (!defined $biotype || $biotype eq '') {
    $biotype = 'cds';
  }
  
  my %biotypes = ( # input     =>  biotype name in Nameserver
		  'cds'        => 'cds',
		  'transcript' => 'transcript',
		  'pseudogene' => 'pseudogene',
		  'transposon' => 'transposable-element-gene',
		  'transposable-element-gene' => 'transposable-element-gene',
		 );
  
  if (!exists $biotypes{lc $biotype}) {die "Unknown biotype specified in merge_genes($from, $into): $biotype\n"}
  $biotype = $biotypes{lc $biotype};


  my $payload = '{"data": {';
  $payload .= '"biotype": "'.$biotype.'"'; #was gene/biotype
  $payload .= '}, ';
  $payload .= '"prov": {';
  $payload .= '"why": "'.$why.'"' if (defined $why); #was provenance/why
  $payload .= '}}';
  
  if ($self->noise()) {print $payload,"\n"}
  
  my $content = $self->curl("POST", "gene/$into/merge/$from", $payload);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}

#======================================================================
# split_genes($gene, $biotype_orig, $sequence_new, $biotype_new, $other_names_new, $why)
# Split a gene, creating a new gene as a product.
# POST /api/gene/{identifier}/split
# Args:
# $gene - string, gene identifier of gene to split
# $biotype_orig - string, biotype to assign to original gene, defaults to 'CDS' if undefined, one of 	'cds'	'transcript'	'pseudogene'	'transposon'
# $sequence_new - string, Sequence-name of the gene to create from the split
# $biotype_new - string, biotype to assign to the newly created gene, defaults to 'CDS' if undefined, one of 	'cds'	'transcript'	'pseudogene'	'transposon'
# $other_names_new - array ref of strings - other_names of the new gene, or undef if not specified
# $why - string, optional reason for the split, ignored if undefined
#
# Returns:
#       {
#          'created' => {
#                         'id' => 'WBGene00304801' #was gene/id
#                       },
#          'updated' => {
#                         'id' => 'WBGene00000001' #was gene/id
#                       }
#       };



sub split_genes {
  my ($self, $gene, $biotype_orig, $sequence_new, $biotype_new, $other_names_new, $why) = @_;
  
  if (!defined $biotype_orig || $biotype_orig eq '') {
    $biotype_orig = 'cds';
  }
  if (!defined $biotype_new || $biotype_new eq '') {
    $biotype_new = 'cds';
  }
  
  my %biotypes = ( # input     =>  biotype name in Nameserver
		  'cds'        => 'cds', #was biotype/cds
		  'transcript' => 'transcript', #was biotype/transcript
		  'pseudogene' => 'pseudogene', #was  biotype/pseudogene
		  'transposon' => 'transposable-element-gene', #was biotype/transposable-element-gene
		  'transposable-element-gene' => 'transposable-element-gene', #was biotype/transposable-element-gene
		 );
  
  if (!exists $biotypes{lc $biotype_orig}) {die "Unknown biotype specified in split_genes($gene): $biotype_orig\n"}
  $biotype_orig = $biotypes{lc $biotype_orig};
  
  if (!exists $biotypes{lc $biotype_new}) {die "Unknown biotype specified in split_genes($gene): $biotype_new\n"}
  $biotype_new = $biotypes{lc $biotype_new};

  my $payload = '{"data": {';
  $payload .= '"biotype": "'.$biotype_orig.'",'; #was gene/biotype
  $payload .= '"product": {';
  $payload .= '"sequence-name": "'.$sequence_new.'",'; #was gene/sequence-name
  $payload .= '"biotype": "'.$biotype_new.'"'; #was gene/biotype
  if (defined $other_names_new && scalar @{$other_names_new}) {
    $payload .= ',';
    $payload .= '"other-names": [';
    $payload .= join(',', map { "\"$_\"" } @{$other_names_new});
    $payload .= ']';
  }
  $payload .= '}},';
  $payload .= '"prov": {';
  $payload .= '"why": "'.$why.'"' if (defined $why) ; #was provenance/why
  $payload .= '}}';
  
  if ($self->noise()) {print $payload,"\n"}
  
  my $content = $self->curl("POST", "gene/$gene/split", $payload);
  if ($self->noise()) {print Dumper $content};

  return $content;
}
#======================================================================
# unmerge_genes()
# unmerge genes

# unmerges one or more sets of genes 
# /api/gene/{identifier}/merge/{from-identifier}
# Args:
# $from - string, gene identifier of gene that is dead
# $into - string, gene identifier of gene that is live
# $why - string, optional reason for the merger, ignored if undefined

sub unmerge_genes {
  my ($self, $from, $into, $why) = @_;

  my $payload = '{"data": {';
  $payload .= '}, ';
  $payload .= '"prov": {';
  $payload .= '"why": "'.$why.'"' if (defined $why); #was provenance/why
  $payload .= '}}';
  
  if ($self->noise()) {print $payload,"\n"}
  
  my $content = $self->curl("DELETE", "gene/$into/merge/$from", $payload);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}


#======================================================================
# kill_gene($identifier, $why)
# $identifier should be a single gene identifier (one of WBGeneID, Sequence-name, CGC-name)
# DELETE /api/gene/{identifier}

# Args:
# $gene - string, gene identifier of the gene to delete
# $why - string, required reason for the deletion


sub kill_gene {
  my ($self, $gene, $why) = @_;
  
  my $payload = '{"prov": {';
  $payload .= '"why": "'.$why.'"'; #was provenance/why
  $payload .= '}}';
  
  if ($self->noise()) {print $payload,"\n"}
  
  my $content = $self->curl("DELETE", "gene/$gene", $payload);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}

#======================================================================
# undo($op_id, $why)
# $op_id should be a value returned from the NameService in response to another operation.
# Operations can be undone when they are the last operation to occur for a given entity type (e.g Gene).
# $why is optional remark text.

#======================================================================
# find_genes($pattern)
# returns a hash ref of any genes whose name matches the input pattern.
# pattern can be any of a regular expression or part of the name of a WBGeneID, CGC name or a Sequence_name
# patterns are case-sensitive
# the pattern is contrained to match at the start of the name, use '.*' at the start of the pattern to match anywhere in the name
# example patterns:
# WBGene00000001
# unc-111
# unc
# AC3.3
# AC3
# WBGene0000000[1-3]
# .*-aat
# Bma-

# Returns hash-ref of list of hashes
# $VAR1 = {
#           'matches' => [
#                         {
#                            'id' => 'WBGene00000003', #was gene/id
#                            'sequence-name' => 'F07C3.7', #was gene/sequence-name
#                            'cgc-name' => 'aat-2', #was gene/cgc-name
#                          },
#                          {
#                            'id' => 'WBGene00000002', #was gene/id
#                            'sequence-name' => 'F27C8.1', #was gene/sequence-name
#                            'cgc-name' => 'aat-1', #was gene/cgc-name
#                          },
#                        ]
#         };


sub find_genes {
  my ($self, $pattern) = @_;

  my $encoded = CGI::escape($pattern);
  
  my $content = $self->curl("GET", "gene?pattern=$encoded");
  if ($self->noise()) {print Dumper $content}
  return $content
}

#======================================================================
# update_gene($identifier)
# updates a single specified gene - it can be specified by WBGeneID, CGC_name or Sequence_name
# it requires an exact match of the specified identifier
# PUT /api/gene/{identifier}
# Args: 
#       $identifier - string, name of the gene
#       $species    - string, full species name (but Matt agrees that this is superfluous for identification purposes and may remove it)
#       $CGC_name   - string to replace existing CGC_name
#       $sequence_name - string, to replace existing sequence_name
#       $biotype    - string, to replace existing biotype
#       $force      - optional -  set this to true to force an invalid CGC name or Sequence name that is not in the normal format and normal validation will be turned off

# Can't delete data here by using a null or blank string - will have to use Delete in /api/batch/gene/cgc-name for example

sub update_gene {
  my ($self, $identifier, $species, $CGC_name, $sequence_name, $biotype, $force) = @_;


  my $encoded = CGI::escape($identifier);
  
  # we currently need to specify the species as well as the gene ID - this is superfluous and should be changed in the Nameserver
  my $data = '"data": {';
  $data .= ' "cgc-name": "'.$CGC_name.'",' if (defined $CGC_name && $CGC_name ne ''); #was gene/cgc-name
  $data .= ' "sequence-name": "'.$sequence_name.'",' if (defined $sequence_name && $sequence_name ne ''); #was gene/sequence-name
  $data .= ' "biotype": "'.$biotype.'",' if (defined $biotype && $biotype ne ''); #was gene/biotype
  $data .= ' "species": "'.$species.'"}'; #was gene/species


  my $payload = '{'.$data.', "prov": { }'; # we don't have to provide a provenenance, but we do have to show we have an empty one
  $payload .= ', "force": "1"' if (defined $force && $force);
  $payload .= '}';

  if ($self->noise()) {print "payload: $payload\n"}

  my $content = $self->curl("PUT", "gene/$encoded", $payload);
  if ($self->noise()) {print Dumper $content}
  return $content
}

#======================================================================
# delete_gene($identifier)
# deletes a single specified gene - it can be specified by WBGeneID, CGC_name or Sequence_name
# it requires an exact match of the specified identifier
# DELETE /api/gene/{identifier}
# Args: 
#       $identifier - string, name of the gene
#       $why        - string, reason for deleting the gene


sub delete_gene {
  my ($self, $identifier, $why) = @_;

  $why =~ s/\"//g; # remove quotes from the payload string
  if (!defined $why || $why eq '') {die "the reason for deleting the gene must be given\n"}

  my $encoded = CGI::escape($identifier);
  
  my $payload = '{ "prov": { "why": "'.$why.'" }}'; #was provenance/why
  if ($self->noise()) {print "payload: $payload\n"}

  my $content = $self->curl("DELETE", "gene/$encoded", $payload);
  if ($self->noise()) {print Dumper $content}
  return $content
}


#======================================================================
# info_gene($gene_identifier)
# get all details of a single named gene
# the identifier can be the name of a WBGeneID, CGC name or a Sequence_name
# The gene name is case-sensitive
# returns hash-ref of details. If the ID doesn't exist then it returns 'undef'.
#
# $VAR1 = {
#           'history' => [
#                          {
#                            'who' => { #was provenance/who
#                                                  'id' => 'WBPerson1971', #was person/id
#                                                  'name' => 'Keith Bradnam' #was person/name
#                                                },
#                            'when' => '2004-04-07T11:29:19Z', #was provenance/when
#                            'changes' => [],
#                            'what' => 'new-gene', #was provenance/what   event/new-gene
#                            'how' => 'importer' #was provenance/how   agent/importer
#                          }
#                        ],
#           'status' => 'live', #was gene/status   gene.status/live
#           'id' => 'WBGene00000001', #was gene/id
#           'sequence-name' => 'Y110A7A.10', #was gene/sequence-name
#           'biotype' => 'cds', #was gene/biotype  biotype/cds 
#           'cgc-name' => 'aap-1', #was gene/cgc-name
#           'species' => 'Caenorhabditis elegans' #was gene/species
#         };

sub info_gene {
  my ($self, $gene) = @_;

  my $encoded = CGI::escape($gene);

  my $payload = undef; # don't want to pass in a payload, but want a dummy parameter so $not_found is set correctly
  my $not_found = 1; # don'y throw an error if the gene is not found
  my $content = $self->curl("GET", "gene/$encoded", $payload, $not_found);
  
  if ($self->noise()) {print Dumper $content}
  return $content

}

#======================================================================
#======================================================================
# PERSON operations
#======================================================================
# new_person
# POST /api/person

# Args:
# $data - hash-ref with values being the new data to update with
#         the hashes contain one or more of the keys ('WBPerson', 'Email', 'Name').
#         the email address must end "@wormbase.org"
# e.g.
# $data = {"Email" => "jbloggs@wormbase.org", "Name" => "Joe Bloggs", "WBPerson" => "WBPerson 0000001"};

#
sub new_person {
  my ($self, $data) = @_;
  
  if ($self->noise()) {print "new person","\n"}
  

#{
#  "active?": true,
#  "email": "matthew.russell@wormbase.org",
#  "id": "WBPerson33035",
#  "name": "Matt Russell", 
#} 
  my $payload = '{"data": {"active?": true,';
  
  my $email = $data->{'Email'};
  my $wbperson = $data->{'WBPerson'};
  my $name = $data->{'Name'};
  
  $payload .= '"email": "'.$email.'",' if (defined $email);
  $payload .= '"id": "'.$wbperson.'",' if (defined $wbperson);
  $payload .= '"name": "'.$name.'",' if (defined $name);
  
  chop $payload; # remove last ','
  $payload .= '}, "prov": {}}';

  if ($self->noise()) {print $payload,"\n"}
  
  my $content = $self->curl("POST", "person", $payload);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}

#======================================================================
# kill_person($identifier)
# $identifier should be a single Person identifier (one of WBPersonID, email, Person's name)
# DELETE /api/person/{identifier}

# Args:
# $person - string, person identifier of the person to delete

sub kill_person {
  my ($self, $person) = @_;

  if ($self->noise()) {print "delete $person","\n"}

  my $encoded = CGI::escape($person);

  my $content = $self->curl("DELETE", "person/$encoded");
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}
#======================================================================
# info_person($identifier)
# summarise information on a single individual
# $identifier should be a single Person identifier (one of WBPersonID, email)
# GET /api/person/{identifier}

# Args:
# $person - string, person identifier of the person to look at
# If the person is not found then it returns 'undef'.

sub info_person {
  my ($self, $person) = @_;
  
  if ($self->noise()) {print "info $person","\n"}

  my $encoded = CGI::escape($person);

  my $payload = undef; # don't want to pass in a payload, but want a dummy parameter so $not_found is set correctly
  my $not_found = 1; # don't throw an error if the person is not found

  my $content = $self->curl("GET", "person/$encoded", $payload, $not_found);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}
#======================================================================
# update_person($identifier)
# $identifier should be a single Person identifier (one of WBPersonID, email, Person's name)
# PUT /api/person/{identifier}

# Args:
# $person - string, person identifier of the person to update
# $data - hash-ref with values being the new data to update with
#         the hashes contain one or more of the keys ('WBPerson', 'Email', 'Name').
#         the email address must end "@wormbase.org"
# e.g.
# $data = {"WBPerson" => "jbloggs@wormbase.org"};

#
sub update_person {
  my ($self, $person, $data) = @_;
  
  if ($self->noise()) {print "update $person","\n"}
  

#{
#  "active?": true,
#  "email": "matthew.russell@wormbase.org",
#  "id": "WBPerson33035",
#  "name": "Matt Russell", 
#} 
  my $payload = '{"data": {"active?": true,';
  
  my $email = $data->{'Email'};
  my $wbperson = $data->{'WBPerson'};
  my $name = $data->{'Name'};
  
  $payload .= '"email": "'.$email.'",' if (defined $email);
  $payload .= '"id": "'.$wbperson.'",' if (defined $wbperson);
  $payload .= '"name": "'.$name.'",' if (defined $name);
  
  chop $payload; # remove last ','
  $payload .= '}, "prov": {}}';

  if ($self->noise()) {print $payload,"\n"}

  my $encoded = CGI::escape($person);

  my $content = $self->curl("PUT", "person/$encoded", $payload);
  if ($self->noise()) {print Dumper $content}
  
  return $content;
}

#======================================================================
#======================================================================
#======================================================================
#======================================================================
# RECENT operations
#======================================================================
#
# These are used to find gene and variation changes that can then be
# read into the geneace ACeDB database
#
#======================================================================

#======================================================================
# recent_gene(from, until, agent)
# get all details of any gene changes between two dates, inclusive

# The dates are in ISO date format "YYYY-MM-DDThh:mm:ssZ" and times are in the UTC time-zone e.g. `date --utc +%Y-%m-%dT%H:%M:%SZ`
# The agent is one of "agent/web" or "agent/console" (the "Agent" type that made the request). This defaults to both agent types.

# If a call to recent/gene/$agent is made twice, then the second call
# will return a '304' error number to indicate that nothing has
# changed.


sub recent_gene {
  my ($self, $from_date, $until_date, $agent) = @_;

  if ($agent ne 'web' && $agent ne 'console') {die "Invalid agent specified: $agent\n";}

  my $encoded_from_date = CGI::escape($from_date);
  my $encoded_until_date = CGI::escape($until_date);
  my $encoded_agent = CGI::escape($agent);

  my $type = "recent/gene/$agent?";
  if (defined $from_date) {$type .= "from=$encoded_from_date"}
  if (defined $until_date) {$type .= "\&until=$encoded_until_date"}

  my $content = $self->curl("GET", $type);

  if ($self->noise()) {print Dumper $content}
  return $content

}

#======================================================================
# recent_variation(from, until, agent)
# get all details of any variation changes between two dates, inclusive

# The dates are in ISO date format "YYYY-MM-DDThh:mm:ssZ" and times are in the UTC time-zone e.g. `date --utc +%Y-%m-%dT%H:%M:%SZ`
# The agent is one of "agent/web" or "agent/console" (the "Agent" type that made the request). This defaults to both agent types.

# If a call to recent/gene/$agent is made twice, then the second call
# will return a '304' error number to indicate that nothing has
# changed.


sub recent_variation {
  my ($self, $from_date, $until_date, $agent) = @_;

  if ($agent ne 'web' && $agent ne 'console') {die "Invalid agent specified: $agent\n";}

  my $encoded_from_date = CGI::escape($from_date);
  my $encoded_until_date = CGI::escape($until_date);
  my $encoded_agent = CGI::escape($agent);

  my $type = "recent/variation/$agent?";
  if (defined $from_date) {$type .= "from=$encoded_from_date"}
  if (defined $until_date) {$type .= "\&until=$encoded_until_date"}

  my $content = $self->curl("GET", $type);

  if ($self->noise()) {print Dumper $content}
  return $content

}



#======================================================================
# recent_strain(from, until, agent)
# get all details of any strain changes between two dates, inclusive

# The dates are in ISO date format "YYYY-MM-DDThh:mm:ssZ" and times are in the UTC time-zone e.g. `date --utc +%Y-%m-%dT%H:%M:%SZ`
# The agent is one of "agent/web" or "agent/console" (the "Agent" type that made the request). This defaults to both agent types.

# If a call to recent/gene/$agent is made twice, then the second call
# will return a '304' error number to indicate that nothing has
# changed.


sub recent_strain {
  my ($self, $from_date, $until_date, $agent) = @_;

  if ($agent ne 'web' && $agent ne 'console') {die "Invalid agent specified: $agent\n";}

  my $encoded_from_date = CGI::escape($from_date);
  my $encoded_until_date = CGI::escape($until_date);
  my $encoded_agent = CGI::escape($agent);

  my $type = "recent/strain/$agent?";
  if (defined $from_date) {$type .= "from=$encoded_from_date"}
  if (defined $until_date) {$type .= "\&until=$encoded_until_date"}

  my $content = $self->curl("GET", $type);

  if ($self->noise()) {print Dumper $content}
  return $content

}






#======================================================================
#======================================================================
#======================================================================
# BATCH GENES operations
#======================================================================
# new_genes(\@data, $species)
# batch creates one or more new uncloned genes (they just have a CGC-name and no Sequence-name) or cloned genes (they have a Sequence name)
# Args:
# $data - array-ref of hashes
#         the hashes contain either the keys ('species', 'gcg-name', 'biotype', 'other-names') for an uncloned gene 
#                            or ('species', 'sequence-name', 'biotype', and optionally 'cgc-name', 'other-names') for a cloned gene.
# $why - string, optional reason for creating the genes
# e.g.
# $data = [{"cgc-name" => "abc-31", "species" => "Caenorhabditis elegans", "biotype" => 'CDS', "other-names" => ["other1", "other2"]}, { ... }];
# or $data = [{"sequence-name" => "ZK1320.13", "species" => "Caenorhabditis elegans", "biotype" => 'CDS', "other-name" => ["other1", "other2"]}, { ... }];
#

sub new_genes {
  my ($self, $data, $why) = @_;
  
  if (!defined $why) {$why = ''}
  
  
#                    input     =>  biotype name in Nameserver
  my %biotypes = ( 
		  'cds'        => 'cds', #was biotype/cds
		  'transcript' => 'transcript', #was biotype/transcript
		  'pseudogene' => 'pseudogene', #was biotype/pseudogene
		  'transposon' => 'transposable-element-gene', #was biotype/transposable-element-gene
		  'transposable-element-gene' => 'transposable-element-gene', #was biotype/transposable-element-gene
		 );
  
  
  my $payload = '{"data": [';
  foreach my $gene_data (@{$data}) {
    my $species = $gene_data->{'species'};
    my $cgc_name = $gene_data->{'cgc-name'};
    my $sequence_name = $gene_data->{'sequence-name'};
    my $biotype;
    if (exists $gene_data->{'biotype'}) {
      $biotype = $biotypes{lc $gene_data->{'biotype'}};
    } 
    my @other_names = @{$gene_data->{'other-names'}};
    
    $payload .= '{';
    $payload .= '"species": "'.$species.'",' if (defined $species); #was gene/species
    $payload .= '"cgc-name": "'.$cgc_name.'",' if (defined $cgc_name); #was gene/cgc-name
    $payload .= '"sequence-name": "'.$sequence_name.'",' if (defined $sequence_name); #was gene/sequence-name
    $payload .= '"biotype": "'.$biotype.'",' if (defined $biotype); #was gene/biotype

    if (scalar @other_names) {
      $payload .= '"other-names": [';
      $payload .= join(',', map { "\"$_\"" } @other_names);
      $payload .= '],';
    }

    chop $payload; # remove last ','

    $payload .= '},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  # returns 
  # {"id":"5dcd84f2-f17a-44c6-8a2c-e096983181f3","ids":[{"id":"WBGene00305168","cgc-name":"abc-1231"}]}
  return $self->batch('POST', 'gene', $payload);
}

#======================================================================
# update_genes($data)
# batch updates one or more genes to add or change CGC name, Sequence name, Biotype
# status - done
#
# Args:
# $data - array-ref of hashes
#         the hashes contain the keys ('id', then any of: 'species', 'gcg-name', 'sequence-name', 'biotype')
#         the WBGene id specifies the gene to update
# $why - string, optional reason for creating the genes
# e.g.
# $data = [{"id" => "WBGene000001", "cgc-name" => "abc-31", "species" => "Caenorhabditis elegans", "biotype" => 'CDS'}, { ... }];
#
# to force an invalid CGC name or Sequence name that is not in the normal fomat, add 'force'=>1 to the $data and normal validation will be turned off
#
# Returns:
# 


sub update_genes {
  my ($self, $data, $why) = @_;
  
  if (!defined $why) {$why = ''}
  
  
  my %biotypes = ( # input     =>  biotype name in Nameserver
		  'cds'        => 'cds', #was biotype/cds
		  'transcript' => 'transcript', #was biotype/transcript 
		  'pseudogene' => 'pseudogene', #was biotype/pseudogene
		  'transposon' => 'transposable-element-gene', #was biotype/transposable-element-gene
		  'transposable-element-gene' => 'transposable-element-gene', #was biotype/transposable-element-gene
		 );
  
  
  my $force;
  my $payload = '{"data": [';
  foreach my $gene_data (@{$data}) {
    my $id = $gene_data->{'id'};
    my $species = $gene_data->{'species'};
    my $cgc_name = $gene_data->{'cgc-name'};
    my $sequence_name = $gene_data->{'sequence-name'};
    my $biotype = $gene_data->{'biotype'};
    if (defined $biotype) {
      if (!exists $biotypes{lc $biotype}) {die "Unknown biotype specified in new_genes(): $biotype\n"}
      $biotype = $biotypes{lc $biotype};
    }
    $force = $gene_data->{'force'}; 

    $payload .= '{';
    $payload .= '"id": "'.$id.'",'; #was gene/id
    $payload .= '"species": "'.$species.'",' if (defined $species); #was gene/species
    $payload .= '"cgc-name": "'.$cgc_name.'",' if (defined $cgc_name); #was gene/cgc-name
    $payload .= '"sequence-name": "'.$sequence_name.'",' if (defined $sequence_name); #was gene/sequence-name
    $payload .= '"biotype": "'.$biotype.'",' if (defined $biotype); #was gene/biotype
    
    chop $payload; # remove last ','
    $payload .= '},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}';
  $payload .= ', "force": "1"' if (defined $force && $force);
  $payload .= '}';
  if ($self->noise()) {print $payload,"\n"}
  
  return $self->batch('PUT', 'gene', $payload);
  
}

#======================================================================
# kill_genes($ids)
# batch deletes one or more genes
#
# Args:
# $ids - array-ref of gene IDs or sequence-names of CGC-names
# $why - string, optional reason for deleting the genes


sub kill_genes {
  my ($self, $ids, $why) = @_;
  
  if (!defined $why) {$why = ''}
  
  
  my $payload = '{"data": [';
  $payload .= join(',', map { ($_ =~ /^WBGene/) ? ("{\"id\": \"$_\"}") : (($_ =~ /-/) ? ("{\"cgc-name\": \"$_\"}") : ("{\"sequence-name\": \"$_\"}"))  } @{$ids});
  $payload .= '], "prov": {';
  $payload .= '"why": "$why"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  return $self->batch('DELETE', 'gene', $payload);
  
}

#======================================================================
# remove_gene_names($data)
# removes gene CGC-names (not the genes themselves)
# status - done
#
# Args:
# @names - list of CGC names to delete


sub remove_gene_names {
  my ($self, @names) = @_;

  my $payload = '{"data": [';
  $payload .= join(',', map { "{\"cgc-name\": \"$_\"}" } @names);
  $payload .= ']}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->batch('DELETE', 'gene/cgc-name', $payload);

}

#======================================================================
# resurrect_genes(@ids)
# Resurrect dead gene.
# Args:
# @ids - array-ref of gene IDs to delete -can be WBGene, CGC, Sequence identifiers.

# input data:
#{
#  "data": [
#    "AAH1.1"
#  ]
#}

sub resurrect_genes {
  my ($self, $ids) = @_;


  my $payload = '{"data": [';
  $payload .= join(',', map { ($_ =~ /^WBGene/) ? ("{\"id\": \"$_\"}") : (($_ =~ /-/) ? ("{\"cgc-name\": \"$_\"}") : ("{\"sequence-name\": \"$_\"}"))  } @{$ids});
  $payload .= ']}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->batch('POST', 'gene/resurrect', $payload);


}

#======================================================================
# suppress_genes($identifiers, $why)


# Args:
# @ids - list of gene IDs to delete -can be WBGene, CGC, Sequence identifiers.

# input data:
#{
#  "data": [
#    "{"sequence-name": AAH1.1"}, ..."
#  ]
#}

sub suppress_genes {
  my ($self, $ids) = @_;


  my $payload = '{"data": [';
  $payload .= join(',', map { ($_ =~ /^WBGene/) ? ("{\"id\": \"$_\"}") : (($_ =~ /-/) ? ("{\"cgc-name\": \"$_\"}") : ("{\"sequence-name\": \"$_\"}"))  } @{$ids});
  $payload .= ']}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->batch('POST', 'gene/suppress', $payload);


}
#======================================================================
# remove_cgc_name_genes($data)
# batch removes CGC-names from one or more genes
#
# This is a near duplicate of the routine remove_gene_names(), but they were both in the original, so I've left them both.
#
# Args:
# $data - array-ref of CGC-names


sub remove_cgc_name_genes {
  my ($self, $data) = @_;
  
  my $payload = '{"data": [';
  $payload .= join(',', map { "{ \"cgc-name\": \"$_\" }" } @{$data});
  $payload .= ']}';
  if ($self->noise()) {print $payload,"\n"}
  
  return $self->batch('DELETE', 'gene/cgc-name', $payload);
  
}
#======================================================================
# remove_other_name_genes($data)
# batch removes other-names from specified genes
#
# Args:
# $data - array-ref of keys of GeneIDs and values of array-ref of their other-names
#  hash-ref of: (WBGene0000001 => ['othername1', 'othername2'])

sub remove_other_name_genes {
  my ($self, $data, $why) = @_;
  
  my $payload = '{"data": [';
  foreach my $id (keys %{$data}) {
    my $value = $data->{$id};
    $payload .= '{"id": "'.$id.'",';

    $payload .= '"other-names": [';
    $payload .= join(',', map { "\"$_\"" } @{$value});
    $payload .= ']';

    $payload .= '},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  my $batch = $self->batch('DELETE', 'gene/update-other-names', $payload);
  
}
#======================================================================
# add_other_name_genes($data)
# batch adds other-names to specified genes
#
# Args:
# $data - array-ref of keys of GeneIDs and values of array-ref of their other-names
#  hash-ref of: (WBGene0000001 => ['othername1', 'othername2'])

sub add_other_name_genes {
  my ($self, $data, $why) = @_;
  
  my $payload = '{"data": [';
  foreach my $id (keys %{$data}) {
    my $value = $data->{$id};
    $payload .= '{"id": "'.$id.'",';

    $payload .= '"other-names": [';
    $payload .= join(',', map { "\"$_\"" } @{$value});
    $payload .= ']';

    $payload .= '},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  my $batch = $self->batch('PUT', 'gene/update-other-names', $payload);
  
}

#======================================================================
# batch_merge_genes($data)
# batch merges genes
#
# Args:
# $data - an array-ref of hash with keys of from-gene (gene to die) and into-gene (gene to live) and into-biotype
# $why - string, optional reason for the merger, ignored if undefined



sub batch_merge_genes {
  my ($self, $data, $why) = @_;
  
  if (!defined $why) {$why = ''}

  my %biotypes = ( # input     =>  biotype name in Nameserver
		  'cds'        => 'cds', #was biotype/cds
		  'transcript' => 'transcript', #was biotype/transcript
		  'pseudogene' => 'pseudogene', #was biotype/pseudogene
		  'transposon' => 'transposable-element-gene', #was biotype/transposable-element-gene
		  'transposable-element-gene' => 'transposable-element-gene', #was biotype/transposable-element-gene
		 );

  my @values = values %biotypes;

  foreach my $gene_data (@{$data}) {
    my $biotype = $gene_data->{'into-biotype'};
    # if the biotype is not already in the datomic style, convert it
    if (!exists $biotypes{lc $biotype}) {die "Unknown biotype specified in batch_merge_genes(): $biotype\n"}
    $gene_data->{'into-biotype'} = $biotypes{lc $biotype};
  }

  my $payload = '{"data": [';
  foreach my $hash (@{$data}) {
    $payload .= '{';
    $payload .= '"from-gene": "'.$hash->{'from-gene'}.'",';
    $payload .= '"into-gene": "'.$hash->{'into-gene'}.'",';
    $payload .= '"into-biotype": "'.$hash->{'into-biotype'}.'"';
    $payload .= '},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  return $self->batch('POST', 'gene/merge', $payload);
  
}
#======================================================================
# batch_split_genes($data)
# batch splits genes
#
# Args:
# $data - an array-ref of hash with keys of from-gene (gene to die) and into-gene (gene to live) and into-biotype
# $why - string, optional reason for the merger, ignored if undefined
# Returns: 
# $result - hash ref with keys of sequence names and values of WBGeneIDs created
# $batchid

# Args:
#{
#  "data": [
#    {
#      "from-id": "WBGene00000421",
#      "new-biotype": "cds", #was biotype/cds
#      "product-sequence-name": "AAH1.1",
#      "product-biotype": "cds" #was biotype/cds
#      "product-other-names": ["name1","name2"]
#    }
#  ],
#  "prov": {
#    "who": { #was provenance/who
#      "id": "WBPerson33035", #was person/id 
#      "email": "matthew.russell@wormbase.org" #was person/email
#    },
#    "how": "web", #was provenance/how   agent/web
#    "what": "new-gene", #was provenance/what   event/new-gene
#    "when": "2019-05-01T16:03:25.715Z[GMT]", #was provenance/when   
#    "why": "<express reason here>" #was provenance/why
#  }
#}

# Returns:
# {'AAH1.1' => 'WBGene01200123', 'B20H19.2' => 'WBGene01200124'}
# $batchid

sub batch_split_genes {
  my ($self, $data, $why) = @_;
  
  if (!defined $why) {$why = ''}

  my %biotypes = ( # input     =>  biotype name in Nameserver
		  'cds'        => 'cds', #was biotype/cds
		  'transcript' => 'transcript', #was biotype/transcript
		  'pseudogene' => 'pseudogene', #was biotype/pseudogene
		  'transposon' => 'transposable-element-gene', #was biotype/transposable-element-gene
		  'transposable-element-gene' => 'transposable-element-gene', #was biotype/transposable-element-gene
		 );

  my @values = values %biotypes;

  foreach my $gene_data (@{$data}) {
    my $biotype = $gene_data->{'new-biotype'};
    # if the biotype is not already in the datomic style, convert it
    if (!exists $biotypes{lc $biotype}) {die "Unknown biotype specified in batch_merge_genes(): $biotype\n"}
    $gene_data->{'new-biotype'} = $biotypes{lc $biotype};

    $biotype = $gene_data->{'product-biotype'};
    # if the biotype is not already in the datomic style, convert it
    if (!exists $biotypes{lc $biotype}) {die "Unknown biotype specified in batch_merge_genes(): $biotype\n"}
    $gene_data->{'product-biotype'} = $biotypes{lc $biotype};
  }

  my $payload = '{"data": [';
  foreach my $hash (@{$data}) {
    $payload .= '{';
    $payload .= '"from-id": "'.$hash->{'from-id'}.'",';
    $payload .= '"new-biotype": "'.$hash->{'new-biotype'}.'",';
    $payload .= '"product-sequence-name": "'.$hash->{'product-sequence-name'}.'",';
    $payload .= '"product-biotype": "'.$hash->{'product-biotype'}.'"';

    if (exists $hash->{'product-other-names'}) {
      $payload .= ',';
      $payload .= '"product-other-names": [';
      $payload .= join(',', map { "\"$_\"" } @{$hash->{'product-other-names'}});
      $payload .= ']';
    }

    $payload .= '}';
  }
  
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}
  
  my $batch = $self->batch('POST', 'gene/split', $payload);

# now find the Gene IDs that have just been created
  my %result;
  foreach  my $gene_data (@{$data}) {
    my $sequence_name = $gene_data->{"product-sequence-name"};
    my $info = $self->info_gene($sequence_name);
    my $id = $info->{'id'};
    $result{$sequence_name} = $id;
  }

  return (\%result, $batch);
}

#======================================================================
# Entities
#======================================================================
#
# Example entity-types are: 'sequence-feature', 'strain', 'variation' (and any others that have been created)
#
# Entities are a simple way to make new simple Object types
#
# They have three attributes
# generic? indicates if the entity is "generic" or not (if it is true it supports only the fields "status", "id" and "name")
#          Only entity/gene currently has generic? set to false and you can ignore the generic? entity for the foreseeable future.
#          No other type will have generic? set to false unless a new "concrete" type is coded into the Names service.
# enabled? indicates if the endpoints are turned on for that entity type (can be disabled via a API call)
# named?   indicates if the entity requires a name when creating/updating

#======================================================================
# Check if a named entity-type is valid and get its attributes
# Args:
#      entity_type - name of entity to check - lowercase string, e.g. "variation"
# Returned:
#      $valid - true if this is a valid entity type
#      $named - true if this entity type can take names
#      $enabled - true if the entity type endpoints are turned on

sub check_entity {
  my ($self, $entity_type) = @_; # name of entity to check
  
  my $content = $self->curl("GET", "entity");
  if ($self->noise()) {print Dumper $content}
  my $entities = $content->{'entity-types'};
  my $enabled;
  my $named;
  my $valid = 0;
  foreach my $entity (@{$entities}) {
    if ($entity->{'entity-type'} eq $entity_type) {
      $valid = 1;
      $enabled = $entity->{'enabled?'};
      $named = $entity->{'named?'};
      last;
    }
  }
  return ($valid, $named, $enabled);
}

#======================================================================
# SINGLE ENTITY operations
#======================================================================

#======================================================================
# find entities by match to a pattern
# find_entity('variation', $pattern)
# returns a hash ref of any entity whose name matches the input pattern.
# pattern can be any of a regular expression or part of the name of a WBStrainID, or name
# patterns are case-sensitive
# the pattern is contrained to match at the start of the name, use '.*' at the start of the pattern to match anywhere in the name
# example patterns:
# .*int

# Args:
# $entity_type - name of type of entity - lowercase string, e.g. 'variation'
# $pattern - string to search for
sub find_entity {
  my ($self, $entity_type, $pattern) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  my $encoded = CGI::escape($pattern);
  
  my $content = $self->curl("GET", "entity/$entity_type/?pattern=$encoded");
  if ($self->noise()) {print Dumper $content}
  return $content;
 
}
#======================================================================
# create new entity object
# new_entity('variation', $names, $why)

# Assign identifiers and associate names, creating new entity-types


# Args:
# $entity_type - string e.g. 'variation'
# $name - name to create entity-type IDs for.
# $why is optional remark text.

# Returns a sequence of primary entity-type identifiers in the same order as the input names,
# and a unique batch code (for undo).

sub new_entity {
  my ($self, $entity_type, $name, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  if (!defined $why) {$why = ''}

  my $payload = '{"data": {"name":' . $name . '';
  $payload .= '}, "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->curl('POST', "entity/$entity_type", $payload);

}
#======================================================================
# delete an entity
# kill_entity('variation', $name, $why)

# Args:
# $entity_type - entity_type e.g. 'sequence-feature'
# $name - name or ID of entity-type to kill
# $why is optional remark text.

sub kill_entity {
  my ($self, $entity_type, $name, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  if (!defined $why) {$why = ''}

  my $payload = '{prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->curl('DELETE', "entity/$entity_type/$name", $payload);

}
#======================================================================
# summarize a single existing entity
# Args:
# $entity_type - 'variation','strain' etc.
# $name - ID of the entity to be described
#
# Returns hash-ref like:
# {"id":"WBVar00296473","name":"th7","status":"live","history":[{"when":"2020-01-24T10:25:36.101Z","batch/id":"5e2ac665-3d0e-4ca4-bd27-8265f72c20dd","t":"2020-01-24T10:26:45Z","what":"import","who":{"name":"Matthew Russell","email":"matthew.russell@wormbase.org","id":"WBPerson33035"},"how":"importer","changes":[{"attr":"id","value":"WBVar00296473","added":true},{"attr":"name","value":"th7","added":true},{"attr":"status","value":"live","added":true}]}]}
# if the entity is not found, then it returns 'undef'

sub info_entity {
  my ($self, $entity_type, $name) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  my $encoded = CGI::escape($name);

  my $payload = undef; # don't want to pass in a payload, but want a dummy parameter so $not_found is set correctly
  my $not_found = 1; # don'y throw an error if the gene is not found

  my $content = $self->curl("GET", "entity/$entity_type/$encoded", $payload, $not_found);
  
  if ($self->noise()) {print Dumper $content}
  return $content;

}

#======================================================================
# update an existing entity
# update_entity('variation', $ID, $new_name, $why)

# Add a new name to an existing entity

# Args:
# $entity_type - 'variation','strain' etc.
# $id - ID of the entity to be updated
# $new_name - the new name of the entity
# $why is optional remark text.

# Upon successful completion, returns a sequence of primary variation
# identifiers, and a unique operation batch code.
sub update_entity {
  my ($self, $entity_type, $ID, $new_name, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}
  if (!$named) {die "$entity_type is not a named entity_type\n"}

  my $encoded = CGI::escape($ID);

  if (!defined $why) {$why = ''}

  my $payload = '{"data": {"name": "' . $new_name . '"}, "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  my $content = $self->curl('PUT', "entity/$entity_type/$encoded", $payload);

  if ($self->noise()) {print Dumper $content}
  return $content;
}
#======================================================================
# resurrect a killed entity
sub resurrect_entity {
  my ($self, $entity_type, $ID) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  my $encoded = CGI::escape($ID);

  my $content = $self->curl("POST", "entity/$entity_type/$encoded/resurrect");
  
  if ($self->noise()) {print Dumper $content}
  return $content;
}
#======================================================================

#======================================================================
# BATCH ENTITY operations
#======================================================================

# delete entities
# kill_entities('variation', $names, $why)

# Args:
# $names - an array-ref of entity-type IDs to kill
# $why is optional remark text.

sub kill_entities {
  my ($self, $entity_type, $names, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  if (!defined $why) {$why = ''}

  my $payload = '{"data": [';
  $payload .= join(',', map { "\"$_\"" } @{$names});
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne ''); #was provenance/why
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  return $self->batch('DELETE', "entity/$entity_type", $payload);

}
#======================================================================
# new
# create new entities
# new_entities('variation', $names, $why)

# Assign identifiers and associate names, creating new entity-types

# $names is a number of variations to create
#     or for non-variation entity-types, ...
# $names should be sequence of one or more names to create,
# $why is optional.

# Args:
# $entity_type - string of name of type of entities to create
# $names - an array-ref of names to create entity-type IDs for. OR a number of variations to create.
# $why is optional remark text.

# Returns a sequence of primary entity-type identifiers in the same order as the input names,
# and a unique batch code (for undo).

sub new_entities {
  my ($self, $entity_type, $names, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  if (!defined $why) {$why = ''}

  my $payload;
  
  if (!$named) {

    $payload = '{"data": {"n" ';
    $payload .= $names;
    $payload .= '}, "prov": {';
    $payload .= '"why": "'.$why.'"' if ($why ne '');
    $payload .= '}}';

  } else {

    $payload = '{"data": [';
    $payload .= join(',', map { "{ \"name\": \"$_\" }" } @{$names});
    $payload .= '], "prov": {';
    $payload .= '"why": "'.$why.'"' if ($why ne '');
    $payload .= '}}';
  }

  if ($self->noise()) {print $payload,"\n"}

  return $self->batch('POST', "entity/$entity_type", $payload);

}
#======================================================================
# update existing entities
# update_entities('variation', $ID, $new_name, $why)

# Add new names to existing entities

# Args:
# $entity_type - 'variation','strain' etc.
# $data - hash-ref with keys = ID and values = new name
# $why is optional remark text.

# Upon successful completion, returns a sequence of primary variation
# identifiers, and a unique operation batch code.
sub update_entities {
  my ($self, $entity_type, $data, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}
  if (!$named) {die "$entity_type is not a named entity_type\n"}


  if (!defined $why) {$why = ''}

  my $payload = '{"data": [';
  foreach my $id (keys %{$data}) {
    my $name = $data->{$id};
    $payload .= '{"id": "'.$id.'",';
    $payload .= '"name": "'.$name.'"},';
  }
  
  chop $payload; # remove last ','
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  my $content = $self->batch('PUT', "entity/$entity_type", $payload);

  if ($self->noise()) {print Dumper $content}
  return $content;
}
#======================================================================

# resurrect killed entities
# Args:
# entity_type - string e.g. strain
# $IDs - array-ref - list of IDs to resurrect
# $why is optional remark text.
sub resurrect_entities {
  my ($self, $entity_type, $IDs) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}

  my $payload = '{"data": [';
  $payload .= join(',', map { "\"$_\"" } @{$IDs});
  $payload .= ']}';

  if ($self->noise()) {print $payload,"\n"}

  my $content = $self->batch("POST", "entity/$entity_type/resurrect", $payload);
  
  if ($self->noise()) {print Dumper $content}
  return $content;
}
#======================================================================

# remove_names from a batch of entities
# Args:
# entity_type - string e.g. strain
# $names - array-ref - list of naames to remove
# $why is optional remark text.
sub remove_names_from_entities {
  my ($self, $entity_type, $names, $why) = @_;
  my ($valid, $named, $enabled) = $self->check_entity($entity_type);
  if (!$valid) {die "$entity_type is not a valid entity_type\n"}
  if (!$named) {die "$entity_type is not a named entity_type\n"}

  if (!defined $why) {$why = ''}

  my $payload = '{"data": [';
  $payload .= join(',', map { "\"$_\"" } @{$names});
  $payload .= '], "prov": {';
  $payload .= '"why": "'.$why.'"' if ($why ne '');
  $payload .= '}}';
  if ($self->noise()) {print $payload,"\n"}

  my $content = $self->batch("DELETE", "entity/$entity_type/name", $payload);
  
  if ($self->noise()) {print Dumper $content}
  return $content;
}


#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================
# BATCH ASSORTED REST operations
#======================================================================

# recent_batch
# GET /api/recent/batch
# List recent batch activity within a time-range.
# Args: 
#      $from_date - start date of query in UTC ISO date format (2019-03-01T11:34:56Z)
#      $until_date - start date of query in UTC ISO date format (2019-05-01)

sub recent_batch {
  my ($self, $from_date, $until_date) = @_;


  my $encoded_from = CGI::escape($from_date);
  my $encoded_until = CGI::escape($until_date);

  my $type = "recent/batch?";
  if (defined $from_date) {$type .= "from=$from_date"}
  if (defined $until_date) {$type .= "&until=$until_date"}

  my $content = $self->curl("GET", $type);

  if ($self->noise()) {print Dumper $content}
  return $content;
}

#======================================================================
#======================================================================
# batch
# The valid http methods and entry points are: 
#   DELETE (kill)
#              where entities are: 'entity/variation', 'entity/sequence-feature', 'entity/variation',  
#                                  'entity/variation/name', 'entity/sequence-feature/name', 'entity/variation/name',  
#                                  'gene', 'gene/cgc-name',
#   POST (create)
#              where entities are: 'entity/variation', 'entity/sequence-feature', 'entity/variation',
#                                  'entity/variation/resurrect', 'entity/sequence-feature/resurrect', 'entity/variation/resurrect',
#                                  'gene', 'gene/resurrect', 'gene/suppress', 'gene/merge', 'gene/split',
#   PUT (update)
#              where entities are: 'entity/variation', 'entity/sequence-feature', 'entity/variation',
#                                  'gene',
#   GET 
#              where entities are: '{batch-id}'

#   
# Args: 
#      method - 'DELETE', 'POST', 'PUT', 'GET'
#      type - The {entity} value for the batch API entry points to complete the URL
#      payload - JSON description of the data to create

sub batch {
  my ($self, $method, $type, $payload) = @_;

  my $content = $self->curl($method, "batch/$type", $payload);
  return $content;

}

#======================================================================
# Curl - generic access to WormBase REST entry-points using the utility 'curl'
#======================================================================
# Args:
#  $method - string - one of GET, POST, PUT, DELETE
#  $type - string - end of URL of REST entry-point
#  $payload - string - JSON data to provide to the REST entry-point
#  $not_found - boolean - if true then 'HTTP/1.1 404 Not Found' is not an error. If the gene is not found, then message => "Resource not found" is returned.
#
# Returns perl data structure of decoded JSON results

sub curl {

  my ($self, $method, $type, $payload, $not_found)=@_;
  
  my $cmd="curl -X $method -H \"Authorization: Token $self->{'TOKEN'}\" -H \"content-type: application/json\"  -v ";
  $cmd.="--data \'$payload\'" if (defined $payload);
  my $test = '';
  if ($self->test()) {$test = 'test-'}
  my $url = "https://${test}names.wormbase.org/api/$type";
  my $res;
  
  # to capture both STDOUT and STDERR, it is easiest and safest to
  # redirect them separately to files and then read from those files.
  #      -- Perl Cookbook
  
  my $out="/tmp/NSAPI.$$.stdout";
  my $err="/tmp/NSAPI.$$.stderr";
  
  if ($self->noise()) {print "$cmd \"$url\" 1>$out 2>$err\n"}
  if (system("$cmd \"$url\" 1>$out 2>$err")) {
    open(ERR, "<$err") || die "Can't open $err\n";
    undef $/;
    my $stderr = <ERR>;
    close(ERR);
    $/ = "\n";
    print "stderr = $stderr\n";
    system("rm -f $out $err");
    die "Program terminated.\n";

  }

  open(OUT, "<$out") || die "Can't open $out\n";
  $res = <OUT>;
  close(OUT);
  
  open(ERR, "<$err") || die "Can't open $err\n";
  undef $/;
  my $stderr = <ERR>;
  close(ERR);
  $/ = "\n";

  system("rm -f $out $err");

  if (index($stderr, 'HTTP/1.1 401 Unauthorized') != -1) {
    print "\n\n";
    $self->print_authentication_instructions;
    die "\n";
  }

  if (!defined $res || $res eq '') {
    print "\n\ncurl command:\n$cmd $url\n";
    print "\n\ncurl STDERR output:\n\n";
    print "$stderr\n";
    die "\n\nERROR: There is nothing returned in the results.\n";
  }

  my $content =  decode_json($res);

  # if having a missing entry is acceptable and the entry is missing, then don't do the error trapping, just return 'undef' as a flag that the entry was not found
  if (defined $not_found && $not_found) {
    if (index($stderr, '< HTTP/1.1 404 Not Found') != -1) {
      return undef;
    }
  }

  if ((index($stderr, '< HTTP/1.1 2') == -1) || 
      exists $content->{'errors'} || 
      exists $content->{'problems'} || 
      exists $content->{'data'}{'problems'} || 
      (exists $content->{'class'} && $content->{'class'} eq 'java.lang.Exception') || 
      (exists $content->{'type'} && $content->{'type'} eq 'conflict')
     ) {
    print "\n\n";
    print "curl command:\n$cmd $url\n";
    print "\n\ncurl STDERR output:\n\n";
    print "\n$stderr\n";
    print "\n\nAPI result data:\n\n";
    print Dumper $content;
    print "\n\n";
    die "ERRORS FOUND.\n";
  }
  
  return $content;
}

#======================================================================
#======================================================================
# Get the Google OAUTH2 authentication token.
# If it is not available in the expected file, 
# tell the user how to get it and where to store it.

sub get_token {
  my ($self) = @_;

  my $tokenfile = glob("~/.nameserver/token.s");
  if (!-e $tokenfile) {
    $self->print_authentication_instructions;
    die "\n";
  }
  my $mode = (stat($tokenfile))[2] & 07777;
  if ($mode != 0400) {die "Permissions on $tokenfile are not user-readonly\n"}
  if (-e "${tokenfile}~") {carp "You should tidy up and delete the file ${tokenfile}~\n";}
  
  open (TOK, "< $tokenfile") || die "Can't open ~/.nameserver/token.s\n";
  my $TOKEN = <TOK>;
  chomp $TOKEN;
  $TOKEN =~ /\"/g;
  close(TOK);
  return $TOKEN;
}

#======================================================================
# prints instructions on how to set up the required authentication
sub print_authentication_instructions {
  my ($self) = @_;

  print <<"END_MESSAGE"

Each NameServer user needs to get a personal Google OAUTH2
Authentication token.

You must store this token in a file which is readable only by you.

You get the token by running a script with 
'client ID' and 'Client secret' data.

To get the required data client_id and client_secret data:
Go to the site:
https://developers.google.com/console
(you should be logged in as your wormbase.org account)
Click on the tab \"Credentials\" on the left side.
Look for the application listed with the Type \"Other\"
This should be on a line like: 
\"WormBase Names Service (Programmatic Access)   date    Other   hex-string\"

If you do not see this line, then check you are in the right project:
    At the top of the page is a title \"Google APIs\" followed by a tab
    with a drop-down triangle. Click on the triangle and select \"WormBase
    Names Service\".

    A window should pop up titled \"Select from\" followed by a tab.

    Click on this tab and select the correct project: \"WORMBASE.ORG\".

    You should now see a line like: 

    \"WormBase Names Service (Programmatic Access)   date    Other   hex-string\"


Click on the link: \"WormBase Names Service (Programmatic Access)\" to get

Client ID	xxxxxxxxxxxxxxxxxxxxxxx.apps.googleusercontent.com
Client secret	xxxxxxxxxxxxxxxxxxxx


Run the bash script:
https://github.com/WormBase/names/blob/develop/scripts/obtaintoken.sh
which takes the two command-line parameters you just found:
  bash obtaintoken.sh  <Client ID>  <Client secret>

This will prompt you with \"Paste Google one-time-code:\" and will
fire up a sign-in page on your web browser.

Select your wormbase.org account on the sign-in page and you will see
an authentication code to be given to the obtaintoken.sh script.

Copy and paste the code by clicking on the rectangle icon to the right
of the code.

Then in your terminal where you are running the script, paste it in by
clicking on the 'Edit' menu at the top of your terminal window and
selecting 'Paste'.  (Or cut and paste it any other way which works for
you).

Then press the RETURN key to continue running the obtaintoken.sh script.

The obtaintoken.sh script will then show some output JSON data.
In the output JSON data, the ID Token you require is on the line
looking something like:

 \"id_token\": \"eyJhbGciOiJSUzI1NiIsImtpZCI6IjdjMzA5ZTNhMWMxOTk5Y2IwNDA0YWI3MTI1ZWU0MGI3Y2RiY2FmN2QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhY2NvdW50cy5nb29nbGUuY29tIiwiYXpwIjoiNTE0ODMwMTk2NzU3LXBkM2dlbDBmNzRwajMyNDNqb2ExdTYzbHZjZHQyZ25kLmFwcHMuZ29vZ2xldXNlcmNvbnRlbnQuY29tIiwiYXVkIjoiNTE0ODMwMTk2NzU3LXBkM2dlbDBmNzRwajMyNDNqb2ExdTYzbHZjZHQyZ25kLmFwcHMuZ26vZ2xldXNlcmNvbnRlbnQuY29tIiwic3ViIjoiMwE0MzEwMjU4MjAyNTU2MDEyMzM0IiwiaGQiOiJ3b3JtYmFyZS5vcmciLCJlbWFpbCI6Imdhcnkud2lsbGlhbXNAd29ybWJhc2Uub3JnIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsImF0X2hhc2giOiJ1ZFR5UGQ0aGRzMzg0UzlSUGhhekpBIiwiZ29vZ2xlIjp7ImdpYyI6IkFMYXczYlNrRTZmSkt0cjl4clUyMW16bTJ5YU9HeGpyZW5OYmQ2SV9McDc0d2hzNTRBIn0sImlhdCI6MTK0OTk4NzM4NywiZXhwIjoxNTQ5OTkwOTg3fQ.CJaqXF1wkM_DnT2gEEy727bhcSqWoBRFft2k_OTeH9C-UBfCEgDFQgo02UFYBjm-FIPubUEWB9VUXoSBLoBZ-_SNNXBSMMpNJHXZaazhYqQAgP6S8ti7nXnjPkeNdhgZBw-dGEhQCF3jd_AWi6E0LH5IF_v8XnUoEli6_2YqHy-Y3ai1n4nWDyEmM0LsgbKCwFqOzS0dkAp20f6tEaPfkTL1-M0F4M5ZSEOPpGCryYLCcGjJmj-Y5tb5aSanPtiGmxBq9IVKWXH35NgSIChRvYpiNMY_ymWeA6K9j3HBZz_Bx5naVTAOPXRnLwX78uOQIwX9gqcA5FUGM-p_0k3Rhw\"

Copy the data string after \"id_token\":  into a file named: 
~/.nameserver/token.s
You should not copy the quote marks.

Then make this file readable only by you:
chmod 400 ~/.nameserver/token.s

It should be necessary to only set this token up once. 
However, it doesn't hurt to repeat this procedure, as required.
END_MESSAGE
}



#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================
# old network access stuff using LWP.
# works fine on EBI, but the Sanger proxy stuff stops it working
# so it has been replaced by calls to 'curl'
#======================================================================
# batch
# The valid http methods and entry points are: 
#   DELETE (kill)
#              where entities are: 'gene' 'variation', 'sequence-feature', 
#                                  'gene/cgc-name', 'variation', 'variation/name', 
#                                  'sequence-feature'
#   POST (create)
#              where entities are: 'gene', 'gene/resurrect', 'gene/suppress', 
#                                  'gene/merge', 'gene/split', 'variation', 
#                                  'variation/resurrect', 'sequence-feature', 
#                                  'sequence-feature/resurrect'
#   PUT (update)
#              where entities are: 'gene', 'variation
#   GET 
#              where entities are: '{batch-id}'
#   
# Args: 
#      method - 'DELETE', 'POST', 'PUT', 'GET'
#      type - The {entity} value for the batch API entry points to complete the URL
#      payload - JSON description of the data to create

sub old_batch {
  my ($self, $method, $type, $payload) = @_;

  
  my $url = "https://names.wormbase.org/api/batch/$type";
  print $url,"\n";
  
  my $header;

  if ($method eq 'GET') {
    $header = ['Accept' => 'application/json',
	       'Authorization' => "Token $self->{'TOKEN'}"
	      ];
  } else {
    $header = [
	       'Content-Type' => 'application/json; charset=utf-8',
	       'Authorization' => "Token $self->{'TOKEN'}"
	      ];
  }


  my $request = HTTP::Request->new($method, $url, $header, $payload);
  my $ua = LWP::UserAgent->new();
#  $ua->proxy('http','http://wwwcache.sanger.ac.uk:3128');
  $ua->env_proxy;
  my $response = $ua->request($request);
  if ($response->is_error) {$self->old_print_errors($response)}
  
  my $content = decode_json($response->decoded_content);
  
  if ($self->noise()) {print Dumper $content}
  return $content;

}

#======================================================================
# Print out the available error details and die

sub old_print_errors {
  my ($self, $response) = @_;

  print "ERROR: ", $response->status_line, "\n";
  if ($response->decoded_content ne '') {
    my $content = decode_json($response->decoded_content);
    print "========== Data Dumper ==========================\n";
    print Dumper $content;
    print "=================================================\n";
    
    print $content->{message},"\n";
    if (ref($content->{errors}) eq 'ARRAY') {
      foreach my $err_msg (@{$content->{errors}}) {
	print "$err_msg\n";
      }
    }
  }

  die ("Processing stopped.\n");
}


1;

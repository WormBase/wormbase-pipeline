#!/software/bin/perl -w

package NameDB_handler;
use Carp;
use NSAPI;
use warnings;
use strict;
use Data::Dumper;

#author gw3

=head1 

NameDB_handler

This is an API layer above the NSAPI.pm module specifically for handling GeneID manipulation

=item Synopsis

my $db = NameDB_handler->new($wormbase, $species);
	
my $gene_id = $db->merge_genes($id_to_keep, $id_to_kill);
my $split   = $db->split_genes($cds_name, $name_type,$gene_id);

New genes created are given the same species name as the Wormbase->{species} name.

=cut

#use base NSAPI;
our @ISA = qw( NSAPI );


# status - done
sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  $self->{'wormbase'} = shift;    

  $self->{'full_name'} = {
		      elegans   => 'Caenorhabditis elegans',
		      briggsae  => 'Caenorhabditis briggsae',
		      remanei   => 'Caenorhabditis remanei',
		      brenneri  => 'Caenorhabditis brenneri',
		      japonica  => 'Caenorhabditis japonica',
		      brugia    => 'Brugia malayi',
		      ovolvulus => 'Onchocerca volvulus',
		      sratti    => 'Strongyloides ratti',
		      tmuris    => 'Trichuris muris',
		     };

  $self->{'short_species'} = $self->{'wormbase'}->species;
  $self->{'long_species'} = $self->{'wormbase'}->long_name;


  $self->{'db'} = NSAPI->connect();
  
  #read in clone list to validate CDS names with
  # *** %clones does not appear to be used now, but it doesn't hurt to have it ***
  my $path = "$ENV{CVS_DIR}/NAMEDB";
  my $clone_file = "$path/clonelist";

  #untaint file
  unless( $clone_file =~ m/^(.+)$/ ) {
    $self->dienice("tainted file\n");
  }
  $clone_file = $1;
  
  #read clones in
  open(CLONE,"<$clone_file") or $self->dienice("cant open $clone_file\n");
  my %clones;
  while (<CLONE>) {
    chomp;
    $self->{clones}->{uc($_)} = 1;
  }
  close CLONE;
  
  return $self;
}

#######################################################################


=head2 

close

status - done

This closes the NSAPI object.

=item Synopsis

$db->close;

=cut

sub close {
  my ($self) = @_;
  $self->{'db'}->disconnect();
}


#######################################################################

=head2 test

status - in progress

This tests some things in this package.

=item Synopsis

$db->test();	

=cut

sub test {
  my ($self) = @_;

# %data - hash keyed by type of entity to be tested, one of: 'gene', 'variation', 'strain', 'feature', 'person', 'species', 'stats'
# $data{$entity} - hash keyed by priority - low number priority is simple stuff that should work
#                                - high numbered stuff relies on lots of simple stuff
# $data{$entity}{$priority} - hash keyed by name of test
# $data{$entity}{$priority}{$test_name} - hash with keys 'setup' (required initial condition), 'test' (what needs to be run), 
#                                                        'expect' (expected results), 'cleanup' (restore to initial conditions), 'description' (what)

my $RESULT;
my @RESULT;
my %data = (
	    'gene' => {
		       0 => {


#			     zexit => {test => sub {exit(0)}, description => 'exit'}

			    },

		       1 => {
##			     idExists_1 => {test => sub {idExists($self, "WBGene00000001")}, description => 'Gene exists'},
##			     idExists_2 => {test => sub {!idExists($self, "Bogus_gene")}, description => "Bogus gene doesn't exist"},
##			     idExists_3 => {test => sub {!idExists($self, "Bogus gene")}, description => 'Gene ID with spaces'},
##			     idExists_4 => {test => sub {!idExists($self, "Bogus_Name")}, description => "Bogus name doesn't exist"},
##			     idExists_5 => {test => sub {idExists($self, "WBGene00000052")}, description => 'Gene exists even when it is Dead'},
##			     idLive_1 => {test => sub {idLive($self, "WBGene00000001")}, description => 'Gene is Live'},
##			     idLive_2 => {test => sub {!idLive($self, "WBGene00000052")}, description => 'Gene is Dead'},
##			     printAllNames => {test => sub {printAllNames($self, 'WBGene00000001')}, description => 'printAllNames'},
##			     validate_name_1 => {test => sub {validate_name($self, "AC3.3", "Sequence", $self->{long_species})}, description => "validate_name('AC3.3')"},
##			     validate_name_2 => {test => sub {validate_name($self, "abc-20", "Sequence", $self->{long_species})}, description => "validate_name('abc-20', 'CGC')"},
##			     validate_name_3 => {test => sub {!validate_name($self, "AC3.3.3", "Sequence", $self->{long_species})}, description => "don't validate_name('AC3.3.3')"},
##			     validate_name_4 => {test => sub {validate_name($self, "abc-2.1", "CGC", $self->{long_species})}, description => "validate_name('abc-2.1', 'CGC')"},
##			     validate_name_5 => {test => sub {!validate_name($self, "abc-2.A", "CGC", $self->{long_species})}, description => "don't validate_name('abc-2.A', 'CGC')"},
##			     validate_name_6 => {test => sub {validate_name($self, "abc-2", "CGC", $self->{long_species})}, description => "validate_name('abc-2', 'CGC')"},
##			     idGetByTypedName_1 => {test => sub {!idGetByTypedName($self, 'Sequence', 'BOGUS.10')},  description => "Bogus Sequence name doesn't exist"},
##			     idGetByTypedName_2 => {test => sub {idGetByTypedName($self, 'CGC', 'aap-1')},  description => "CGC name exists"},
##			     idGetByTypedName_3 => {test => sub {idGetByTypedName($self, 'Sequence', 'Y110A7A.10')},  description => "Sequence name exists"},
##			     check_pre_exists_1 => {test => sub {!check_pre_exists($self, "aat-1", "CGC")}, description => 'CGC name already exists'},
##			     check_pre_exists_2 => {test => sub {check_pre_exists($self, "xyz-999", "CGC")}, description => "CGC name doesn't already exist"},

##			     addName_1 => {test => sub {addName($self, 'WBGene00000263', 'CGC', 'xxx-1')},  cleanup => sub {delName($self, 'xxx-1')}, description => 'Set the CGC name'},
##			     addName_2 => {test => sub {addName($self, 'WBGene00000263', 'Sequence', 'BOG.1')},  cleanup => sub {addName($self, 'WBGene00000263', 'Sequence', 'F23H11.5')}, description => 'Set the Sequence name'},
##			     addName_3 => {test => sub {addName($self, 'WBGene00000263', 'Biotype', 'Transposon')},  cleanup => sub {addName($self, 'WBGene00000263', 'Biotype', 'CDS')}, description => 'Set the Biotype'},


##			     make_new_obj_1 => {test => sub {$RESULT = make_new_obj($self, 'xxxx-5', 'CGC')},  cleanup => sub {kill_gene($self, $RESULT->{'id'}); delName($self, $RESULT->{'cgc-name'})}, description => 'Make a new uncloned gene'},
##			     make_new_obj_2 => {test => sub {$RESULT = make_new_obj($self, 'AYC3.303', 'Sequence', 'Pseudogene')},  cleanup => sub {kill_gene($self, $RESULT->{'id'});}, description => 'Make new cloned gene'},

##			     delName => {setup => sub {addName($self, 'WBGene00000263', 'CGC', 'xxx-1')}, test => sub {delName($self, 'xxx-1')}, description => 'Delete a CGC name from a Gene'},
##			     idKill => {test => sub {idKill($self, 'WBGene00000263')},  cleanup => sub {idResurrect($self, 'WBGene00000263')}, description => 'Set Gene status to Dead'},
##			     idResurrect => {setup => sub {idKill($self, 'WBGene00000263')}, test => sub {idResurrect($self, 'WBGene00000263')}, description => 'Set Gene status to Live'},

##			     idSplit => {test => sub {idSplit($self, 'WBGene00000263', 'CDS', 'DDDD.1', 'CDS', 'Testing idSplit')},  cleanup => sub {idKill($self, 'DDDD.1')}, description => 'Split a gene'},
##			     idMerge => {test => sub {idMerge($self, 'WBGene00000263', 'WBGene00000264', 'CDS')},  cleanup => sub {idUnmerge($self, 'WBGene00000263', 'WBGene00000264')}, description => 'Merge two genes'},
##			     idUnmerge => {setup => sub {idMerge($self, 'WBGene00000263', 'WBGene00000264', 'CDS')}, test => sub {idUnmerge($self,'WBGene00000263', 'WBGene00000264')},  description => 'Unmerge two genes'},

##			     new_genes_1 => {test => sub {$RESULT = new_genes($self, [{'cgc-name' => 'gary-1', "species" => "Caenorhabditis elegans"}, {'cgc-name' => 'gary-2', "species" => "Caenorhabditis elegans"}])},  cleanup => sub {my $id = $RESULT->[0]{'id'}; kill_gene($self, $id); delName($self, 'gary-1'); delName($self, 'gary-2')}, description => 'Make two new uncloned gene'},
##			     new_genes_2 => {test => sub {$RESULT = new_genes($self, [{'sequence-name' => 'HOGWART.1', 'biotype' => 'CDS', "species" => "Caenorhabditis elegans"}])},  cleanup => sub {my $id = $RESULT->[0]{'id'}; kill_gene($self, $id)}, description => 'Make a new cloned gene'},


##			     idGetHistory => {test => sub {idGetHistory($self, 'console', '2019-01-01', '2019-11-30')}, description => 'History of events'},

##			     idAllNames_1 => {test => sub {idAllNames($self, 'sequence-name', 'WBGene00000263')}, description => 'Get sequence name of WBGene'},
##			     idAllNames_2 => {test => sub {idAllNames($self, 'cgc-name', 'WBGene00000265')}, description => 'Get cgc name of WBGene'},
##			     idAllNames_3 => {test => sub {idAllNames($self, 'sequence-name', 'aat-1')}, description => 'Get sequence name of cgc id'},
##			     idAllNames_4 => {test => sub {idAllNames($self, 'cgc-name', 'aat-1')}, description => 'Get cgc name of cgc id'},
##			     idAllNames_5 => {test => sub {idAllNames($self, 'sequence-name', 'AC3.3')}, description => 'Get sequence name of sequence id'},
##			     idAllNames_6 => {test => sub {idAllNames($self, 'cgc-name', 'AC3.3')}, description => 'Get cgc name of sequence id'},


##			     idTypedNames_1 => {test => sub {idTypedNames($self, 'CGC', 'WBGene00000265')}, description => 'Get cgc name of WBGene'},
##			     idTypedNames_2 => {test => sub {idTypedNames($self, 'Sequence', 'WBGene00000265')}, description => 'Get sequence name of WBGene'},
##			     idTypedNames_3 => {test => sub {idTypedNames($self, 'Biotype', 'WBGene00000265')}, description => 'Get biotype name of WBGene'},

##			     idPublicName_1 => {test => sub {idPublicName($self, 'WBGene00000265')}, description => 'Get public name of WBGene00000265'},
##			     idPublicName_2 => {test => sub {idPublicName($self, 'WBGene00000266')}, description => 'Get public name of WBGene00000266'},
##			     idPublicName_3 => {test => sub {idPublicName($self, 'WBGene00000267')}, description => 'Get public name of WBGene00000267'},

##			     idSearch_1 => {test => sub {idSearch($self, 'CGC', 'brd')},  description => 'Find all gene with CGC "brd*"'},
##			     idSearch_2 => {test => sub {idSearch($self, 'Sequence', 'AC3')},  description => 'Find all gene with Sequence name "AC3*"'},
##			     idSearch_3 => {test => sub {idSearch($self, 'AC')},  description => 'Find all gene with CGC or Sequence "AC*"'},

##			     idGetByAnyName_1 => {test => sub {idGetByAnyName($self, 'brd-1')}, description => 'Get any name'},

			    },
		       2 => {
##			     validate_id => {test => sub {validate_id($self, "WBGene00000001")}, description => "Gene Exists and is Live"},
##			     remove_all_names => {setup => sub {add_name($self, 'WBGene00000263', 'xxx-1', 'CGC')}, test => sub {remove_all_names($self, 'xxx-1')}, description => 'Remove CGC name'},


#			     recent_gene => {setup => sub  {}, test => sub {}, cleanup => sub sub {}, description => ''},
			     find_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     getNameTypes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     new_gene => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
			    },
		       3 => {

##			     add_name_1 => {test => sub {add_name($self, 'WBGene00000263', 'xxx-1', 'CGC')},  cleanup => sub {delName($self, 'xxx-1')}, description => 'Set the CGC name'},
##			     add_name_2 => {test => sub {add_name($self, 'WBGene00000263', 'BOG.1', 'Sequence')},  cleanup => sub {addName($self, 'WBGene00000263', 'Sequence', 'F23H11.5')}, description => 'Set the Sequence name'},
##			     add_name_3 => {test => sub {add_name($self, 'WBGene00000263', 'Transposon', 'Biotype')},  cleanup => sub {addName($self, 'WBGene00000263', 'Biotype', 'CDS')}, description => 'Set the Biotype'},
##			     force_name_1 => {test => sub {force_name($self, 'WBGene00000263', 'invalidCGC', 'CGC')},  cleanup => sub {delName($self, 'invalidCGC')}, description => 'Force set an invalid CGC name'},
##			     force_name_2 => {test => sub {force_name($self, 'WBGene00000263', 'invalidSequence', 'Sequence')},  cleanup => sub {addName($self, 'WBGene00000263', 'Sequence', 'F23H11.5')}, description => 'Force set an invalid Sequence name'},

#			     merge_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},


#			     batch_split_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     batch_merge_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     remove_cgc_name_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     suppress_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     resurrect_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     kill_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     update_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     new_genes => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     remove_name => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     kill_gene => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     split_gene => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#			     print_history => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
			    },
		      },

	    'variation' => {
			    1 => {
##				  find_variations_1 => { test => sub {find_variations($self, 'a83')}, description => 'Find variations by pattern a83'},
##				  find_variations_2 => { test => sub {find_variations($self, 'bogus')}, description => 'Find variations by pattern bogus'},
#				  resurrect_variations => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#				  update_variations => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#				  new_variations => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
				 },
			   },
#	    'feature' => {
#			    1 => {
#				  resurrect_features => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#				  kill_features => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#				  new_features => {setup => sub {}, test => sub {},  cleanup => sub {}, description => ''},
#				 },
#			 },
#	    'debug' => {
#			1 => {
#			      debug => {setup => sub {print "OK"}, description => 'debug'},
#			     },
#		       },
	   );


  my $errors=0;
  foreach my $entity (sort {$a cmp $b} keys %data) {
    foreach my $priority (sort {$a <=> $b} keys %{$data{$entity}}) {
      foreach my $test_name (sort {$a cmp $b} keys %{$data{$entity}{$priority}}) {
	my $setup = $data{$entity}{$priority}{$test_name}{setup} if (exists $data{$entity}{$priority}{$test_name}{setup});
	my $test = $data{$entity}{$priority}{$test_name}{test} if (exists $data{$entity}{$priority}{$test_name}{test});
	my $cleanup = $data{$entity}{$priority}{$test_name}{cleanup} if (exists $data{$entity}{$priority}{$test_name}{cleanup});
	my $description = $data{$entity}{$priority}{$test_name}{description};
	print "$entity\t$priority\t$test_name\t$description\t";
	if (defined $setup) {$setup->()}
	my $result = '';
	if (defined $test) {$result = $test->()} 
	if (!defined $result || !$result) {print "ERROR"; $errors++} elsif ($result eq '1') {print "OK"} else {print Dumper $result}
	if (defined $cleanup) {$cleanup->()}
	print "\n";
      }
    }
  }
  print "\n\n$errors ERRORS\n";
  exit(0);

# the commented out tests have generally already been run and have
# added their data to the database, so will thow an error if run again
# - change the data to be added if you wish to test these

  print "test create variation:             ", $self->new_variations(['nix8', 'nix9']) ? "OK"  : "ERROR", "\n";
  print "test for a gene split              ", $self->split_gene("WBGene00000001", 'AC3.132', 'CDS', 'Caenorhabditis elegans', 'Testing') ?  "OK"  : "ERROR", "\n";
#  print "test if a bogus gene exists:       ", $self->{'db'}->info_gene("Bogus_Name") ? "ERROR"  : "OK", "\n";
  print "test update variation:         ", $self->update_variations({'WBVar01000844' => 'nixon1'}) ? "OK"  : "ERROR", "\n";

  print "test if NSAPI object is working:   ", ($self->{'db'}->ping eq 'Bonjour!') ? "OK" : "ERROR", "\n";


  print "test of dienice: \n"; $self->dienice("This message should be displayed.")

}

#######################################################################
=head2 noise

status - done

This turns on/off debugging output in the curl call.

=item Synopsis


$db->noise(1);

=cut

sub noise {

  my $self = shift;
  my $arg = shift;
  $self->{'db'}->noise($arg);

}

#######################################################################
=head2 print_authentication_instructions

status - done

This displays the instructions to set up the authentication to use the Nameserver

=item Synopsis


$db->print_authentication_instructions();

=cut

sub print_authentication_instructions {
  my $self = shift;
  $self->{'db'}->print_authentication_instructions();

}



#######################################################################


=head2 info_gene

status - done

This calls NSAPI::info_gene
The data structure returned is:

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



=item Synopsis


$db->info_gene($id);	


=cut


sub info_gene {
  my $self = shift;
  my $id = shift;
 
  return $self->{'db'}->info_gene($id);
}

#######################################################################

=head2 printAllNames

status - done

This prints the Sequence-name and/or CGC name (if they exist) for a gene.
This has been simplified from the original because now CGC name and Sequence name are unique for a gene.

=item Synopsis


$db->printAllNames($id);	

=cut

sub printAllNames
  {
    my $self = shift;
    my $id = shift;

    my $info = $self->{'db'}->info_gene($id);
    if (exists $info->{'message'}) {print "Not found\n";}
    print "Current gene names for $id\n";
    print "cgc-name:      ", $info->{'cgc-name'}, "\n" if (exists $info->{'cgc-name'});
    print "sequence-name: ", $info->{'sequence-name'}, "\n" if (exists $info->{'sequence-name'});
    return 1;
  }

#######################################################################

=head2 validate_name

status - done

This checks whether the input sequence_name conforms to the accepted format.

Args:
      name - name to test
      type - one of 'CGC' or 'Sequence' or 'Biotype'
      binomial species name - e.g. 'Caenorhabditis elegans'

=item Synopsis


$db-> ();	

=cut

sub validate_name {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $long_species = shift;

  if (exists $self->{'full_name'}{$long_species}) {$long_species = $self->{'full_name'}{$long_species}}

  #is this a valid name type?
  my @types = $self->getNameTypes;
  if ( grep {$_ eq $type} @types) {    		

    #check name structure matches format eg CDS = clone.no
    my $name_checks = {
		       'Caenorhabditis elegans'  => { 
						     "CGC" => '^[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
						     "Sequence" => '(^([A-Z0-9_cel]+)\.\d+$)|(^([A-Z0-9_cel]+)\.t\d+$)',
						    },
		       'Caenorhabditis briggsae' => {
						     "CGC" => '(^Cbr-[a-z21]{3,4}-[1-9]\d*(\.\d+)?)|(^Cbr-[a-z21]{3,4}\([a-z]+\d+\)?$)',
						     "Sequence" => '^CBG\d{5}$',
						    },
		       
		       'Caenorhabditis remanei ' => {
						     "CGC" => '^Cre-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
						     "Sequence" => '^CRE\d{5}$',
						    },
		       'Caenorhabditis brenneri' => {
						     "CGC" => '^Cbn-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
						     "Sequence" => '^CBN\d{5}$',
						    },
		       'Pristionchus pacificus' => {
						    "CGC" => '^Ppa-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
						    "Sequence" => '^PPA\d{5}$',
						   },
                       'Caenorhabditisjaponica' => {
						    "CGC" => '^Cjp-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',	
						    "Sequence" => '^CJA\d{5}$',
						   },
		       'Brugia malayi' => {
					   'CGC' => '^Bma-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
					   'Sequence' => '^Bm\d+$',
					  },
		       'Onchocerca volvulus' => {
						 'CGC' => '^Ovo-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
						 'Sequence' => 'OVOC\d+$',
						},
		       'Strongyloides ratti' => {
						 'CGC' => '^Sra-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
						 'Sequence' => 'SRAE_[\dXM]\d+$',
						},
		       'Trichuris muris' =>  {
					      'CGC' => '^Tmu-[a-z21]{3,4}-[1-9]\d*(\.\d+)?$',
					      'Sequence' => '^TMUE_[0123M]\d{9}$',
					     }

		      };

    if ( $type ne 'Biotype' && $name !~ /$name_checks->{$long_species}->{$type}/ ) {
      $self->dienice("$name is the incorrect format for $long_species $type match ".$name_checks->{$long_species}->{$type});
      return undef;
    }
#    if (($species eq 'elegans') and ($type eq 'Sequence') and !(defined $self->{'clones'}->{uc($2)}) ) {
#      $self->dienice("$name isnt on a valid clone $2");
#      return undef;
#    }
  } else {
    $self->dienice("$type is a invalid typename: @types");
    return undef;
  }
  return $name;
}

#######################################################################

=head2 validate_id

status - done

This checks whether the input gene ID exists and is Live.

=item Synopsis


$db->validate_id($id);

=cut

sub validate_id
  {
    my $self = shift;
    my $id = shift;
    unless ( $self->idExists($id) ) {
      $self->dienice("$id does not exist");
      return undef ;
    }
    unless( $self->idLive($id) == 1) {
      $self->dienice("$id is not live");
      return undef;
    }
    return $id;
  }

#######################################################################

=head2 check_pre_exists

status - done

This checks for a pre-existing gene with the specified Sequence or CGC name.
Returns undef if the name is already used else returns the input name.

=item Synopsis


$db->check_pre_exists('abc-1', 'CGC');	

=cut

sub check_pre_exists
  {
    my ($self, $name, $type) = @_;
    my $id = $self->idGetByTypedName($type => $name); #must be typed name as CDS and Sequence often have same name
    if (defined $id ) {
      return undef;
    }
    return $name;
  }

#######################################################################

=head2 

status - done

Add a name to the identifier.  $nametype must be one of the valid name
types for the domain.  If the identifier already has a name of this
type, the name will be replaced.

Valid nametypes: CGC, Sequence, Biotype


=item Synopsis


$db->add_name($gene_id, 'abc-20', 'CGC');	

=cut

sub add_name {
  my ($self, $gene_id, $name, $type) = @_;
  unless ($self and $gene_id and $name and $type) {
    $self->dienice("bad parameters");
    return undef;
  }
  
  #confirm GeneID and name are valid
  $self->validate_name($name, $type, $self->{'long_species'})    or return undef;
  $self->validate_id($gene_id)          or return undef; # gene should exists and be Live
  $self->check_pre_exists($name, $type) or return undef;
  if ( $self->addName($gene_id, $type => $name) ) {
    return $name;  
  } else {
    $self->dienice("$name not added to $gene_id");
    return undef;
  }
}

#######################################################################

=head2 

status - done

This is like add_name, but it is for use by gene name curator to add a
non-standard GCG or Sequence name that would fail validation.

Add a name to the identifier which does not have to conform to the standard ID format.  

$nametype must be one of the valid name types for the domain.  If the
identifier already has a name of this type, the name will be replaced.

It is likely that the database will not allow invalid names to be
entered, even if force_name() does not reject it.

Valid nametypes: CGC, Sequence, Biotype

=item Synopsis

$db->force_name($gene_id, 'abc-Balpha.1.1', 'CGC');	

=cut

# adding this for use by gene name curator to add lists of non-standard gene names.
sub force_name {
  my ($self, $gene_id, $name, $type) = @_;
  unless ($self and $gene_id and $name and $type) {
    $self->dienice("bad parameters");
    return undef;
  }
  
  #confirm GeneID is valid
  $self->validate_id($gene_id)          or return undef;
  $self->check_pre_exists($name, $type) or return undef;
  my $force = 1; # turn off validation 
  if ($self->addName($gene_id, $type, $name, $force)) {
    return $name;  
  } else {
    $self->dienice("$name not added to $gene_id");
    return undef;
  }
}




#######################################################################

=head2 

status - done

This creates a gene.
It takes a type for the gene which can be one of:
CGC or CDS to give the gene a GCC name (uncloned gene) or a Sequence name (cloned gene).

The biotype is assumed to be 'CDS', but can optionally be specified - valid biotypes: CDS, Pseudogene, Transcript, Transposon

New genes created are given the same species name as the $self->{long_species} name.

Returns hash-ref with key 'id' and value of created WBGeneID

e.g.
{
         'id' => 'WBGene00305174'
         'cgc-name' => 'xxxx-5'
}

or
{
         'id' => 'WBGene00305175'
         'sequence-name' => 'AYC3.303'
}

=item Synopsis


$id = $db->make_new_obj('abc-1', 'CGC');	
$id = $db->make_new_obj('H07G02.3', 'Sequence', 'CDS');	

=cut

sub make_new_obj
  {
    my $self = shift;
    my $name = shift;
    my $type = shift;
    my $biotype = shift; # optional - default is 'CDS'

    my %hash;
    $hash{'species'} = $self->{'long_species'};
    $hash{'cgc-name'} = $name if ($type eq 'CGC');
    $hash{'sequence-name'} = $name if ($type eq 'CDS' || $type eq 'Sequence'); # the script GENEACE/newgene.pl expects to set the Sequence_name by passing the type 'CDS'
    $hash{'sequence-name'} = $name if ($type eq 'Sequence'); # but I prefer to use the type Sequence to set the Sequence_name 
    if ($type ne 'CGC') { # uncloned gene are not permitted to have a biotype
      $hash{'biotype'} = 'CDS'; # default biotype
      $hash{'biotype'} = $biotype if ($biotype); # the optional biotype has been specified
    }
    my @array = (\%hash);
    
    my $info = $self->{'db'}->new_genes(\@array);

    return $info->{'ids'}[0];
  }

#######################################################################

=head2 dienice

Just does a croak().

status - done

This 

=item Synopsis


$db->dienice("Message.");	

=cut
# status - TBD
sub dienice {
  my $self = shift;
  my($errmsg) = @_;
  croak "$errmsg\n";
  return $errmsg;
}

#######################################################################

=head2 

status - TBD

This prints out the history of a gene - it is only used by idstats.pl
It has not been implemented.

=item Synopsis


$db-> ();	

=cut
# status - TBD
sub print_history {
  my $self = shift;
  my $id = shift;

  $self->validate_id($id);

  my $history = $self->idGetHistory($id);
  print "$id\n";
  foreach my $event (@{$history}) {
    print $event->{version}," ";
    print $event->{date}," ";
    print $event->{event}," ";
    print $event->{name_type}," " if (defined $event->{name_type});
    print $event->{name_value}," " if (defined $event->{name_value});
    print $event->{user},"\n";
  }
  return;
}

#######################################################################

=head2 

status - done

This removes all CGC names from the specified gene.

=item Synopsis
$db->remove_all_names($id);	


=cut

sub remove_all_names   {
  my $self = shift;
  my $id = shift;

  # a gene only has one CGC name
  my $name = $self->idAllNames('cgc-name', $id);
  # pass the CGC name as a list, as delName takes an array of names to batch delete.
  if (defined $name) {$self->delName(($name))};

  return $name;
}

#######################################################################

=head2 

status - done

This 

=item Synopsis


$db->merge_genes($live_gene, $dead_gene);	

=cut

sub merge_genes {
  my ($self, $gene_gene, $merge_gene) = @_;
  
  my ($gene_id, $merge_id);
  $gene_id  = ($self->idGetByAnyName($gene_gene)->[0]  or $gene_gene); #allow use of any name or gene_id
  $merge_id = ($self->idGetByAnyName($merge_gene)->[0] or $merge_gene); #allow use of any name or gene_id
  $self->validate_id($gene_id) or return undef;
  $self->validate_id($merge_id) or return undef;
  
  if ( $gene_id eq $merge_id) {
    $self->dienice("FAILED : $gene_gene and $merge_gene are the same!");
    return undef;
  }
  
  #enforce retention of CGC named gene
  my $info1 = $self->{'db'}->info_gene($gene_id); # info_gene returns {'message' => 'Resource not found'}
  if (exists $info1->{'message'}) {$info1 = undef}
  my $info2 = $self->{'db'}->info_gene($merge_id); # info_gene returns {'message' => 'Resource not found'}
  if (exists $info2->{'message'}) {$info2 = undef}
  my $name1 = $info1->{'cgc-name'};
  my $name2 = $info2->{'cgc-name'};
  
  unless( ($ENV{'USER'} eq 'pad') or ($ENV{'USER'} eq 'gw3' )){
    if ( defined $name1 ) {
      if (defined $name2 ) {
	#both genes have CGC name - confirm with geneace 
	$self->dienice("FAILED: Both genes have CGC names ".$name1." and ".$name2.".  The correct course of action should be determined by the CGC admin (pad)\nPlease contact Geneace curator to resolve this\n");
      } else {
	#gene being eaten has a CGC name and eater doesnt
	$self->dienice("FAILED: $merge_gene has a CGC name ".$name2." and should probably be retained");
      }
      return undef;
    }
  }
  # warn that a gene with a CGC name has been killed
  if (defined $name2 ) { 
    print "$merge_id had a CGC name : $name2\n";
  }
  #always remove names from merged id
  $self->remove_all_names($merge_id);
  #if this is a merger between a CGC gene and a sequence we need to transfer the Seq & CDS names too
  # the CGC gene id.
  my $sequence1 = $info1->{'sequence-name'};
  my $sequence2 = $info2->{'sequence-name'};
  unless ( defined $sequence1 ) {
    $self->addName($gene_id, 'Sequence', $sequence2) if ($sequence2);
  }
  #do the merge
  my $biotype = $info2->{'biotype'};
  if ($self->idMerge($merge_id, $gene_id, $biotype)) {
    return ([$gene_id, $merge_id]);
  } else {
    $self->dienice("FAILED : merge failed");
    return undef;
  }
}

#######################################################################

=head2 

status - done

This takes a gene to split and the sequence name and biotype of a gene
to create and splits the gene.

=item Synopsis


$db->split_gene($gene_id, $sequence_new, $biotype_new, $long_species, $why);	

=cut
# status - done
sub split_gene {
  my ($self, $gene_id, $sequence_new, $biotype_new, $long_species, $why) = @_;
#$gene_id, $type, $gene, $long_species) = @_;

  if (defined $long_species && exists $self->{'full_name'}{$long_species}) {$long_species = $self->{'full_name'}{$long_species}}


  unless ($sequence_new && $gene_id){
    $self->dienice("bad parameters");
    return undef;
  }

  my $info = $self->{'db'}->info_gene($gene_id);
  if (exists $info->{'message'}) {$info = undef}
  if (!defined $long_species) {
    $long_species = $info->{'species'};
  }

  # we assume the the original gene that is being split will not have its biotype changed.
  # If this is not the case, then the $biotype_orig would have to be explicitly specified.
  # For now, just use the existing biotype.
  my $biotype_orig =  $info->{'biotype'};


  $self->validate_id($gene_id) or return undef;
  $self->validate_name($sequence_new, 'Sequence', $long_species) or return undef;
  if (!$self->check_pre_exists($sequence_new, 'Sequence')) {
    $self->dienice("FAILED: split gene $gene_id failed because $sequence_new already exists."); #error msg
    return undef;
  }
  my $id = $self->idSplit($gene_id, $biotype_orig, $sequence_new, $biotype_new, $why);

  if ( $id =~ /WBGene/ ) {
    return "$id";
  } else {
    $self->dienice("FAILED: split gene $gene_id failed"); #error msg
    return undef;
  }
}

#######################################################################

=head2 

status - done

This takes a gene ID and a reason why the gene is being killed.
It kills the gene.

=item Synopsis


$db->kill_gene($gene_id, $why);	

=cut

sub kill_gene {
  my ($self, $gene_id, $why) = @_;
  $self->validate_id($gene_id) or return undef;
  if (!defined $why) {$why = "No evidence for this gene."}
  if ($self->idKill($gene_id, $why)) {
    return $gene_id;
  } else {	
    $self->dienice("cant kill GeneID $gene_id");
    return undef;
  }
}

#######################################################################

=head2 

status - done

This 

=item Synopsis


$db-> ();	

=cut
# status - done
sub remove_name {
  my ($self, $gene_id, $type, $name, $long_species) = @_;
  unless ($self and $gene_id and $type and $name) {
    $self->dienice("bad parameters");
    return undef;
  }
		
  $self->validate_id($gene_id) or return undef;

  $self->validate_name($name, $type, $self->{long_species}) or return undef;
  my $exist_id = $self->idGetByTypedName($type,$name);
  if ( !$exist_id or ("$exist_id" ne "$gene_id") ) {
    $self->dienice("$name is not a name of $gene_id");
    return undef;
  }
  print "\nGeneID = $gene_id\nType = $type\nName = $name\nSpecies = $long_species\n";
  if ( $self->delName($name) ) {
    return $name;
  }
}

#######################################################################

=head2 

status - done

This 

=item Synopsis
$db->new_gene('H07G02.3', 'Sequence', $species);	


The second parameter is always 'Sequence' so that we create a mapped gene.

=cut
# status - done
sub new_gene {
  my ($self, $name, $type, $long_species) = @_;
	
  $self->validate_name($name, $type, $long_species)  or return undef;
  $self->check_pre_exists($name, $type) or return undef;

  my $id = $self->make_new_obj($name, $type, 'CDS');
  return $id ? $id : undef;
}
 

######################################################################
# routines from NameDB.pm
######################################################################

=head2 getNameTypes

=over 4

=item @types = $db->getNameTypes() 

status - done

Return a list of the valid types of name.

Example:

 @types = $db->getNameTypes()
  

=back

=cut

sub getNameTypes {
  my $self   = shift;
  return ('CGC', 'Sequence', 'Biotype');
}


#-------------------- public ID creation-----------------------

=head2 Identifier Creation

=over 4

=item $new_id = $db->idCreate([$domain]) # status - retired

Create a new unique identifier in the indicated domain and returns the
ID.  The nametype and name are assigned to the identifier.

Example:

  $new_id = $db->idCreate;

=back

=cut

sub idCreate {
  my $self   = shift;
  my $domain = shift;

  die "idCreate is no longer used\nYou should call something like new_gene() instead";
}


#-------------------- public_id retrieval-----------------------

#=head2 Identifier Lookup

=cut

#######################################################################

=item idExists

 # status - done

test if a gene ID exists

Examples:

$test = $db->idExists("WBGene00000001")

=cut

sub idExists {
  my $self = shift;
  my $id = shift;

  if ($id =~ /\s+/) {return 0} # spaces in the ID are invalid

  my $info = $self->{'db'}->info_gene($id); # info_gene returns hashref {message' => 'Resource not found'} if the ID does not exist
  if (exists $info->{'message'}) {$info = undef}

  return (defined $info) ? 1 : 0; 
}


#######################################################################

=item $ids = $db->idGetByTypedName($nametype, $name) 

# status - done

Look up identifier(s) using an exact match to a typed name.

Name should be one of 'CGC' or 'Sequence'

Examples:

  $ids = $db->idGetByTypedName('CGC', 'aap-1');
  $ids = $db->idGetByTypedName('Sequence', 'Y110A7A.10');

=cut

sub idGetByTypedName {
  my $self = shift;
  my ($type, $name, $domain) = @_;

  my $info = $self->{'db'}->info_gene($name);
  if (exists $info->{'message'}) {$info = undef}  # info_gene returns {message => 'Unable to find any entity for given identifier.'} if the ID does not exist
  if (defined $info && $type eq 'CGC') {if (exists $info->{'cgc-name'}) {return  $info->{'id'}}}
  if (defined $info && $type eq 'Sequence') {if (exists $info->{'sequence-name'}) {return  $info->{'id'}}}
  if ($type eq 'Public_name' && defined $info) { # match either CGC-name or Sequence-name
    if (exists $info->{'sequence-name'}) {return  $info->{'id'}}
    if (exists $info->{'cgc-name'}) {return  $info->{'id'}}
  }

  return undef;
}


=item $ids = $db->idGetByAnyName($name) # status - done

An exact match for a name.
As for the previous method, except that the search will return a match
on any of the nametypes defined for the domain:

  @ids = $db->idGetByAnyName('A12H8.1');
  @ids = $db->idGetByAnyName('brd-1');

=cut

sub idGetByAnyName {
  my $self = shift;
  my ($name) = @_;

  my @array;
  my $info = $self->{'db'}->info_gene($name);
  if (exists $info->{'message'}) {$info->{'id'} = undef}  # info_gene returns {message => 'Unable to find any entity for given identifier.'} if the ID does not exist
  push @array,  $info->{'id'};
  return wantarray ? @array : \@array;

}

=item @ids = $db->idSearch($type => $wildcardpattern) # status - done

=item @ids = $db->idSearch($wildcardpattern)

Perform a wildcard search, optionally constrained by name type, where $type is 'GCG' or 'Sequence'.

=cut

sub idSearch {
  my $self = shift;
  my ($type,$pattern);
  if (@_ == 1) {
    $pattern = shift;
  } else {
    ($type, $pattern) = @_;
  }

  my @array;
  $pattern = '.*' unless defined $pattern;

  my $info = $self->{'db'}->find_genes($pattern);

  foreach my $i (@{$info->{'matches'}}) {
    if (defined $type) {
      if ($type eq 'CGC' || $type eq 'cgc-name') {
	push @array,  $i->{'cgc-name'};
	
      } elsif ($type eq 'Sequence' || $type eq 'sequence-name') {
	push @array,  $i->{'sequence-name'};

      }
    } else {
	push @array,  $i->{'cgc-name'} if (exists $i->{'cgc-name'});
	push @array,  $i->{'sequence-name'} if (exists $i->{'sequence-name'});
    }

  }

  return wantarray ? @array : \@array;

}


=head2 Methods on Individual Identifiers

=over 4

=item $db->idPublicName($id) # status - make as simple call returning CGC or Sequence name

Given an identifier's ID, return its "best" name.  The best name is
the 'CGC' name or if that doesn't exist, the 'Sequence' name.

=cut

sub idPublicName {
  my $self = shift;
  my ($id) = @_;

  my $info = $self->{'db'}->info_gene($id);
  if (exists $info->{'message'}) {$info = undef}  # info_gene returns {message => 'Resource not found'} if the ID does not exist

  my $name = $info->{'cgc-name'} if (exists $info->{'cgc-name'});
  if (! defined $name && exists $info->{'sequence-name'}) {
    $name = $info->{'sequence-name'};
  }

  return $name;
}

#######################################################################

=item $flag = $db->idLive($id) 

# status - done

Return true if the indicated ID is alive, false otherwise.

Examples:

$test = $db->idLive("WBGene00000001")

=cut

sub idLive {
  my $self = shift;
  my ($id) = @_;


  my $info = $self->{'db'}->info_gene($id);
  if (exists $info->{'message'}) {$info = undef}  # info_gene returns {message => 'Resource not found'} if the ID does not exist

  return ($info->{'status'} eq 'live') ? 1 : 0; 
}


=item @names = $db->idTypedNames($nametype, $id) # status - done

Returns the names of the indicated type that are assigned to this
identifier.  In a scalar context will return an array reference
containing zero or more assigned names.  In a list context will return
a list of zero or more names.

For unique names, the easiest way to work with this is to call in a
list context:

  ($genbank_name) = $db->idTypedNames('CGC', $id);
Returns an array of the CGC name (or Sequence, or Biotype) name assigned to this gene.
e.g. (K04C2.4'), ('brd-1') or ('cds')

Args:
     $type - field of gene information to return, e.g. 'CGC', 'cgc-name',  'Sequence', 'sequence-name', 'Biotype' or 'biotype' (See NSAPI::info_gene()))
     $gene_id - name of gene . Can be WBGeneID, Sequence name or CGC name


=cut

# return array of names of given type
sub idTypedNames {
  my $self = shift;
  my ($type, $public_id) = @_;

  my @array;
  my $info = $self->{'db'}->info_gene($public_id);
  if (exists $info->{'message'}) {$info = undef}  # info_gene returns {message => 'Resource not found'} if the ID does not exist
  if (defined $info) {

    if (($type eq 'CGC' || $type eq 'cgc-name') && exists $info->{'cgc-name'}) {push @array,  $info->{'cgc-name'}}
    if (($type eq 'Sequence' || $type eq 'cgc-sequence-name') && exists $info->{'sequence-name'})  {push @array,  $info->{'sequence-name'}}
    if (($type eq 'Biotype' || $type eq 'biotype') && exists $info->{'biotype'})  {push @array,  $info->{'biotype'}}
  }

  return wantarray ? @array : \@array;

}

=item @names = $db->idAllNames($type, $id) 
# status - done

Returns the CGC name or Sequence name assigned to this gene.

Args:
     $type - field of gene information to return, e.g. 'cgc-name',  'sequence-name' (See NSAPI::info_gene()))
     $gene_id - name of gene . Can be WBGeneID, Sequence name or CGC name

=cut

sub idAllNames {
  my $self = shift;
  my ($type, $public_id) = @_;

  my $info = $self->{'db'}->info_gene($public_id);
  if (exists $info->{'message'}) {$info = undef}  # info_gene returns {message => 'Resource not found'} if the ID does not exist

  return $info->{$type};

}

=item @history = $db->idGetHistory($id,[$after, [,$domain]]) 

# status - used in idstat.pl
This returns  the history of all gene between dates - it is only used by idstats.pl

The dates are in the ISO date format
YYYY-MM-DD HH:MM:SS
or
date --utc +%Y-%m-%dT%H:%M:%SZ`; # '2019-06-07T15:04:15Z'

The agent is one of "web" or "console" (the "Agent" type that made the request).

=cut

sub idGetHistory {
  my $self = shift;
  my ($agent, $from, $until) = @_;

  my $info = $self->{'db'}->recent_gene($from, $until, $agent);

  return $info;
}

#######################################################################

=item $db->addName($id, $nametype, $name) 

# status - done

Add a name to the gene identifier.  

$nametype must be one of the valid name types.  If the
identifier already has a name of this type, the name will be replaced
if it is unique.

Valid nametypes: CGC, Sequence, Biotype

 Args: 
       $public_id - string, name of the gene
       $nametypes - string,  one of 'CGC', 'Sequence', 'Biotype'
       $name - string, name of CGC, Sequence, or Biotype to be set e.g. 'AC3.3' or 'abc-1' or 'CDS'
       $force      - optional integer -  set this to true to force an invalid CGC name or Sequence name that is not in the normal format and normal validation will be turned off

 Can't delete data here by using a null or blank string - will have to use Delete in /api/batch/gene/cgc-name for example

=cut


#-------------------- add a typed name -----------------------
sub addName {
  my $self = shift;
  my ($public_id, $nametype, $name, $force) = @_;
  
  
  my %hash;
  $hash{'id'} = $public_id;
  $hash{'species'} = $self->{'long_species'};
  $hash{'cgc-name'} = $name if ($nametype eq 'CGC');
  $hash{'biotype'} = $name if ($nametype eq 'Biotype');
  $hash{'sequence-name'} = $name if ($nametype eq 'Sequence');
  $hash{'force'} = 1 if (defined $force && $force);
  my @array = (\%hash);
  
  my $info = $self->{'db'}->update_genes(\@array);


}

#######################################################################

=item $db->delName($name) # status - done

Removes the indicated CGC name.

It assumes you know that the name exists.

Fails with an error if the name does not exist.

=cut

#-------------------- delete a CGC name -----------------------
sub delName {
  my $self = shift;
  my ($name) = @_;

  

  my $info = $self->{'db'}->remove_gene_names(($name));


}

=head2 idKill

=item $db->idKill($id, $why) 

# status - done

Kill the gene identifier, making its "live" flag false.

=cut

sub idKill {
  my $self        = shift;
  my ($id, $why) = @_;

  my $info = $self->{'db'}->kill_genes([$id], $why);

}

=head2 idResurrect

=item $db->idResurrect($id) # status - done

Resurrect the identifier.

=cut

sub idResurrect {
  my $self        = shift;
  my ($public_name) = @_;

  my $info = $self->{'db'}->resurrect_genes([$public_name]);
}

=head2 idMerge

=item $db->idMerge($from_id, $into_id, $biotype, $why) 

# status - done

Merge the identifier referred to by $from_id into the identifier
referred to by $into_id.  $from_id is killed, becoming
inactive.  All names that are not globally unique are copied from
$from_id to $into_id.

=cut

# merge $from gene_id into $into gene_id
# $from identifier is killed, and $into inherits all of $from's names
#
# Args:
#       $from - string, gene identifier of gene to die
#       $into - string, gene identifier of gene to live
#       $biotype - string, biotype of resulting merged gene, defaults to 'CDS' if undefined, one of 	'cds'	'transcript'	'pseudogene'	'transposon'
#       $why - string, optional reason for the merger, ignored if undefined

sub idMerge {
  my $self = shift;
  my ($from, $into, $biotype, $why) = @_;


  my $info = $self->{'db'}->merge_genes( $from, $into, $biotype, $why);

}

=head2 idSplit

=item $new_id = $db->idSplit( $source_gene_id, $biotype_orig, $sequence_new, $biotype_new, $why)
# status - done

Split the identifier given by $source_id into two, returning the new
identifier.  The new identifier is given the biotype and sequence_name
indicated, but does not inherit any other names automatically.

The new gene_id is returned.

=cut

# split out a new object from current one
# and assign indicated sequence_name
sub idSplit {
  my $self = shift;
  my ( $source_gene_id, $biotype_orig, $sequence_new, $biotype_new, $why) = @_;

  my $info = $self->{'db'}->split_genes( $source_gene_id, $biotype_orig, $sequence_new, $biotype_new, $why);
  # info contains:
  # {'created' => {id => 'WBGene00305173'}, 'updated' => {'id' => 'WBGene00000263'}}
  return $info->{'created'}{'id'};

}

=Head2 idUnmerge

=item $new_id = $db->idUnmerge($from_id, $into_id, $why) 

# status - done

UnMerge the identifier referred to by $from_id from the identifier
referred to by $into_id.  $from_id is resurrected, becoming
active.

# $from identifier is resurrected, and split from $into
#
# Args:
#       $from - string, gene identifier of dead gene to resurrect
#       $into - string, gene identifier of live gene to split from
#       $why - string, optional reason for the unmerger, ignored if undefined

=cut

sub idUnmerge {
  my $self = shift;
  my ($from, $into, $why) = @_;


  my $info = $self->{'db'}->unmerge_genes( $from, $into, $why);
  return $info;


}

#############################################################################

#======================================================================
# GENES
#======================================================================

#############################################################################
=head2 new_genes

This creates a batch of genes.

It takes a type for the gene which can be one of:
CGC or CDS to give the gene a GCC name (uncloned gene) or a Sequence name (cloned gene).

The biotype is one of: CDS, Pseudogene, Transcript, Transposon

New genes created are given the same species name as the $self->{long_species} name.

Args:
     
 $data - array-ref of hashes
         the hashes contain either the keys ('species', 'gcg-name', 'biotype') for an uncloned gene 
                            or ('species', 'sequence-name', 'biotype', and optionally 'cgc-name') for a cloned gene.
 $why - string, optional reason for creating the genes


 e.g.
 $data = [{"cgc-name" => "abc-31", "species" => "Caenorhabditis elegans", "biotype" => 'CDS'}, { ... }];
 or $data = [{"sequence-name" => "ZK1320.13", "species" => "Caenorhabditis elegans", "biotype" => 'CDS'}, { ... }];

Returns:
 $data - array-ref of hashes with keys 'id' and 'cgc-name' / 'sequence-name' (the 'id' keys have the new WBGeneID values)
 $batch - the batch ID for this database insertion.

=item Synopsis


$db->new_genes([{"cgc-name" => 'abc-1', "species" => "Caenorhabditis elegans"}, {cgc-name => 'abc-2', "species" => "Caenorhabditis elegans"}]);	
$db->new_genes([{"sequence-name" => 'H07G02.3', 'biotype' => 'CDS', "species" => "Caenorhabditis elegans"}]);	

=cut

sub new_genes {
  
  my ($self, $data, $why) = @_;
  
  my $info = $self->{'db'}->new_genes($data, $why);
  # $info contains
  # {"id":"5dcd84f2-f17a-44c6-8a2c-e096983181f3","ids":[{"id":"WBGene00305168","cgc-name":"abc-1231"} , ... ]}

  my @ids = @{$info->{'ids'}};
  return (\@ids, $info->{'id'}); # array-ref of hashes with keys 'id' and 'cgc-name' / 'sequence-name', and batch-id
  

}
#############################################################################
=head2 update_genes

This updates a batch of genes.
It is a batch version of addName().

The gene is specified by its WBGeneID.
The type of things to be updated are given by the type: one of  ('species', 'gcg-name', 'biotype', 'sequence-name')


Args:
     
 $data - array-ref of hashes
         the hashes contain the keys ('species', 'gcg-name', 'biotype', 'sequence-name')
 $why - string, optional reason for creating the genes
 $force - optional, if set then the normal checking for a correct CGC name format is turned off, any format is allowed.

 e.g.
 $data = [{"id" => 'WBGene00304791', "cgc-name" => "abc-31", "biotype" => 'CDS'}, { ... }];
 or $data = [{"id" => 'WBGene00304792', "sequence-name" => "ZK1320.13", "biotype" => 'CDS'}, { ... }];

Returns:
 $batch - the batch ID for this database insertion.

=cut

sub update_genes {
  
  my ($self, $data, $why, $force) = @_;

  # convert short species names to binomial names
  foreach my $hash (@{$data}) {
    if (exists $hash->{'species'}) {
      if (exists $self->{'full_name'}{$hash->{'species'}}) {
	$hash->{'species'} = $self->{'full_name'}{$hash->{'species'}};
      }
    }
    if (defined $force) {
      $hash->{'force'} = $force;
    }
  }

  my $info = $self->{'db'}->update_genes($data, $why);

  return ($info->{'updated'}{'id'});

}
#############################################################################
=head2 kill_genes

This deletes a batch of genes.

The gene is specified by its WBGeneID.


Args:
     
 $data - array-ref of gene IDs
 $why - string, optional reason for creating the genes


Returns:
 $batch - the batch ID for this database insertion.

=cut

sub kill_genes {
  
  my ($self, $data, $why) = @_;

  my $info = $self->{'db'}->kill_genes($data, $why);

  return ($info->{'delete'}{'id'});

}

#############################################################################
# Resurrect a batch of dead genes.
#
# Args:
# $names - an array-ref of gene IDs to resurrect


sub resurrect_genes {
  my ($self, $names) = @_;

  my $info = $self->{'db'}->resurrect_genes($names);
  return $info;

}
#############################################################################
# Suppress a batch of genes.
#
# Args:
# $names - an array-ref of gene IDs to suppress


sub suppress_genes {
  my ($self, $names) = @_;

  my $info = $self->{'db'}->suppress_genes($names);
  return $info;

}
#############################################################################
# Remove a batch of CGC-names from their genes
#
# Args:
# $names - an array-ref of CGC-names to remove from their genes


sub remove_cgc_name_genes {
  my ($self, $names) = @_;

  my $info = $self->{'db'}->remove_cgc_name_genes($names);
  return $info;

}
#############################################################################
# Merge a batch of pairs of genes
#
# Args:
# $names - an array-ref of hash with keys of from-gene (gene to die) and into-gene (gene to live) and into-biotype
# $why - string, optional reason for creating the genes


sub batch_merge_genes {
  my ($self, $names, $why) = @_;

  my $info = $self->{'db'}->batch_merge_genes($names, $why);
  return $info;

}
#############################################################################
# Split a batch of genes
#
# Args:
# $names - an array-ref of hash with keys of "from-id" existing gene, "new-biotype" existing gene's new biotype, "product-sequence-name" Sequence-name of gene to create, "product-biotype" biotype of gene to create
# $why - string, optional reason for creating the genes

# Returns:
# 
sub batch_split_genes {
  
  my ($self, $data, $why) = @_;
 
  return $self->{'db'}->batch_split_genes($data, $why);

}

#############################################################################

=head2 find_genes

status - done

This gets the Sequence-name and/or CGC name (if they exist) for a gene.

=item Synopsis


$db->find_genes($pattern);	

 returns a hash ref of any genes whose name matches the input pattern.
 pattern can be any of a regular expression or part of the name of a WBGeneID, CGC name or a Sequence_name
 patterns are case-sensitive
 the pattern is contrained to match at the start of the name, use '.*' at the start of the pattern to match anywhere in the name
 example patterns:
 WBGene00000001
 unc-111
 unc
 AC3.3
 AC3
 WBGene0000000[1-3]
 .*-aat
 Bma-

 Returns array-ref of hashes
      {
        'id' => 'WBGene00000003',
        'sequence-name' => 'F07C3.7',
        'cgc-name' => 'aat-2'
      },
      {
        'id' => 'WBGene00000002',
        'sequence-name' => 'F27C8.1',
        'cgc-name' => 'aat-1'
      },


=cut



sub find_genes {
  my $self = shift;
  my $pattern = shift;
  
  my $info = $self->{'db'}->find_genes($pattern);
  return $info->{'matches'};
  
}

#############################################################################

#======================================================================
# VARIATIONS
#======================================================================

#############################################################################
# new_variations($names, $why)
# creates new Variation IDs
#
# Args:
# $names - array ref of names to give the new variations. Names have a format of letters then digits (e.g. 'pex2')
# $why - string provenence reason for creating the variations (optional)
#
# Returns a hash ref keyed by the new WBVar IDs with the names as values and the batch ID

sub new_variations {
  my ($self, $names, $why) = @_;
  my %new_ids;

  my $info = $self->{'db'}->new_varations($names, $why);

  my @ids = @{$info->{'ids'}};
  
  for my $i (0 .. $#ids) {
    my $id  = $ids[$i]{'id'};
    my $name = $ids[$i]{'name'};
    $new_ids{$id} = $name;
  }
  
  return (\%new_ids, $info->{'id'});

}
#############################################################################
# update_variations($new_names, $why)
# updates Variation IDs with new names
#
# Args:
# $names - hash ref of (key) WBVar IDs and (values) names to give the new variations
# $why - string provenence reason for creating the variations (optional)
#
# Returns an array ref of $info->{'updated'}->{'id'}

sub update_variations {
  my ($self, $new_names, $why) = @_;

  my $info = $self->{'db'}->update_variations($new_names, $why);
  return $info;

}

#############################################################################
# kill a set of variation names
#
# Args:
# $names - an array-ref of variation IDs to kill
# $why is optional remark text.


sub kill_variations {
  my ($self, $names, $why) = @_;

  my $info = $self->{'db'}->kill_variations($names, $why);
  return $info;

}

#############################################################################
#Resurrect a batch of dead variations.
#
# Args:
# $names - an array-ref of variation IDs to resurrect


sub resurrect_variations {
  my ($self, $names) = @_;

  my $info = $self->{'db'}->resurrect_variations($names);
  return $info;

}
#############################################################################
# find_variations($pattern)
# returns an array-ref of hash-ref of any variations whose name matches the input pattern.
# pattern can be any of a regular expression or part of the name of a WBVarID, or name
# patterns are case-sensitive
# the pattern is contrained to match at the start of the name, use '.*' at the start of the pattern to match anywhere in the name
# example patterns:
# .*int

# Args:
# $pattern - string
#
# Returns an array-ref of hash of 
#      'id' => 'WBVar01000645'
#      'name' => 'gk427998'
#
# If nothing is found, an empty array-ref is returned.

sub find_variations {
  my ($self, $pattern) = @_;

  my $info = $self->{'db'}->find_variations($pattern);
  return $info->{'matches'};

}

#############################################################################

#======================================================================
# FEATURES
#======================================================================

#############################################################################
# new_features($number)
# creates new Feature IDs
#
# Args:
# $number - number of Features to create
#
# Returns an array-ref of new Feature IDs with the batch ID

sub new_features {
  my ($self, $number) = @_;
  my %new_ids;

  my $info = $self->{'db'}->new_features($number);
  my @ids = @{$info->{'ids'}};
  
  return (\@ids, $info->{'id'});

}
#############################################################################
# kill a set of feature names
#
# Args:
# $names - an array-ref of feature IDs to kill
# $why is optional remark text.


sub kill_features {
  my ($self, $names, $why) = @_;

  my $info = $self->{'db'}->kill_features($names, $why);
  return $info;

}

#############################################################################
#Resurrect a batch of dead features.
#
# Args:
# $ids - an array-ref of feature IDs to resurrect


sub resurrect_features {
  my ($self, $ids) = @_;

  my $info = $self->{'db'}->resurrect_features($ids);
  return $info;

}

#############################################################################

#############################################################################

#======================================================================
# STRAINS
#======================================================================

#############################################################################
# new_strains($number)
# creates new Feature IDs
#
# Args:
# $number - number of Strains to create
#
# Returns an array-ref of new Feature IDs with the batch ID

sub new_strains {
  my ($self, $number) = @_;
  my %new_ids;

  my $info = $self->{'db'}->new_strains($number);
  my @ids = @{$info->{'ids'}};
  
  return (\@ids, $info->{'id'});

}
#############################################################################
# kill a set of feature names
#
# Args:
# $names - an array-ref of feature IDs to kill
# $why is optional remark text.


sub kill_strains {
  my ($self, $names, $why) = @_;

  my $info = $self->{'db'}->kill_strains($names, $why);
  return $info;

}

#############################################################################
#Resurrect a batch of dead strains.
#
# Args:
# $ids - an array-ref of feature IDs to resurrect


sub resurrect_strains {
  my ($self, $ids) = @_;

  my $info = $self->{'db'}->resurrect_strains($ids);
  return $info;

}

#############################################################################

#############################################################################
# RECENT CHANGES
#############################################################################

# The date is in ISO format UTC zone
# date --utc +%Y-%m-%dT%H:%M:%SZ`; # '2019-06-07T15:04:15Z'

sub recent_gene {
  my ($self, $from_date, $until_date, $agent) = @_;

  my $info = $self->{'db'}->recent_gene($from_date, $until_date, $agent);
  return $info;

}

#############################################################################

# The date is in ISO format UTC zone
# date --utc +%Y-%m-%dT%H:%M:%SZ`; # '2019-06-07T15:04:15Z'

sub recent_variation {
  my ($self, $from_date, $until_date, $agent) = @_;

  my $info = $self->{'db'}->recent_variation($from_date, $until_date, $agent);
  return $info;

}

#############################################################################

# The date is in ISO format UTC zone
# date --utc +%Y-%m-%dT%H:%M:%SZ`; # '2019-06-07T15:04:15Z'

sub recent_strain {
  my ($self, $from_date, $until_date, $agent) = @_;

  my $info = $self->{'db'}->recent_strain($from_date, $until_date, $agent);
  return $info;

}

#############################################################################

1;

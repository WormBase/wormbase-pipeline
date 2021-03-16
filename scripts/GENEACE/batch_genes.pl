#!/software/bin/perl -w
#
# batch_genes.pl
# 
# by Gary Williams                  
#
# This is for managing Gene IDs in the new (datomic) NameServer system
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/perl';
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use FileHandle;          # Or `IO::Handle' or `IO::'-anything-else used for autoflush

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;


=pod

=head batch_gene.pl

=item Options:

  -action    one of "new", "update", "kill", "resurrect", "suppress", "remove-cgc", "add-other-name", "remove-other-name", "merge", "split", "find", "help"
  -file      TAB or comma delimited file containing input IDs and old/new names  <Mandatory>
  -test      use the test nameserver

    for action "new":
    column 1 - either 'CGC' if this is an uncloned gene, or the biotype for a cloned gene, (one of 'cds','transcript','pseudogene','transposon' in upper or lower-case)
    column 2 - CGC name or Sequence name to give to the newly created WBGene ID

    example:
    CGC	bwa-1
    CGC	bwa-2
    Pseudogene	AC3.22

    

    for action "update":
    It should be noted that when updating CGC names, the command-line option -species is not used. The species of each of the given Gene IDs is checked against the new CGC name to ensure that the prefix is correct.
    If the format of a new CGC name is incorrect, it will be reported as an error, unless the command-line option -force is given.

    column 1 - WBGene ID to update
    column 2 - thing to update. One of 'CGC', 'Sequence', 'Biotype', 'Species'
    column 3 - new name to give to this WBGene ID e.g. 'bat-1' or 'AC3.24' or 'Pseudogene' or 'elegans'

    example:
    WBGene00304791	CGC	bat-1
    WBGene00304790	Sequence	AC3.24
    WBGene00304791	Biotype	Pseudogene
    WBGene00304792	Species	elegans


    for action "kill":
    column 1 - Gene ID / name to kill
    
    example:
    WBGene00304791
    WBGene00304792
    AAH1.1

    for action "resurrect":
    column 1 - Gene ID to resurrect
    
    example:
    WBGene00304791
    WBGene00304792


    for action "suppress":
    column 1 - Gene ID to suppress
    
    example:
    WBGene00304791
    WBGene00304792


    for action "remove-cgc"
    column 1 - CGC name to remove

    example:
    tag-100
    bat-1

    for action "add-other-name"
    column 1 - Gene ID
    column 2 - names of Other_name to add to Gene, separated by commas and/or spaces

    example:
    WBGene00304791	ehk09
    WBGene00304792	alpha-darwin, beta-wallace

    for action "remove-other-name"
    column 1 - Gene ID
    column 2 - names of Other_name to add to Gene, separated by commas and/or spaces

    example:
    WBGene00304791	ehk09
    WBGene00304792	alpha-darwin, beta-wallace

    for action "merge"

    example:
    column 1 - 'from' Gene ID - this is the gene that will die
    column 2 - 'into' Gene ID - this is the gene that will live
    column 3 - biotype of the resulting gene. One of 'cds','transcript','pseudogene','transposon'

    example:
    WBGene00304791	WBGene00304792	CDS
    WBGene00304793	WBGene00304794	pseudogene
    

    for action "split"

    example:
    column 1 - 'from' Gene ID - this is the existing gene that is to be split
    column 2 - new biotype to set on existing gene to be split. One of 'cds','transcript','pseudogene','transposon'
    column 3 - Sequence-name of the gene to be created.
    column 4 - biotype for the gene to be created. One of 'cds','transcript','pseudogene','transposon'

    example:
    WBGene00304791	CDS	H50C1.3	Pseudogene
    WBGene00304665	Pseudogene	AC3.4	CDS


    for action "find":
    column 1 -  patterns to find
     pattern can be any of a regular expression or part of the name of a WBVarID, or name
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

    for action "help" - no input is required, instructions on how to set up authentication to use the Nameserver are output


  -output    output file holding ace results
  -why       optional string describing the reason for performing the action
  -debug     limits reports to specified user <Optional>
  -species   used to specify non elegans genes
  

e.g. perl batch_genes.pl -species elegans -action new -file gene_name_data -output results_file.ace

=cut




######################################
# variables and command-line options # 
######################################

my ($test, $help, $debug, $verbose, $store, $wormbase, $logfile);
my ($species, $file, $output, $action, $why, $force);
my $BATCH_SIZE = 500; # maximum entries to put into any one batch API call

GetOptions (
	    "test"       => \$test,
	    "help"       => \$help,
            "debug=s"    => \$debug,
            "logfile:s"  => \$logfile,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "file:s"     => \$file,
	    "output:s"   => \$output,
	    "action:s"   => \$action,
	    "why:s"      => \$why,
	    "force"      => \$force, # if set, then accept any CGC name even it it does not conform to our format conventions 
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -organism => $species,
			     -test => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = $logfile ? Log_files->make_log($logfile, $debug) : Log_files->make_build_associated_log($wormbase);


if (!defined $file) {die "-file file not specified\n"}
if (!-e $file) {die "-file file doesn't exist\n"}

if (!defined $output) {die "-output file not specified\n"}

if (!defined $action) {die "-action not specified\n"}


open (IN, "<$file") || $log->log_and_die("Can't open file $file");
open (OUT, ">$output") || $log->log_and_die("Can't open file $output");
OUT->autoflush(1);       # empty the buffer after every line when writing to OUT

my $db = NameDB_handler->new($wormbase, $test);


if ($action eq 'new') {
  new_gene();
} elsif ($action eq 'update') {
  update_gene();
} elsif ($action eq 'kill') {
  kill_gene();
} elsif ($action eq 'resurrect') {
  resurrect_gene();
} elsif ($action eq 'suppress') {
  suppress_gene();
} elsif ($action eq 'remove-cgc') {
  remove_cgc_gene();
} elsif ($action eq 'add-other-name') {
  add_other_name();
} elsif ($action eq 'remove-other-name') {
  remove_other_name();
} elsif ($action eq 'merge') {
  merge_gene();
} elsif ($action eq 'split') {
  split_gene();
} elsif ($action eq 'find') {
  find_gene();
} elsif ($action eq 'help') {
  help_authentication();
} else {
  die "-action action '$action' not recognised\n"
}



close(OUT);
close(IN);
$db->close;


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

sub new_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my ($type, $name) = split /[\t,]/, $line;
    $type = lc $type;
    
    
    if ($type eq 'cgc') {
      push @names, {"cgc-name" => $name};
    } else {
      push @names, {"biotype" => $type, "sequence-name" => $name};
    }

    $count++;
    print OUT "//\tgene name '$line' queued for creation\n";

    if ($count == $BATCH_SIZE) {
      my ($new_ids, $batch) = $db->new_genes(\@names, $why);
      foreach my $hash (@{$new_ids}) {
	my $id = $hash->{'id'};
	my $cgc_name = $hash->{'cgc-name'};
	my $biotype = $hash->{'biotype'};
	my $seq_name = $hash->{'sequence-name'};
	my $species = $hash->{'species'};
	print OUT "\nGene : $id\n";
	print OUT "CGC_name $cgc_name\n" if (defined $cgc_name);
	print OUT "Public_name $cgc_name\n" if (defined $cgc_name);
	print OUT "Sequence_name $seq_name\n" if (defined $seq_name);
	print OUT "Public_name $seq_name\n" if (defined $seq_name);
	print OUT "Species $species\n" if (defined $species);
	print OUT "Version 1\n";
	print OUT "Version_change 1 now WBPerson1983 Event Imported \"Initial conversion from geneace\"\n";
	print OUT "Live\n";
	print OUT "\n\n";
      }

      $count = 0;
      @names = ();
      print OUT "// batch '$batch' created\n";
    }
  }

  if ($count) {
    my ($new_ids, $batch) = $db->new_genes(\@names, $why);
    foreach my $hash (@{$new_ids}) {
      my $id = $hash->{'id'};
      my $cgc_name = $hash->{'cgc-name'};
      my $biotype = $hash->{'biotype'};
      my $seq_name = $hash->{'sequence-name'};
      my $species = $hash->{'species'};
      print OUT "\nGene : $id\n";
      print OUT "CGC_name $cgc_name\n" if (defined $cgc_name);
      print OUT "Public_name $cgc_name\n" if (defined $cgc_name);
      print OUT "Sequence_name $seq_name\n" if (defined $seq_name);
      print OUT "Public_name $seq_name\n" if (defined $seq_name);
      print OUT "Species $species\n" if (defined $species);
      print OUT "Version 1\n";
      print OUT "Version_change 1 now WBPerson1983 Event Imported \"Initial conversion from geneace\"\n";
      print OUT "Live\n";
      print OUT "\n\n";
    }
    print OUT "// batch '$batch' created\n";
  }


}


##########################################

sub update_gene {
  my @names;
  my %names;
  my $batch;
  my $count = 0;

# should merge all data from a single ID together to make all update
# changes at the same time - this allows you to change both the
# species and the sequence-name at the same time so that the format of
# the new sequence-name is checked against the regex for the new
# species.

  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my ($id, $type, $name) = split /[\t,]/, $line;
    $type = lc $type;
    
    # -force on the command-line makes any CGC name acceptable
    # otherwise check that the CGC name is valid
    if (!$force && $type eq 'cgc' && !defined &check_cgc($id, $name)) {
      $log->write_to("ERROR, Invalid CGC name - not changed - use -force to force a change : $line\n");
      $log->error;

    } else {

      $names{$id}{'id'} = $id; 
      if ($type eq 'cgc') {
	$names{$id}{"cgc-name"} = $name;
      } elsif ($type eq 'sequence') {
	$names{$id}{"sequence-name"} = $name;
      } elsif ($type eq 'biotype') {
	$names{$id}{"biotype"} = $name;
      } elsif ($type eq 'species') {
	$names{$id}{"species"} = $name;
      }
    }
  }

  foreach my $id (keys %names) {
    $count++;
    print OUT "//\tgene name '$id' queued for creation\n";
    
    push @names, $names{$id};

    print OUT "\nGene : $id\n";
#    print OUT "Public_name $name\n\n";


    if ($count == $BATCH_SIZE) {
	my $batch = $db->update_genes(\@names, $why, $force);
      $count = 0;
      @names = ();
      print OUT "// batch '$batch' updated\n";
    }
  }

  if ($count) {
    my $batch = $db->update_genes(\@names, $why, $force);
    print OUT "// batch '$batch' updated\n";
  }


}
##########################################
# returns undef if the CGC name is not valid for the species of the $id gene
sub check_cgc {
  my ($id, $name) = @_;

  my $info = $db->info_gene($id);
  my $species = $info->{'species'};
  return $db->validate_name($name, 'CGC', $species);

}
##########################################

sub kill_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    print OUT "\nGene $line\n";
    print OUT "Dead\n\n";
    $count++;
    print OUT "//\tgene '$line' queued for being killed\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->kill_genes(\@names, $why);
      $count = 0;
      @names = ();
      $batch = $info->{'dead'}{'id'}; # was batch/id
      print OUT "// batch '$batch' killed\n";
    }
  }

  if ($count) {
    my $info = $db->kill_genes(\@names, $why);
    $batch = $info->{'dead'}{'id'}; # was batch/id
    print OUT "// batch '$batch' killed\n";
  }

}

##########################################

sub resurrect_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    print OUT "\nGene : $line\n";
    print OUT "Live\n\n";
    $count++;
    print OUT "//\tgene '$line' queued for being resurrected\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->resurrect_genes(\@names, $why);
      $count = 0;
      @names = ();
      $batch = $info->{live}{'id'}; # was batch/id
      print OUT "// batch '$batch' resurrected\n";
    }
  }

  if ($count) {
    my $info = $db->resurrect_genes(\@names, $why);
    $batch = $info->{live}{'id'}; # was batch/id
    print OUT "// batch '$batch' resurrected\n";
  }

}

##########################################

sub suppress_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    print OUT "\nGene : $line\n";
    print OUT "Suppressed\n\n";
    $count++;
    print OUT "//\tgene '$line' queued for being suppressed\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->suppress_genes(\@names, $why);
      $count = 0;
      @names = ();
      $batch = $info->{suppressed}{'id'}; # was batch/id
      print OUT "// batch '$batch' resurrected\n";
    }
  }

  if ($count) {
    my $info = $db->suppress_genes(\@names, $why);
    $batch = $info->{suppressed}{'id'}; # was batch/id
    print OUT "// batch '$batch' resurrected\n";
  }

}

##########################################

sub remove_cgc_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    $count++;
    print OUT "//\tCGC-name '$line' queued for being suppressed\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->remove_cgc_name_genes(\@names);
      $count =  0;
      @names = ();
      $batch = $info->{retracted}{'id'}; # was batch/id
      print OUT "// batch '$batch' resurrected\n";
    }
  }

  if ($count) {
    my $info = $db->remove_cgc_name_genes(\@names);
    $batch = $info->{retracted}{'id'}; # was batch/id
    print OUT "// batch '$batch' resurrected\n";
  }

}
##########################################

sub add_other_name {
  my %names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my @f = split /[\s,]+/, $line;
    my $geneid = shift @f;
    $names{$geneid} = [@f];
    $count++;
    print OUT "//\tOther-name '$line' queued for being added\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->add_other_name_genes(\%names, $why);
      $count =  0;
      %names = ();
      $batch = $info->{'id'}; # was batch/id
      print OUT "// batch '$batch'\n";
    }
  }

  if ($count) {
    my $info = $db->add_other_name_genes(\%names, $why);
    $batch = $info->{'id'}; # was batch/id
    print OUT "// batch '$batch'\n";
  }

}
##########################################

sub remove_other_name {
  my %names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my @f = split /[\s,]+/, $line;
    my $geneid = shift @f;
    $names{$geneid} = [@f];
    $count++;
    print OUT "//\tOther-name '$line' queued for being removed\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->remove_other_name_genes(\%names, $why);
      $count =  0;
      %names = ();
      $batch = $info->{retracted}{'id'}; # was batch/id
      print OUT "// batch '$batch' removed\n";
    }
  }

  if ($count) {
    my $info = $db->remove_other_name_genes(\%names, $why);
    $batch = $info->{retracted}{'id'}; # was batch/id
    print OUT "// batch '$batch' removed\n";
  }

}

##########################################

sub merge_gene {
  my @data;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}

    my ($from, $into, $biotype) = split /[\t,]/, $line;

    push @data, {"from-gene" => $from, "into-gene" => $into, "into-biotype" => $biotype};

    $count++;
    print OUT "//\tgenes '$line' queued for merge\n";
    print OUT "\nGene : $into\n";
    print OUT "\n\n";


    if ($count == $BATCH_SIZE) {
      my ($info) = $db->batch_merge_genes(\@data, $why);
      $batch = $info->{'id'}; # was batch/id

      $count = 0;
      @data = ();
      print OUT "// batch '$batch' created\n";
    }
  }

  if ($count) {
    my ($info) = $db->batch_merge_genes(\@data, $why);
    $batch = $info->{'id'}; # was batch/id
    print OUT "// batch '$batch' created\n";
  }

}

##########################################

sub split_gene {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my ($from, $new_biotype, $product_sequence_name, $product_biotype) = split /[\t,]/, $line;
    
    push @names, {"from-id" => $from, "new-biotype" => $new_biotype, "product-sequence-name" => $product_sequence_name, "product-biotype" => $product_biotype};


    $count++;
    print OUT "//\tgene name '$line' queued for creation\n";

    if ($count == $BATCH_SIZE) {
      my ($new_ids, $batch) = $db->batch_split_genes(\@names, $why);
      foreach my $hash (@{$new_ids}) {
	my $from_id = $hash->{'from-id'};
	my $from_biotype = $hash->{'new-biotype'};
	my $product_id = $hash->{'product-id'};
	my $product_seq_name = $hash->{'product-sequence-name'};
	my $product_biotype = $hash->{'product-biotype'};
	print OUT "\nGene : $product_id\n";
	print OUT "Sequence_name $product_seq_name\n";
	print OUT "Public_name $product_seq_name\n";
	print OUT "Version 1\n";
	print OUT "Version_change 1 now WBPerson1983 Event Imported \"Initial conversion from geneace\"\n";
	print OUT "Live\n";
	print OUT "\n\n";
      }

      $count = 0;
      @names = ();
      print OUT "// batch '$batch' split\n";
    }
  }

  if ($count) {
    my ($new_ids, $batch) = $db->batch_split_genes(\@names, $why);
    foreach my $hash (@{$new_ids}) {
      my $from_id = $hash->{'from-id'};
      my $from_biotype = $hash->{'new-biotype'};
      my $product_id = $hash->{'product-id'};
      my $product_seq_name = $hash->{'product-sequence-name'};
      my $product_biotype = $hash->{'product-biotype'};
      print OUT "\nGene : $product_id\n";
      print OUT "Sequence_name $product_seq_name\n";
      print OUT "Public_name $product_seq_name\n";
      print OUT "Version 1\n";
      print OUT "Version_change 1 now WBPerson1983 Event Imported \"Initial conversion from geneace\"\n";
      print OUT "Live\n";
      print OUT "\n\n";
    }
    print OUT "// batch '$batch' split\n";
  }

}

##########################################

sub find_gene {
  my $pattern;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}

    my ($pattern) = split /[\t,]/, $line;
    my $info = $db->find_genes($pattern);

    print OUT "// Genes matching '$pattern'\n";
    foreach my $result (@{$info}) {
      my $id = $result->{'id'}; # was gene/id
      my $seqname = $result->{'sequence-name'}; # was gene/sequence-name
      if (!defined $seqname) {$seqname = '.'}
      my $cgcname = $result->{'cgc-name'}; # was gene/cgc-name
      if (!defined $cgcname) {$cgcname = '.'}
      print OUT "Gene $id // $seqname  $cgcname\n";
    }
  }

}

##########################################
sub help_authentication {
  $db->print_authentication_instructions();
}

##########################################
##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut

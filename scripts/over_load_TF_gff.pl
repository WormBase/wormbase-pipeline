#!/software/bin/perl -w                 
#
# This is to add information of the Transcription Factor name and ID to TF_binding_site Features
# These GFF lines have a Type column of 'TF_binding_site'
# It also adds the Feature Public_name to lines where this exists.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-12-04 10:15:46 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase,$database);
my ($species, $gff_file);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database,
	    "file:s"     => \$gff_file
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$wormbase->{autoace} = $database if $database;
$species = $wormbase->species;

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################

my %TF;
my %PN;


#load TF and Public_name details from table maker query
my $table = $wormbase->table_maker_query($wormbase->autoace, &write_def_file);
while(<$table>) {
  s/\"//g; #"
  next if (/acedb/ or /\/\//);
  chomp;
  my ($feature, $tf_id, $tf_name, $public_name) = split(/\t/,$_);
  if (defined $feature && $tf_id ne '') {
    $TF{$feature}->{'id'} = $tf_id;
    $TF{$feature}->{'name'} = $tf_name;
  }
  if (defined $feature && $public_name ne '') {
    $PN{$feature} = $public_name;
  }
}

my $stat = 0;

my @gff_files;
if ($gff_file) {
  if (not -e $gff_file or -z $gff_file) {
    $log->log_and_die("Non-existent or zero length GFF file");
  }
  @gff_files = ($gff_file);
} else {
  if ($wormbase->assembly_type eq 'contig'){
    @gff_files = ($wormbase->species);
  } else {
    @gff_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
  }
  for(my $i=0; $i < @gff_files; $i++) {
    $gff_files[$i] = sprintf("%s/%s.gff", $wormbase->gff, $gff_files[$i]);
    if (not -e $gff_files[$i] or -z $gff_files[$i]) {
      $log->log_and_die("Non-existent or zero-length GFF file $gff_files[$i]");
    }
  }
}

foreach my $file (@gff_files) {
  
  open(GFF,"<$file") or $log->log_and_die("cant open $file");
  open(NEW,">$file.tmp") or $log->log_and_die("cant open $file tmp file\n");
  while( <GFF> ) {
    chomp;	
    print NEW "$_";
    #CHROMOSOME_V    TF_binding_site  TF_binding_site     155950  155951  .       +       .       Feature "WBsf046537"
    my ($feature) = (/Feature \"(\S+)\"/);
    if (defined $feature && exists $TF{$feature}) {
      my $id = $TF{$feature}{'id'};
      my $name = $TF{$feature}{'name'};
      print NEW " ; TF_ID \"$id\" ; TF_name \"$name\"";
    }
    if (defined $feature && exists $PN{$feature}) {
      my $public_name = $PN{$feature};
      print NEW " ; Public_name \"$public_name\"";
    }
    print NEW "\n";
  }
  $wormbase->run_command("mv -f $file.tmp $file", $log);
}

##################
# Check the files
##################

if($gff_file){
    $log->write_to("Not checking ad hoc file\n");
}
else { 
  foreach my $file (@gff_files) {
    my $minsize = ($file=~/random|un/)?170000:1500000;
    $wormbase->check_file($file, $log,
			      minsize => $minsize,
			      lines => ['^##',
					"^\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			      );
    }
}


# Close log files and exit
$log->write_to("\n\nChanged $stat lines\n");
$log->write_to("----------\n\n");

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

sub write_def_file {
	my $def = '/tmp/overload_TF_GFF.def';
	open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
	my $txt = <<END;
// Spread sheet definition for the ACeDB software 
// User: gw3
// Date: 2012-07-13_14:29:33

// %n (%%n in the graphic) are parameter to be given on the command line in tace
// or by default by the Parameters given in this file
// \%n (%n in the graphic) are substituted by the value of column n at run time
// Line starting with // are ignored, starting with # are comments

Sortcolumn 1

Colonne 1 
Subtitle Feature 
Width 12 
Optional 
Visible 
Class 
Class Feature 
From 1 
 
Colonne 2 
Subtitle TF 
Width 40 
Optional 
Visible 
Class 
Class Transcription_factor 
From 1 
Tag Transcription_factor 
 
Colonne 3 
Subtitle Name 
Width 40 
Optional 
Visible 
Text 
From 2 
Tag Name 
 
Colonne 4 
Width 80 
Optional 
Visible 
Class 
Class Text 
From 1 
Tag Public_name 

// End of these definitions
END

	print TMP $txt;
	close TMP;
	return $def;
}




__END__
# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


=pod

=head2 NAME - over_load_TF_gff.pl

=head1 USAGE

=over 4

=item over_load_TF_gff.pl  [-options]

=back

This script adds names and TranscriptionFactor object IDs to TF_binding_site lines

=over 4

=item None at present.

=back

$0.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

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

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut

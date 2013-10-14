#!/usr/local/bin/perl
#
# map_feature2gene.pl
#
# make the connections between TSL features and CDS/Gene
#
# Dan Lawson
#
# Usage : map_feature2gene.pl [-options]
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-10-14 09:39:32 $
#################################################
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;

#####################################
# variables and command-line options #
######################################
my $maintainers = 'All';

my $test;       # In test build mode
my $help;       # Help/Usage page
my $verbose;    # Verbose mode
my $noload;     # do not upload results to autoace
my $debug;      # Debug mode
my $store;      # configuration file
my $species;

GetOptions(
    'help'    => \$help,
    'noload'  => \$noload,
    'verbose' => \$verbose,
    'test'    => \$test,
    'debug=s' => \$debug,
    'store=s' => \$store,
    'species=s'=>\$species,
);

# Help pod if needed
&usage(0) if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user
my $tace   = $wb->tace;                                        # tace executable path
my $dbdir  = $wb->autoace;                                     # Database path
my $output = $wb->acefiles."/TSL_feature_connections.ace";    # output file path

# create log
my $log = Log_files->make_build_log($wb);

#######################
# MAIN BODY OF SCRIPT #
#######################

open( OUTPUT, ">$output" ) || die "Can't open output file $output\n";

#####################
# open a connection #
#####################

print "// Connect to database $dbdir\n" if ($verbose);
my $db = Ace->connect(
    -path    => $dbdir,
    -program => $tace
  )
  || do { print "Connection failure: ", Ace->error; die};

print "// done\n" if ($verbose);

my @TSL = $db->fetch( -class => 'Feature', -name => '*' );

# Loop through each Feature object

my $method;
my @sequence_name;
my $sequence;
my @CDS;
my $seq;
my $i;
my $cds;
my %match;
my $sum;

foreach my $feature (@TSL) {

    $method = $feature->Method(1);
    next unless ( ( $method eq "SL1" ) || ( $method eq "SL2" ) );
    @sequence_name = $feature->Defined_by_sequence(1);
    undef %match;    # reset hash

    foreach $i (@sequence_name) {
        print "// defined by sequence $i\n" if ($verbose);

        $sequence = $db->fetch( Sequence => $i );
        @CDS      = $sequence->Matching_CDS(1);

        foreach $cds (@CDS) {
            next if ( $cds eq "" );
            $match{$cds} = 1;
            print "// matches CDS '$cds'\n" if ($verbose);
        }
    }

    # write output and error tracking
    unless ( %match ) {
        print "\n// No connection to a CDS found for Feature : $feature\n";
    }
    else {
        print OUTPUT "Feature : \"$feature\"\n";
        foreach my $j ( keys %match ) {
            print OUTPUT "Associated_with_CDS $j\n";
        }
    }
    print OUTPUT "\n";

    $sequence->DESTROY() if $sequence;
    $feature->DESTROY()  if $feature;
}

close OUTPUT;
$db->close;

# Upload file to autoace (if you have been asked to)
unless ($noload) {
  $log->write_to("Loading file to autoace\n");
  $wb->load_to_database( $wb->autoace, $output, 'TSL_CDS_connect', $log );
}

$log->mail();
exit(0);

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}
############################################

__END__

=pod

=head2 NAME - map_feature2gene.pl

=head1 USAGE

=over 4

=item map_feature2gene.pl [-options]

=back

map_feature2gene.pl connects TSL Feature objects to CDS objects
via the defined_by tag in the database.

mandatory arguments:

=over 4

=item none

=back

mmap_feature2gene.pl optional arguments:

=over 4

=item -debug, Debug mode

=item -verbose, Verbose mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=item -outfile specify location for the output Ace (if you want to use it outside of the build)

=item -store specifiy a configuration file

=back

=cut


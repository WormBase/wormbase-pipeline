# Author: jkh1

=head1 NAME

 Treefam::Config

=head1 DESCRIPTION

 This modules contains default configuration settings,
 mostly database connection details. It only exports
 variables.

=head1 CONTACT

 jkh1@sanger.ac.uk

=cut

package Treefam::Config;

use strict;
use Exporter;

our @ISA = ('Exporter');

our @EXPORT = qw($APIVERSION $TFDBHOST $TFDBNAME $TFDBUSER $TFDBPASS $TFDBPORT $TFSEQDIR $NJTREE);

our $APIVERSION=4;

# Treefam database
our $TFDBNAME='treefam_7'; # was 4 - mh6
our $TFDBUSER='anonymous';
our $TFDBHOST='vegasrv.sanger.ac.uk';
our $TFDBPASS='';
our $TFDBPORT=3308;

# path to TFSEQ
our $TFSEQDIR='/nfs/treefam/treeFam/data/release-4.0/TFSEQ';

# njtree program
our $NJTREE='/nfs/faculty/jkh1/programs/njtree';

# HMMER
our $HMMER_HOST = 'cbi1';
our $HMMER = '';

1;

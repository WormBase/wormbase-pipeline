#!/usr/local/bin/perl -w
######################### Overlap.pl ################################
#
# Copyright (c) 1999  Richard Bruskiewich (rbsk@sanger.ac.uk)
# Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation
#
# Name: Overlap.pl
#
# Function: 
#    This script extracts feature sets from 
#    an input GFF file, 
#
# History:
#   (Re)Created: 19/11/99 - adapted from old Features.pl
#
####################################################################
BEGIN {
 unshift (@INC,"/nfs/disk80/rbsk/perl/PerlModules");
}
use strict ;
use Cwd ;
use MyUtilLib ;
use GFF ;
use GFF::Analysis qw(normalize_mRNA) ;

use vars qw( $opt_T $opt_D $opt_V $opt_h 
	     $opt_p $opt_i $opt_x $opt_y  
	     $opt_s $opt_r $opt_S $opt_f $opt_g 
	     $opt_G $opt_H $opt_t $opt_u $opt_M $opt_m );

######## Forward Function Declarations ####### 

sub help ;
sub trace { print STDERR @_ if ($opt_D) ; } # debugging function in -verbose mode
sub dbug { print STDERR @_ if ($opt_T or $opt_D) ; }  # debugging function in -verbose mode
sub admonish { print STDERR "\n@_\n\n" ; help ; }

######## File Scope Initializations ##########

$| = 1 ; # flush stdout every time
my $def_p = Cwd::cwd ;

my $def_x = 'gff' ;
my $def_s = '\w+' ;
my $def_r = '\w+' ;
my $def_f = '\w+' ;
my $def_g = '\w+' ;
my $def_t = '\w+' ;
my $def_u = '\w+' ;

######## Get Command Line Parameters #########

my $commandline = "$0 @ARGV" ; # echo command line before getopts

use Getopt::Std;
getopts('hTDVi:p:x:y:s:r:S:f:g:G:H:t:u:M:m:');

if($opt_h){
    help;
} elsif(!defined($opt_i)) {
    admonish "*** Error: You must provide a GFF input file name!\n" ;
}

trace("$commandline\n*** Overlap.pl Processing Started ***\n\n") ;

my $path  = defined($opt_p) ? $opt_p: $def_p ;
trace("GFF input file path:\t\t$path\n") ;

my $input_suffix  = defined($opt_x) ? $opt_x : $def_x ;
trace("GFF input file 1 suffix:\t\t$input_suffix\n") ;

my $inputfile1 = filAbsPath("$opt_i.$input_suffix",$path) ;
open(GFFIN1,"<$inputfile1")
    or admonish "*** Error: Cannot open input file $inputfile1: $!!" ;
my $GFFIN1 = \*GFFIN1 ;

my $GFFIN2 ;
my $inputfile2 ;
if(defined($opt_y)) {
    trace("GFF input file 2 suffix:\t\t$opt_y\n")  ;
    $inputfile2 = filAbsPath("$opt_i.$opt_y",$path) ;
    open(GFFIN2,"<$inputfile2")
	or admonish "*** Error: Cannot open input file $inputfile2: $!!" ;
    $GFFIN2 = \*GFFIN2 ;
} else {
    $inputfile2 = $inputfile1 ;
    $GFFIN2 = \*GFFIN1 ; # default to same input file
}

my $source_pattern1  = defined($opt_s) ? $opt_s : $def_s ;
trace("Source field pattern 1:\t\t$source_pattern1\n") ;

my $feature_pattern1 = defined($opt_f) ? $opt_f : $def_f ;
trace("Feature field pattern 1:\t\t$feature_pattern1\n") ;

my $source_pattern2  = defined($opt_r) ? $opt_r : $def_r ;
trace("Source field pattern 2:\t\t$source_pattern2\n") ;

my $feature_pattern2 = defined($opt_g) ? $opt_g : $def_g ;
trace("Feature field pattern 2:\t\t$feature_pattern2\n") ;

my $tag1 = defined($opt_G) ? $opt_G : undef ;
trace("Group tag name:\t\t$tag1\n") ;

my $tag2 = defined($opt_H) ? $opt_H : undef ;
trace("Group tag name:\t\t$tag2\n") ;

my $tag_value_filter1 = defined($opt_t) ? $opt_t : $def_t ;
trace("Group tag value filter pattern:\t$tag_value_filter1\n") ;

my $tag_value_filter2 = defined($opt_u) ? $opt_u : $def_u ;
trace("Group tag value filter pattern:\t$tag_value_filter2\n") ;

my $clustertag = defined($opt_M) ? $opt_M : undef ;
trace("Cluster tag name:\t\t$clustertag\n") if defined($clustertag) ;

my $cluster_tag_value_filter = defined($opt_m) ? $opt_m : $tag_value_filter2 ;
trace("Cluster tag value filter pattern:\t$cluster_tag_value_filter\n") ;

############### Read in GFF subset # 1 ###############

trace("Reading GFF set #1 from input file:\t$inputfile1\n") ;

my $feature_filter1 = sub {
    my $self = shift ;
    my $match = eval { $self->feature() =~ /$feature_pattern1/i } ; # match the feature field
    if($@) { 
	die "*** Error: Improper feature pattern!" ;
    }
    if($match) {
	if(defined($tag1)) {
	    my $tagvalue = $self->group_value($tag1) ;
	    return 0 if !defined($tagvalue) ;
	    $match = eval { $tagvalue =~ /$tag_value_filter1/i } ; # match the feature field
	    if($@) { 
		die "*** Error: Improper feature pattern!" ;
	    }
	    return 0 if(!$match) ;
	}
	return $self ;
    } else {
        return 0 ;
    }
} ;

my $source_filter1 = sub {
    my $line = shift ;
    my $match = eval { $line =~ /^\S+\s+$source_pattern1\s+/i } ; # match the source field
    if($@) { 
	die "*** Error: Improper source pattern!" ;
    }
    if($match) {
	return $line ;
    } else {
        return 0 ;
    }
} ;

my $gff1 = new GFF::GeneFeatureSet(2) ;
GFF::trace(1) if $opt_T or $opt_D ;
$gff1->read( $GFFIN1, $feature_filter1, $source_filter1 ) ;
GFF::trace(0) if $opt_T or $opt_D ;

my $count = $gff1->count() ;
if(!$count) {
    trace("No $source_pattern1($feature_pattern1) features found in file $inputfile1. Exiting Overlap.pl...\n") ;
    exit ;
}

if($opt_V) {
    print STDERR "\nDiagnostic dump of GFF input set 2:\n\n" ;
    $gff1->dump_header(\*STDERR) ;
    $gff1->dump(\*STDERR) ;
}
trace("$count GFF records read in for subset #1\n") ;

############### Read in GFF subset # 2 ###############

trace("Reading GFF set #2 from input file:\t$inputfile2\n") ;

my $feature_filter2 = sub {
    my $self = shift ;
    my $match = eval { $self->feature() =~ /$feature_pattern2/i } ; # match the feature field
    if($@) { 
	die "*** Error: Improper feature pattern!" ;
    }
    if($match) {
	if(defined($tag2)) {
	    my $tagvalue = $self->group_value($tag2) ;
	    return 0 if !defined($tagvalue) ;
	    $match = eval { $tagvalue =~ /$tag_value_filter2/i } ; # match the feature field
	    if($@) { 
		die "*** Error: Improper feature pattern!" ;
	    }
	    return 0 if(!$match) ;
	}
	return $self ;
    } else {
        return 0 ;
    }
} ;


my $source_filter2 = sub {
    my $line = shift ;
    my $match = eval { $line =~ /^\S+\s+$source_pattern2\s+/i } ; # match the source field
    if($@) { 
	die "*** Error: Improper source pattern!" ;
    }
    if($match) {
	return $line ;
    } else {
        return 0 ;
    }
} ;

my $gff2 = new GFF::GeneFeatureSet(2) ;
GFF::trace(1) if $opt_T or $opt_D ;
$gff2->read( $GFFIN2, $feature_filter2, $source_filter2 ) ;
GFF::trace(0) if $opt_T or $opt_D ;

unless($GFFIN2 == $GFFIN1) {
    close($GFFIN2) ;
}
close($GFFIN1) ;

$count = $gff2->count() ;
if(!$count) {
    trace("No $source_pattern2($feature_pattern2) features found in file $inputfile2. Exiting Overlap.pl...\n") ;
    exit ;
}

if($opt_V) {
    print STDERR "\nDiagnostic dump of GFF input set 2:\n\n" ;
    $gff2->dump_header(\*STDERR) ;
    $gff2->dump(\*STDERR) ;
}
trace("$count GFF records read in for subset #2\n") ;

############### Look for overlaps ###############

trace("Performing intersect overlap match comparison\n") ;
GFF::trace(1) if $opt_T or $opt_D ;
my $gffo = $gff1->intersect_overlap_matches( $gff2,-1,0,0,1,\*STDERR) ; # any overlap, strand insensitive, report overlaps to STDERR
GFF::trace(0) if $opt_T or $opt_D ;

########### Identify the matching ESTs ##########

my $set_filter = sub {
    my $gf1 = shift ;
    my $gf2 = shift ;
    my ($name1) = ($gf1->group_value($clustertag) =~ /$cluster_tag_value_filter/) ;
    my ($name2) = ($gf2->group_value($clustertag) =~ /$cluster_tag_value_filter/) ;
    if( defined($name1) and  
	defined($name2) and 
	$name1 eq $name2 ) {
	return $name1 ;
    } else {
	0;
    }
} ;

my @gffsets ;
if(defined($opt_M)) {
    @gffsets = $gff2->cluster($set_filter) ;
} else {
    push(@gffsets,$gff2) ; # just one set?
}

################  Output the results ############

trace("Dumping refined GFF matches to STDOUT\n") ;
print STDOUT "Report of $source_pattern2($feature_pattern2) GFF set hits against $source_pattern1($feature_pattern1) set:\n\n" ;
foreach my $gffset (@gffsets) {
    my ($name) = $gffset->region() ;
    print STDOUT "$name ","#"x20,"\n" ;
    $gffset->dump_matches(\*STDOUT) ;
}

trace("\n*** Overlap.pl Processing Completed ***\n\n") ;

######## Subroutine Declarations #############

sub help {
    print <<ENDHELP;
Function:

    This script performs a intersect overlap comparison between two
    defined gff subsets taken from a single given input file.
    The GFF match output is reported to STDOUT.

Usage:

   Overlap.pl <input spec> [output spec] [options]
or
   Overlap.pl -h

Input:

  -i  filename  # Input file

  -x word       [$def_x]
                # GFF input file 1 suffix; input file is <INPUTFILE>.<x_suffix>

  -y word       # Optional) GFF input file 2 suffix; input file is <INPUTFILE>.<y_suffix>
                # if not defined, defaults to taking the second GFF subset from GFF file 1. 

Options:

  -p  pathname  [$def_p]
                # Path for input file (default: current working directory)

  -s pattern    [$def_s]
                # simple identifier or Perl regex pattern to match
                # the GFF input <source> field for the first file (default: all)

  -r pattern    [$def_r]
                # simple identifier or Perl regex pattern to match
                # the input <source> field  for the first file (default: all)

  -f pattern    [$def_f]
                # simple identifier or Perl regex pattern to match
                # the input <feature> field  for the second file (default: all)

  -g pattern    [$def_g]
                # simple identifier or Perl regex pattern to match
                # the input <feature> field  for the second file (default: all)

  -G pattern    # optional simple identifier or Perl regex pattern to match [group] tag name
                # in GFF subset # 1

  -H pattern    # optional simple identifier or Perl regex pattern to match [group] tag name
                # in GFF subset # 2

  -t pattern    # optional simple identifier or Perl regex pattern to 
                # match [group] tag value in GFF subset # 1

  -u pattern    # optional simple identifier or Perl regex pattern to 
                # match [group] tag value in GFF subset # 2

  -M pattern    # optional simple identifier or Perl regex pattern to match [group] tag name
                # in GFF subset # 2 for clustering operations
                # if defined, the script performs a filtered comparison driven
                # clustering operation on GFF set # 2, based
                # upon the -M (and -m) values, and reports the matches accordingly

  -m pattern    # optional simple identifier or Perl regex pattern to 
                # filter out [group] tag values from two features for clustering;

  -S word       # <source> field identifier for <source> relabelling of 
                # output GFF (default: no relabelling)

  -h            # for help

  -T            # run trace output

  -D            # verbose debug output

  -V            # very verbose debug output (dumps intermediate GFF sets)

ENDHELP

    exit 0;
}


__END__
# Need  POD docs here

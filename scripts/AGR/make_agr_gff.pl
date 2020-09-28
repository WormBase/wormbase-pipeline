#!/usr/bin/env perl
#
# Script to create GFF according to AGR 1.0.1.3 schema as specified in
# https://docs.google.com/document/d/1yAECtOs1VCEs3mhplJMXqg1akBQJyJbjPnvG2RrE5Aw/edit


use strict;
use Getopt::Long;
use URI::Escape;

use lib $ENV{CVS_DIR};
use Modules::AGR;

my ($bgi_json, $gff_in, $gff_out, $ws_version);

GetOptions (
    "gffout:s"    => \$gff_out,
    "gffin=s"     => \$gff_in,
    "bgijson=s"   => \$bgi_json,
    "wsversion=s" => \$ws_version,
    )||die();

$gff_in   || die "You must supply an input GFF file\n";
$gff_out  || die "You must supply an output file name\n";
$bgi_json || die "You must supply an BGI JSON file\n";

my %bgi_genes = %{AGR::get_bgi_genes($bgi_json)};
my $date      = AGR::get_rfc_date();

my $in_fh = &open_gff_file($gff_in);
open(my $out_fh, ">$gff_out") or die "Could not open $gff_out for writing\n";
print $out_fh "# WormBase release $ws_version\n" if defined $ws_version;

my %headers; # do not print redundant headers

while(<$in_fh>) {
    if (/^\#/){
	next if $headers{$_};
	print $out_fh $_;
	$headers{$_}=1;
	if (/gff-version 3/){
	    print $out_fh 
		"#!date-produced $date\n".
		"#!data-source WormBase\n".
		"#!assembly WBcel235\n";
	}
    next;
    }
    
    if (/^\S+\s+WormBase\s+gene\s+/) {
	chomp;
	my @l = split(/\t/, $_);
	my %attr;
	foreach my $kv (split(/;/, $l[8])) {
	    my ($k, $v) = split(/\=/, $kv);
	    $attr{$k} = $v;
	}
	my $gid = $attr{curie};

	$attr{gene_id} = $attr{curie}; # Duplication seems strange but requested in AGR doc
	
	# use symbol as Name
	$attr{Name} = $bgi_genes{$gid}->{symbol};
	if (exists $bgi_genes{$gid}->{name}) {
	    $attr{long_name} = $bgi_genes{$gid}->{name};
	}
	if (exists $bgi_genes{$gid}->{geneSynopsis}) {
	    $attr{description} = uri_escape($bgi_genes{$gid}->{geneSynopsis},"\t=%;,&");
	}
	$attr{Ontology_term}='SO:0000704';
	
	$l[8] = join(";", map { $_ . "=" . $attr{$_} } keys %attr);
	print $out_fh join("\t", @l), "\n";
	
    } elsif(/^\S+\s+WormBase\s+mRNA/){
	change_transcript($_,$out_fh,'SO:0000234')
    } elsif(/^\S+\s+WormBase\s+pseudogenic_transcript/){
	change_transcript($_,$out_fh,'SO:0000516')
    } elsif(/^\S+\s+WormBase\s+ncRNA/){
	change_transcript($_,$out_fh,'SO:0000655')
    } elsif(/^\S+\s+WormBase\s+CDS\s/){
	change_cds($_,$out_fh)
    } elsif (/^\S+\s+WormBase\s+/){
	print $out_fh $_;
    }
}
close($out_fh) or die "Could not close $gff_out after writing (probably file system full)\n";

exit(0);

sub change_cds {
    my ($line, $out_fh) = @_;

    # Changes protein_id attribute to WB curie version of wormpep ID, replacing ENA ID
    $line =~ s/;protein_id=/;ena_id=/;
    $line =~ s/;wormpep=/;protein_id=WB:/;

    print $out_fh $line;
}

sub change_transcript{
    my($line,$outf,$soID)=@_;
    my @l = split(/\t/, $line);
    my $transcriptID="$1" if $line =~ /Name=([^;\n]+)/;
    my $geneID="$1" if $line =~/(WBGene\d+)/;
    # $line=~s/;Parent/;curie=WB:$transcriptID;curie=WB:$geneID;Ontology_term=$soID;Parent/;
    $line=~s/;Parent/;transcript_id=WB:$transcriptID;curie=WB:$transcriptID;Ontology_term=$soID;Parent/;
    print $outf $line;
}

sub open_gff_file {
    my ($file) = @_;
    
    my $fh;
    if ($file =~ /\.gz$/) {
    open( $fh, "gunzip -c $file |") or die "Could not open gunzip stream to $file\n";
    } else {
	open( $fh, $file) or die "Could not open $file for reading\n";
    }
    
    return $fh;
}


#!/usr/bin/env perl
use strict;
use Const::Fast;
use Path::Class;
use Text::Tabs;
use Getopt::Long;
use Smart::Comments;

const my $WB_MODEL_SPEC => $ENV{WORMPUB} . "/wormbase-pipeline/wspec/models.wrm";

GetOptions(
    "out|o=s" => \my $out_file
    );

my $model_lines = get_and_parse_model_lines();
my ($class_slots, $hash_slots, $comments, $controlled_vocabs) = extract_slot_details($model_lines);

$controlled_vocabs = filter_controlled_vocabs($controlled_vocabs);
$class_slots = merge_controlled_vocabs($class_slots, $controlled_vocabs);
$hash_slots = merge_controlled_vocabs($hash_slots, $controlled_vocabs);

my ($classes, $associations, $slots);
for my $class (keys %$class_slots) {
    $classes->{$class}{'description'} = "$class class";
    for my $slot_value (@{$class_slots->{$class}}) {
	($classes, $associations, $slots) = process_slot_value($class, $slot_value->{'slot'},
							       $slot_value->{'value'}, $classes,
							       $associations, $slots, $hash_slots,
							       $comments, '', $controlled_vocabs);		       
    }
}

($classes, $slots) = avoid_duplicate_class_and_slot_names($classes, $slots);
print_yaml($out_file, $classes, $associations, $slots);


sub avoid_duplicate_class_and_slot_names {
    my ($classes, $slots) = @_;

    my (%final_slots, %changed_slots);
    for my $slot (keys %$slots) {
	if (exists $classes->{$slot} or $slot eq 'date') {
	    $final_slots{$slot . '_slot'} = $slots->{$slot};
	    $changed_slots{$slot}++;
	}
	else {
	    $final_slots{$slot} = $slots->{$slot};
	}
    }
    
    for my $class (keys %$classes) {
	for my $ix (0 .. scalar @{$classes->{$class}{'slots'}} - 1) {
	    if (exists $changed_slots{$classes->{$class}{'slots'}[$ix]}) {
		$classes->{$class}{'slots'}[$ix] .= '_slot';
	    }
	}
    }
    
    return ($classes, \%final_slots);
}


sub print_yaml {
    my ($out_file, $classes, $associations, $slots) = @_;

    my $header = <<'YAML_HEADER';

id: https://github.com/WormBase/wormbase-pipeline/tree/master/wspec/wormbase.yaml
name: WormBase-Prototype-BioLinkML-Schema
title: WormBase Prototype BioLinkML Schema
description: >-
  WormBase Schema Prototype with LinkML
license: https://creativecommons.org/publicdomain/zero/1.0/

version: 1.0.0

prefixes:
  WB: 'https://wormbase.org/'
  Alliance: 'http://alliancegenome.org/'
  biolinkml: 'https://w3id.org/biolink/biolinkml/'

default_prefix: WB
default_range: string

imports:
  - biolinkml:types

classes:

  named thing:
    is_a: entity
    description: "a databased entity or concept/class"

  entity:
    description: >-
      Root Biolink Model class for all things and informational relationships, real or imagined.
    abstract: true

  association:
    is_a: entity
    description: >-
      A typed association between two entities, supported by evidence
    comments:
      - This is roughly the model used by biolink and ontobio at the moment
    slots:
      - subject
      - predicate
      - object
    slot_usage:
      type:
        description: rdf:type of biolink:Association should be fixed at rdf:Statement
      category:
        range: association
    exact_mappings:
      - biolink:Association

YAML_HEADER


    my $footer = <<'YAML_FOOTER';

  node property:
    description: >-
      A group for any property that holds between a node and a value

  related to:
    description: >-
      A relationship that is asserted between two named things
    domain: named thing
    range: named thing
    multivalued: true
    inherited: true
    symmetric: true

  subject:
    is_a: association slot
    local_names:
      ga4gh: annotation subject
      neo4j: node with outgoing relationship
    description: >-
      Connects an association to the subject of the association.
      For example, in a gene-to-phenotype association, the gene is subject and phenotype is object.
    required: true
    range: named thing
    exact_mappings:
      - biolink:subject

  object:
    is_a: association slot
    description: >-
      Connects an association to the object of the association.
      For example, in a gene-to-phenotype association, the gene is subject and phenotype is object.
    required: true
    range: named thing
    local_names:
      ga4gh: descriptor
      neo4j: node with incoming relationship
    exact_mappings:
      - biolink:object

  predicate:
    is_a: association slot
    description: >-
      A high-level grouping for the relationship type. AKA minimal predicate.
      This is analogous to category for nodes.
    domain: association
    notes: >-
      Has a value from the Biolink related_to hierarchy. In RDF,  this
      corresponds to rdf:predicate and in Neo4j this corresponds to the
      relationship type. The convention is for an edge label in snake_case
      form. For example, biolink:related_to, biolink:causes, biolink:treats
    range: predicate type
    required: true
    local_names:
      ga4gh: annotation predicate
      translator: predicate
    exact_mappings:
      - biolink:predicate

  association slot:
    abstract: true
    domain: association
    aliases: ['edge property', 'statement property']
    description: >-
      any slot that relates an association to another entity
    exact_mappings:
      - biolink:association_slot

types:

  predicate type:
    typeof: uriorcurie
    description: >-
      A CURIE from the biolink related_to hierarchy.
      For example, biolink:related_to, biolink:causes, biolink:treats.

YAML_FOOTER

    my $out_fh = file($out_file)->openw;
    $out_fh->print($header);
    for my $class (sort keys %$classes) {
	replace_pipes_and_print($out_fh, "  $class:\n");
        replace_pipes_and_print($out_fh, "    description: " . $classes->{$class}{'description'} . "\n");
	replace_pipes_and_print($out_fh, "    is_a: named thing\n");
	replace_pipes_and_print($out_fh, "    slots:\n");
	for my $slot (@{$classes->{$class}{'slots'}}) {
	    replace_pipes_and_print($out_fh, "      - $slot\n");
	}
	replace_pipes_and_print($out_fh, "\n");
    }
    for my $assoc (sort keys %$associations) {
	replace_pipes_and_print($out_fh, "  $assoc:\n    is_a: association\n");
	replace_pipes_and_print($out_fh, "    defining_slots:\n");
	for my $ds (qw(subject predicate object)) {
	    replace_pipes_and_print($out_fh, "      - $ds\n");
	}
	replace_pipes_and_print($out_fh, "    slot_usage:\n");
	replace_pipes_and_print($out_fh, "      subject:\n        range: " . $associations->{$assoc}{'subject_range'} . "\n");
	replace_pipes_and_print($out_fh, "      predicate:\n        subproperty_of: " . $associations->{$assoc}{'predicate_subproperty_of'} . "\n");
	replace_pipes_and_print($out_fh, "      object:\n        range: " . $associations->{$assoc}{'object_range'} . "\n");
	replace_pipes_and_print($out_fh, "\n");
    }
    replace_pipes_and_print($out_fh, "\nslots:\n\n");
    for my $slot (sort keys %$slots) {
	print $slot . "\n";
	replace_pipes_and_print($out_fh, "  $slot:\n");
	replace_pipes_and_print($out_fh, "    is_a: " . $slots->{$slot}{'is_a'} . "\n");
	replace_pipes_and_print($out_fh, "    description: >-\n      " . $slots->{$slot}{'description'} . "\n") if exists $slots->{$slot}{'description'};
	my @domains = keys %{$slots->{$slot}{'domain'}};
	my $domain = @domains > 1 ? 'named thing' : $domains[0];
	replace_pipes_and_print($out_fh, "    domain: $domain\n");
	my @ranges = keys %{$slots->{$slot}{'range'}};
	my $range = @ranges > 1 ? 'named thing' : $ranges[0];
	replace_pipes_and_print($out_fh, "    range: $range\n");
	replace_pipes_and_print($out_fh, "    values_from: " . $slots->{$slot}{'values_from'} . "\n") if exists $slots->{$slot}{'values_from'};
	replace_pipes_and_print($out_fh, "\n");
    }
    $out_fh->print($footer);

    return;
}


sub replace_pipes_and_print{
    my ($out_fh, $string) = @_;

    $string =~ s/\|/_/g;
    $out_fh->print($string);

    return;
}
 
   
sub process_slot_value {
    my ($class, $slot, $value, $classes, $associations, $slots, $hash_slots, $comments, $prefix, $controlled_vocabs) = @_;

    my $slot_key = substr($slot, length($prefix));

    if ($value =~ /#(\w+)$/) {
	for my $hash_slot_value (@{$hash_slots->{$1}}) {
	    ($classes, $associations, $slots) = process_slot_value($class, $slot . '|' . $hash_slot_value->{'slot'},
								   $value . ' ' . $hash_slot_value->{'value'},
								   $classes, $associations, $slots, $hash_slots,
								   $comments, $slot . '|', $controlled_vocabs);
	}
	$value =~ s/\s*#\w+$//;
    }
    
    my $us_slot = $slot =~ s/\|/_/gr;  # Need to replace pipes now to avoid duplicate slots from e.g. has_gene|contains and has_gene_contains
    if ($value =~ /\?(\w+)\s.*\?(\w+)/) {
	my $us1 = $1 =~ s/\|/_/gr;
	my $us2 = $2 =~ s/\|/_/gr;
	$classes->{"${1}_${2}"}{'description'} = "Linking class for $1 and $2 classes";
	push @{$classes->{"$1_$2"}{'slots'}}, "has_$1";
	push @{$classes->{"$1_$2"}{'slots'}}, "has_$2";
	push @{$classes->{$class}{'slots'}}, "has_$1_$2";
	$slots->{"has_${us1}_${us2}"}{'is_a'} = 'related to';
	$slots->{"has_${us1}_${us2}"}{'domain'}{$class}++;
	$slots->{"has_${us1}_${us2}"}{'range'}{"$1_$2"}++;
	$slots->{"has_$us1"}{'domain'}{"$1_$2"}++;
	$slots->{"has_$us2"}{'domain'}{"$1_$2"}++;
	$slots->{"has_$us1"}{'range'}{$1}++;
	$slots->{"has_$us2"}{'range'}{$2}++;
	$slots->{"has_$us1"}{'is_a'} = 'related to';
	$slots->{"has_$us2"}{'is_a'} = 'related to';
	my $assoc_name = join('_', $class, 'has', $slot, $1, $2, 'association');
	$associations->{$assoc_name}{'subject_range'} = $class;
	$associations->{$assoc_name}{'predicate_subproperty_of'} = "has_$1_$2";
	$associations->{$assoc_name}{'object_range'} = "$1_$2";
	$assoc_name = join('_', $1, $2, 'has', $1, 'association');
	$associations->{$assoc_name}{'subject_range'} = "$1_$2";
	$associations->{$assoc_name}{'predicate_subproperty_of'} = "has_$1";
	$associations->{$assoc_name}{'object_range'} = $1;
	$assoc_name = join('_', $1, $2, 'has', $2, 'association');
	$associations->{$assoc_name}{'subject_range'} = "$1_$2";
	$associations->{$assoc_name}{'predicate_subproperty_of'} = "has_$2";
	$associations->{$assoc_name}{'object_range'} = $2;
    }
    elsif ($value =~ /\?(\w+)/) {
	my $assoc_name = join('_', $class, 'has', $slot, $1, 'association');
	$associations->{$assoc_name}{'subject_range'} = $class;
	$associations->{$assoc_name}{'predicate_subproperty_of'} = "has_$1";
	$associations->{$assoc_name}{'object_range'} = $1;
	push @{$classes->{$class}{'slots'}}, "has_$1";
	$slots->{"has_$1"}{'domain'}{$class}++;
	$slots->{"has_$1"}{'range'}{$1}++;
	$slots->{"has_$1"}{'is_a'} = 'related to';
	$slots->{"has_$1"}{'description'} = $comments->{$slot_key} if exists $comments->{$slot_key};
    }
    else {
	$slots->{$us_slot}{'is_a'} = 'node property';
	$slots->{$us_slot}{'domain'}{$class}++;
	$slots->{$us_slot}{'range'}{get_range($value)}++;
	$slots->{$us_slot}{'description'} = $comments->{$slot_key} if exists $comments->{$slot_key};
	$slots->{$us_slot}{'values_from'} = '[' . join(', ', @{$controlled_vocabs->{$class}{$slot_key}}) . ']' 
	    if $value eq 'controlled_vocabulary';
	push @{$classes->{$class}{'slots'}}, $slot;
    }
   	    
    return ($classes, $associations, $slots);
}


sub get_range {
    my $value = shift;

    return 'boolean' if $value eq '';
    return 'float' if $value eq 'Float';
    return 'int' if $value eq 'Int';
    return 'string';
}

sub extract_slot_details {
    my $model_lines = shift;

    my ($class, $class_or_hash, $previous_indent, $possible_cv);
    my @cv;
    my (%classes, %hashes, %comments, %slot_parts, %controlled_vocabs);

    for my $ix (0 .. scalar(@$model_lines) - 1) {
	my $line = $model_lines->[$ix];

	# Capture and remove comments
	my $comment;
	$line =~ s/\/\/$//;
	if ($line =~ /\/\/\s*(.+)$/) {
	    $comment = $1;
	    $line =~ s/\s*\/\/.+$//;
	}

	# Start of new class/hash definition
	if ($line =~ /^([\?#])(\w+)\s(.+)$/) {
	    $class_or_hash = $1 eq '?' ? 'class' : 'hash';
	    $class = $2;
	    undef %slot_parts;
	    $previous_indent = 0;
	    # Replace class name with spaces so line can be treated like others
	    $line = ' ' x ((length $class) + 2);
	    $line .= $3;
	}
	
	my ($leading_spaces, $text) = $line =~ /^(\s+)(\S.+)$/;
	$text =~ s/\s+/ /g;
	
	my $indent = length($leading_spaces);
	my @parts = split(' ', $text);
	                    
	# Check indentation of next line
	my $next_line = $model_lines->[$ix + 1];
	my ($nl_leading_spaces) = $next_line =~ /^(\s+)\S+/;
	my $next_indent = length $nl_leading_spaces;
	
	# If reduction in indentation, remove slot title parts as required    
	if ($indent < $previous_indent) {
	    for my $slot_indent (sort {$b<=>$a} keys %slot_parts) {
		last if $slot_indent < $indent;
		delete $slot_parts{$slot_indent};
	    }
	}
	 
	# Combine multiple words before next level of indentation
	while (length $parts[0] < $next_indent - $indent) {
	    my $first_part = shift @parts;
	    $slot_parts{$indent} = $first_part;
	    $indent += (length $first_part) + 1;
	    
	}
	$slot_parts{$indent} = shift @parts; 
  
	# If next level increases indent then incorporate first 2 parts into slot tag
	if ($next_indent > $indent) {
	    $slot_parts{$next_indent} = shift @parts;
	}
	
	my @slot_title_parts;
	for my $slot_indent(sort {$a<=>$b} keys %slot_parts) {
	    push @slot_title_parts, $slot_parts{$slot_indent} if defined $slot_parts{$slot_indent};
	}

	# Join with pipes that are later switch to underscores to differentiate between tags
	my $slot_title = join('|', @slot_title_parts);
	
	if ($slot_title =~ /^(.+)_unique\|([^\|]+)$/) {
	    $controlled_vocabs{$class}{$1}{$2} = @parts;
	}
   
	$slot_title =~ s/_unique//g;

	for my $ix (0 .. @parts - 1) {
	    $parts[$ix] =~ s/_unique//g;
	}
	
	if ($class_or_hash eq 'class') {
	    push @{$classes{$class}}, {slot => $slot_title, value => join(' ', @parts)};
	}
	else {
	    push @{$hashes{$class}}, {slot => $slot_title, value => join(' ', @parts)};
	}
	$comments{$slot_title} = $comment if $comment;
 	
	$previous_indent = $indent;
    }

    return (\%classes, \%hashes, \%comments, \%controlled_vocabs);
}


sub merge_controlled_vocabs {
    my ($slot_values, $controlled_vocabs) = @_;

    my %cvs_added;
    my %merged_slot_values;
    for my $class (keys %$slot_values) {
	for my $slot_value (@{$slot_values->{$class}}) {
	    my ($slot_stem) = $slot_value->{'slot'} =~ /^(.+)\|[^\|]+$/;
	    if (exists $controlled_vocabs->{$class}{$slot_stem}) {
		unless (exists $cvs_added{$class}{$slot_stem}) {
		    push @{$merged_slot_values{$class}}, {slot => $slot_stem, 
							  value => 'controlled_vocabulary'};
		    $cvs_added{$class}{$slot_stem} = 1;
		}
	    }
	    else {
		push @{$merged_slot_values{$class}}, $slot_value;
	    }
	}
    }
    
    return \%merged_slot_values;
}


sub filter_controlled_vocabs {
    my $possible_cvs = shift;

    my %valid_cvs;
    for my $class (keys %$possible_cvs) {
	for my $group (keys %{$possible_cvs->{$class}}) {
	    my $is_valid = 1;
	    for my $term (keys %{$possible_cvs->{$class}{$group}}) {
		if ($possible_cvs->{$class}{$group}{$term} > 0) {
		    $is_valid = 0;
		    last;
		}
	    }
	    push @{$valid_cvs{$class}{$group}}, keys %{$possible_cvs->{$class}{$group}} if $is_valid;
	}
    }

    return \%valid_cvs;
}

    
sub get_and_parse_model_lines {
    my @model_lines = file($WB_MODEL_SPEC)->slurp(chomp => 1);
    my @parsed_lines;
    for my $line (@model_lines) {
	$line = expand($line);
	next if $line =~ /^\s*\/\//;
	next if $line =~ /^\s*$/;
	
	# Remove evidence hashes
	$line =~ s/Evidence\s#Evidence/Evidence Text/;
	$line =~ s/\s#Evidence//g;

	# Convert text classes to plain text
	$line =~ s/\?Text/Text/g;
	$line =~ s/\?LongText/Text/g;
	$line =~ s/\?Keyword/Text/g;
	$line =~ s/\?WWW_server/Text/g;
	$line =~ s/\?DNA/Text/g;
	$line =~ s/\?Peptide/Text/g;

	# Remove XREFs
	$line =~ s/\sXREF\s\w+//g;

	$line =~ s/ UNIQUE/_UNIQUE/g;
	push @parsed_lines, lc $line;
    }

    return \@parsed_lines;
}



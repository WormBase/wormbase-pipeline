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
my ($slot_details, $comments, $controlled_vocabs, $unique_slots) = extract_slot_details($model_lines);

$controlled_vocabs = filter_controlled_vocabs($controlled_vocabs);
$slot_details = merge_controlled_vocabs($slot_details, $controlled_vocabs);

my ($classes, $slots);
for my $class (keys %$slot_details) {
    $classes->{$class}{'description'} = "$class class";
    for my $slot_value (@{$slot_details->{$class}}) {
	($classes, $slots) = process_slot_value($class, $slot_value->{'slot'},
						$slot_value->{'value'}, $classes,
						$slots, $comments, $controlled_vocabs,
						$unique_slots);		       
    }
}

($classes, $slots) = avoid_duplicate_class_and_slot_names($classes, $slots);
print_yaml($out_file, $classes, $slots);


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
	next unless exists $classes->{$class}{'non_multivalued_slots'};
	for my $ix (0 .. scalar @{$classes->{$class}{'non_multivalued_slots'}} - 1) {
	    if (exists $changed_slots{$classes->{$class}{'non_multivalued_slots'}[$ix]}) {
		$classes->{$class}{'non_multivalued_slots'}[$ix] .= '_slot';
	    }
	}
    }
    
    return ($classes, \%final_slots);
}


sub extract_slot_details {
    my $model_lines = shift;

    my ($class, $previous_indent, $possible_cv);
    my @cv;
    my (%slots, %comments, %slot_parts, %controlled_vocabs, %unique_slots);

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
	    $class = $2;
	    undef %slot_parts;
	    $previous_indent = 0;
	    # Replace class name with spaces so line can be treated like others
	    $line = ' ' x ((length $class) + 2);
	    $line .= $3;
	    $class =~ s/_unique//;
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
   
	my $has_unique_tag = $slot_title =~ /_unique/ ? 1 : 0;
	if ($has_unique_tag) {
	    $slot_title =~ s/_unique//g;
	    $unique_slots{$class}{$slot_title} = 1;
	}

	for my $ix (0 .. @parts - 1) {
	    $parts[$ix] =~ s/_unique//g;
	}
	
	push @{$slots{$class}}, {slot => $slot_title, value => join(' ', @parts)};
	$comments{$slot_title} = $comment if $comment;
 	
	$previous_indent = $indent;
    }

    return (\%slots, \%comments, \%controlled_vocabs, \%unique_slots);
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


sub get_range {
    my $value = shift;

    return 'boolean' if $value eq '';
    return 'float' if $value eq 'Float';
    return 'int' if $value eq 'Int';
    return 'string';
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


sub print_yaml {
    my ($out_file, $classes, $slots) = @_;

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
	my %printed_class_slots;
	for my $slot (@{$classes->{$class}{'slots'}}) {
	    next if exists $printed_class_slots{$slot};
	    replace_pipes_and_print($out_fh, "      - $slot\n");
	    $printed_class_slots{$slot} = 1;
	}
	my %printed_class_nmv_slots;
	if (exists $classes->{$class}{'non_multivalued_slots'}) {
	    replace_pipes_and_print($out_fh, "    slot_usage:\n");
	    for my $nmv_slot (@{$classes->{$class}{'non_multivalued_slots'}}) {
		next if exists $printed_class_nmv_slots{$nmv_slot};
		replace_pipes_and_print($out_fh, "      ${nmv_slot}:\n");
		replace_pipes_and_print($out_fh, "        multivalued: false\n");
		$printed_class_nmv_slots{$nmv_slot} = 1;
	    }
	}
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


sub process_slot_value {
    my ($class, $slot, $value, $classes, $slots, $comments, $controlled_vocabs, $unique_slots) = @_;

    # Need to replace pipes now to avoid duplicate slots (e.g. has_gene|contains and has_gene_contains)
    my $us_slot = $slot =~ s/\|/_/gr;

    # Deal with rows that link to two different object types

    # Check number of classes and hashes
    my $nr_classes = $value =~ tr/\?#//;
    die "Code cannot currently deal with rows linking to > 3 classes or hashes\n$value\n" if $nr_classes > 3;
    if ($value =~ /\?(\w+)\s.*\?(\w+)\s.*[\?#](\w+)/) {
	$classes->{"${us_slot}_${1}_${2}_${3}"}{'description'} = "Linking class for $1, $2, and $3 classes";
	push @{$classes->{"${us_slot}_${1}_${2}_${3}"}{'slots'}}, "has_${1}";
	push @{$classes->{"${us_slot}_${1}_${2}_${3}"}{'slots'}}, "has_${2}";
	push @{$classes->{"${us_slot}_${1}_${2}_${3}"}{'slots'}}, "has_${3}";
	$slots->{"${us_slot}_${1}_${2}_${3}"}{'is_a'} = 'related to';
	$slots->{"${us_slot}_${1}_${2}_${3}"}{'domain'}{$class}++;
	$slots->{"${us_slot}_${1}_${2}_${3}"}{'range'}{"${us_slot}_${1}_${2}_${3}"}++;
	$slots->{"has_${1}"}{'domain'}{"${us_slot}_${1}_${2}_${3}"}++;
	$slots->{"has_${2}"}{'domain'}{"${us_slot}_${1}_${2}_${3}"}++;
	$slots->{"has_${3}"}{'domain'}{"${us_slot}_${1}_${2}_${3}"}++;
	$slots->{"has_${1}"}{'range'}{$1}++;
	$slots->{"has_${2}"}{'range'}{$2}++;
	$slots->{"has_${3}"}{'range'}{$3}++;
	$slots->{"has_${1}"}{'is_a'} = 'related to';
	$slots->{"has_${2}"}{'is_a'} = 'related to';
	$slots->{"has_${3}"}{'is_a'} = 'related to';
    }
    elsif ($value =~ /\?(\w+)\s.*[\?#](\w+)/) {
	$classes->{"${us_slot}_${1}_${2}"}{'description'} = "Linking class for $1 and $2 classes";
	push @{$classes->{"${us_slot}_${1}_${2}"}{'slots'}}, "has_${1}";
	push @{$classes->{"${us_slot}_${1}_${2}"}{'slots'}}, "has_${2}";
	push @{$classes->{$class}{'slots'}}, "${us_slot}_${1}_${2}";
	$slots->{"${us_slot}_${1}_${2}"}{'is_a'} = 'related to';
	$slots->{"${us_slot}_${1}_${2}"}{'domain'}{$class}++;
	$slots->{"${us_slot}_${1}_${2}"}{'range'}{"${us_slot}_${1}_${2}"}++;
	$slots->{"has_${1}"}{'domain'}{"${us_slot}_${1}_${2}"}++;
	$slots->{"has_${2}"}{'domain'}{"${us_slot}_${1}_${2}"}++;
	$slots->{"has_${1}"}{'range'}{$1}++;
	$slots->{"has_${2}"}{'range'}{$2}++;
	$slots->{"has_${1}"}{'is_a'} = 'related to';
	$slots->{"has_${2}"}{'is_a'} = 'related to';	
    }
    elsif ($value =~ /\?(\w+)/) {
	push @{$classes->{$class}{'slots'}}, "${us_slot}_${1}";
	$slots->{"${us_slot}_${1}"}{'domain'}{$class}++;
	$slots->{"${us_slot}_${1}"}{'range'}{$1}++;
	$slots->{"${us_slot}_${1}"}{'is_a'} = 'related to';
	$slots->{"${us_slot}_${1}"}{'description'} = $comments->{$slot} if exists $comments->{$slot};
    }
    else {
	$slots->{$us_slot}{'is_a'} = 'node property';
	$slots->{$us_slot}{'domain'}{$class}++;
	$slots->{$us_slot}{'range'}{get_range($value)}++;
	$slots->{$us_slot}{'description'} = $comments->{$slot} if exists $comments->{$slot};
	if ($value eq 'controlled_vocabulary') {
	    $slots->{$us_slot}{'values_from'} = '[' . join(', ', @{$controlled_vocabs->{$class}{$slot}}) . ']';
	}
	elsif (exists $unique_slots->{$class}{$slot}) {
	    push @{$classes->{$class}{'non_multivalued_slots'}}, $us_slot;
	}
	push @{$classes->{$class}{'slots'}}, $slot;
    }
   	    
    return ($classes, $slots);
}


sub replace_pipes_and_print{
    my ($out_fh, $string) = @_;

    $string =~ s/\|/_/g;
    $out_fh->print($string);

    return;
}
 
   



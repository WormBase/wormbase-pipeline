#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Text::Tabs;
use Getopt::Long;
use Const::Fast;

const my @RESTRICTED_NAMES => qw(type class slot domain range inherited symmetric description name version);

GetOptions(
    "model|m=s" => \my $model_file, # annotated models file models.wrm.annot
    "out|o=s"   => \my $out_file
    );

my $model_lines = read_model($model_file);
my $expanded_enum_lines = expand_enum_lines($model_lines);
my ($slots, $comments, $enums, $descriptions, $class_name_map) = extract_slot_details($expanded_enum_lines);
$class_name_map = avoid_duplicate_class_names($class_name_map);
($slots, $enums) = avoid_class_slot_name_overlap($slots, $class_name_map, $enums);
my $combined_slots = combine_slots_with_shared_name($slots);
print_yaml($slots, $comments, $enums, $descriptions, $class_name_map, $combined_slots, $out_file);


sub avoid_class_slot_name_overlap {
    my ($slots, $class_name_map, $enums) = @_;

    my %class_names = map {biolinkify_name($_) => 1} values %$class_name_map;
    my %restricted = map {$_ => 1} @RESTRICTED_NAMES;
    my (%renamed_slots, %renamed_enums);
    
    for my $class (keys %$slots) {
	for my $slot (keys %{$slots->{$class}}) {
	    my $renamed_slot = (exists $class_names{biolinkify_name($slot)} or 
				exists $restricted{biolinkify_name($slot)}) ? $slot . '_slot' : $slot;
	    $renamed_slots{$class}{$renamed_slot} = $slots->{$class}{$slot};
	}
	for my $enum (keys %$enums) {
	    if ($enums->{$enum}{'attributes'}) {
		$renamed_enums{$enum}{'values'} = $enums->{$enum}{'values'};
		for my $attr (keys %{$enums->{$enum}{'attributes'}}) {
		    my $renamed_attr =  (exists $class_names{biolinkify_name($attr)} or
			exists $restricted{biolinkify_name($attr)}) ? $attr . '_slot' : $attr;
		    $renamed_enums{$enum}{'attributes'}{$renamed_attr} = $enums->{$enum}{'attributes'}{$attr};
		}
	    }
	    else {
		$renamed_enums{$enum} = $enums->{$enum};
	    }
	}	    
    }

    return (\%renamed_slots, \%renamed_enums);
}


sub avoid_duplicate_class_names {
    # Class names can end up as duplicated once biolinkified without this
    my $class_name_map = shift;

    my %biolinkified_names = map {$_ => biolinkify_name($class_name_map->{$_})} keys %$class_name_map;
    my %reversed;
    for my $original_name (keys %biolinkified_names) {
	$reversed{$biolinkified_names{$original_name}}++;
	my $nr = $reversed{$biolinkified_names{$original_name}};
	$class_name_map->{$original_name} = $class_name_map->{$original_name} . "_v$nr" if $nr > 1;
    }

    return $class_name_map;
}


sub biolinkify_name {
    my $name = shift;

    $name =~ s/[_\-]/ /g;
    $name =~ s/^[?#]//;

    return lc $name;
}


sub combine_slots_with_shared_name {
    my $slots = shift;

    my %combined_slots;
    for my $class (keys %$slots) {
	for my $slot (keys %{$slots->{$class}}) {
	    $combined_slots{biolinkify_name($slot)}{'domain'}{$class} = 1;
	    $combined_slots{biolinkify_name($slot)}{'range'}{$slots->{$class}{$slot}} = 1;
	}
    }

    return \%combined_slots;
}


sub expand_enum_lines {
    my $lines = shift;

    my @expanded_enum_lines;
    my ($is_enum, $enum_string);
    for my $ix (0 .. scalar @$lines - 1) {
	if ($lines->[$ix] =~ /^[^\/]+ENUM/) {
	    $is_enum = 1;
	    ($enum_string) = $lines->[$ix] =~ /^(.+\sENUM(\s*\^\w+)?\s)/;
	    push @expanded_enum_lines, $lines->[$ix];
	    next;
	}
	
	my ($indent) = $lines->[$ix] =~ /^(\s*)\S/;
	if ($is_enum and length $indent >= length $enum_string) {
	    push @expanded_enum_lines, $enum_string . substr($lines->[$ix], length $enum_string);
	}
	else {
	    $is_enum = 0;
	    push @expanded_enum_lines, $lines->[$ix];
	}
    }

    return \@expanded_enum_lines;
}


sub extract_slot_details {
    my $lines = shift;

    my ($class, $previous_indent, $branch);
    my (%slots, %comments, %enums, %class_name_map, %descriptions);

    $class_name_map{'Text'} = 'string';
    $class_name_map{'Float'} = 'float';
    $class_name_map{'Int'} = 'integer';
    $class_name_map{'boolean'} = 'boolean';
    $class_name_map{'named thing'} = 'named thing';
    $class_name_map{'related to'} = 'related to';
    $class_name_map{'node property'} = 'node property';
    $class_name_map{'DateType'} = 'datetime';

    for my $ix (0 .. scalar(@$lines) - 1) {
	my $line = $lines->[$ix];

	# Capture and remove comments
	my $comment;
	$line =~ s/\/\/\s*$//;
	if ($line =~ /\/\/\s*(.+)$/) {
	    $comment = $1;
	    $comment =~ s/"/'/g;
	    $line =~ s/\s*\/\/.+$//;
	}	
		
        # Start of new class/hash definition
	if ($line =~ /^[\?#](\w+)\s(.+)$/) {
	    $class = $1;
	    $previous_indent = 0;

	    # Replace class name with spaces so line can be treated like others
	    $line = ' ' x ((length $class) + 2);
	    $line .= $2;

	    my $old_class;
	    ($old_class, $class) = process_name($class);
	    $class_name_map{$old_class} = $class;
	    $descriptions{$class} = ucfirst(biolinkify_name($class)) . ' class';
	}

	my ($leading_spaces) = $line =~ /^(\s*)\S.+$/;
	my $indent = length $leading_spaces;

	my $next_indent = next_line_indent($lines, $ix);
	
	# If start of a branch
	if ($next_indent > $indent) {
	    $line = substr($line, $next_indent);
	    $line = (' ' x $next_indent) . $line;
	}

	#Remove leading and trailing spaces and ensure all spaces are single
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	$line =~ s/\s+/ /g;

	# Process ENUM lines
	if ($line =~ /ENUM/) {
	    $line =~ s/ ENUM / /;
	    $line =~ s/\s\^/^/g;
	    my @parts = split(' ', $line);

	    my $enum_name = shift @parts;
	    if ($enum_name =~ /^.+\^(.+)$/) {
		$enum_name = $1;
	    }
	    $enum_name = $class . '_' . $enum_name;

	    my $enum_value = pop @parts;
	    my @enum_attributes;
	    while ($enum_value =~ /^[#?]/ or $enum_value eq 'Text' or
		   $enum_value eq 'Int' or $enum_value eq 'Float' or
		   $enum_value eq 'DateType'
		) {
		$enums{$enum_name . '_enum'}{'attributes'}{$enum_value}++;
		$enum_value = pop @parts;
	    }

	    if ($enum_value =~ /^.+\^(.+)$/) {
		$enum_value = $1;
	    }
	    push @{$enums{$enum_name . '_enum'}{'values'}}, $enum_value;
	    
	    $slots{$class}{$enum_name} = $enum_name . '_enum';
	    $comments{$class}{$enum_name} = $comment if $comment;
	    $class_name_map{$enum_name . '_enum'} = $enum_name . '_enum';
	    next;
	}

	# Process non-ENUM lines
	$line =~ s/\s\^/^/g;
	
	my @parts = split(' ', $line);
	if (scalar @parts == 1) {
	    $slots{$class}{$parts[0]} = 'boolean';
	    $comments{$class}{$parts[0]} = $comment if $comment;
	    next;
	}

	my @labels;
	my %attributes;
	for my $part (@parts) {
	    if ($part =~ /^Text$/ or $part =~ /^Float$/ or $part =~ /^Int$/ or 
		$part =~ /^Text\^/ or $part =~ /^Float\^/ or $part =~ /^Int\^/ or
		$part =~ /^[\?#]/ or $part =~ /^DateType\^/ or $part =~ /^DateType$/) {
		$part =~ s/^[?#]//;
		my ($range, $attribute) = process_name($part);
		$attributes{$attribute} = $range;
	    }
	    else {
		my ($old_label, $label) = process_name($part);
		push @labels, $label;
	    }
	}
	
	# Deal with simple cases with single attribute per line
	if (scalar keys %attributes == 1) {
	    for my $attribute (keys %attributes) {
		my $range = $attributes{$attribute};
		$range =~ s/^[?#]//;
		$slots{$class}{join('_', @labels)} = $range;
		$comments{$class}{join('_', @labels)} = $comment if $comment;
	    }
	    next;
	}
	
	# Create new classes for cases with multiple attributes
	my $new_class;
	if (@labels) {
	    $new_class = join('_', $class, @labels);
	    $descriptions{$new_class} = 'Linking class for ' . biolinkify_name(join('_', @labels)) .
		' properties of ' . biolinkify_name($class) . ' class';
	}
	else {
	    my @attr_names = map {biolinkify_name($_)} keys %attributes;
	    $new_class = join('_', $class, @attr_names);
	    $descriptions{$new_class} = 'Linking class for ' . biolinkify_name(join('_', @attr_names)) .
		' properties of ' . biolinkify_name($class) . ' class';
	}
	$slots{$class}{$new_class} = $new_class;
	for my $attribute (keys %attributes) {
	    my $range = $attributes{$attribute};
	    $range =~ s/^[?#]//;
	    $slots{$new_class}{$attribute} = $range;
	}
	$class_name_map{$new_class} = $new_class;
    }

    return (\%slots, \%comments, \%enums, \%descriptions, \%class_name_map);
}


sub next_line_indent {
    my ($lines, $ix) = @_;

    my $next_indent = 0;
    if ($ix + 1 != scalar @$lines) {
	my ($next_line_leading_spaces) = $lines->[$ix + 1] =~ /^(\s*)\S/;
	$next_indent = length $next_line_leading_spaces;
    }

    return $next_indent;
}


sub print_yaml {
    my ($slots, $comments, $enums, $descriptions, $cnmap, $combined_slots, $out_file) = @_;

    my ($header, $footer) = yaml_header_and_footer();

    my $out_fh = file($out_file)->openw;
    $out_fh->print($header);
    
    # Print classes
    for my $class (keys %$slots) {
	$out_fh->print('  ' . biolinkify_name($cnmap->{$class}) . ":\n" .
		       '    description: "' . $descriptions->{$class} . "\"\n" .
		       "    is_a: named thing\n" .
		       "    slots:\n"
	    );
	for my $class_slot (keys %{$slots->{$class}}) {
	    $out_fh->print('      - ' . biolinkify_name($class_slot) . "\n");
	}
	my $slot_usage_printed = 0;
	for my $class_slot (keys %{$slots->{$class}}) {
	    my $slot_printed = 0;
	    if (scalar keys %{$combined_slots->{biolinkify_name($class_slot)}{'range'}} > 1) {
		if (!$slot_usage_printed) {$out_fh->print("    slot_usage:\n"); $slot_usage_printed = 1;}
		if (!$slot_printed) {$out_fh->print('      ' . biolinkify_name($class_slot) . ":\n"); $slot_printed = 1;}
		$out_fh->print('        range: ' . biolinkify_name($cnmap->{$slots->{$class}{$class_slot}}) . "\n");
	    }
	    if ($comments->{$class}{$class_slot} and scalar keys %{$combined_slots->{biolinkify_name($class_slot)}{'domain'}} > 1) {
		if (!$slot_usage_printed) {$out_fh->print("    slot_usage:\n"); $slot_usage_printed = 1;}
		if (!$slot_printed) {$out_fh->print('      ' . biolinkify_name($class_slot) . ":\n"); $slot_printed = 1;}
		$out_fh->print('        description: "' . $comments->{$class}{$class_slot} . "\"\n");
	    }
	}
	$out_fh->print("\n");
    }

    # Print linking classes for enums with additional attributes
    for my $enum (keys %$enums) {
	next unless exists $enums->{$enum}{'attributes'};
	my $enum_name = biolinkify_name($enum); 
	$out_fh->print("  $enum_name:\n");
	$out_fh->print("    description: Linking class for $enum_name enum and associated attributes\n");
	$out_fh->print("    is_a: named thing\n" .
		       "    slots:\n");
	$out_fh->print('      - ' . substr(biolinkify_name($enum), 0, -4) . "values slot\n");
	for my $attribute (keys %{$enums->{$enum}{'attributes'}}) {
	    $out_fh->print('      - ' . biolinkify_name($attribute) . "\n");
	}
	$out_fh->print("\n");
    }
    $out_fh->print("\n");

    # Print enums
    $out_fh->print("enums:\n\n");
    for my $enum (keys %$enums) {
	my $enum_name = exists $enums->{$enum}{'attributes'} ? $enum . '_values' : $enum;
	$enum_name = biolinkify_name($enum_name);
	$out_fh->print("  $enum_name:\n");
	$out_fh->print("    permissible_values:\n");
	for my $value (@{$enums->{$enum}{'values'}}) {
	    $out_fh->print("      $value:\n");
	    $out_fh->print("        text: $value\n");
	}
	$out_fh->print("\n");
    }
    $out_fh->print("\n");

    # Print slots
    $out_fh->print("slots:\n\n");
    for my $slot (keys %$combined_slots) {
	my ($type, $domain, $range) = slot_details($combined_slots, $slot, $enums);
	$domain = $cnmap->{$domain};
	$range = $cnmap->{$range};
	$out_fh->print('  ' . biolinkify_name($slot) . ":\n");
	$out_fh->print("    is_a: $type\n");
	$out_fh->print('    domain: ' . biolinkify_name($domain) . "\n");
	$out_fh->print('    range: ' . biolinkify_name($range) . "\n");
	$out_fh->print("\n");
    }

    # Print slots for linking classes create for enums with additional attributes
    for my $enum (keys %$enums) {
	next unless exists $enums->{$enum}{'attributes'};
	$out_fh->print('  ' . substr(biolinkify_name($enum), 0, -4) . 'values slot' . ":\n" .
		       "    is_a: node property\n" .
		       '    domain: ' . biolinkify_name($enum) . "\n" .
		       '    range: ' . biolinkify_name($enum . '_values') . "\n" .
		       "\n" );
    }

    $out_fh->print($footer);

    return;
}

	
sub process_name {
    my $name = shift;

    return ($name, $name) unless $name =~ /\^/;

    if ($name =~ /^(\S+)\s*\^(.+)$/) {
	return ($1, $2);
    }

    die "Cannot process alternate name: $name\n";
}


sub read_model {
    my @model_lines = file($model_file)->slurp(chomp => 1);
    my @parsed_lines;
    for my $line (@model_lines) {
	$line = expand($line);

	# Skip commented and empty lines
	next if $line =~ /^\s*\/\//;
	next if $line =~ /^\s*$/;

	# Convert models not in schema to plain text
	$line =~ s/\?Text/Text/g;
	$line =~ s/\?LongText/Text/g;
	$line =~ s/\?Keyword/Text/g;
	$line =~ s/\?WWW_server/Text/g;
	$line =~ s/\?DNA/Text/g;
	$line =~ s/\?Peptide/Text/g;
	$line =~ s/#Colour/Text/g;
	$line =~ s/Free_dup/Text/g;
	$line =~ s/Phenotype2GO/Text/g;
	$line =~ s/\sNull\s/ Null_value /g;

	# Remove XREFs and #Ordered
	$line =~ s/\s(NO|IN|OUT)XREF\s(\w+)?/ /g;
	$line =~ s/\s#Ordered//g;

	# Remove UNIQUE labels preserving indentation
	$line =~ s/UNIQUE/      /g;
	
	push @parsed_lines, $line;
    }

    return \@parsed_lines;
}


sub slot_details {
    my ($combined_slots, $slot, $enums) = @_;

    my ($slot_domain, $slot_range);
    my $type = 'related to';
    for my $domain (keys %{$combined_slots->{$slot}{'domain'}}) {
	$slot_domain = defined $slot_domain ? 'named thing' : $domain;
    }
    for my $range (keys %{$combined_slots->{$slot}{'range'}}) {
	$slot_range = defined $slot_range ? 'named thing' : $range;
	$type = 'node property' if $range eq 'boolean' or
	    $range eq 'integer' or $range eq 'float' or
	    $range eq 'string' or exists $enums->{$range};
    }

    return ($type, $slot_domain, $slot_range);
}


sub yaml_header_and_footer {
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


YAML_FOOTER

    return ($header, $footer);
}

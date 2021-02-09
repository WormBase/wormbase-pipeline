#!/usr/bin/env perl
use strict;
use Const::Fast;
use Path::Class;
use Text::Tabs;
use Getopt::Long;

const my $WB_MODEL_SPEC => $ENV{WORMPUB} . "/wormbase-pipeline/wspec/models.wrm";

GetOptions(
    "out|o=s" => \my $out_file
    );

my $model_lines = get_and_parse_model_lines();
my ($class_slots, $hash_slots, $comments) = extract_slot_details($model_lines);

my ($classes, $associations, $slots);
for my $class (keys %$class_slots) {
    $classes->{$class}{'description'} = "$class class";
    for my $slot_value (@{$class_slots->{$class}}) {
	($classes, $associations, $slots) = process_slot_value($class, $slot_value->{'slot'},
							       $slot_value->{'value'}, $classes,
							       $associations, $slots, $hash_slots,
							       $comments, '');		       
    }
}

print_yaml($out_file);


sub print_yaml {
    my $out_file = shift;

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

YAML_HEADER

    my $out_fh = file($out_file)->openw;
    $out_fh->print("$header\n\nclasses:\n\n");
    for my $class (keys %$classes) {
	my %slots_printed = ();
	$out_fh->print("  $class:\n");
        $out_fh->print("    description: " . $classes->{$class}{'description'} . "\n");
	$out_fh->print("    is_a: named_thing\n");
	$out_fh->print("    slots:\n");
	for my $slot (@{$classes->{$class}{'slots'}}) {
	    $out_fh->print("      - $slot\n") unless exists $slots_printed{$slot};
	    $slots_printed{$slot}++;
	}
	$out_fh->print("\n");
    }
    for my $assoc (keys %$associations) {
	$out_fh->print("  $assoc\n    is_a: association\n");
	$out_fh->print("    defining_slots:\n");
	for my $ds (qw(subject predicate object)) {
	    $out_fh->print("      - $ds\n");
	}
	$out_fh->print("    slot_usage:\n");
	$out_fh->print("      subject:\n        range: " . $associations->{$assoc}{'subject_range'} . "\n");
	$out_fh->print("      predicate:\n        subproperty_of: " . $associations->{$assoc}{'predicate_subpropert_of'} . "\n");
	$out_fh->print("      object:\n        range: " . $associations->{$assoc}{'object_range'} . "\n");
	$out_fh->print("\n");
    }
    $out_fh->print("\nslots:\n\n");
    for my $slot (keys %$slots) {
	$out_fh->print("  $slot:\n");
	$out_fh->print("    is_a: " . $slots->{$slot}{'is_a'} . "\n");
	$out_fh->print("    description: >-\n      " . $slots->{$slot}{'description'}) if exists $slots->{$slot}{'description'};
	my @domains = keys %{$slots->{$slot}{'domain'}};
	my $domain = @domains > 1 ? 'named_thing' : $domains[0];
	$out_fh->print("    domain: $domain\n");
	$out_fh->print("    range: " . $slots->{$slot}{'range'} . "\n");
    }

    return;
}
 
   
sub process_slot_value {
    my ($class, $slot, $value, $classes, $associations, $slots, $hash_slots, $comments, $prefix) = @_;

    my $comment_key = substr($slot, length($prefix));

    if ($value =~ /#(\w+)$/) {
	$value =~ s/\s*#\w+$//;
	for my $hash_slot_value (@{$hash_slots->{$1}}) {
	    ($classes, $associations, $slots) = process_slot_value($slot . '_' . $hash_slot_value->{'slot'},
								   $value . ' ' . $hash_slot_value->{'value'},
								   $classes, $associations, $slots, $hash_slots,
								   $comments, $slot . '_');
	}
    }
    elsif ($value =~ /\?(\w+)\s.*\?(\w+)/) {
	$classes->{"${1}_${2}"}{'description'} = "Linking class for $1 and $2 classes";
	push @{$classes->{"$1_$2"}{'slots'}}, "has_$1";
	push @{$classes->{"$1_$2"}{'slots'}}, "has_$2";
	push @{$classes->{$class}{'slots'}}, "has_$1_$2";
	$slots->{"has_$1_$2"}{'domain'}{$class}++;
	$slots->{"has_$1"}{'domain'}{"$1_$2"}++;
	$slots->{"has_$2"}{'domain'}{"$1_$2"}++;
	$slots->{"has_$1"}{'range'} = $1;
	$slots->{"has_$2"}{'range'} = $2;
	$slots->{"has_$1"}{'is_a'} = 'related_to';
	$slots->{"has_$2"}{'is_a'} = 'related_to';
	my $assoc_name = join('_', $class, 'has', $slot, $1, $2, 'association');
	$associations->{$assoc_name}{'subject_range'} = $class;
	$associations->{$assoc_name}{'predicate_subproperty_of'} = $slot;
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
	$associations->{$assoc_name}{'predicate_subproperty_of'} = "has_$slot";
	$associations->{$assoc_name}{'object_range'} = $1;
	$slots->{"has_$slot"}{'domain'}{$class}++;
	$slots->{"has_$slot"}{'range'} = $1;
	$slots->{"has_$slot"}{'is_a'} = 'related_to';
	$slots->{"has_$slot"}{'description'} = $comments->{$comment_key} if exists $comments->{$comment_key};
    }
    else{
	$slots->{$slot}{'is_a'} = 'node_property';
	$slots->{$slot}{'domain'}{$class}++;
	$slots->{$slot}{'range'} = $value eq '' ? 'boolean' : 'string';
	$slots->{$slot}{'description'} = $comments->{$comment_key} if exists $comments->{$comment_key};
	push @{$classes->{$class}{'slots'}}, $slot;
    }
   	    
    return ($classes, $associations, $slots);
}


sub extract_slot_details {
    my $model_lines = shift;

    my ($class, $class_or_hash, $indent, $previous_indent, $prefix, $previous_text);
    my (%classes, %hashes, %comments);
    my (@slot_title_parts, @indent_lengths);

    for my $ix (0 .. scalar(@$model_lines) - 1) {
	my $line = $model_lines->[$ix];

	# Capture and remove comments
	my $comment;
	$line =~ s/\/\/$//;
	if ($line =~ /\/\/\s*(.+)$/) {
	    $comment = $1;
	    $line =~ s/\s*\/\/.+$//;
	}

	if ($line =~ /^([\?#])(\w+)\s(.+)$/) {
	    $class_or_hash = $1 eq '?' ? 'class' : 'hash';
	    $class = $2;
	    $previous_indent = (length $class) + 2;
	    $prefix = '';
	    @indent_lengths = ($previous_indent);
	    @slot_title_parts = ('');
	    
	    # Replace class name with spaces so line can be treated like others
	    $line = ' ' x ((length $class) + 2);
	    $line .= $3;
	}
	
	my ($leading_spaces, $text) = $line =~ /^(\s+)(\S.+)$/;
	$text =~ s/\s+/ /g;
	
	$indent = length($leading_spaces);
	my @parts = split(' ', $text);
	
	pop @slot_title_parts;
	if ($indent == $previous_indent) {
	    push @slot_title_parts, $parts[0];
	    my $next_line = $model_lines->[$ix + 1];
	    my ($nl_leading_spaces) = $next_line =~ /^(\s+)\S+/;
	    if (length $nl_leading_spaces > $indent) {
		shift @parts;
		push @slot_title_parts, $parts[0];
	    }
	}
	elsif ($indent > $previous_indent) {
	    push @slot_title_parts, $parts[0];
	    push @indent_lengths, $indent;
	}
	else {
	    my $x = -1;
	    while ($indent_lengths[$x] > $indent) {
		pop @indent_lengths;
		pop @slot_title_parts;
		$x--;
	    }
	    push @slot_title_parts, $parts[0];
	    
	    my $next_line = $model_lines->[$ix + 1];
	    my ($nl_leading_spaces) = $next_line =~ /^(\s+)\S+/;
	    if (length $nl_leading_spaces > $indent) {
		shift @parts;
		push @slot_title_parts, $parts[0];
	    }
	}
	
	shift @parts;
	
	my $slot_title = join('_', @slot_title_parts);
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

    return (\%classes, \%hashes, \%comments);
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

	# Convert Text class to plain text
	$line =~ s/\?Text/Text/g;

	# Remove XREFs
	$line =~ s/\sXREF\s\w+//g;

	$line =~ s/ UNIQUE/_UNIQUE/g;
	push @parsed_lines, lc $line;
    }

    return \@parsed_lines;
}



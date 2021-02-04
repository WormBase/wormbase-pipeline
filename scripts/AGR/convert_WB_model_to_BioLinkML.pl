#!/usr/bin/env perl
use strict;
use Const::Fast;
use Path::Class;

const my $WB_MODEL_SPEC => $ENV{WORMPUB} . "/wormbase-pipeline/wspec/models.wrm";

my $model_lines = get_and_parse_model_lines();

my ($class, $class_or_hash, $indent, $previous_indent, $prefix, $previous_text);
my (@slot_titles, @indent_lengths);

for my $ix (0 .. scalar(@$model_lines) - 1) {
    my $line = $model_lines->[$ix];

    # Capture and remove comments
    my $comment;
    $line =~ s/\/\/$//;
    if ($line =~ /\/\/\s*(.+)$/) {
	$comment = $1;
	$line =~ s/\s*\/\/.+$//;
    }

    if ($line =~ /^\?(\w+)\s(.+)$/) {
	$class = $1;
	$class_or_hash = $class;
	$previous_indent = (length $class) + 2;
	$prefix = '';
	@indent_lengths = ($previous_indent);
	@slot_titles;

	# Replace class name with spaces so line can be treated like others
	$line = ' ' x ((length $class) + 2);
	$line .= $2;
    }

    my ($leading_spaces, $text) = $line =~ /^(\s+)(\S.+)$/;
    $text =~ s/\s+/ /g;

    $indent = length($leading_spaces);
    my @parts = split(' ', $text);
    if ($indent == $previous_indent) {
	pop @slot_titles;
	push @slot_titles, $parts[0];
	my $next_line = $model_lines->[$ix + 1];
	my ($nl_leading_spaces) = $next_line =~ /^(\s+)\S+/;
	if (length $nl_leading_spaces > $indent) {
	    shift @parts;
	    push @slot_titles, $parts[0];
	}
    }
    elsif ($indent > $previous_indent) {
	pop @slot_titles;
	push @slot_titles, $parts[0];
	push @indent_lengths, $indent;
    }
    else {
	for my $ix (1 .. @indent_lengths) {
	    pop @slot_titles;
	    if ($indent_lengths[-$ix] > $indent) {
		pop @indent_lengths;
	    }
	    else {
		last;
	    }
	}
	push @slot_titles, @parts[0];

	my $next_line = $model_lines->[$ix + 1];
	my ($nl_leading_spaces) = $next_line =~ /^(\s+)\S+/;
	if (length $nl_leading_spaces > $indent) {
	    shift @parts;
	    push @slot_titles, $parts[0];
	}
    }

    shift @parts;

    for my $ix (0 .. @slot_titles - 1) {
	$slot_titles[$ix] =~ s/_UNIQUE//;
    }
    	
    $previous_indent = $indent;
}


sub get_and_parse_model_lines {
    my @model_lines = file($WB_MODEL_SPEC)->slurp(chomp => 1);
    my @parsed_lines;
    for my $line (@model_lines) {
	next if $line =~ /^\s*\/\//;
	next if $line =~ /^\s*$/;
	
	# Remove evidence hashes
	$line =~ s/\s#Evidence//g;

	# Remove XREFs
	$line =~ s/\sXREF\s\w+//g;

	$line =~ s/ UNIQUE/_UNIQUE/g;
	push @parsed_lines, $line;
    }

    return \@parsed_lines;
}

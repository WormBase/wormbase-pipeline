=head1 NAME

ProteinConsequence::RunPartialVep - attempt to complete failed VEP analyses

=cut

=head1 DESCRIPTION

Attempts to complete failed VEP analyses.
Removes lines successfully processed from VEP input VCF.
Runs VEP on remaining variations.
Conjoins output with original failed output.

=cut

package ProteinConsequence::RunPartialVep;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    $self->param('vep_failure', 1);
    my $batch_id = $self->required_param('batch_id');
    my $failed_file = $self->required_param('output_dir') . "/$batch_id/$batch_id";
    my $input_file = $self->create_input_file($failed_file);
    my $output_file = $input_file . '.vep.txt';

    my $species = $self->required_param('species');
    my $output_dir = $self->required_param('output_dir');
    my $fasta_file = "$output_dir/$species.fa.gz";
    my $gff_file = "$output_dir/$species.gff.gz";

    my $plugin_str = 'ProtFuncAnnot,mod=WB,pass=' .
	$self->required_param('password');


    # Run VEP on new file
    $self->dbc->disconnect_when_inactive(1);
    my $vep_dir = $self->required_param('vep_dir');
    my $cmd = "$vep_dir/vep -i $input_file -gff $gff_file -fasta $fasta_file --hgvs --hgvsg -shift_hgvs=0 --symbol --distance 0 --output_file $output_file --numbers --force_overwrite --plugin $plugin_str";
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }

    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	$self->warning("ERROR: $flat_cmd - exit code: $exit_code: $stderr") if $exit_code != 0;
    }
    else {
	$self->join_results($failed_file, $input_file);
	$self->remove_header($failed_file) unless $failed_file =~ /\w0+1$/;
	$self->param('vep_failure', 0);
    }
}


=head2 create_input_file

    Title:    create_input_file
    Function: create VCF file for VEP input that only include unprocessed variations
    Args:     filename (string) of VCF input file for failed VEP analysis
    Returns:  filename (string) of file created for re-analysis

=cut

sub create_input_file {
    my ($self, $failed_file) = @_;

    my ($last_id, $last_pos, $last_allele) = $self->last_vep_result_printed($failed_file);

    my $input_file = $failed_file . '_part';
    for my $file ($input_file, $input_file . '.vep.txt') {
	system("rm $file") if -e $file;
    }

    my $cmd;
    if ($last_pos == 0) {
	$cmd = "cp $failed_file $input_file";
    }
    else {
	my $grep_cmd = 'grep -n $' . "'^" . $last_id . '\t' . $last_pos .
	    '\t' . $last_allele . '\t' . "' " . $failed_file . ' |';
	open (GREP, $grep_cmd) 
	    or die "Couldn't grep file $failed_file for $last_id, $last_pos, $last_allele";
	my $lines_printed;
	while (<GREP>) {
	    ($lines_printed) = $_ =~ /^(\d+):/;
	}
	close (GREP);

	my $input_lines = `wc -l < $failed_file`;
	die "Line count for $failed_file failed: $?" if $?;
	my $lines_remaining = $input_lines - $lines_printed;
	$cmd = "tail -n $lines_remaining $failed_file > $input_file";
    }

    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    die "Couldn't create input file: $flat_cmd: $exit_code: $stderr" if $exit_code != 0;

    return $input_file;
}


=head2 join_results

    Title:    join_results
    Function: combines original and new VEP output files
    Args:     filename (string) of original VCF input file for failed VEP analysis,
              filename (string) of VCF input file created for re-analysis by VEP
    Returns:  n/a

=cut

sub join_results {
    my ($self, $failed_file, $input_file) = @_;

    my $failed_results_file = $failed_file . '.vep.txt';
    my $input_results_file = $input_file . '.vep.txt';

    system("grep -v '^#' $input_results_file > $input_results_file.copy");
		
    my $cmd = "cat $input_results_file.copy >> $failed_results_file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    die "Couldn't concatenate $failed_results_file and $input_results_file: $flat_cmd - $stderr"
	unless $exit_code == 0;
		
    system("rm $input_file*");
    
    return;
}


=head2 last_vep_result_printed
    
    Title:    last_vep_result_printed
    Function: retrieves details of last result successfully printed in output of failed VEP analysis
    Args:     filename (string) of VCF input for failed VEP analysis
    Returns:  variation_id, position, allele

=cut

sub last_vep_result_printed {
    my ($self, $failed_file) = @_;

    my $last_pos = 0;
    my $last_allele = '';
    my $last_id = '';
    my $last_line = '';

    my $vep_file = $failed_file . '.vep.txt';
    return ($last_id, $last_pos, $last_allele) unless -e $vep_file;

    open(VEP, "<$vep_file") or die "Couldn't open $vep_file";
    open(TMP, ">$vep_file.tmp");
    while (<VEP>) {
	last if $_ !~ /\n$/;
	print TMP $_;
	next if $_ =~ /^#/;
	my @columns = split("\t", $_);
	$last_id = $columns[0];
	$last_pos = $columns[1];
	$last_allele = $columns[2];
    }
    close(VEP);
    close(TMP);

    my ($exit_code, $stderr, $flat_cmd) =
	$self->run_system_command("mv $vep_file.tmp $vep_file");
    die "$flat_cmd: $exit_code: $stderr" unless $exit_code == 0;

    return ($last_id, $last_pos, $last_allele);
}


sub write_output {
    my $self = shift;

    if ($self->param('vep_failure')) {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 4);
    }
    else {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 2);
    }
}
 

1;

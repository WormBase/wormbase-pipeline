#!/usr/bin/perl -w
use Tk;
use strict;
use Getopt::Long;
require Tk::Dialog;


my ($species, $nocheck, $noload);
GetOptions (
	    "species:s" => \$species, # defaults to 'elegans'
	    "nocheck"   => \$nocheck, # by default it asks if you have saved your session, this stops it
	    "noload"    => \$noload,  # by default it loads the ACE file into the database, this stops it
	   );
if (!defined $species) {$species = 'elegans'}

my @cmds = qw/Replace Isoform Split Merge To_CDS To_Pseudogene To_Transcript Create_Gene Delete Last_reviewed Just_Make_History Check_Gene_name/;


my $mw = MainWindow->new;
$mw->title( "Curate $species" );
$mw->configure(-background => 'LightSteelBlue2');

my $tf = Top_frame($mw);
my $bf = Bottom_frame($mw);
my $mf = Middle_frame($mw);
my $status = Status($mw);


# need the old and new text widgets to exist before the cmd menu as the menu tries to reconfigure the text widgets
my $old = Old($mf); 
my $new = New($mf);

my $cmd_variable;
my $cmd = Cmd($tf, \$cmd_variable);

my $force_variable;
my $force = Force($tf, \$force_variable);

my $class_variable;
my $class = Class($tf, \$class_variable);

my $submit = Submit($bf);

MainLoop;

#############################
sub Top_frame {
  my ($f) = @_;
  my $rf = $f->Frame( -background => "LightSteelBlue2", 
#		      -label      => "New",
		    )->pack(-pady => '5',
			    -padx => '5',
			    -side => 'top');
  return $rf;
}
#############################
sub Middle_frame {
  my ($f) = @_;
  my $rf = $f->Frame( -background => "LightSteelBlue2", 
#		      -label      => "New",
		    )->pack(-pady => '5',
			    -padx => '5',
			    -side => 'bottom');
  return $rf;
}
#############################
sub Bottom_frame {
  my ($f) = @_;
  my $rf = $f->Frame( -background => "LightSteelBlue2", 
#		      -label      => "New",
		    )->pack(-pady => '5',
			    -padx => '5',
			    -side => 'bottom');
  return $rf;
}
#############################
sub Cmd {
  my ($f, $cmd_var_ref) = @_;
  
  my $cmd = $f->Optionmenu(
			   -variable => $cmd_var_ref,
			   -options  => [@cmds],
			   -background => 'LightSteelBlue',
			   -command  => [\&command_chosen],
			  )->pack(-pady => '5',
				  -padx => '5',
				  -side => 'left',
				  -anchor => 'w',
				 );
  return $cmd;
}
#############################
sub Force {
  my ($f, $force_var_ref) = @_;
  my $force = $f->Checkbutton(
			      -text => 'Force', 
			      -variable => $force_var_ref,
			      -background => 'LightSteelBlue',
			     )->pack(-pady => '5',
				     -padx => '5',
				     -side => 'right',
				     -expand => 1,
				    );

  return $force;
}
#############################
sub reset_Force {
  my ($f) = @_;
  $f->deselect;
}
#############################
sub Old {
  my ($f) = @_;
  my $old = $f->Text(-width  => 20,
		     -height => 5,
		     -wrap   => 'none',
		     -background => "white",
		    )->pack(-pady => '5',
			    -padx => '5',
			    -side => 'left');
  return $old;
}
#############################
sub New {
  my ($f) = @_;
  my $new = $f->Text(-width  => 20,
		     -height => 5,
		     -wrap   => 'none',
		     -state  => 'normal',
		     -background => "white",
		    )->pack(-pady => '5',
			    -padx => '5',
			    -side => 'right');
  return $new;
}
#############################
sub Status {
  my ($f) = @_;
  my $status = $f->Text(-width  => 60,
			-height => 1,
			-wrap   => 'none',
			-background => "LightCyan",
		       )->pack(-pady => '5',
			       -padx => '5',
			       -side => 'bottom');
  return $status;
}
#############################
sub Class {
  my ($f, $class_var_ref) = @_;
  my @classes = qw/
CDS Pseudogene Transcript Transposon_CDS
/;
  
  my $class = $f->Optionmenu(
			    -variable => $class_var_ref,
			    -options  => [@classes],
			    -background => 'LightSteelBlue',
			   )->pack(-pady => '5',
				   -padx => '5',
				   -side => 'left'
				  );
  return $class;
}
#############################
sub Submit {
  my ($f) = @_;
  my $submit = $f->Button( -text => "Save your session and Make It So",
			   -command => [\&submit],
			   -background => 'LightSteelBlue',
                              )->pack(
                                      -pady => '5',
                                      -padx => '5',
                                      -side => 'right',
                                     );

}
#############################
# delete the contents of a text widget
sub delete_text {
  my ($f) = @_;
  $f->delete('0.0', 'end');
}
#############################
# append to contents of a text widget
sub add_text {
  my ($f, $text) = @_;
  $f->insert('end', $text);
}
#############################
# get contents of a text widget
sub get_text {
  my ($f) = @_;
  return $f->get('0.0', 'end');
}
#############################
# inactivate text widget
sub inactivate_text {
  my ($f) = @_;
  $f->configure(-state  => 'disabled',
		-background => "grey",
	       );
}
#############################
# activate text widget
sub activate_text {
  my ($f) = @_;
  $f->configure(-state  => 'normal',
		-background => "white",
	       );
}
#############################
# activate text widget
sub one_line_text {
  my ($f) = @_;
  $f->configure(-height => 1,
	       );
}
#############################
# activate text widget
sub many_lines_text {
  my ($f) = @_;
  $f->configure(-height => 5,
	       );
}
#############################
# SUBMIT COMMAND
sub submit {

#  print "Submit command\n";

  if (!$nocheck) {confirm_message('Session', 'Have you saved your session?')}

  my $old_text = get_text($old);
  if (!defined $old_text && $cmd_variable ne 'Create_Gene') {
    error_warning('ERROR', "No existing sequence name was given");
    return;
  }
  my @old_text = split /\n/, $old_text;
  my $o = join ',', @old_text;

  my $new_text = get_text($new);
  my @new_text = split /\n/, $new_text;
  my $n = join ',', @new_text;

  my $cmd_var = lc $cmd_variable;

  # run script
  my $c = "curate.pl -species $species -$cmd_var -class $class_variable ";
  if ($o ne '') {$c .= "-old '$o' "};
  if ($n ne '') {$c .= "-new '$n' "};

  if ($noload) {$c .= "-noload "}

  if ($force_variable) {$c .= "-force "}
  reset_Force($force);

  print "$c\n";
  my $output = qx(/software/bin/perl \$CVS_DIR/$c  2>&1);

  # display output text
  if ($output =~ /ERROR/) {
    error_warning('ERROR', $output);
  } elsif ($output =~ /-force/) {
    error_warning('WARNING', $output);
  } elsif ($output =~ /Nameserver/) {
    wide_confirm_message('Nameserver Action Required', "Nameserver Action Required\n$output");
  } else {
    if ($noload) {
      wide_confirm_message("OK -noload so no change.", "OK -noload so no change.\n\n$output");
    } else {
      wide_confirm_message("OK - that's been done.", "OK - that's been done.\n\n$output");
    }
  }

}
#############################
# the command menu has been changed.
# set up other widgets as required
sub command_chosen {
  my ($selection) = @_;

#  print "selection = $selection\n";

  # clear text and set default
  delete_text($old);
  delete_text($new);
  activate_text($old);
  activate_text($new);
  if ($selection eq 'Replace') {
    one_line_text($old);
    one_line_text($new);
    add_text($new, "temp_gene_1");
    delete_text($status); add_text($status, "Replace old structure by the new structure")
  } elsif ($selection eq 'Isoform') {
    one_line_text($old);
    one_line_text($new);
    add_text($new, "temp_gene_1");
    delete_text($status); add_text($status, "Add the new structure as an isoform of the old one")
  } elsif ($selection eq 'Split') {
    one_line_text($old);
    many_lines_text($new);
    add_text($new, "temp_gene_1");
    delete_text($status); add_text($status, "Split old locus, the first new locus inherits the GeneID")
  } elsif ($selection eq 'Merge') {
    many_lines_text($old);
    one_line_text($new);
    add_text($new, "temp_gene_1");
    delete_text($status); add_text($status, "Merge old loci forming new locus, Gene to Live is first name")
  } elsif ($selection eq 'To_CDS') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Convert the old name into a CDS")
  } elsif ($selection eq 'To_Pseudogene') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Convert the old name into a Pseudogene")
  } elsif ($selection eq 'To_Transcript') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Convert the old name into a Transcript (Method/Type = ncRNA)")
  } elsif ($selection eq 'Create_Gene') {
    inactivate_text($old);
    one_line_text($old);
    one_line_text($new);
    add_text($new, "temp_gene_1");
    delete_text($status); add_text($status, "Create a new Gene")
  } elsif ($selection eq 'Delete') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Delete the old structure")
  } elsif ($selection eq 'Last_reviewed') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Set the 'Last_reviewed' tag")
  } elsif ($selection eq 'Just_Make_History') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Just make a History object")
  } elsif ($selection eq 'Check_Gene_name') {
    inactivate_text($new);
    one_line_text($old);
    one_line_text($new);
    delete_text($status); add_text($status, "Check for a CGC Gene name")
  } else {
    error_warning('ERROR', "Unknown command from command option menu: '$selection'");
    die ("Unknown command from command option menu: '$selection'\n");
  }

}
#############################

#############################

sub error_warning {
  my ($title, $text) = @_;
  
  my $D = $mw->Dialog(
		      -width  => 100,
		      -title => "$title",
		      -text  => "$text",
		      -default_button => 'ok',
		      -buttons        => ['ok']
		     );
  $D->configure(-background => 'Red',
		-foreground => 'black');
  
  my $choice = $D->Show;
}

#############################
sub confirm_message {
  my ($title, $text) = @_;
  
  my $D = $mw->Dialog(
		      -title => "$title",
		      -text  => "$text",
		      -default_button => 'ok',
		      -buttons        => ['ok']
		     );
  $D->configure(-background => 'DarkSeaGreen', # was cyan PaleGreen
		-foreground => 'black');
  
  my $choice = $D->Show;
}
#############################
sub wide_confirm_message {
  my ($title, $text) = @_;
  
  my $D = $mw->Dialog(
		      -width  => 100,
		      -title => "$title",
		      -text  => "$text",
		      -default_button => 'ok',
		      -buttons        => ['ok']
		     );
  $D->configure(-background => 'DarkSeaGreen', # was cyan PaleGreen
		-foreground => 'black');
  
  my $choice = $D->Show;
}


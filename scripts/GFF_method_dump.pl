#!/usr/local/bin/perl5.8.0 -w
#
# GFF_method_dump.pl
#
# by Anthony Rogers
#
# Selectively dump GFF for certain acedb methods
#
# dumps the method through sace / tace and concatenates them.
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2014-10-22 13:14:48 $


use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Storable;
use File::stat;

my ($help, $debug, $test, $quicktest, $database, $species, @methods, @chromosomes, $dump_dir, @clones, $list,$host, $giface, $giface_client, $gff3,$port,$fprefix );
my @sequences;
my $store;
GetOptions (
  "help"           => \$help,
  "debug=s"        => \$debug,
  "test"           => \$test,
  "store:s"        => \$store,
  "species:s"      => \$species,
  "quicktest"      => \$quicktest,
  "database:s"     => \$database,
  "dump_dir:s"     => \$dump_dir,
  'host:s'         => \$host,
  'port:s'         => \$port,
  'giface:s'       => \$giface,
  'gifaceclient:s' => \$giface_client,
  'gff3'           => \$gff3,

  # ive added method and methods for convenience
  "method:s"      => \@methods,
  "methods:s"     => \@methods,
  "chromosomes:s" => \@chromosomes,
  "clone:s"       => \@clones,
  "clones:s"      => \@clones,
  "list:s"        => \$list,
  "fprefix:s"     => \$fprefix,
    );


my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism=> $species
			   );
}
$species = $wormbase->species; #incase defaulting to elegans
my $log = Log_files->make_build_log($wormbase);

@methods     = split(/,/,join(',',@methods));
@chromosomes = split(/,/,join(',',@chromosomes));
@sequences = split(/,/,join(',',@clones)) if @clones;

if ($host) {
  $log->log_and_die("You must supply a valid saceclient binary to be used\n")
      if not defined $giface_client or not -x $giface_client;
  $log->log_and_die("You must supply port when dumping using server\n") 
      if not $port;
} else {
  $log->log_and_die("You must supply a valid giface binary to be used\n")
      if not defined $giface or not -x $giface;
}

if ($wormbase->assembly_type eq 'contig' and not $host) {
  $log->log_and_die("no host passed for contig assembly");
}

$fprefix = "" if not defined $fprefix;
my $file_suffix = ($gff3) ? "gff3" : "gff";

$database = $wormbase->autoace unless $database;
$dump_dir = "/tmp/GFF_CLASS" unless $dump_dir;

# ensure that there is no '/' at the end of $dump_dir because 'gif
# seqfeatures' removes everything after a '//' as it looks like a
# comment in acedb!
$dump_dir =~ s#/$##g;

&check_options;
`mkdir -p $dump_dir/tmp` unless -e "$dump_dir/tmp";

#make sure dump_dir is writable
system("touch $dump_dir/tmp_file.$$") and $log->log_and_die ("cant write to $dump_dir\n");
unlink "$dump_dir/tmp_file.$$";

unless($host) {
  open (NONSERVERWRITE,"| $giface $database") or $log->log_and_die ("failed to open giface connection to $database\n");
}

$log->write_to("dumping methods:".join(",",@methods)."\n");

my $count=0; # debug hack

if ( @methods ) {
  foreach my $method ( @methods ) {
    if($host) {
      my $out_file = "$dump_dir/${fprefix}${method}.${file_suffix}";
      if (-e $out_file) {
        unlink $out_file;
      }
      $wormbase->run_command("touch $out_file", 'no_log');
      
      foreach my $sequence ( @sequences ) {
        my $file = "$dump_dir/tmp/${sequence}_${method}.$$";
        $file .= ($gff3) ? ".gff3" : ".gff2";
        my $cmd = sprintf("gif seqget %s +method %s; seqfeatures -version %s -file %s", 
                          $sequence, 
                          $method,
                          ($gff3) ? "3" : "2",
                          $file);
        
        #
        # For reasons I cannot get to the bottom of, the client sometimes fails on larger 
        # sequences, even though the server is still happy. So if we get a fail here, we
        # try once again
        #

        foreach my $try (1,2) {
          open (my $gif_out,"echo '$cmd' | $giface_client $host -port $port -userid wormpub -pass blablub |") or $log->log_and_die("$!\n");
          while (my $line = <$gif_out>) {
            if ($line =~ 'ERROR') {
              $log->error;
              $log->write_to("ERROR detected while GFF dumping $sequence:\n\t$line\n\n");
              print "ERROR detected while GFF dumping $sequence:\n\t$line\n\n";
            }
          }
          my $success = close($gif_out);
          if ($success) {
            last;
          } else {
            if ($try == 1) {
              unlink $file if -e $file;
              continue;
            } else {
              $log->log_and_die("Failed twice to run command '$cmd'\n");
            }
          }
        }
          
        my $check = &check_the_file($file, $sequence);
        if ($check) {
          $log->log_and_die("Something went wrong with dumping of $sequence: $check\n");
        }
        
        $wormbase->run_command("cat $file >> $out_file",'no_log');
        unlink $file;
      }
    } else {
      foreach my $sequence (@sequences) {
        my $file = "$dump_dir/${fprefix}${sequence}_${method}.${file_suffix}";
        
        my $command = sprintf("gif seqget %s +method %s; seqactions -hide_header; seqfeatures -version %s -file %s\n", 
                              $sequence, 
                              $method,
                              ($gff3) ? "3" : "2",
                              $file);
        print NONSERVERWRITE $command;
      }
    }
  }
} else {
  if($host) {
    my $out_file = "$dump_dir/${fprefix}${species}.${file_suffix}";

    if (-e $out_file) {
      unlink $out_file;
    }
    $wormbase->run_command("touch $out_file", 'no_log');

    foreach my $sequence (@sequences) {
      my $file = "$dump_dir/tmp/${sequence}.gff_dump.$$"; 
      $file .= ($gff3) ? ".gff3" : ".gff2";
      my $cmd = sprintf("gif seqget %s; seqfeatures -version %s -file %s", 
                        $sequence, 
                        ($gff3) ? "3" : "2",
                        $file);
      $log->write_to("Running '$cmd'\n");
      open (WRITEDB,"echo '$cmd' | $giface_client $host -port $port -userid wormpub -pass blablub |") or $log->log_and_die("$!\n");
      while (my $line = <WRITEDB>) {
	if ($line =~ 'ERROR') {
	  $log->error;
	  $log->write_to("ERROR detected while GFF dumping $sequence:\n\t$line\n\n");
	  print "ERROR detected while GFF dumping $sequence:\n\t$line\n\n";
	}
      }
      close(WRITEDB) or $log->log_and_die("Did not successfully close command: $cmd\n");
	
      my $check = &check_the_file($file, $sequence);
      if ($check) {
        $log->log_and_die("Something went wrong with dumping of $sequence: $check\n");
      }

      $wormbase->run_command("cat $file >> $out_file", 'no_log');
      unlink $file;
    }
  } else {
    foreach my $sequence (@sequences) {
      my $out_file = "$dump_dir/${fprefix}${sequence}.${file_suffix}";
      my $command = sprintf("gif seqget %s; seqactions -hide_header; seqfeatures -version %s -file %s\n",
                            $sequence, 
                            ($gff3) ? "3" : "2",
                            $out_file);
      print "$command";
      print NONSERVERWRITE $command;
    }
  }
  $count++;
}

unless($host) {
  close(NONSERVERWRITE) or $log->log_and_die ("failed to close giface connection to $database\n");
}

$log->write_to("dumped $count sequences\n");

##################
# Check the files
##################
if ($host) { # contig sequences are contenated
  if ( @methods ) {
    foreach my $method ( @methods ) {
      
      my $file = "$dump_dir/${method}.${file_suffix}";
      $wormbase->check_file($file, $log,
			    lines => ['^##', 
				      "^\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			    gff => 1,
			   );   
    }


  } else { 
    my $file = "$dump_dir/${species}.${file_suffix}";
    $wormbase->check_file($file, $log,
			  lines => ['^##', 
				    "^\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			  gff => 1,
                          maxsize => 3771054952,
			 );   

  }

} else { # chromosome sequences are kept separate
  if ( @methods ) {
    foreach my $method ( @methods ) {

      foreach my $sequence ( @sequences ) {
	my $file = "$dump_dir/${sequence}_${method}.${file_suffix}";
	$wormbase->check_file($file, $log,
			      lines => ['^##', 
					"^${sequence}\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			      gff => 1,
			     );   
	
	
      }
    }

  } else { 
    foreach my $sequence ( @sequences ) {
      my $file = "$dump_dir/${sequence}.${file_suffix}";

      if($wormbase->species eq 'elegans') {
      
	my %sizes = (
		     'CHROMOSOME_I'       => 230000000,
		     'CHROMOSOME_II'      => 250000000,
		     'CHROMOSOME_III'     => 230000000,
		     'CHROMOSOME_IV'      => 290000000,
		     'CHROMOSOME_V'       => 330000000,
		     'CHROMOSOME_X'       => 220000000,
		     'CHROMOSOME_MtDNA'   =>   2000000,
		    );
	$wormbase->check_file($file, $log,
			      minsize => $sizes{$sequence},
			      lines => ['^##', 
					"^${sequence}\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			      gff => 1,
			   );   
      } elsif ($wormbase->species eq 'briggsae') {

	my %sizes = (
		     'chr_I'          => 150000000,
		     'chr_II'         => 200000000,
		     'chr_III'        => 400000000,
		     'chr_IV'         => 200000000,
		     'chr_V'          => 250000000,
		     'chr_X'          => 250000000,
		    );
	if (exists $sizes{$sequence}) {
	  $wormbase->check_file($file, $log,
				minsize => $sizes{$sequence},
				lines => ['^##', 
					  "^${sequence}\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
				gff => 1,
			       );   
	}
      } else {
	$wormbase->check_file($file, $log,
			      lines => ['^##', 
					"^${sequence}\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			      gff => 1,
			     );   
	
      }

    }
  }
}


$log->mail();
exit(0);

###############################
sub check_the_file {
  my ($f, $seq) = @_;

  eval {
    my $found_seq_line = 0;

    open(my $fh, $f) or die "Could not open file $f\n";
    while(<$fh>) {
      /^\#\#sequence\-region\s+$seq/ and $found_seq_line = 1;
    }
    close($fh) or die "Could not close file $f\n";

    if (not $found_seq_line) {
      die "Could not find sequence-region line in $f\n";
    }
  };
  return $@;
}


######################
sub check_options {

  unless($list or @clones ) {
    
    # -chromosomes
    
    my @chrom =  $wormbase->get_chromosome_names(-mito => 1, -prefix => 1);
    
    my %chroms = map {$_ => 1} $wormbase->get_chromosome_names(-mito => 1, -prefix => 1);
    
    unless (@chromosomes ) {
      @sequences= @chrom;
      print "Dumping for all chromosomes\n";
    } 
    else {
      foreach (@chromosomes) {
    	if ( $chroms{$_} ) {
          push( @sequences,$_);
        }
        else {
          $log->log_and_die ("$_ is not a valid chromosome\n");
        }
      }
    }
  }
  
  &process_list if $list;
  
  # -database
  if ( $database ){
    if( -e "$database" ) {
      if( -e "$database/wspec/models.wrm" ) {
        print "$database OK\nDumping @methods for chromosomes @chromosomes\n";
        return;
      }
    }
  }
  else {
    $log->log_and_die ("You must enter a valid database\n");
  }
  $log->log_and_die ("$database is not a valid acedb database\n");	
}

#####################
sub process_list {
  open(LIST,"<$list") or $log->log_and_die ("bad list $list\n");
  while(<LIST>) {
    chomp;
    push(@sequences,$_);
  }
}

=pod 

=head1 NAME - GFF_method_dump.pl

=head2 USAGE

=over 4

This script will GFF dump specified methods from a database

It is use by dump_gff_batch.pl so if you change it make sure it is still compatible !

=back

=item MANDATORY ARGS:

=over 4

-methods     Methods to dump eg curated,history (Comma separated list)

=back

=item OPTIONAL ARGS:

=over 4

-database    Database to dump from ( default /wormsrv2/autoace )

-chromosomes Chromosomes to dump as comma separated list eg I,II,X ( defaults to all )

=back

=item OUTPUT:

=over 4

A separate file is written for each method for each chromosome and is named 

CHROMOSOME_($chrom)_($method).gff

=back

=item EXAMPLES:

=over 4

GFF_method_dump.pl -database wormsrv2/autoace -chromosomes II,V -method curated,TRANSCRIPT

will GFF dump separate curated and TRANSCRIPT files for both chromosomes II and V ie

  CHROMOSOME_II_curated.gff
  CHROMOSOME_II_TRANSCRIPT.gff
  CHROMOSOME_V_curated.gff
  CHROMOSOME_V_TRANSCRIPT.gff

=back

=item REQUIREMENTS:

=over 4

Access to ~acedb for giface binary

=back

=item WARNING:

=over 4

At time of writing the version of acedb (4.9y) adds extra data to the GFF output. Eg if you ask for method curated you also get a load of ALLELES that have no method, so the GFF files pr
oduced should be treated with caution and may need post-processing.

=back

=cut

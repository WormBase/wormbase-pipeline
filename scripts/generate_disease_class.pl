#!/software/bin/perl -w
#
# generate_disease_class.pl
# 

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $omimfile, $acefile, $noload);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noload"     => \$noload,
            "omimfile"   => \$omimfile,
            "acefile"    => \$acefile,
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

$acefile = $wormbase->acefiles . "/Human_disease.ace" if not defined $acefile;
$omimfile = $wormbase->misc_static . "/omim.txt" if not defined $omimfile;

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my ($entry_hash, $old2newmap) = &parse_omim_flatfile($omimfile);

my $omim2wbgene = &get_omim_to_wormbase_gene_mapping();

open(my $outfh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

foreach my $num (sort { $a <=> $b } keys %$entry_hash) {
  my $e = $entry_hash->{$num};
  print $outfh "\nDisease : \"$num\"\n";
  print $outfh "Type Gene\n" if $e->{is_gene};
  print $outfh "Type Disease\n" if $e->{is_disease};
  print $outfh "Name \"$e->{label}\"\n";
  foreach my $syn (@{$e->{synonms}}) {
    print $outfh "Synonym \"$syn\"\n";
  }
  
  if (exists $omim2wbgene->{$num}) {
    foreach my $g (sort keys %{$omim2wbgene->{$num}}) {
      print $outfh "Gene $g\n";
    }
  }
}

$log->mail;
exit(0);




#############################
sub parse_omim_flatfile {
  my ($file, $omim2gene) = @_;

  my %old_to_new;
  my %removed;
  my %entries;

  local $/ = "*RECORD*";

  open(my $mim_io, $file) or $log->log_and_die("Could not open $file for reading\n");

  my $dummy = <$mim_io>; # first record is empty with *RECORD* as the
                         # record seperator
  while (<$mim_io> ) {
    #get the MIM number
    my $number = 0;
    my $label = undef;
    my $is_morbid = 0;
    my $type =undef;

    if(/\*FIELD\*\s+NO\n(\d+)/){
      $number = $1;

      if(/\*FIELD\*\sTI\n([\^\#\%\+\*]*)\d+\s+([\S\s]*)\*FIELD\*\sTX\n/m){
        #print;
        $type = $1;
        my $id_lines = $2;

        my ($first_line, @others) = split(/\n/, $id_lines);
        my $syns = join(" ", @others);

        $label = $first_line; # taken from description as acc is meaning less
        $label =~ s/\;\s[A-Z0-9]+$//; # strip gene name at end

        my @syn = map { $_ =~ s/^\s*// and $_ } grep { $_ } split(/;/, $syns);

        if($type eq "*"){ # gene only
          $entries{$number} = { 
            acc        => $number,
            label      => $label,
            is_gene    => 1,
            is_disease => 0,
            synonyms   => \@syn,
          };

        }
        elsif((!defined $type) or ($type eq "") or ($type eq "#") or ($type eq "%")){ #phenotype only
          $entries{$number} = { 
            acc        => $number,
            label      => $label,
            is_disease => 1,
            is_gene    => 0,
            synonyms   => \@syn,
          };
        }
        elsif($type eq "+"){ # both
          
          $entries{$number} = { 
            acc        => $number,
            label      => $label,
            is_gene    => 1,
            is_disease => 1,
            synonyms   => \@syn,
          };

        }
        elsif($type eq "^"){
          if(/\*FIELD\*\sTI\n[\^]\d+ MOVED TO (\d+)/){
            $old_to_new{$number} = $1;
          }
          else{
            $removed{$number} = 1;
          }        
        }
      }
    }
  }

  foreach my $mim (keys %old_to_new){
    my $old= $mim;
    my $new= $old_to_new{$old};
    while(defined($old_to_new{$new})){
      $new = $old_to_new{$new};
    }
    
    if(!defined($removed{$new})){
      push @{$entries{$new}->{synonyms}}, $old;
    }
  }

  return (\%entries, \%old_to_new);
}


##########################
sub get_omim_to_wormbase_gene_mapping {
  
  my %omim2gene;

  my $tb_file = &get_omim_tm_def;
  my $tace = $wormbase->tace;
  my $db = $wormbase->autoace;

  my $tb_cmd = "Table-maker -p \"$tb_file\"\nquit\n";
  open(my $tace_fh, "echo '$tb_cmd' | $tace $db |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"/ and do {
      my ($gid, $type, $omim_id) = ($1, $2, $3);
      $omim2gene{$omim_id}->{$gid} = 1;
    };
  }
  close($tace_fh);
  unlink $tb_file;

  return \%omim2gene;
}


##########################
sub get_omim_tm_def {

  my $tmdef = "/tmp/omim2wbgene.def";

   open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition "OMIM"
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Database_field 
Right_of 2 
Tag  HERE  
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Text
Right_of 3 
Tag  HERE  

EOF

  print $qfh $query;
  close($qfh);
  return $tmdef;
}


__END__

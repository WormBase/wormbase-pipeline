#!/usr/local/perl -w


=pod

=head1 NAME

  make_common_data.pl

=head1 DESCRIPTION

  Converts a tab separated text file in to a COMMON_DATA dat file

=head1 OPTIONS

  -file    The tab delimited file containing the data
  -name    Name to call the common data file 
  -key     The field in the line of text to use as the hash key ( default 0 )
  -value   The field in the line of text to use as the hash value ( default 1 )
  -array   The value field can be an array. Constructed from all fields after the key ( see example )
  -recip   Construct a reciprocal hash and dat file eg accession2 clone and clone2accession.  Supply name of reciprocal dat file
  -append  Where keys may have multipe values they can be appended to the current value separated with a space

=head1 EXAMPLES

=item  FILE file.gff lines as such . . 

CHROMOSOME_I    curated CDS     11641   16585   .       +       .       CDS "Y74C9A.2"

  make_common_data.pl -file file.gff -name CDS2start -key 9 -value 3 -recip start2CDS

  would make 2 dat files 

  CDS2start.dat
  $VAR1 = {
          'Y74C9A.2' => '11641'
	  };

  start2CDS.dat
  $VAR1 = {
          '11641' => 'Y74C9A.2'
	  };

=item  FILE orthologs.txt

  2L52.1  CBG07041        5.5e-33 0.359   200     seg-on

  make_common_data.pl -file orthologs.txt -array -name elegans2ortho.dat

  would make (something like this )

  $VAR1 = {
	   '2L52.1' => arrayref:([CBG07041, 5.5e-33, 0.359, 200, seg-on])
	  }

=head1 NOTE

  any '"' characters in the data will be removed.

  will write to /wormsrv2/autoace/COMMON_DATA if available, other wise ~wormpub/TEST_BUILD/autoace/COMMON_DATA

=cut

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my ($file, $name, $key, $value, $array, $recip, $append);
GetOptions( 'file:s' => \$file,
	    'name:s' => \$name,
	    'key:s'  => \$key,
	    'value:s'=> \$value,
	    'array'  => \$array,
	    'recip:s'=> \$recip,
	    'append' => \$append
	  );

die "cant have value and array set\n" if ( $value and $array );

my %store;
$key = 0 unless $key;
$value = 1 unless $value;

open( DATA ,"<$file") or die "cant open $file\t$!\n";
while( <DATA> ) {
  s/\"//g;
  my @in = split;
  my $key_data = $in[$key];
  my $value_data;
  
  if( $array ) {
      $value_data = splice( @in,$key+1);
    }
    else {
      $value_data = $in[$value];
    }
    if ( $store{$key_data} ) {
      if( $append ) {
	$store{$key_data} .= " $value_data";
      }
      else {
	warn "duplicate key : $key_data will be overwritten\n";
      }
    }
    else {
      $store{$key_data} = $value_data;
    }
}
close DATA;

my $common_data_dir = (-e "/wormsrv2/autoace/COMMON_DATA") ? "/wormsrv2/autoace/COMMON_DATA" : glob("~wormpub/TEST_BUILD/autoace/COMMON_DATA");


my ($base,$path,$type) = fileparse( $file );
$base = $name if $name;
my $dat_file = "$common_data_dir"."/$base.dat";
open (DAT,">$dat_file") or die "cant open $dat_file\t$!\n";
print DAT Data::Dumper->Dump([\%store]);
close DAT;

print "commnon data dat file $dat_file created\n";

if( $recip ) {

  my %recip_store;
  foreach my $storekey ( keys %store ) {
    foreach my $storevalue ( split(/\s+/, $store{"$storekey"} ) ) {
      if ( $recip_store{ $storevalue } ) {
	if( $append ) {
	  $recip_store{$storevalue} .= "$storekey";
	}
	else {
	  warn "duplicate key in recip : $storevalue\n";
	}
      }
      else {
	$recip_store{$storevalue} = "$storekey";
      }
    }
  }
  my $dat_file = "$common_data_dir"."/$recip.dat";
  open (DAT,">$dat_file") or die "cant open $dat_file\t$!\n";
  print DAT Data::Dumper->Dump([\%recip_store]);
  close DAT;
  print "commnon data dat file $dat_file created\n";
}

exit(0);

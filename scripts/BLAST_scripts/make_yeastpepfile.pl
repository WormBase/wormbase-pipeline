#!/usr/local/bin/perl5.6.1 -w

# krb 020829

# converts yeastX.pep file with yeast peptide entries into ace file
# Puts SGDID as Accssion in Database field.

my $dir = glob("~wormpipe/BlastDB");
my $acefile = "$dir/yeast.ace";

my $source_file = shift;
my $pepfile = shift;
use strict;


open (SOURCE,"<$source_file");
open (PEP,">$pepfile");
open (ACE,">$acefile");

while (<SOURCE>) {
  if( />/ ) { 
    if (/ORFP:(\S+)\s+(\S+)\s+SGDID:(\w+)/) {
      my $ID = $1;
      print ACE "\nProtein : \"SGD:$ID\"\n";
      print ACE "Peptide \"SGD:$ID\"\n";
      print ACE "Species \"Saccharomyces cerevisiae\"\n";
      print ACE "DB_remark \"SGD gene name is $2\"\n";
      print ACE "Database \"SGD\" \"$ID\" \"$3\"\n\n";
      print ACE "Peptide : \"SGD:$ID\"\n"; 	
      
      print PEP "\>$ID\n";
    }
    else {
      print $_;
    }
  }
  else { 
    print PEP $_;
    print ACE $_;
  }
}

close(SOURCE);
close(PEP);
close(ACE);

print "\n\nabout to copy (scp) $acefile to /wormsrv2/wormbase/ensembl_dumps/\n";
system ("scp -r $acefile wormpub\@wormsrv2:/wormsrv2/wormbase/ensembl_dumps/") and warn "copy $acefile failed\n";

print "\n$source_file is now converted to $pepfile...now removing $source_file\n";
system ("rm -f $source_file") and warn "removing $source_file failed\n";
__END__


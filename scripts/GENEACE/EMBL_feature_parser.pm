#!/usr/local/bin/perl5.8.0 -w

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-06-11 15:47:18 $

package EMBL_feature_parser;

use strict;

# initialize a hash to bless
sub init {
  my $class = shift;
  my $this = { version  => $_[0],
               flatfile => undef,
  };
  bless($this, $class);
}

# download EMBL release automatically, output to a file the entries without protein and DNA sequences
sub download_EMBL_release {

  my $this = shift;
  my $version = $this->{version};

  # directory, files used for EMBL update
  my $output_dir = "/wormsrv1/chaokung/EMBL";
  my $download = "$output_dir/embl";
  my $flatfile= "$output_dir/$version"."_EMBL_entry_no_seqs";

  print "\nDownloading EMBL Release $version . . .\n";

  $this->{flatfile}=$flatfile;

  # query parameters

  `getz -e "(([emblrelease-Division:inv] & [emblrelease-Organism:caenorhabditis*]) & (([emblrelease-FtKey:polya_site] | [emblrelease-FtKey:polya_signal]|[emblrelease-FtQualifier:gene]) > parent )) " > $download`;

  print "\nDeleting protein and DNA sequences from EMBL flatfile...\n";

  open(IN, $download) || die "Can't read in file $download!";

  open(EMBL, ">$flatfile") || die "Can't write to file!";

  ## fetch only ID, AC, SV, KW, De, OS, FT (no DNA/protein seq)
  while(<IN>){
    chomp;
    if ($_ =~ /(^ID[\w\W.]+)|(^AC[\w\W.]+)|(^SV[\w\W.]+)|(^KW[\w\W.]+)|(^DE[\w\W.]+)|(^OS[\w\W.]+)|(^\/\/)/){
      print EMBL $_,"\n";
    }
    if ( ($_ =~ /^FT[\w\W.]+/) && ($_ !~ /^FT\s+\/translation[\w\W.]+/) && ($_ !~ /^FT\s+[A-Z]+/) ){
      print EMBL $_,"\n";
    }
  }
  close EMBL;
  close IN;
  system("rm -f $download");

}

# generic EMBL entry parser to retrieve any feature and its corresponding qualifier(s)
sub get_feature_info {

  my $this = shift;

  # warn message if get_feature_info subroutine is called without parameter
  if (!@_){my @info = caller(); warn_msg(@info)}

  my $flatfile = $this->{flatfile};
  my @Features_qualifiers = @_;

  my (%FEATURE, $feature, @qualifiers, $qualifier);
  my ($AC, $seq_version, @ALL_ACs);#, $f_count);
  my $AC_line = 0; my $chrom =(); my $coord =();

  # loop through each feature and qualifier passed into this routine
  for (my $i=0; $i< scalar @Features_qualifiers; $i=$i+2){

    $feature = $Features_qualifiers[$i];
    my $f_count = 0;

    @qualifiers = @{$Features_qualifiers[$i+1]};
    foreach $qualifier (@qualifiers){

      # initialization
      my $one_line = 0;
      my $multi_line = 0;
      my $count_q = 0;
      my $q_info = ();

      open(IN, $flatfile) || die "Can't open the file $flatfile!";

      while(<IN>){
        chomp;
        my $each_line = $_;

        ################################################################
        #  get AC:  get only 1st accession of each AC line if multiple
        ################################################################
        if ($each_line =~ /^AC\s+(\w+);$/ || $each_line =~ /^AC\s+(\w+);.+$/){
          $AC_line++;
          $AC = $1 if $AC_line == 1; # some AC may have > 1 AC lines, read in only the first AC line
	  print $AC, "\n";
          next;
        }

        ###################
        # get seq version
        ###################
        if ($each_line =~ /^SV.+\.(\d)+/){
          $seq_version = $1;
	  print $seq_version, "\n";
          next;	
        }

        ###################
        #  get chromosome
        ###################
        if ($each_line =~ /^FT\s+\/chromosome=\"(.+)\"/){
          $chrom = $1;
          my %CHR = (1=>"I", 2=>"II", 3=>"III", 4=>"IV", 5=>"V"); # some appear as 1,2,3..
          if (exists $CHR{$chrom}){$chrom = $CHR{$chrom}};
          $chrom =~ s/[^IVX]//g;	# some appear as LG II, linkage group I, chromosome I ..
          next;
        }

        ############################################
        #  generic EMBL Features retrieval lines
        ############################################

        if ($each_line =~ /^FT\s+$feature\s+(\d+\.\.\d+)$/ ||
            $each_line =~ /^FT\s+$feature\s+(\d+)$/ ){  # get coords like 1100 or 100..120

          $coord = $1;
          $f_count = 1;
          $chrom = "NA" if !defined $chrom;
          next;
        }

        # fetch info of one-liner qualifier of =text OR ="text" format when feature counter $f_count is 1
        if ($f_count == 1 && ( $each_line =~ /^FT\s{19,19}\/$qualifier=\"(.+)\"$/ ||
                               $each_line =~ /^FT\s{19,19}\/$qualifier=(\w+)$/ )){

	  push(@{$FEATURE{$AC}}, $seq_version, $chrom, $feature, $coord, $qualifier, $1);
          $one_line = 1;
          next;
        }

	if ( $each_line =~ /^FT\s{19,19}\/$qualifier=\"(.+)\"$/ ||
             $each_line =~ /^FT\s{19,19}\/$qualifier=(\w+)$/ ){

	  push(@{$FEATURE{$AC}}, $seq_version, $qualifier, $1);
          $one_line = 1;
          next;
        }

        # fetch qualifier that has multiple lines
        if ($f_count == 1 && $_ =~ /^FT\s+\/$qualifier=\"(.+)$/){

          $q_info .= $1." ";
          $count_q++;
          next; 
        }
        # fetch lines in between start and end of a qualifier 
        if ($count_q == 1 && $_ !~ /^FT\s+\/$qualifier.+/ && $_ =~ /^FT\s+(.+)[^\"]$/){
          $q_info .= $1." ";
          next;
        }

        # fetch the line which marks the end of an qualifier
        if ($count_q == 1 && $_ !~ /^FT\s+\/$qualifier.+/  && $_ =~ /^FT\s+(.+)\"$/){

          $q_info .= $1." ";

          push(@{$FEATURE{$AC}}, $seq_version, $chrom, $feature, $coord, $qualifier, $q_info);

          $q_info = ();
          $count_q = 0;
          $multi_line = 1;
          next;
        }
        # reinitialize feature counter $f_count to zero when next line is unwanted feature key
	# assign NA to qualifier and its info
        if ( $each_line =~ /^FT\s{3,3}(\w+)\s+.+/ && $each_line !~ /^FT\s+$feature\s+.+/ && $f_count == 1){

          $f_count = 0;
          next;
        }

        #######################
        #  at end of an entry
        #######################
        if ($each_line =~ /^\/\//){

          # assign NA to qualifier and its info to feature which has no qualifier info
          if ($multi_line == 0 && $one_line == 0 && defined $coord){
            push(@{$FEATURE{$AC}}, $seq_version, $chrom, $feature, $coord, "NA", "NA");
          }

          # reinitialization at end of each accession or screwed up
          $chrom=(); $coord =(); $AC_line = 0; $f_count = 0; #$uf_count = 0;
          $count_q = 0; $one_line = 0;
        }
      }
    }
  }
  close IN;
  return %FEATURE;
}

sub protID_2_seq_name {

  ################################################################
  # map EMBL protein_id to WormBase sequence name by wormpep.table
  ################################################################

  my $this = shift;

  # warn if protID_2_seq_name routine is called without parameter
  if (!@_){my @info = caller(); warn_msg(@info)}

  my $prot_id = shift;

  my $temp = "/tmp/seq_protid.tmp"; 
  `cut -f 1,6,7 /wormsrv2/WORMPEP/wormpep111/wormpep.table111 > $temp`;
  system("chmod 777 $temp");

  my $seq_name = `grep $prot_id $temp | cut -f 1`;
  if (!$seq_name){
    $prot_id = `getz -qDes '[swall-prd:$prot_id]'`;
    chomp $prot_id;
    $prot_id =~ s/.+://;
    $seq_name = `grep $prot_id $temp | cut -f 1`; 
    chomp $seq_name;
  }
    return $seq_name;
}	

sub warn_msg {
  ############################################################################
  # general WARN message when calling methods of this package w/o parameter(s)
  ############################################################################

  my $script = $_[1];
  my $line   = $_[2];
  warn "\n--- ERROR: No parameters specified from $script at line $line ---\n\n";
  exit(0);
}

1;


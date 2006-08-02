#!/usr/local/bin/perl5.8.0 -w
#
# overlapcheck.pl
#
# by Kerstin Jekosch
#
# checks whether genes overlap, ESTs match introns and repeats match exons                                   
# sorts output for stl and cam clones
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-08-02 16:22:45 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use IO::Handle;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);



##########################
# MAIN BODY OF SCRIPT
##########################


my @chrom = qw(I II III IV V X);
my (%exon, %est, %genes, %repeat, %intron, %camace, %stlace);

my %estorientation = $wormbase->FetchData('estorientation');     # EST accession => orientation [5|3]

my $tace      = $wormbase->tace;        # TACE PATH
my $ace_dir   =  $wormbase->autoace;

# set parameters for ESTs and repeats (score threshold, overlap threshold)
my $ESTth      = '98.0';
my $ESTlength  = '10';
my $repth      = '15';
my $repolstart = '10';
my $repolmid   = '23';
my $repolend   = '10';

#################################
# I. get clone out of databases #
#################################
    
# NB at some time this should be changed to extract this information out of database('camace')
# instead of the primary('camace') and primary('stlace'), because someone might be changing the
# primary databases at the end of the build.

my $camace = $wormbase->primary('camace');
my $camdb  = Ace->connect(-path => $camace) || die "Couldn't connect to camace\n", Ace->error;

my @camclones = $camdb->fetch(-query => 'FIND Genome_Sequence');
foreach my $camclone (@camclones) {
  my $string = $camclone->Confidential_remark(1);
  if ((defined $string) && (($string =~ /Louis/) || ($string =~ /not in Cambridge/))) {
    next;
  }
  else {$camace{$camclone} = 1;}
}

my $stlace = $wormbase->primary('stlace');
my $stldb    = Ace->connect(-path => $stlace) || die "Couldn't connect to stlace\n", Ace->error;
my @stlclones = $stldb->fetch(-query => 'FIND Genome_Sequence');
foreach my $stlclone (@stlclones) {
  $stlace{$stlclone} = 1;
}



#############################################
# II. get data for one chromosome at a time #
#############################################

my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES

foreach my $chrom (@chrom) {
    print "\nProcessing chromosome $chrom\n";
    open (GFF, "$chromosomes_dir/CHROMOSOME_$chrom.gff") or die "cant open $chromosomes_dir/CHROMOSOME_$chrom.gff\n"; 
    #my $currentdb = $wormbase->database('currentdb');
    #open (GFF, "$currentdb/CHROMOSOMES/CHROMOSOME_$chrom.gff"); 
    open (CAMOL, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.overlapping_genes_cam") 
	|| die "Cannot open output file $chrom $!\n"; 
    open (CAMEST, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.EST_in_intron_cam")    
	|| die "Cannot open output file $chrom $!\n"; 
    open (CAMREP, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.repeat_in_exon_cam")  
	|| die "Cannot open output file $chrom $!\n"; 
    open (STLOL, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.overlapping_genes_stl") 
	|| die "Cannot open output file $chrom $!\n"; 
    open (STLEST, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.EST_in_intron_stl")    
	|| die "Cannot open output file $chrom $!\n"; 
    open (STLREP, ">$ace_dir/CHECKS/CHROMOSOME_$chrom.repeat_in_exon_stl")   
	|| die "Cannot open output file $chrom $!\n"; 
    
  %exon  = %est = %genes = %intron = %repeat = (); 
  my (%exoncount, %introncount, %estcount, %repeatcount, $name);
    my (%exon_orient, %intron_orient, %est_orient);

  while (<GFF>) {
    if (/^CHROM/) {
      my @fields = ();
      @fields = split(/\t/,$_);

      #################
      # get the exons #
      #################
      # %exon{AH6.1}        = [34454,36677]
      # %exon_orient{AH6.1} = +
      
      
      if ($fields[2] =~ /exon/ and $fields[1] =~ /curated/) {
	$fields[8] =~ /^CDS \"(\S+)\"/;
	$name = $1;
	$exoncount{$name}++; 
	my $exonname = $name.".".$exoncount{$name};
	$exon{$exonname} = [$fields[3],$fields[4]]; 
	$exon_orient{$name} = $fields[6];
    }   
      
      ###################
      # get the introns #
      ###################
      
      if ($fields[2] =~ /intron/ and $fields[1] =~ /curated/) {
	$fields[8] =~ /^CDS \"(\S+)\"/;
	$name = $1;
	$introncount{$name}++; 
	my $intronname = $name.".".$introncount{$name};
	$intron{$intronname} = [$fields[3],$fields[4]]; 
	$intron_orient{$name} = $fields[6];
      }
      
      ############################################################################
      # get the ESTs above a certain thresholdthreshold and certain match length #
      ############################################################################
      
      elsif (($fields[1] =~ /^BLAT_EST_BEST/ && $fields[5] > $ESTth) && (($fields[4] - $fields[3]) > $ESTlength)) {
	$fields[8] =~ /Sequence:(\S+)\"/;
	$name = $1;
	$estcount{$name}++;
	my $estname = $name.".".$estcount{$name};
	my @names = split (/ /, $fields[8]);
	$est{$estname} = [$fields[3],$fields[4],$names[2],$names[3]];  
	$est_orient{$name} = $fields[6];
      } 
      
      ###################
      # get the repeats #
      ###################

      
      elsif (($fields[5] ne ".") && ($fields[1] eq"RepeatMasker") && ($fields[5] > $repth)) {
	my @descr = split (/ /);
	$descr[1] =~ /\"Motif:(\S+)\"/; 
	$name = $1;
	$repeatcount{$name}++;
	my $repeatname = $name.".".$repeatcount{$name};
	$repeat{$repeatname} = [$fields[3],$fields[4]];  
      } 
    }
  }        
  
  #########################
  # make exons into genes #
  #########################
  
  foreach $name (sort keys %exoncount) {
    my $v = $exoncount{$name};
    my $w = "$name.$v";
    $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
  }
  
  #################
  # hunt isoforms #
  #################
  
  my %iso = ();
  my @alphabet = ("a" .. "z");
  foreach $name (sort keys %genes) {
    if ($name =~ /\da$/) {
      my ($newname) = ($name =~ /(\S+\.\d+)a$/);
      foreach my $letter (sort @alphabet) {
	my $k = $newname.$letter;
	if (exists $genes{$k}) {
	  push @{$iso{$name}}, $k; 
	}
	else {
	  last;
	}
      }    
    }
  }
  
  
  

  #####################
  # get "index" lists #
  #####################
  
  my @exonlist    = sort { ${$exon{$a}}[0]     <=> ${$exon{$b}}[0]   || $a cmp $b  } keys %exon;
  my @estlist        = sort { ${$est{$a}}[0]        <=> ${$est{$b}}[0]    || $a cmp $b  } keys %est;
  my @genelist    = sort { ${$genes{$a}}[0]    <=> ${$genes{$b}}[0]  || $a cmp $b  } keys %genes;
  my @repeatlist = sort { ${$repeat{$a}}[0]   <=> ${$repeat{$b}}[0] || $a cmp $b  } keys %repeat; 
  my @intronlist  = sort { ${$intron{$a}}[0]   <=> ${$intron{$b}}[0] || $a cmp $b  } keys %intron; 


  ###############################
  # III. find overlapping genes #
  ###############################
    
    my @geneoutput = ();
    for (my $x = 0; $x < @exonlist; $x++) {
	for (my $y = $x + 1; $y < @exonlist; $y++) {
	    my $name  = $exonlist[$x];
	    my $other = $exonlist[$y];
	    undef my $one; undef my $two;
	    ($one) = ($name  =~ /(\S.*?\.\d+)[a-z]\./);
	    ($two) = ($other =~ /(\S.*?\.\d+)[a-z]\./);
	    if  (!$one)  {$one = "n/a";}
	    if  (!$two)  {$two = "n/b";}
	    if  (($name eq $other) || ($one eq $two)) {                          # this should shortcut self matches!
		next;
	    }
	    elsif (${$exon{$name}}[1] >= ${$exon{$other}}[0])  {
		(my $onename)    = ($name  =~ /(\S.*?\.\S.*?)\./);
		(my $othername)  = ($other =~ /(\S.*?\.\S.*?)\./);
		push (@geneoutput, [sort ($onename,$othername)] );
	    }   
	    else {last;}
	}
    }
    
 # output for camace
    my @camgenes         = find_database(\@geneoutput,\%camace);   
    my %finalgenescamace = sort_by_gene(\@camgenes);
    warn "problem with genematch  reference\n" if ref(\@camgenes) ne 'ARRAY';
    foreach my $pair (sort keys %finalgenescamace) {
        my @single = split (/:/, $pair); 
        print CAMOL "Gene $single[1]\toverlaps with gene $single[0]\n" if ($single[0] ne $single[1]);
    } 
	
    # output for stlace    
    my @stlgenes         = find_database(\@geneoutput,\%stlace);
    my %finalgenesstlace = sort_by_gene(\@stlgenes);
    warn "problem with genematch  reference\n" if ref(\@stlgenes) ne 'ARRAY';
    foreach my $pair (sort keys %finalgenesstlace) {
	my @single = split (/:/, $pair); 
	print STLOL "Gene $single[1]\toverlaps with gene $single[0]\n" if ($single[0] ne $single[1]);;
    } 



	  
##################################
# IV. find ESTs matching introns #
##################################
	
  my $lastfail  = 0;
  my @intronoutput = ();
  
  for (my $x = 0; $x < @estlist; $x++) {
      my $testest = $estlist[$x];
      
      for (my $y = $lastfail; $y < @genelist; $y++) {
	  my $testgene  = $genelist[$y];
	  
	  if (!defined $introncount{$testgene}) {
	      next;
	  }
	  elsif ($est{$testest}->[0] > $genes{$testgene}->[1]){
	      $lastfail = $y;
	      next;
	  }
	  elsif ($est{$testest}->[1] < $genes{$testgene}->[0]){
	      last;
	  }
	  else {
	      for (my $z = 1; $z <= $introncount{$testgene}; $z++) {
		  my $intronstart = $intron{"$testgene.$z"}->[0];
		  my $intronend   = $intron{"$testgene.$z"}->[1];
		  my $eststart  = $est{$testest}->[0];
		  my $estend    = $est{$testest}->[1];
		  if ( not (($eststart > $intronend) || ($estend < $intronstart))) {
		      my ($finalest)  = ($testest =~ /^(\S+)\.\d+/);  	

		      # report to log file
		      print "\nLogging error for $finalest which extends into intron $intronstart - $intronend for $testgene\n";
		      print "EST runs $est_orient{$finalest}\tCDS runs $exon_orient{$testgene}\n";
		      print "EST is stored as a $estorientation{$finalest} read\n";
		      
		      # throw away matches on the opposite strand
		      next if (($exon_orient{$testgene} eq "+") && ($estorientation{$finalest} eq "5") && ($est_orient{$finalest} ne "+"));
		      next if (($exon_orient{$testgene} eq "+") && ($estorientation{$finalest} eq "3") && ($est_orient{$finalest} ne "-"));
		      next if (($exon_orient{$testgene} eq "-") && ($estorientation{$finalest} eq "5") && ($est_orient{$finalest} ne "-"));
		      next if (($exon_orient{$testgene} eq "-") && ($estorientation{$finalest} eq "3") && ($est_orient{$finalest} ne "+"));
			  
		      print "data ok - push to array\n";

		      # push to array
		      push (@intronoutput, [$finalest,$testgene,$intronstart,$intronend]);
		  }


	      }
	  }
      }   
  }

    # output for camace
    my @camintrons        = find_database(\@intronoutput,\%camace);
    my %finalintroncamace = sort_by_gene(\@camintrons); 
    my %camest = ();
    foreach my $pair (sort keys %finalintroncamace) {
        my @single = split (/:/, $pair);
        push @{$camest{$single[1]}}, "$single[0] $single[2] - $single[3]";
    }   
    my %camout = erase_isoforms(\%camest,\%iso);
    foreach my $x (sort keys %camout) {
        @{$camout{$x}} = sort @{$camout{$x}};
        print CAMEST "Introns of gene $x match ESTs @{$camout{$x}}\n"; 
    }           

    # output for stlace
    my @stlintrons        = find_database(\@intronoutput,\%stlace);
    my %finalintronstlace = sort_by_gene(\@stlintrons); 
    my %stlest = ();
    foreach my $pair (sort keys %finalintronstlace) {
        my @single = split (/:/, $pair);
        push @{$stlest{$single[1]}}, "$single[0] $single[2] - $single[3]";
    }   
    my %stlout = erase_isoforms(\%stlest,\%iso);
    foreach my $x (sort keys %stlout) {
        @{$stlout{$x}} = sort @{$stlout{$x}};
        print STLEST "Introns of gene $x match ESTs @{$stlout{$x}}\n"; 
    }           



	
##################################
# V. find repeats matching exons #
##################################
	
    my $lastrepeatfail  = 0;
    my @repeatoutput    = ();
        
    for (my $x = 0; $x < @repeatlist; $x++) {
        my $testrepeat = $repeatlist[$x];
        for (my $y = $lastrepeatfail; $y < @genelist; $y++) {
            my $testgene  = $genelist[$y];
            if (!defined $exoncount{$testgene}) {
                next;
            }
            elsif ($repeat{$testrepeat}->[0] > $genes{$testgene}->[1]){
                $lastrepeatfail = $y;
                next;
            }
            elsif ($repeat{$testrepeat}->[1] < $genes{$testgene}->[0]){
                last;
            }
            else {
                for (my $z = 1; $z <= $exoncount{$testgene}; $z++) {
                    my $exonstart = $exon{"$testgene.$z"}->[0];
		    my $exonend   = $exon{"$testgene.$z"}->[1];
		    my $repeatstart  = $repeat{$testrepeat}->[0];
		    my $repeatend    = $repeat{$testrepeat}->[1];
		    if (($repeatstart < $exonstart) && ($repeatend > $exonstart)) {
                        my $overlap = ($repeatend - $exonstart + 1);
                        my ($finalrepeat)  = ($testrepeat =~ /^(\S+)\.\d+/);
                        push (@repeatoutput, [$finalrepeat,$testgene]) if ($overlap >= $repolstart);
                    }
                    elsif (($repeatstart > $exonstart) && ($repeatend < $exonend)) {
                        my $overlap = ($repeatend - $repeatstart + 1);
                        my ($finalrepeat)  = ($testrepeat =~ /^(\S+)\.\d+/);
                        push (@repeatoutput, [$finalrepeat,$testgene]) if ($overlap >= $repolmid);
                    }
                    elsif (($repeatstart < $exonend) && ($repeatend > $exonend)) {
                        my $overlap = ($exonend - $repeatstart + 1);
                        my ($finalrepeat)  = ($testrepeat =~ /^(\S+)\.\d+/);
                        push (@repeatoutput, [$finalrepeat,$testgene]) if ($overlap >= $repolend);
                    }
                }
            }
        }   
    }
    
    # output for camace
    my @camrepeats      = find_database(\@repeatoutput,\%camace);
    my %finalrepcamace  = sort_by_gene(\@camrepeats);
    my %camrepout       = ();
    foreach my $pair (sort keys %finalrepcamace) {
        my @single = split (/:/, $pair); 
        push @{$camrepout{$single[0]}}, $single[1];
    }
    foreach my $gene (sort keys %camrepout) {
        print CAMREP "$gene matches @{$camrepout{$gene}}\n";
    }
    
    # output for stlace
    my @stlrepeats      = find_database(\@repeatoutput,\%stlace);
    my %finalrepstlace  = sort_by_gene(\@stlrepeats);
    my %stlrepout          = ();
    foreach my $pair (sort keys %finalrepstlace) {
        my @single = split (/:/, $pair); 
        push @{$stlrepout{$single[0]}}, $single[1];
    }
    foreach my $gene (sort keys %stlrepout) {
        print STLREP "$gene matches @{$stlrepout{$gene}}\n";
    }
}


# Tidy up
$camdb->close;
$stldb->close;

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


###################
# VI. subroutines #
###################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################



sub find_database {
    
    my @messy = @{$_[0]};
    my %ace= %{$_[1]};     
    my @output;
    carp "find_database not called with references\n" if (ref($_[0]) ne 'ARRAY' || ref($_[1]) ne 'HASH');
    
    foreach my $testpair (@messy) {
        my ($gene) = ($testpair->[1] =~ /(\S+)\./);
        carp "Second element of pair submitted to find_database looks like $gene\n" if !defined $gene;
        push (@output, $testpair) if (exists ($ace{$gene}));
    }
    return @output;   
}	

#################

sub erase_isoforms {
    my %dbout = ();
    my %dbest = %{$_[0]};
    my %iso   = %{$_[1]};    
    foreach my $testest (sort keys %dbest) {
        for (my $x = 0; $x < @{$dbest{$testest}}; $x++) {
            my $testgene = $dbest{$testest}->[$x];
            if ($testgene =~ /\da$/) {
                my $isocount = 0;
                for (my $y = 0; $y < @{$iso{$testgene}}; $y++) {
                    for (my $z = 0; $z < @{$dbest{$testest}}; $z++) {
                        if ($iso{$testgene}->[$y] eq $dbest{$testest}->[$z]) {
                            $isocount++;
                        }
                    }     
                }
                if ($isocount == @{$iso{$testgene}}) {
                    my $number = length @{$iso{$testgene}};
                    push @{$dbout{$testgene}}, $testest; 
                }
            }
            elsif ($testgene =~ /\d$/) {
                push @{$dbout{$testgene}}, $testest;    
            }
        }
    }
    return %dbout;
}

#################

sub sort_by_gene {  

    my @output = @{$_[0]};
    my %final_output;
	
    foreach my $out (@output) {
    	my $both  = $out->[1].":".$out->[0].":".$out->[2].":".$out->[3] if ($out->[0] && $out->[1] && $out->[2]); # dirty fix to skip noninitialized references
        $final_output{$both}++,
    }
    return %final_output;
}



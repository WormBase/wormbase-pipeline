#!/usr/local/bin/perl5.6.0 
#
# checks gene models for EST and cDNA evidence and adds a tag to the sequence object if the whole model is confirmed
# 010905 by Kerstin Jekosch

use Carp;
use IO::Handle;
use Getopt::Std;
use Wormbase;
$|=1;

# variables 
$opt_e="";   # check ESTs
$opt_m="";   # check cDNAs
$opt_h="";   # help
$opt_v="";   # verbose mode
my @chrom = qw(I);
my (%ESTs,%ESTlist,%genes,%genelist,$todo,%recall);
#my $WS = &get_wormbase_version();
my $wormbase = 'WS58'; 
getopts ('emhv');
&printhelp() unless ($opt_e || $opt_m);
&printhelp() if ($opt_e && $opt_m);
&printhelp() if ($opt_h);
$todo = 'BLAT_EST_BEST' if ($opt_e);
$todo = 'BLAT_mRNA_BEST' if ($opt_m); 

# directories
my $gffdir = "/wormsrv2/autoace/GFF_SPLITS/$wormbase";
my $outdir = "/nfs/griffin/kj2/testruns";

#############
# get genes #
#############

foreach my $chrom (@chrom) {

    my @refgen = &get_gec('CDS_exon',$chrom);
    print "Get exons out of $gffdir/CHROMOSOME_$chrom.CDS_exon.gff\n" if ($opt_v);
    %genes     = %{$refgen[0]}; 
    @genelist  = @{$refgen[1]};
    print "Got the genes for chrom $chrom\n" if ($opt_v);

    ########################
    # get EST and/or cDNAs #
    ########################

    my @refest = &get_gec($todo,$chrom);
    print "Get $todo out of $gffdir/CHROMOSOME_$chrom.$todo.gff\n" if ($opt_v);
    %ESTs      = %{$refest[0]};
    @ESTlist   = @{$refest[1]};
    print "Got the ESTs for chrom $chrom\n" if (($todo =~ /EST/) && $opt_v);
    print "Got the cDNAs for chrom $chrom\n" if (($todo =~ /mRNA/) && $opt_v);

    #######################
    # compare coordinates #
    #######################

    my $lastfail  = 0;
    my @output = ();
    print "Start searching" if ($opt_v);

    ##################
    # loop over ESTs #
    ##################

    for (my $x = 0; $x < @ESTlist; $x++) {

        # progress indicator
        print ".";

        my $testest   = $ESTlist[$x];
        my $estblocks = (scalar @{$ESTs{$testest}})-1;  # -1 for array pos last exon
        my $eststart  = $ESTs{$testest}->[0][1];
        my $estend    = $ESTs{$testest}->[$estblocks][2];

        ###################
        # loop over genes #
        ###################

        for (my $y = $lastfail; $y < @genelist; $y++) {
            my $testgene   = $genelist[$y];
            my $geneblocks = (scalar @{$genes{$testgene}})-1; # -1 for array pos last exon
            my $genestart = $genes{$testgene}->[0][1];
            my $geneend   = $genes{$testgene}->[$geneblocks][2];

            # next gene if gene is left of EST 
            if ($eststart > $geneend) {
                $lastfail = $y;
                next;
            }

            # next EST if gene is right of EST
            elsif ($estend < $genestart) {
                last;
            }

            ############################
            # EST matches gene somehow #
            ############################

            else {

                # loop over EST exons
                for (my $v = 0; $v <= $estblocks; $v++) {

                    my $eexstart  = $ESTs{$testest}->[$v][1];
                    my $eexend    = $ESTs{$testest}->[$v][2];

                    # loop over gene exons
                    for (my $w = 1; $w < $geneblocks; $w++) {

                        $gexstart = $genes{$testgene}->[$w][1];
                        $gexend   = $genes{$testgene}->[$w][2];

                        ####################
                        # WHAT WE ALL WANT #
                        ####################

                        if (($gexstart = $eexstart) && ($gexend = $eexend)) {
                            &confirm($testgene,$eexstart,$eexend);
                        }

                        #################
                        # SPECIAL CASES #
                        #################

                        # $v = 0 (eststart)
                        elsif ($v == 0) {
                            if (($eexstart >= $gexstart) && ($eexstart < $gexend) && ($eexend == $gexend)) {
                                &confirm($testgene,$eexstart,$eexend);
                            }
                            # $v = 0 AND $v = $estblocks (single exon EST)
                            elsif ($v == $estblocks) {
                                if (($eexstart >= $gexstart) && ($eexend <= $gexend)) {
                                    &confirm($testgene,$eexstart,$eexend);
                                }
                            }
                        }
                        # $v = $estblocks (estend)
                        elsif ($v == $estblocks) {
                            if (($eexstart == $gexstart) && ($eexend > $gexstart) && ($eexend <= $gexend)) {
                                &confirm($testgene,$eexstart,$eexend);
                            }
                        }
                        # $w = 0 (genestart)
                        elsif ($w == 0) {
                            if (($gexstart >= $eexstart) && ($gexstart < $eexend) && ($gexend == $eexend)) {
                                &confirm($testgene,$eexstart,$eexend);
                            }
                            # $w = 0 AND $w = $geneblocks (single exon gene)
                            elsif ($w = $geneblocks) {
                                if (($gexstart >= $eexstart) && ($gexend <= $eexend)) {
                                   &confirm($testgene,$eexstart,$eexend); 
                                }
                            }
                        }
                        # $w = $geneblocks (geneend) 
                        elsif ($w == $geneblocks) {
                            if (($gexstart == $eexstart) && ($gexend > $eexstart) && ($gexend <= $eexend)) {
                                &confirm($testgene,$eexstart,$eexend);
                            }
                        }   
                    }    
                }
            }        
        }
    }
    print "\n" if ($opt_v);

    ##############################################
    # check for confirmation and produce ouutput #
    ##############################################

    open(ACE,">$outdir/chrom_$chrom.confirmed_$todo.ace");
    print "Writing to outfile\n" if ($opt_v);
    foreach my $gene (@genelist) {
        my $true = 1;
        my $exons = scalar @{$genes{$gene}};
        for (my $c = 0; $c < $exons; $c++) {
            my $first = $genes{$gene}->[$c][1];
            my $last  = $genes{$gene}->[$c][2];
            for (my $d = $first; $d <= $last; $d++) {
                if (!exists $recall{$gene}->{$d}) {
                    $true = 0;
                    last;
                }
            }
        }
        if ($true == 1) {
            print ACE "Sequence : \"$gene\"\n";
            print ACE "Confirmed_by EST\n\n" if ($todo =~ /EST/);
            print ACE "Confirmed_by cDNA\n\n" if ($todo =~ /mRNA/);
        }    
    }
}



################
#     subs     #
################

sub get_gec {
    my $filename = shift;
    my $chrom    = shift; 
    my $name     = "";
    my %hash     = ();
    open (GFF, "$gffdir/CHROMOSOME_$chrom.$filename.gff") or die "Cannot open $gffdir/CHROMOSOME_$chrom.$filename.gff $!\n";
    while (<GFF>) {
        next if /\#/;
        my @f = split /\t/;
        my $name;
        ($name) = ($f[8] =~ /\"(\S+)\"/) if ($filename eq 'CDS_exon');            
        ($name) = ($f[8] =~ /\"Sequence:(\S+)\"/) if ($filename ne 'CDS_exon');
        my $first  = 1;
        my $second = 2;
        my @g = split /\s/, $f[8];
        $first  = $g[2] if ($g[2]);
        $second = $g[3] if ($g[3]);
        my $strand = '+';
        if ($first > $second) {$strand = '-';}
        push @{$hash{$name}} , [$strand,$f[3],$f[4],$first,$second]; 
    }
    @list = sort { $hash{$a}->[0][1] <=> $hash{$b}->[0][1] || $a cmp $b } keys %hash;  
    return (\%hash,\@list); 
}

################

sub confirm  {
    my $gene  = shift;
    my $start = shift;
    my $end   = shift;
    for (my $n = $start; $n <= $end; $n++) {
        $recall{$gene}->{$n} = 1;
    } 
}  

################

sub printhelp {
    exec ('perldoc',$0);
}    

__END__

=pod

=head2   NAME - confirm_genes.pl

=head1 USAGE

=over 4

=item  confirm_genes.pl -options

=back

confirm_genes.pl parses GFF files and looks for gene models that are completely supported by 
EST or cDNA matches, depending on your choice. These gene models get an appropriate tag.

mandatory arguments:

=over 4

=item -e for checking for EST matches

=item or

=item -m for checking for cDNA matches 

=back

optional arguments

=over4

=item -v for being verbose

=item -h for help

=back

=cut




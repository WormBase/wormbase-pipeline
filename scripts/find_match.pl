#!/usr/local/bin/perl -w
#
# find_match.pl
# checks whether a chosen feature overlaps gene introns/exons                                   
# sorts output for stl and cam clones
#
# by Kerstin Jekosch
# 10/08/01


use Carp;
use Ace;
use IO::Handle;
use Getopt::Long;
$|=1;

my @chrom = qw(III);
my (%structure, %feature, %genes, %camace, %stlace);

my $structure="";
my $feature="";
my $dir="";
my $h="";

GetOptions (
	    "structure=s"   => \$structure,
	    "feature=s"     => \$feature,
	    "dir=s"    => \$dir,
	    "h"       => \$h,
);

length ($structure)==0 && &PrintHelp;
length ($feature)==0 && &PrintHelp;
length ($dir)==0 && &PrintHelp;

print $dir, "\n";
print $feature, "\n";
print $structure, "\n";

# get clone out of databases 
my $camdb     = Ace->connect(-path => '/wormsrv2/camace/') || die "Couldn't connect to camace\n", Ace->error;
my @camclones = $camdb->fetch(-query => 'FIND Genome_Sequence');
foreach my $camclone (@camclones) {
	my $string = $camclone->Confidential_remark(1);
	if ((defined $string) && (($string =~ /not in Cambridge LINK/) || ($string =~ /Louis/))) {
		next;
	}
	else {$camace{$camclone} = 1;}
}

my $stldb     = Ace->connect(-path => '/wormsrv2/stlace/') || die "Couldn't connect to stlace\n", Ace->error;
my @stlclones = $stldb->fetch(-query => 'FIND Genome_Sequence');
foreach my $stlclone (@stlclones) {
    $stlace{$stlclone} = 1;
}



# get data for one chromosome at a time 
foreach my $chrom (@chrom) {

    open (GFF, "/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_$chrom.gff"); 
    open (CAMOUT, ">$dir/CHROMOSOME_$chrom.$feature_cam") 
	|| die "Cannot open output file CAMOUT $chrom $!\n"; 
    open (STLOUT, ">$dir/CHROMOSOME_$chrom.$feature_stl") 
	|| die "Cannot open output file STLOUT $chrom $!\n"; 
		
    %exon = %intron = %feature = %genes = (); 
    my (%exoncount, %introncunt, %featurecount, $name);
    
    while (<GFF>) {
        if (/^CHROM/) {
            my @fields = ();
            @fields = split(/\t/,$_);

            # get the exons/introns 
            if ($fields[2] =~ /exon/ and $fields[1] =~ /curated/) {
                $fields[8] =~ /\"(\S+)\"/;
                $name = $1;
                $exoncount{$name}++; 
                my $exonname = $name.".".$exoncount{$name};
                $exon{$exonname} = [$fields[3],$fields[4]]; 
            }   
            
            # get the introns 
	    if ($fields[2] =~ /intron/ and $fields[1] =~ /curated/) {
                $fields[8] =~ /^Sequence \"(\S+)\"/;
                $name = $1;
                $introncount{$name}++; 
                my $intronname = $name.".".$introncount{$name};
                $intron{$intronname} = [$fields[3],$fields[4]]; 
            }
    
            # get the features 
            elsif ($fields[1] =~ /$feature/) {
                $fields[8] =~ /\"(\S+)\"/;
                $name = $1;
                $featurecount{$name}++;
                my $featurename = $name.".".$featurecount{$name};
                my @names = split (/ /, $fields[8]);
                $feature{$featurename} = [$fields[3],$fields[4],$names[2],$names[3]];  
            } 
        }        
    }    

    # make exons into genes 
    foreach $name (sort keys %exoncount) {
	my $v = $exoncount{$name};
	my $w = "$name.$v";
        $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
    }

    # get "index" lists
    my @exonlist   = sort { ${$exon{$a}}[0]     <=> ${$exon{$b}}[0]   || $a cmp $b  } keys %exon;
    my @intronlist = sort { ${$intron{$a}}[0]   <=> ${$intron{$b}}[0] || $a cmp $b  } keys %intron;    
    my @featurelist= sort { ${$feature{$a}}[0]  <=> ${$feature{$b}}[0]|| $a cmp $b  } keys %feature;
    my @genelist   = sort { ${$genes{$a}}[0]    <=> ${$genes{$b}}[0]  || $a cmp $b  } keys %genes;
    my @repeatlist = sort { ${$repeat{$a}}[0]   <=> ${$repeat{$b}}[0] || $a cmp $b  } keys %repeat; 


	

    # find feature matching structure 
    my $lastfail      = 0;
    my @featureoutput = ();
    
    if ($structure =~ /exon/) {
        @structurelist = @exonlist;
        %structure     = %exon;
        %structurecount= %exoncount;
    }
    elsif ($structure =~ /intron/) {
        @structurelist = @intronlist;
        %structure     = %intron;
        %structurecount= %introncount;
    }
    else { 
        die "Structure definition missing!\n";
    } 
        
    for (my $x = 0; $x < @featurelist; $x++) {
        my $testfeature = $featurelist[$x];
        for (my $y = $lastfail; $y < @genelist; $y++) {
            my $testgene  = $genelist[$y];
            if (!defined $structurecount{$testgene}) {
                next;
            }
            elsif ($feature{$testfeature}->[0] > $genes{$testgene}->[1]){
                $lastfail = $y;
                next;
            }
            elsif ($feature{$testfeature}->[1] < $genes{$testgene}->[0]){
                last;
            }
            else {
                for (my $z = 1; $z <= $structurecount{$testgene}; $z++) {
                    my $genestart = $structure{"$testgene.$z"}->[0];
		    my $geneend   = $structure{"$testgene.$z"}->[1];
		    my $featurestart  = $feature{$testfeature}->[0];
		    my $featureend    = $feature{$testfeature}->[1];
		    if ( not (($featurestart > $geneend) || ($featureend < $genestart))) {
                       my ($finalfeature)  = ($testfeature =~ /^(\S+)\.\d+/);  	
		        push (@featureoutput, [$finalfeature,$testgene]);
	            }
                }
            }
        }   
    }
    
    # print output
    my @camfeature         = &find_database(\@featureoutput,\%camace);
    my %finalfeaturecamace = &sort_by_gene(\@camfeature);
    
    foreach my $pair (sort keys %finalfeaturecamace) {
        my @single = split (/:/, $pair); 
        print CAMOUT "$single[1]\tmatches gene $single[0]\n";
    }

    my @stlfeature         = &find_database(\@featureoutput,\%stlace);
    my %finalfeaturestlace = &sort_by_gene(\@stlfeature);
    
    foreach my $pair (sort keys %finalfeaturestlace) {
        my @single = split (/:/, $pair); 
        print STLOUT "$single[1]\tmatches gene $single[0]\n";
    }
}



    
#################

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

sub sort_by_gene {  

    my @output = @{$_[0]};
    my %final_output;
	
    foreach my $out (@output) {
    	my $both  = $out->[1].":".$out->[0];
        $final_output{$both}++,
    }
    return %final_output;
}

#################

sub PrintHelp {
    exec ('perldoc',$0);
}



__END__

=pod

=head2   NAME - find_match.pl

=head1 USAGE

=over 4

=item  find_match.pl -f <feature> -s <intron/exon> -d <dir>

=back

find_match parses GFF files, compares the given feature to the selected structure (introns or exons) 
and reports overlaps.

autoace_minder mandatory arguments:

=over 4

=item -f feature (takes second field GFF definitions, e.g. EST_GENOME for ESTs) 

=item -s structure to compare to (exon or intron)

=item -d output directory

=back

=cut

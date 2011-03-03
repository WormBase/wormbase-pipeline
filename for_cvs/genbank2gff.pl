#!/usr/bin/perl -w

# minimally changed from scripts/Bio-DB-GFF/genbank2gff.pl
# it is a GFF2/3 like output format:
#   * tags are not separated by semicolons
#   * no parent / child IDs

# get bioperl from here: https://github.com/bioperl/bioperl-live.git

use lib '/tmp/bioperl-live';
use Bio::SeqIO;

@ARGV = qw(-) unless @ARGV;

# user-specifiable params
my $use_source = 0;
my $pr_comment = 0;
my $pr_trans = 0;
my $no_white = 0;
my $no_tags = 0;

# fixed params
my $width = 50;
my $default_source = "GenBank";

my $usage = "";
$usage .= "Usage: $0\n";
$usage .= "       [-usesource]  (use bioperl 'source' tag, rather than '$default_source')\n";
$usage .= "       [-comment]    (also print feature table as indented comments)\n";
$usage .= "       [-trans]      (print translation tags)\n";
$usage .= "       [-nowhite]    (suppress text-with-whitespace tags)\n";
$usage .= "       [-notags]     (suppress ALL tags)\n";
$usage .= "       <GenBank files...>\n";

my @argv;
while (@ARGV) {
    my $opt = shift;
    unless ($opt =~ /^-./) { push @argv, $opt; next }
    if ($opt eq '-usesource') { $use_source = 1 }
    elsif ($opt eq '-comment') { $pr_comment = 1 }
    elsif ($opt eq '-trans') { $pr_trans = 1 }
    elsif ($opt eq '-nowhite') { $no_white = 1 }
    elsif ($opt eq '-notags') { $no_tags = 1 }
    else { die $usage . "Unknown option: $opt\n" }
}
die $usage unless @argv == 1;

foreach my $filename (@argv) {
    my $stream = Bio::SeqIO->new (-file => $filename, -format => 'GenBank');
    while (my $seq = $stream->next_seq) {
    	if ($pr_comment) {
    	    print "\# Features for ", $seq->display_id, "\n";
    	    print "\# ", $seq->desc, "\n";
    	    comment_sf ("", $seq->top_SeqFeatures);
    	}
    	foreach my $sf ($seq->all_SeqFeatures) {
    	    my @tags = $sf->all_tags;
    	    @tags = grep (!/translation/i, @tags) unless $pr_trans;
    	    my @tagval;
    	    my $description;
    	    unless ($no_tags) {
        		@tagval = map { $tag = $_; map ("$tag=$_", $sf->each_tag_value ($tag)) } @tags;
        		if ($no_white) {
        		    @tagval = grep (!/\s/, @tagval);
        		} else {
        		    grep (s/=(.*\s.*)/=\'$1\'/, @tagval);  # quote whitespace
        		}
        		
#        		map {s/=/ /}@tagval; #GFF2ify it
        		
        	    # if codon_start != 1 change the phase
        	    $description=join (" ", @tagval);
        	    if ($description=~/codon_start=(\d)/){
    	            $sf->frame($1-1);
        	    }
    	    }
    	    
    	    
    	    my $source = $use_source ? $sf->source_tag : $default_source;
    	    my $score = $sf->score || '.';
            my $phase = $sf->frame;
            my $strand = ($sf->strand > 0)?'+':'-';
    	    my @gff = ($seq->display_id, $source, $sf->primary_tag, $sf->start, $sf->end, $score, $strand, $phase, $description);
    	    # iron out a few bioperl wrinkles
    	    $gff[5] = '+' unless defined $gff[5];  # strand
    	    @gff = map (defined($_) ? $_ : ".", @gff);  # everything else

    	    # display
    	    print join ("\t", @gff), "\n";

            # coding_exon bit
            if ($sf->location->isa('Bio::Location::SplitLocationI') && $sf->primary_tag eq 'CDS'){
               my @subFeature = sort{if($sf->strand <0){$b->start <=> $a->start}else{$a->start <=> $b->start}} $sf->location->sub_Location;
               $description=~s/codon_start=(\d)\s*//;
               for my $loc (@subFeature){
     	         my @gff = ($seq->display_id, $source,'coding_exon', $loc->start, $loc->end, $score, $strand,$phase, $description);
    	         # iron out a few bioperl wrinkles
    	         $gff[5] = '+' unless defined $gff[5];  # strand
    	         @gff = map (defined($_) ? $_ : ".", @gff);  # everything else
                 print join ("\t",@gff),"\n";
                 $phase= ($loc->end - $loc->start - $phase +1)%3;
               }
            }
            elsif ($sf->primary_tag eq 'CDS'){
                 $description=~s/codon_start=(\d)\s*//;
                 my @gff = ($seq->display_id, $source,'coding_exon', $sf->start, $sf->end, $score, $strand,$phase, $description);
    	         # iron out a few bioperl wrinkles
    	         $gff[5] = '+' unless defined $gff[5];  # strand
    	         @gff = map (defined($_) ? $_ : ".", @gff);  # everything else
                 print join ("\t",@gff),"\n";
            }
    	}
    }
}

sub comment_sf {
    my ($indent, @sf) = @_;
    foreach my $sf (@sf) {
	print "\# ", $indent, "\"", $sf->primary_tag, "\" start=", $sf->start, " end=", $sf->end, "\n";
	    my @subsf = $sf->sub_SeqFeature;
	    if (@subsf) { comment_sf($indent."  ", @subsf) }
    }
}

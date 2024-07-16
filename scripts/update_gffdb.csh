#!/bin/csh
# small nasty thing to submit into lsf
# will update the GFF-database to current_db

    cd $CVS_DIR/Modules
    
if ($1 == '-autoace') then
foreach f(~wormpub/BUILD/autoace/CHROMOSOMES/CHROMOSOME*.gff)
	echo "... adding new tags from autoace..."
		perl -mGFF_sql -e 'my $a=GFF_sql->new({-build => 1});$a->generate_tags($ARGV[0])' $f
		perl $CVS_DIR/slurm_submit.pl perl -mGFF_sql -e 'my $a=GFF_sql->new({-build => 1}); foreach $i (@ARGV){$i=~/(CHROMOSOME_\w+)\.gff/;my $name=$1;$a->load_gff($i,$name,1)}' $f
	   
else
foreach f(~wormpub/DATABASES/current_DB/CHROMOSOMES/CHROMOSOME*.gff)
	echo "... adding new tags from current_DB..."
		perl -mGFF_sql -e 'my $a=GFF_sql->new({});$a->generate_tags($ARGV[0])' $f
		perl $CVS_DIR/slurm_submit.pl perl -mGFF_sql -e 'my $a=GFF_sql->new({}); foreach $i (@ARGV){$i=~/(CHROMOSOME_\w+)\.gff/;my $name=$1;$a->load_gff($i,$name,1)}' $f

	endif
end 

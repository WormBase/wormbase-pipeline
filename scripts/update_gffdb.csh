#!/bin/csh
# small nasty thing to submit into lsf
# will update the GFF-database to current_db

if (`uname` == "OSF1" ) then
	pushd ./
	cd $CVS_DIR/Modules
	foreach f(~wormpub/DATABASES/current_DB/CHROMOSOMES/CHROMOSOME*.gff)
		echo "... adding new tags ..."
		if ($1 == '-build') then
			perl -mGFF_sql -e 'my $a=GFF_sql->new({-build => 1});$a->generate_tags($ARGV[0])' $f
			bsub -J loaddb_build -o /dev/null -e /dev/null perl -mGFF_sql -e 'my $a=GFF_sql->new({-build => 1}); foreach $i (@ARGV){$i=~/(CHROMOSOME_\w+)\.gff/;my $name=$1;$a->load_gff($i,$name,1)}' $f
		else
			perl -mGFF_sql -e 'my $a=GFF_sql->new({});$a->generate_tags($ARGV[0])' $f
			bsub -J loaddb -o /dev/null -e /dev/null perl -mGFF_sql -e 'my $a=GFF_sql->new({}); foreach $i (@ARGV){$i=~/(CHROMOSOME_\w+)\.gff/;my $name=$1;$a->load_gff($i,$name,1)}' $f
		endif
	end
	cd =1
else
	echo "needs to be run from an alpha on LSF"
	exit 1
endif

/^DNA/d
/^Sequence/d
/^Assembly_Tags/d
/^Clone/d
/^Allele/d
/^;/d
/^>/d
s/<.*>//
s/<[^>]*$//
s/^ [^<]*>/ /
y/ACGT-N/acgtnn/
s/*//g
/^$/d

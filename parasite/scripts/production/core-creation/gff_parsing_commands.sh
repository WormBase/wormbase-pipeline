#Remove comment lines
function gff_nc {
  cat $1 | grep "^[^#]"
}

#Only get lines with genes
function gff_gene {
  gff_nc $1 | awk '$3~/gene/'
}

#Only get lines with transcripts
function gff_tr {
  gff_nc $1 | awk '$3~/([a-z]+_)?[a-z]RNA|([a-z]+)?_?transcript/'
}

#Get GFF scaffolds
function gff_scafs {
  gff_nc $1 | cut -f1 | sort | uniq
}

#Get GFF entries types.
function gff_types {
  gff_nc $1 | cut -f3 | sort | uniq
}

#Get attributes
function gff_attr {
  gff_nc $1 | cut -f9
}

#Get object IDs
function gff_id {
  perl -nle 'print for /(?:^|.+)ID=([^;]+)/;'
}

#Get object Name
function gff_name {
  perl -nle 'print for /(?:^|.+)Name=([^;]+)/;'
}

#Get longest
function gff_longest {
  awk 'length > max_length { max_length = length; longest_line = $0 } END { print longest_line }' $1
}

#Get specific attribute
function gff_defattr {
  perl -snle 'print for /(?:^|.+)$attrname=([^;]+)/;' -- -attrname=$1;
}

#Get all unique attribute flags
function gff_flagsq {
  perl -snle 'print for /(?:^|.+)$attrname=([^;]+)/;' -- -attrname=$1;
}






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

function fasta_sequence_from_header {
  header=$1
  fasta_file=$2
  awk -v header="$header" '
    BEGIN { found = 0; }
    {
        if ($0 ~ header) {
            found = 1;
            print $0;
        } else if (found) {
            if (substr($0, 1, 1) == ">") {
                exit;
            } else {
                print $0;
            }
        }
    }
' $fasta_file
}

function fasta_sequence_length_from_header {
  header=$1
  fasta_file=$2
  awk -v header="$header" '
      BEGIN { found = 0; sequence = ""; }
      {
          if ($0 ~ header) {
              found = 1;
              sequence = "";
              print $0;
          } else if (found) {
              if (substr($0, 1, 1) == ">") {
                  exit;
              } else {
                  sequence = sequence $0;
              }
          }
      }
      END {
          if (sequence != "") {
              seq_length = length(sequence);
              print "Sequence Length: " seq_length;
          }
      }
  ' $fasta_file
}

function rename_fasta_headers_with_mapping {
  input_fasta = $1
  mapping_file = $2
  while read p;
    do
      old=$(echo $p | cut -d' ' -f1);
      new=$(echo $p | cut -d' ' -f2);
      sed -i "s/^>${old}/>${new}/" ${input_fasta};
    done < ${mapping_file}
}

function rename_seq_region_file_with_mapping {
  seq_region_file=$1
  mapping_file=$2
  join -t $'\t' $seq_region_file $mapping_file -a1 -1 2 -2 1 \
    | awk '{print $2"\t"$5"\t"$3"\t"$4}'
}
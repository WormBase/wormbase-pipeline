#!/bin/csh -f

unlimit

switch ($#argv)
  case "1":
        breaksw
  default:
        echo 'Usage: dump_links ACEDB_directory' ; exit 1
endsw


cd $1

/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/giface . >! /tmp/chrom_dump.log <<EOF
find sequence SUPERLINK_CB_I
dna -f /tmp/SUPERLINK_CB_I.dna
find sequence SUPERLINK_CB_IR
dna -f /tmp/SUPERLINK_CB_IR.dna
find sequence SUPERLINK_CB_II
dna -f /tmp/SUPERLINK_CB_II.dna
find sequence SUPERLINK_CB_IIIL
dna -f /tmp/SUPERLINK_CB_IIIL.dna
find sequence SUPERLINK_CB_IIIR
dna -f /tmp/SUPERLINK_CB_IIIR.dna
find sequence SUPERLINK_CB_IV
dna -f /tmp/SUPERLINK_CB_IV.dna
find sequence SUPERLINK_CB_V
dna -f /tmp/SUPERLINK_CB_V.dna
find sequence SUPERLINK_CB_X
dna -f /tmp/SUPERLINK_CB_X.dna
EOF

cd /tmp
/bin/cat *.dna | /nfs/disk100/wormpub/bin.ALPHA/composition >! composition.all
/usr/local/bin/perl -ne 'if (/(\d+) total/) { $t = $1 ; print $t ; } if (/- (\d+)/) { $t -= $1 ; print " $t\n" ; }' composition.all >! totals


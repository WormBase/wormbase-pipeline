#!/bin/sh

# Ortholog assignment pipeline
#  Todd Harris (harris@cshl.org)
#  DATE : 05 May 2003

# Change these paths as appropriate to point to the blast output files
# Input files
CE_CB_SEG_ON=blastp-current/ce_vs_cb-seg-on.out
CB_CE_SEG_ON=blastp-current/cb_vs_ce-seg-on.out
CE_CB_SEG_OFF=blastp-current/ce_vs_cb-seg-off.out
CB_CE_SEG_OFF=blastp-current/cb_vs_ce-seg-off.out

# Output files
BEST_MUTUALS_SEG_ON=output/orthologs.best_mutuals.seg-on
BEST_MUTUALS_SEG_OFF=output/orthologs.best_mutuals.seg-off
BEST_MUTUALS=output/orthologs.best_mutuals
ORTHOLOG_SPANS=output/orthologs.ortholog.spans
ORTHOLOGS_100K=output/orthologs.100k_spans
ORTHOLOGS_100K=output/orthologs.200k_spans
ORTHOLOGS_ALL=output/orthologs.all

# Cutoffs, thresholds, etc
MAX_EVAL=1e-10
MEX_EVAL_SEG_OFF=1e-20
CUTOFF=100000

function get_start {
echo "Started at: "
date
}

function get_end {
echo "Finished at: "
date
}

# 1. Assign best mutuals from the seg-on blastp
echo "Finding best mutuals from the seg-on blastp data"
echo $(get_start)
./best_mutuals \
--blast1 ${CE_CB_SEG_ON} \
--blast2 ${CB_CE_SEG_ON} \
--max_eval ${MAX_EVAL} \
--cutoff ${CUTOFF} \
--method seg-on \
> ${BEST_MUTUALS_SEG_ON}
echo $(get_end)

# 1. Assign best mutuals from the seg-on blastp, including searching for orphans
#echo "Finding best mutuals from the seg-on blastp data"
#echo $(get_start)
#./best_mutuals \
#--blast1 ${CE_CB_SEG_ON} \
#--blast2 ${CB_CE_SEG_ON} \
#--orphans 1 \
#--fasta1 /Users/todd/projects/briggsae/data/gene_set-current/longest/ws77.pep2.00_longest_splices \
#--fasta2 /Users/todd/projects/briggsae/data/gene_set-current/cb25.hybrid.pep2.00 \
#--max_eval ${MAX_EVAL} \
#--cutoff ${CUTOFF} \
#--method seg-on \
#> ${BEST_MUTUALS_SEG_ON}
#echo $(get_end)

## 2. Assign best mutuals from the seg-off blastp
echo "Finding best mutuals from the seg-off blastp data"
echo $(get_start)
./best_mutuals \
--blast1 ${CE_CB_SEG_OFF} \
--blast2 ${CB_CE_SEG_OFF} \
--max_eval ${MAX_EVAL_SEG_OFF} \
--cutoff ${CUTOFF} \
--method seg-off \
--orthologs ${BEST_MUTUALS_SEG_ON} \
> ${BEST_MUTUALS_SEG_OFF}
cat ${BEST_MUTUALS_SEG_ON} ${BEST_MUTUALS_SEG_OFF} >> ${BEST_MUTUALS}
echo $(get_end)

# 4. Assign orthologs by ortholog-based spans
echo "Assigning orthologs by ortholog span"
echo $(get_start)
./find_orthos_by_synteny \
--method 'synteny-ortholog_spans' \
--cutoff ${CUTOFF} \
--max_eval ${MAX_EVAL} \
--orthologs ${BEST_MUTUALS} \
> ${ORTHOLOG_SPANS}
cat ${BEST_MUTUALS} ${ORTHOLOG_SPANS} >> ${ORTHOLOGS_ALL}
echo $(get_end)

# 5. Assign orthologs by 100kb spans
echo "Assigning orthologs by range span: 100kb"
echo $(get_start)
./find_orthos_by_synteny \
--method 'synteny-100kb_spans' \
--cutoff ${CUTOFF} \
--max_eval ${MAX_EVAL} \
--range 100000 \
--orthologs ${ORTHOLOGS_ALL} \
> ${ORTHOLOGS_100K}
cat ${ORTHOLOGS_100K_SPANS} >> ${ORTHOLOGS_ALL}
echo $(get_end)

# 6. Assign orthologs by 200kb spans
echo "Assigning orthologs by range span: 200kb"
echo $(get_start)
./find_orthos_by_synteny \
--method 'synteny-200kb_spans' \
--cutoff ${CUTOFF} \
--max_eval ${MAX_EVAL} \
--range 200000 \
--orthologs ${ORTHOLOGS_ALL} \
> ${ORTHOLOGS_200K}
cat ${ORTHOLOGS_200K} >> ${ORTHOLOGS_ALL}
echo $(get_end)

echo "finished..."

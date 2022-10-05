args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two arguments must be supplied (study_dir and study).n", call.=FALSE)
}
study_dir = args[1]
study = args[2]

brc4_unique = read.csv(paste0(study_dir,"/genes.htseq-union.stranded.sum.counts.tpm"), sep="\t", header=0)
rnaseqer = read.csv(paste0(study_dir,"/",study,".pe.genes.tpm.htseq2.irap.tsv"), sep="\t", header=0)

colnames(brc4_unique) = c("gene","value")
colnames(rnaseqer) = c("gene","rnaseqer")

old_cdf = merge(brc4_unique, rnaseqer, by.x="gene", by.y="gene", all.x=F, all.y=F)
cdf = transform(old_cdf, rnaseqer = as.numeric(rnaseqer))

head(cdf)

# png(
#   paste0(study_dir,"/",study,".comparison.png"),
#   width     = 3.25,
#   height    = 5,
#   units     = "in",
#   res       = 1200,
#   pointsize = 6
# )
pdf(
  paste0(study_dir,"/",study,".tpm_comparison_log2.pdf"),
  width     = 3.25,
  height    = 5,
  pointsize = 6
)
par(mfrow=c(1,1))
#plot(log(cdf$value_brc4), log(cdf$value_brc4_nu), pch = 21, main = study, xlab = "Metazoa unique", ylab = "Metazoa non-unique")
#plot(log(cdf$value_brc4_nu), log(cdf$rnaseqer), pch = 21, xlab = "Metazoa non-unique", ylab = "RNASEQer")
plot(log(cdf$value), log(cdf$rnaseqer), pch = 21, main = study, xlab = "e! Metazoa log(Gene TPM)", ylab = "RNASEQer log(Gene TPM)")
dev.off()
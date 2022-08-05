args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two arguments must be supplied (study_dir and study).n", call.=FALSE)
}
study_dir = args[1]
study = args[2]

brc4_unique = read.csv(paste0(study_dir,"/genes.htseq-union.unstranded.counts"), sep="\t", header=0)
brc4_nonunique = read.csv(paste0(study_dir,"/genes.htseq-union.unstranded.nonunique.counts"), sep="\t", header=0)
rnaseqer = read.csv(paste0(study_dir,"/",study,".pe.genes.raw.htseq2.tsv"), sep="\t", header=0)

colnames(brc4_unique) = c("gene","value")
colnames(brc4_nonunique) = c("gene","value")
colnames(rnaseqer) = c("gene","rnaseqer")

brc4 = merge(brc4_unique, brc4_nonunique, by.x="gene", by.y="gene", suffixes = c("_brc4","_brc4_nu"))
cdf = merge(brc4, rnaseqer, by.x="gene", by.y="gene", all.x=F, all.y=F)

# png(
#   paste0(study_dir,"/",study,".comparison.png"),
#   width     = 3.25,
#   height    = 5,
#   units     = "in",
#   res       = 1200,
#   pointsize = 6
# )
pdf(
  paste0(study_dir,"/",study,".comparison_log2.pdf"),
  width     = 3.25,
  height    = 5,
  pointsize = 6
)
par(mfrow=c(1,1))
#plot(log(cdf$value_brc4), log(cdf$value_brc4_nu), pch = 21, main = study, xlab = "Metazoa unique", ylab = "Metazoa non-unique")
#plot(log(cdf$value_brc4_nu), log(cdf$rnaseqer), pch = 21, xlab = "Metazoa non-unique", ylab = "RNASEQer")
plot(log(cdf$value_brc4), log(cdf$rnaseqer), pch = 21, main = study, xlab = "e! Metazoa log(raw gene counts)", ylab = "RNASEQer log(raw gene counts)")
dev.off()
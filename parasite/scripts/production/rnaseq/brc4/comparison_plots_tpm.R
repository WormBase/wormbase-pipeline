args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Three arguments must be supplied (brc4 tmp file, rnaseqer timp file and the run accession).\n", call.=FALSE)
}
brc4_tpm = args[1]
rnaseqer_tpm = args[2]
run = args[3]

brc4_unique = read.csv(brc4_tpm, sep="\t", header=0)
rnaseqer = read.csv(rnaseqer_tpm, sep="\t", header=0)

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
  paste0(brc4_tpm,".tpm_comparison_log2.pdf"),
  width     = 3.25,
  height    = 5,
  pointsize = 6
)
par(mfrow=c(1,1))
#plot(log(cdf$value_brc4), log(cdf$value_brc4_nu), pch = 21, main = study, xlab = "Metazoa unique", ylab = "Metazoa non-unique")
#plot(log(cdf$value_brc4_nu), log(cdf$rnaseqer), pch = 21, xlab = "Metazoa non-unique", ylab = "RNASEQer")
plot(log(cdf$value), log(cdf$rnaseqer), pch = 21, main = run, xlab = "e! Metazoa log(Gene TPM)", ylab = "RNASEQer log(Gene TPM)")
dev.off()
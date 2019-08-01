dir <- "airway2/inst/extdata"
sra <- read.delim(file.path(dir,"SraRunTable.txt"))
files <- file.path(dir,"quants",sra$Run,"quant.sf")
lvls <- c("Untreated", "Dexamethasone")
coldata <- data.frame(
  names=sra$Run,
  files=files,
  donor=sra$cell_line,
  condition=factor(sra$treatment, lvls)
)

library(tximeta)
library(SummarizedExperiment)
library(fishpond)

se <- tximeta(coldata)
y <- se
y <- scaleInfReps(y)
y <- labelKeep(y)
y <- y[ mcols(y)$keep, ]
set.seed(1)
y <- swish(y, x="condition", pair="donor", nperms=16)
hist(mcols(y)$pvalue, col="grey")
with(mcols(y),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
     )
sig <- mcols(y)$qvalue < .05
lo <- order(mcols(y)$log2FC * sig)
hi <- order(-mcols(y)$log2FC * sig)
plotInfReps(y, idx=hi[1], x="condition", cov="donor", xaxis=FALSE)
plotInfReps(y, idx=lo[1], x="condition", cov="donor", xaxis=FALSE)

library(DESeq2)
dds <- DESeqDataSet(se, ~donor + condition)
dds <- dds[rownames(y),]
dds <- DESeq(dds, fitType="local")
res <- results(dds)

table(deseq2=res$padj < .05, swish=sig)

plot(mcols(y)$qvalue, -log10(res$padj), ylim=c(0,30), cex=.3)
idx <- which(res$padj < 1e-5 & mcols(y)$qvalue > .4)

par(mfrow=c(2,5))
for (i in 1:5) 
  plotInfReps(y, idx=idx[i], x="condition", cov="donor", xaxis=FALSE)
for (i in 1:5) 
  plotCounts(dds, gene=idx[i], "condition", transform=FALSE)

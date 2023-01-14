args <- commandArgs(T)

## Loading required pakcages
library(ChIPseeker)
library(clusterProfiler)
library("org.Mm.eg.db")
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

peak <- readPeakFile(args[1])
covplot(peak)   #whole genome

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoBar(peakAnno) 
plotAnnoPie(peakAnno)
upsetplot(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")

dev.off()

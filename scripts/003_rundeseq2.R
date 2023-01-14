#! /usr/bin/Rscript

###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  USAGE: R --vanilla < RunDESeq2.R inputCounts OutputFile [group for each sample]

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./DESeq2_normalize_only.R" }

if (length(commandArgs()) < (7 - offset))
{
	message("\nAnalyze a matrix of raw counts with DESeq2.")
	message(paste("USAGE:   ", this.script, "inputCounts OutputFile group1 group2 [...]"))
	message(paste("Example: ", this.script, "InputDEseq2.txt DESeq2_output.2016.Rscript.txt UHR UHR brain brain\n"))
	q()
}


library(DESeq2)

# First argument is input file
input.filename = commandArgs()[6 - offset]
# Second argument is output file
output.filename = commandArgs()[7 - offset]
# Third and fourth arguments are groups (in same order as input file)
groups = commandArgs()[(8 - offset):length(commandArgs())]

paste("Input filename is", input.filename)
paste("Output filename is", output.filename)
# pca.figure.filename = paste(input.filename, "PCA.pdf", sep=".")
# paste("PDF file of PCA figure is", pca.figure.filename)

# Load columns of expression counts
counts = read.delim(input.filename, row.names=1)

print("Groups are")
groups
control.group = groups[1]
exp.group = groups[length(groups)]
paste("Control group is", control.group)
paste("Experimental group is", exp.group)

########################################################################################################################################################################################################
# Make a DESeqDataSet
# Need to relevel because of inconsistent strangeness
# dds = DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition=factor(groups)), design = ~ condition)
dds = DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition=relevel(factor(groups), ref=control.group)), design = ~ condition)
factor(groups)

###  MAJOR COMMAND:
# Estimate size factors and dispersions and fit generalized linear model
dds = DESeq(dds)

# Print size factors
sizeFactors(dds)

# Do stats 
res = results(dds) #padj<0.1
summary(res)
resOrdered <- res[order(res$padj),]

#up_DEGs <- subset(resOrdered,log2FoldChange > 0)
#down_DEGs <- subset(resOrdered,log2FoldChange < 0)
#up_2fold <- subset(res,log2FoldChange > 1)
#down_2fold <- subset(res,log2FoldChange < -1)
#DEGs_2fold <- subset(res,abs(log2FoldChange)>1 )

## More strictly statistic significance
#resSig <- subset(res, padj < 0.05)
#summary(resSig)
#resSigOrdered <- res[order(resSig$padj),]

# change header for log2FC column
colnames(res)[2] = paste("log2(", paste(exp.group, control.group, sep="/"), ")", sep="")

# Add normalized counts for this version (GB - 19 Mar 2013)
counts.normalized = round(t(t(counts(dds))/sizeFactors(dds)), 2)
colnames(counts.normalized) = paste(colnames(counts), "norm", sep=".")

# Print output (including norm counts)
output.table = cbind(rownames(dds), as.matrix(res), counts, counts.normalized)
colnames(output.table)[1] = "Feature.ID"

######## if the ref samples name hava a priority, you should pay attention to the counts&norm_counts of the outputfile
write.table(output.table, file=output.filename, sep="\t", quote=F, row.names=F)
########################################################################################################################################################################################################

# We may need to modify the function plotPCA() because it makes a max of 12 colors (so >12 samples cause problems)
plotPCA_BaRC = function (x, intgroup = "condition", ntop = 500)
{
	library(genefilter)
	library(RColorBrewer)
    rv = rowVars(exprs(x))
    select = order(rv, decreasing = TRUE)[seq_len(ntop)]
    pca = prcomp(t(exprs(x)[select, ]))
    fac = factor(apply(pData(x)[, intgroup, drop = FALSE], 1,
        paste, collapse = " : "))
	if (nlevels(fac) >= 3) {
		colours = brewer.pal(nlevels(fac), "Paired")
	} else { colours = c("green", "blue") }
    xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
        pch = 16, cex = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
        text = list(levels(fac)), rep = FALSE)))
}

### Do VST transformation (recommended for PCA plot)
# dds.vst = varianceStabilizingTransformation(dds, blind=TRUE)
# rld = rlogTransformation(dds, blind=TRUE)
### Run Principal Components Analysis
# pdf(pca.figure.filename)
# plotPCA(rld, intgroup=c("condition"))
# plotPCA_BaRC(rld)
# dev.off()

#paste("See", output.filename, "for output")

#### MA-plot for all DEGs with padj<0.1
plotMA(res, ylim=c(-2,2), main='DESeq2')
dev.copy(pdf,'deseq2_DEGs_MAplot.pdf')
dev.off()

## MA plot for DEGs with padj<0.05 & LFC>1
resSig <- subset(subset(res, padj < 0.05), abs(log2FoldChange)>1)
plotMA(resSig, ylim=c(-2,2), main='DESeq2_DEGs_padj0.05_LFC1')
dev.copy(pdf,'deseq2_significant_DEGs_MAplot.pdf')
dev.off()

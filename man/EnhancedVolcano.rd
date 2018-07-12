\name{EnhancedVolcano}
\alias{EnhancedVolcano}
\title{Publication-ready volcano plots with enhanced colouring and labeling}
\description{Publication-ready volcano plots with enhanced colouring and labeling}
\usage{
EnhancedVolcano(
	toptable,
	lab,
	x,
	y,
	selectLab = NULL,
	xlim = c(min(toptable[,x], na.rm=TRUE), max(toptable[,x], na.rm=TRUE)),
	ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5),
	xlab = bquote(~Log[2]~ "fold change"),
	ylab = bquote(~-Log[10]~italic(P)),
	axisLabSize = 16,
	pCutoff = 0.05,
	pLabellingCutoff = pCutoff,
	FCcutoff = 2.0,
	title = "",
	titleLabSize = 16,
	transcriptPointSize = 0.8,
	transcriptLabSize = 2.0,
	col=c("grey30", "forestgreen", "royalblue", "red2"),
	colAlpha = 1/2,
	legend=c("NS", "Log2 FC", "P", "P & Log2 FC"),
	legendPosition = "top",
	legendLabSize = 10,
	legendIconSize = 3.0,
	DrawConnectors = FALSE,
	widthConnectors = 0.5,
	colConnectors = "black",
	cutoffLineType = "longdash",
	cutoffLineCol = "black",
	cutoffLineWidth = 0.4)
}
\arguments{
	\item{toptable}{A data-frame of test statistics (if not a data frame, an attempt will be made to convert it to one). Requires at least the following: transcript names as rownames; a column for log2 fold changes; a column for nominal or adjusted p-value}
	\item{AdjustedCutoff}{Adjusted p-value cut-off for statistical significance}
	\item{LabellingCutoff}{Adjusted p-value cut-off for statistical significance for labeling of transcripts}
	\item{FCCutoff}{absolute log2FoldChange cut-off for statistical significance}
	\item{main}{Plot title}
	\item{col}{Colour shading of points for: log2FoldChange <= FCCutoff && padj >= AdjustedCutoff; log2FoldChange > FCCutoff && padj >= AdjustedCutoff; log2FoldChange <= FCCutoff && padj < AdjustedCutoff, log2FoldChange > FCCutoff && padj < AdjustedCutoff}
	\item{DrawConnectors}{Spread out labels and connect to points by lines (TRUE/FALSE)}
}
\details{
...
}
\value{
A \code{\link{ggplot2}} object.
}
\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}
\examples{
library("pasilla")
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
dds <- DESeq(dds)
res <- results(dds)

source("R/EnhancedVolcano.R")
EnhancedVolcano(res,
	lab = rownames(res),
	x = "log2FoldChange",
	y = "pvalue",
	pCutoff = 10e-9,
	FCcutoff = 2.5,
	transcriptLabSize = 3.0,
	title = "DESeq2 results",
	legendPosition = "right",
	legendLabSize = 14,
	col = c("grey30", "forestgreen", "royalblue", "red2"),
	selectLab = c("FBgn0039155","FBgn0003360","FBgn0034434"),
	DrawConnectors = TRUE)
}

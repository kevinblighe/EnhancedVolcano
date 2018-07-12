\name{EnhancedVolcano}
\alias{EnhancedVolcano}
\title{Publication-ready volcano plots with enhanced colouring and labeling}
\description{Publication-ready volcano plots with enhanced colouring and labeling}
\usage{
EnhancedVolcano(toptable,
	AdjustedCutoff = 0.05,
	LabellingCutoff = 0.05,
	FCCutoff = 2.0,
	main = "DESeq2 results",
	col = c("grey30", "forestgreen", "royalblue", "red2"),
	DrawConnectors = FALSE)
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

EnhancedVolcano(res,
	x = "log2FoldChange",
	y = "padj",
	main = "DESeq2 results",
	col = c("grey30", "forestgreen", "royalblue", "red2"),
	DrawConnectors = FALSE)
}

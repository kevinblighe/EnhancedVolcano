\name{EnhancedVolcanoEdgeR}
\alias{EnhancedVolcanoEdgeR}
\title{Publication-ready volcano plots with enhanced colouring and labeling}
\description{Publication-ready volcano plots with enhanced colouring and labeling}
\usage{
EnhancedVolcanoDESeq2(topTable=res,
	AdjustedCutoff=0.05,
	LabellingCutoff0.05,
	FCCutoff=2.0,
	main="DESeq2 results",
	col=c("grey30", "forestgreen", "royalblue", "red2"))
}
\arguments{
	\item{topTable}{A data-frame of test statistics (if not a data frame, an atempt will be made to convert it to one). Requires at least the following: transcript names as rownames; a column named 'log2FoldChange' for DESeq2 or 'logFC' for EdgeR; a column named 'padj' for DESeq2 or 'FDR' for EdgeR}
	\item{AdjustedCutoff}{Adjusted p-value cut-off for statistical significance}
	\item{LabellingCutoff}{Adjusted p-value cut-off for statistical significance for labeling of transcripts}
	\item{FCCutoff}{absolute fold change cut-off for statistical significance (assumes log base 2)}
	\item{main}{Plot title}
	\item{col}{Colour shading of poins for Not Significant, >FCCutoff, <AdjustedCutoff, >FCCutoff && <AdjustedCutoff}
}
\details{
...
}
\value{
A \code{\link{ggplot2}} object.
}
\author{
Kevin Blighe <k.blighe@ucl.ac.uk>
}
\examples{
}

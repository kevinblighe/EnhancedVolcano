\name{EnhancedVolcano}

\alias{EnhancedVolcano}

\title{Publication-ready volcano plots with enhanced colouring and labeling.}

\description{Volcano plots represent a useful way to visualise the results of differential
expression analyses. Here, we present a highly-configurable function that
produces publication-ready volcano plots. EnhancedVolcano [@EnhancedVolcano]
will attempt to fit as many labels in the plot window as possible,
thus avoiding 'clogging' up the plot with labels that could not otherwise
have been read. Other functionality allows the user to identify up to 3
different types of attributes in the same plot space via colour, shape, size, and
shade parameter configurations.}

\usage{
EnhancedVolcano(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1,
    max(toptable[[x]], na.rm=TRUE) + 1),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = 'EnhancedVolcano',
  caption = paste0('Total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 3.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  colCustom = NULL,
  colAlpha = 1/2,
  .legend = c("NS","Log2 FC","P","P & Log2 FC"),
  legendLabels = c('NS', expression(Log[2]~FC),
    "p-value", expression(p-value~and~log[2]~FC)),
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 4.0,
  shade = NULL,
  shadeLabel = NULL,
  shadeAlpha = 1/2,
  shadeFill = "grey",
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = "partial",
  borderWidth = 0.8,
  borderColour = "black")
}

\arguments{
  \item{toptable}{A data-frame of test statistics (if not, a data frame,
  an attempt will be made to convert it to one). Requires at least
  the following: column for transcript names (can be rownames); a column
  for log2 fold changes; a column for nominal or adjusted p-value.
  REQUIRED.}
  \item{lab}{A column name in toptable containing transcript names. Can be
  rownames(toptable). REQUIRED.}
  \item{x}{A column name in toptable containing log2 fold changes. REQUIRED.}
  \item{y}{A column name in toptable containing nominal or adjusted p-values.
  REQUIRED.}
  \item{selectLab}{A vector containing a subset of lab. DEFAULT = NULL.
  OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = c(min(toptable[,x], na.rm=TRUE),
  max(toptable[,x], na.rm=TRUE)). OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = c(0, max(-log10(toptable[,y]),
  na.rm=TRUE) + 5). OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = bquote(~Log[2]~ "fold change").
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = bquote(~-Log[10]~italic(P)).
  OPTIONAL.}
  \item{axisLabSize}{Size of x- and y-axis labels. DEFAULT = 18. OPTIONAL.}
  \item{title}{Plot title. DEFAULT = 'Volcano plot'. OPTIONAL.}
  \item{subtitle}{Plot subtitle. DEFAULT = 'EnhancedVolcano'. OPTIONAL.}
  \item{caption}{Plot caption. DEFAULT =
  paste0('Total = ', nrow(toptable), ' variables'). OPTIONAL.}
  \item{titleLabSize}{Size of plot title. DEFAULT = 18. OPTIONAL.}
  \item{subtitleLabSize}{Size of plot subtitle. DEFAULT = 14. OPTIONAL.}
  \item{captionLabSize}{Size of plot caption. DEFAULT = 14. OPTIONAL.}
  \item{pCutoff}{Cut-off for statistical significance. A horizontal line
  will be drawn at -log10(pCutoff). DEFAULT = 10e-6. OPTIONAL.}
  \item{FCcutoff}{Cut-off for absolute log2 fold-change. Vertical lines will
  be drawn at the negative and positive values of log2FCcutoff. DEFAULT =
  1.0. OPTIONAL.}
  \item{cutoffLineType}{Line type for FCcutoff and pCutoff ("blank",
  "solid", "dashed", "dotted", "dotdash", "longdash", "twodash").
  DEFAULT = "longdash". OPTIONAL.}
  \item{cutoffLineCol}{Line colour for FCcutoff and pCutoff. DEFAULT =
  "black". OPTIONAL.}
  \item{cutoffLineWidth}{Line width for FCcutoff and pCutoff. DEFAULT =
  0.4. OPTIONAL.}
  \item{pointSize}{Size of plotted points for each transcript. Can be
  a single value or a vector of sizes. DEFAULT = 2.0. OPTIONAL.}
  \item{labSize}{Size of labels for each transcript. DEFAULT =
  3.0. OPTIONAL.}
  \item{labCol}{Colour of labels for each transcript. DEFAULT =
  'black'. OPTIONAL.}
  \item{labFace}{Font face of labels for each transcript. DEFAULT
  = 'plain'. OPTIONAL.}
  \item{labhjust}{Horizontal adjustment of label for each
  transcript. DEFAULT = 0.5. OPTIONAL.}
  \item{labvjust}{Vertical adjustment of label for each
  transcript. DEFAULT = 1.5. OPTIONAL.}
  \item{boxedLabels}{Logical, indicating whether or not to draw labels in
  boxes. DEFAULT = FALSE. OPTIONAL.}
  \item{shape}{Shape of the plotted points. Either a single value for
  all points, or 4 values corresponding to < abs(FCcutoff) && > pCutoff,
  > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff. DEFAULT = 19.
  OPTIONAL.}
  \item{shapeCustom}{Named vector / key-value pairs that will over-ride the
  default shape scheme. The order must match that of toptable. Names / keys
  relate to groups / categories; values relate to shape encodings. DEFAULT
  = NULL. OPTIONAL.}
  \item{col}{Colour shading for plotted points, corresponding to
  < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff,
  > abs(FCcutoff) && < pCutoff. DEFAULT = c("grey30", "forestgreen",
  "royalblue", "red2"). OPTIONAL.}
  \item{colCustom}{Named vector / key-value pairs that will over-ride the
  default colour scheme. The order must match that of toptable. Names / keys
  relate to groups / categories; values relate to colour. DEFAULT = NULL.
  OPTIONAL.}
  \item{colAlpha}{Alpha for purposes of controlling colour transparency of
  transcript points. DEFAULT = 1/2. OPTIONAL.}
  \item{.legend}{Plot legend key. DEFAULT = c("NS", "Log2 FC", "P",
  "P & Log2 FC"). OPTIONAL.}
  \item{legendLabels}{Plot legend text labels. DEFAULT = c('NS', expression(Log[2]~FC),
    "p-value", expression(p-value~and~log[2]~FC). OPTIONAL}
  \item{legendPosition}{Position of legend ("top", "bottom", "left",
  "right"). DEFAULT = "top". OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 14. OPTIONAL.}
  \item{legendIconSize}{Size of plot legend icons / symbols. DEFAULT = 4.0.
  OPTIONAL.}
  \item{shade}{A vector of transcript names to shade. DEFAULT = NULL.
  OPTIONAL.}
  \item{shadeLabel}{Label for the transcrips to shade. DEFAULT = NULL.
  OPTIONAL.}
  \item{shadeAlpha}{Alpha for purposes of controlling colour transparency of
  shaded regions. DEFAULT = 1/2. OPTIONAL.}
  \item{shadeFill}{Colour of shaded regions. DEFAULT = "grey". OPTIONAL.}
  \item{shadeSize}{Size of the shade contour lines. DEFAULT = 0.01.
  OPTIONAL.}
  \item{shadeBins}{Number of bins for the density of the shade. DEFAULT = 2.
  OPTIONAL.}
  \item{drawConnectors}{Logical, indicating whether or not to connect plot
  labels to their corresponding points by line connectors. DEFAULT = FALSE.
  OPTIONAL.}
  \item{widthConnectors}{Line width of connectors. DEFAULT = 0.5. OPTIONAL.}
  \item{typeConnectors}{Have the arrow head open or filled ('closed')?
  ('open', 'closed'). DEFAULT = 'closed'. OPTIONAL.}
  \item{endsConnectors}{Which end of connectors to draw arrow head? ('last',
  'first', 'both'). DEFAULT = 'first'. OPTIONAL.}
  \item{lengthConnectors}{Length of the connectors. DEFAULT =
  unit(0.01, 'npc'). OPTIONAL}
  \item{colConnectors}{Line colour of connectors. DEFAULT = 'grey10'. OPTIONAL.}
  \item{hline}{Draw one or more horizontal lines passing through this/these
  values on y-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = NULL. OPTIONAL.}
  \item{hlineType}{Line type for hline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{hlineCol}{Colour of hline. DEFAULT = 'black'. OPTIONAL.}
  \item{hlineWidth}{Width of hline. DEFAULT = 0.4. OPTIONAL.}
  \item{vline}{Draw one or more vertical lines passing through this/these
  values on x-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = NULL. OPTIONAL.}
  \item{vlineType}{Line type for vline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{vlineCol}{Colour of vline. DEFAULT = 'black'. OPTIONAL.}
  \item{vlineWidth}{Width of vline. DEFAULT = 0.4. OPTIONAL.}
  \item{gridlines.major}{Logical, indicating whether or not to draw major
  gridlines. DEFAULT = TRUE. OPTIONAL}
  \item{gridlines.minor}{Logical, indicating whether or not to draw minor
  gridlines. DEFAULT = TRUE. OPTIONAL}
  \item{border}{Add a border for just the x and y axes ('partial') or the
  entire plot grid ('full')? DEFAULT = 'partial'. OPTIONAL.}
  \item{borderWidth}{Width of the border on the x and y axes. DEFAULT = 0.8.
  OPTIONAL.}
  \item{borderColour}{Colour of the border on the x and y axes. DEFAULT =
  "black". OPTIONAL.}
}

\details{
Volcano plots represent a useful way to visualise the results of differential
expression analyses. Here, we present a highly-configurable function that
produces publication-ready volcano plots [@EnhancedVolcano]. EnhancedVolcano
will attempt to fit as many transcript names in the plot window as possible,
thus avoiding 'clogging' up the plot with labels that could not otherwise
have been read.
}

\value{
A \code{\link{ggplot2}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
library("pasilla")
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
  package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
  package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
  colData = coldata,
  design = ~ condition)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
dds <- DESeq(dds)
res <- results(dds)

EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 10e-4,
  FCcutoff = 1.333,
  xlim = c(-5.5, 5.5),
  ylim = c(0, -log10(10e-12)),
  pointSize = 1.5,
  labSize = 2.5,
  shape = c(6, 6, 19, 16),
  title = "DESeq2 results",
  subtitle = "Differential expression",
  caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
  legendPosition = "right",
  legendLabSize = 14,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  colAlpha = 0.9,
  drawConnectors = TRUE,
  hline = c(10e-8),
  widthConnectors = 0.5)
}

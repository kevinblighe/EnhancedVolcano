#' EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling.
#'
#' Volcano plots represent a useful way to visualise the results
#' of differential expression analyses. Here, we present a highly-configurable
#' function that produces publication-ready volcano plots. EnhancedVolcano
#' will attempt to fit as many transcript names in the plot window as possible,
#' thus avoiding 'clogging' up the plot with labels that could not otherwise
#' have been read.
#' 
#' @section EnhancedVolcano functions:
#' EnhancedVolcano(
#'     toptable,
#'     lab,
#'     x,
#'     y,
#'     selectLab = NULL,
#'     xlim = c(min(toptable[,x], na.rm=TRUE), max(toptable[,x], na.rm=TRUE)),
#'     ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5),
#'     xlab = bquote(~Log[2]~ "fold change"),
#'     ylab = bquote(~-Log[10]~italic(P)),
#'     axisLabSize = 16,
#'     pCutoff = 0.05,
#'     pLabellingCutoff = pCutoff,
#'     FCcutoff = 2.0,
#'     title = "",
#'     titleLabSize = 16,
#'     transcriptPointSize = 0.8,
#'     transcriptLabSize = 2.0,
#'     col = c("grey30", "forestgreen", "royalblue", "red2"),
#'     colAlpha = 1/2,
#'     legend = c("NS", "Log2 FC", "P", "P & Log2 FC"),
#'     legendPosition = "top",
#'     legendLabSize = 10,
#'     legendIconSize = 3.0,
#'     DrawConnectors = FALSE,
#'     widthConnectors = 0.5,
#'     colConnectors = "black",
#'     cutoffLineType = "longdash",
#'     cutoffLineCol = "black",
#'     cutoffLineWidth = 0.4)
#' }
#'
#' @docType package
#' @name EnhancedVolcano
NULL

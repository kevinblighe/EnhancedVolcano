#' Publication-ready volcano plots with enhanced colouring and labeling.
#'
#' Volcano plots represent a useful way to visualise the results of
#' differential expression analyses. Here, we present a highly-configurable
#' function that produces publication-ready volcano plots [@EnhancedVolcano].
#' \code{EnhancedVolcano} will attempt to fit as many variable names in
#' the plot window as possible, thus avoiding 'clogging' up the plot with
#' labels that could not otherwise have been read.
#'
#' @param toptable A data-frame of test statistics (if not, a data frame,
#'   an attempt will be made to convert it to one). Requires at least
#'   the following: column for variable names (can be rownames); a column
#'   for log2 fold changes; a column for nominal or adjusted p-value.
#' @param lab A column name in \code{toptable} containing variable names.
#'    Can be \code{rownames(toptable)}.
#' @param x A column name in \code{toptable} containing log2 fold changes.
#' @param y A column name in \code{toptable} containing nominal or adjusted
#'    p-values.
#' @param selectLab A vector containing a subset of lab.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param pCutoff Cut-off for statistical significance. A horizontal line
#'   will be drawn at -log10(pCutoff).
#' @param pCutoffCol Column name of statistical significance values to be used as
#'   the cut-off. A typical usage situation would be to pass nominal [un-adjusted]
#'   p-values as 'y', but adjusted p-values as pCutoffCol. In this way, a
#'   plot is generated via -log10(unadjusted p-value), but cut-offs based on
#'   adjusted p-values.
#' @param FCcutoff Cut-off for absolute log2 fold-change. Vertical lines will
#'   be drawn at the negative and positive values of log2FCcutoff.
#' @param cutoffLineType Line type for \code{FCcutoff} and \code{pCutoff}
#'   ('blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash').
#' @param cutoffLineCol Line colour for \code{FCcutoff} and \code{pCutoff}.
#' @param cutoffLineWidth Line width for \code{FCcutoff} and \code{pCutoff}.
#' @param pointSize Size of plotted points for each variable. Can be
#'   a single value or a vector of sizes.
#' @param labSize Size of labels for each variable.
#' @param labCol Colour of labels for each variable.
#' @param labFace Font face of labels for each variable.
#' @param boxedLabels Logical, indicating whether or not to draw labels in
#'   boxes.
#' @param parseLabels Logical, indicating whether or not to parse expressions
#'   in labels
#' @param shape Shape of the plotted points. Either a single value for
#'   all points, or 4 values corresponding to the default 4 legend labels
#'   specified by \code{legendLabels}.
#' @param shapeCustom Named vector / key-value pairs that will over-ride the
#'   default shape scheme. The order must match that of \code{toptable}.
#'   Names / keys relate to groups / categories; values relate to shape encodings.
#' @param col Colour shading for plotted points, corresponding to
#'   the default 4 legend labels specified by \code{legendLabels}.
#' @param colCustom Named vector / key-value pairs that will over-ride the
#'   default colour scheme. The order must match that of \code{toptable}.
#'   Names / keys relate to groups / categories; values relate to colour.
#' @param colAlpha Alpha for purposes of controlling colour transparency of
#'   variable points.
#' @param colGradient If activated, over-rides the default discrete colour scheme
#'   and replaces it with a continous scheme that shades based on nominal or 
#'   adjusted p-value specified by \code{y}. For example, c('red2', 'blue2').
#' @param colGradientBreaks Break-points for the two colours specified by
#'   colGradient.
#' @param colGradientLabels Labels for the break-points specified by
#'   colGradientBreaks.
#' @param colGradientLimits Limits of the colour scheme specified by
#'   colGradient, i.e., max and min possible p-values.
#' @param legendLabels Plot legend text labels.
#' @param legendPosition Position of legend ('top', 'bottom', 'left',
#'   'right').
#' @param legendLabSize Size of plot legend text.
#' @param legendIconSize Size of plot legend icons / symbols.
#' @param legendDropLevels Logical, drop unused factor levels from legend.
#' @param encircle A vector of variable names to encircle. Requires installation
#'   of package \code{\link[ggalt:geom_encircle]{ggalt}}.
#' @param encircleCol Colour of the encircled line.
#' @param encircleFill Colour fill of the encircled region.
#' @param encircleAlpha Alpha for purposes of controlling colour transparency of
#'   encircled region.
#' @param encircleSize Line width of the encircled line.
#' @param shade A vector of variable names to shade.
#' @param shadeFill Colour of shaded regions.
#' @param shadeAlpha Alpha for purposes of controlling colour transparency of
#'   shaded region.
#' @param shadeSize Size of the shade contour lines.
#' @param shadeBins Number of bins for the density of the shade.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding points by line connectors.
#' @param widthConnectors Line width of connectors.
#' @param typeConnectors Have the arrow head open ('open') or filled ('closed')?
#' @param endsConnectors Which end of connectors to draw arrow head? ('last',
#'   'first', 'both').
#' @param lengthConnectors Length (size) of the connector arrowheads.
#' @param colConnectors Line colour of connectors and line segments.
#' @param max.overlaps Equivalent of max.overlaps in ggrepel. Set to
#'   'Inf' to always display all labels when drawConnectors = TRUE.
#' @param maxoverlapsConnectors See max.overlaps.
#' @param min.segment.length When drawConnectors = TRUE, specifies the minimum
#'   length of the connector line segments.
#' @param directionConnectors direction in which to draw connectors.
#'   'both', 'x', or 'y'.
#' @param arrowheads Logical, indicating whether or not to draw arrow heads or
#'   or just have straight lines.
#' @param hline Draw one or more horizontal lines passing through this/these
#'   values on y-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param hlineType Line type for \code{hline} ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param hlineCol Colour of \code{hline}.
#' @param hlineWidth Width of \code{hline}.
#' @param vline Draw one or more vertical lines passing through this/these
#'   values on x-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param vlineType Line type for \code{vline} ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param vlineCol Colour of \code{vline}.
#' @param vlineWidth Width of \code{vline}.
#' @param gridlines.major Logical, indicating whether or not to draw major
#'   gridlines.
#' @param gridlines.minor Logical, indicating whether or not to draw minor
#'   gridlines.
#' @param border Add a border for just the x and y axes ('partial') or the
#'   entire plot grid ('full')?
#' @param borderWidth Width of the border on the x and y axes.
#' @param borderColour Colour of the border on the x and y axes. 
#' @param raster Logical, indicating whether to rasterize the geom_point layer. 
#'   Requires installation of \code{\link[ggrastr:geom_point_rast]{ggrastr}}.
#'
#' @details
#' Volcano plots represent a useful way to visualise the results of differential expression analyses. Here, we present a highly-configurable function that produces publication-ready volcano plots [@EnhancedVolcano]. \code{EnhancedVolcano} will attempt to fit as many variable names in the plot window as possible, thus avoiding 'clogging' up the plot with labels that could not otherwise have been read.
#'
#' @return A \code{\link{ggplot2}} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' library('pasilla')
#' pasCts <- system.file('extdata', 'pasilla_gene_counts.tsv',
#'   package='pasilla', mustWork=TRUE)
#' pasAnno <- system.file('extdata', 'pasilla_sample_annotation.csv',
#'   package='pasilla', mustWork=TRUE)
#' cts <- as.matrix(read.csv(pasCts,sep='\t',row.names='gene_id'))
#' coldata <- read.csv(pasAnno, row.names=1)
#' coldata <- coldata[,c('condition','type')]
#' rownames(coldata) <- sub('fb', '', rownames(coldata))
#' cts <- cts[, rownames(coldata)]
#' library('DESeq2')
#' dds <- DESeqDataSetFromMatrix(countData = cts,
#'   colData = coldata,
#'   design = ~ condition)
#' 
#' featureData <- data.frame(gene=rownames(cts))
#' mcols(dds) <- DataFrame(mcols(dds), featureData)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' 
#' EnhancedVolcano(res,
#'   lab = rownames(res),
#'   x = 'log2FoldChange',
#'   y = 'pvalue',
#'   pCutoff = 10e-4,
#'   FCcutoff = 1.333,
#'   xlim = c(-5.5, 5.5),
#'   ylim = c(0, -log10(10e-12)),
#'   pointSize = 1.5,
#'   labSize = 2.5,
#'   title = 'DESeq2 results',
#'   subtitle = 'Differential expression',
#'   caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
#'   legendPosition = "right",
#'   legendLabSize = 14,
#'   col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
#'   colAlpha = 0.9,
#'   drawConnectors = TRUE,
#'   hline = c(10e-8),
#'   widthConnectors = 0.5)
#'
#' @import ggplot2
#' @import ggrepel
#' @importFrom methods is
#' 
#' @export
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1.5,
    max(toptable[[x]], na.rm=TRUE) + 1.5),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = bquote(italic(EnhancedVolcano)),
  caption = paste0('total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  pCutoffCol = y,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 5.0,
  labCol = 'black',
  labFace = 'plain',
  boxedLabels = FALSE,
  parseLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = 'top',
  legendLabSize = 14,
  legendIconSize = 5.0,
  legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = 'black',
  encircleFill = 'pink',
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = 'grey',
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  max.overlaps = 15,
  maxoverlapsConnectors = NULL,
  min.segment.length = 0,
  directionConnectors = 'both',
  arrowheads = TRUE,
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
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black', 
  raster = FALSE)
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))
  }

  if(!is.numeric(toptable[[pCutoffCol]])) {
    stop(paste(y, ' is not numeric!', sep=''))
  }
  
  if (raster) {

    has_ggrastr <- ! is(try(find.package("ggrastr"), silent=TRUE), "try-error")

    if (has_ggrastr) {
      geom_point <- ggrastr::geom_point_rast
    } else {
      warning("raster disabled, required package \"ggrastr\" not installed")
    }
  }

  if (!is.null(maxoverlapsConnectors)) {
    max.overlaps <- maxoverlapsConnectors
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'

  toptable$Sig[(toptable[[pCutoffCol]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[pCutoffCol]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))
  # reset pCutoff to corresponding value on y
  # allowing to draw hline at the correct
  # threshold
  if (pCutoffCol != y) {
    pCutoff = max(
      toptable[which(
        toptable[pCutoffCol] <= pCutoff), y]
      )
  }
  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 0.5),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes to geom_point (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = Sig),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE,
        drop = legendDropLevels)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes = list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = Sig),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = Sig,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    }
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == 'partial') {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop('Unrecognised value passed to \'border\'. Must be \'full\' or \'partial\'')
  }

  # Gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (!boxedLabels) {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        xlim = c(NA, NA),
        ylim = c(NA, NA),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE,
        direction = directionConnectors,
        max.overlaps = max.overlaps,
        min.segment.length = min.segment.length)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        xlim = c(NA, NA),
        ylim = c(NA, NA),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE,
        direction = directionConnectors,
        max.overlaps = max.overlaps,
        min.segment.length = min.segment.length)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE)
    }

  } else {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        xlim = c(NA, NA),
        ylim = c(NA, NA),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE,
        direction = directionConnectors,
        max.overlaps = max.overlaps,
        min.segment.length = min.segment.length)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        xlim = c(NA, NA),
        ylim = c(NA, NA),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE,
        direction = directionConnectors,
        max.overlaps = max.overlaps,
        min.segment.length = min.segment.length)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        colour = labCol,
        fontface = labFace,
        parse = parseLabels,
        na.rm = TRUE)

    }
  }

  # encircle
  if (!is.null(encircle)) {

    if (is(try(find.package("ggalt"), silent=TRUE), "try-error")) {
      stop("Please install package \"ggalt\" to access the \"encircle\" features")
    }

    plot <- plot + 
      ggalt::geom_encircle(
        data = subset(toptable,
          rownames(toptable) %in% encircle),
        colour = encircleCol,
        fill = encircleFill,
        alpha = encircleAlpha,
        size = encircleSize,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  # shade
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  plot <- plot + coord_cartesian(clip = 'off')

  return(plot)
}

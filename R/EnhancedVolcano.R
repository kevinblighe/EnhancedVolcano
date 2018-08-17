EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[,x], na.rm=TRUE),
    max(toptable[,x], na.rm=TRUE)),
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
  transcriptLabSize = 3.0,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  colAlpha = 1/2,
  legend = c("NS","Log2 FC","P","P & Log2 FC"),
  legendPosition = "top",
  legendLabSize = 10,
  legendIconSize = 3.0,
  DrawConnectors = FALSE,
  widthConnectors = 0.5,
  colConnectors = "black",
  cutoffLineType = "longdash",
  cutoffLineCol = "black",
  cutoffLineWidth = 0.4)
{
  if(!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call.=FALSE)
  }

  if(!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call.=FALSE)
  }

  if(!is.numeric(toptable[,x])) {
    stop(paste(x[i], " is not numeric!", sep=""))
  }

  if(!is.numeric(toptable[,y])) {
    stop(paste(x[i], " is not numeric!", sep=""))
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[,x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[,y]<pCutoff)] <- "P"
  toptable$Sig[(toptable[,y]<pCutoff) &
    (abs(toptable[,x])>FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig,
    levels=c("NS","FC","P","FC_P"))

  if (min(toptable[,y], na.rm=TRUE) == 0) {
    warning(paste("One or more P values is 0.",
      "Converting to minimum possible value..."),
      call. = FALSE)
    toptable[which(toptable[,y] == 0), y] <- .Machine$double.xmin
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[,x]
  toptable$yvals <- toptable[,y]

  if (!is.null(selectLab)) {
    names.new <- rep("", length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) +

    geom_point(aes(color=factor(Sig)),
      alpha=colAlpha,
      size=transcriptPointSize) +

    scale_color_manual(values=c(NS=col[1],
      FC=col[2],
      P=col[3],
      FC_P=col[4]),
      labels=c(NS=legend[1],
      FC=paste(legend[2], sep=""),
      P=paste(legend[3], sep=""),
      FC_P=paste(legend[4], sep=""))) +

    theme_bw(base_size=24) +

    theme(
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=titleLabSize, face="bold", vjust=1),

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),

      axis.text.x=element_text(angle=0, size=axisLabSize, vjust=1),
      axis.text.y=element_text(angle=0, size=axisLabSize, vjust=1),
      axis.title=element_text(size=axisLabSize),

      legend.position=legendPosition,
      legend.key=element_blank(),
      legend.key.size=unit(0.5, "cm"),
      legend.text=element_text(size=legendLabSize),

      title=element_text(size=legendLabSize),
      legend.title=element_blank()
    ) +

    guides(colour = guide_legend(
      override.aes=list(size=legendIconSize))) +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    ggtitle(title) +

    geom_vline(xintercept=c(-FCcutoff, FCcutoff),
      linetype=cutoffLineType,
      colour=cutoffLineCol,
      size=cutoffLineWidth) +

    geom_hline(yintercept=-log10(pCutoff),
      linetype=cutoffLineType,
      colour=cutoffLineCol,
      size=cutoffLineWidth)

  if (DrawConnectors == TRUE) {
    plot <- plot + geom_text_repel(
      data=subset(toptable,
        toptable[,y]<pLabellingCutoff &
          abs(toptable[,x])>FCcutoff),
            aes(label=subset(toptable,
              toptable[,y]<pLabellingCutoff &
                abs(toptable[,x])>FCcutoff)[,"lab"]),
              size = transcriptLabSize,
              segment.color = colConnectors,
              segment.size = widthConnectors,
              vjust = 1.5)
  } else if (DrawConnectors == FALSE && !is.null(selectLab)) {
    plot <- plot + geom_text(
      data=subset(toptable,
        toptable[,y]<pLabellingCutoff &
          abs(toptable[,x])>FCcutoff),
            aes(label=subset(toptable,
              toptable[,y]<pLabellingCutoff &
                abs(toptable[,x])>FCcutoff)[,"lab"]),
              size = transcriptLabSize,
              check_overlap = FALSE,
              vjust = 1.5)
  } else if (DrawConnectors == FALSE && is.null(selectLab)) {
    plot <- plot + geom_text(
      data=subset(toptable,
        toptable[,y]<pLabellingCutoff &
          abs(toptable[,x])>FCcutoff),
            aes(label=subset(toptable,
              toptable[,y]<pLabellingCutoff &
                abs(toptable[,x])>FCcutoff)[,"lab"]),
              size = transcriptLabSize,
              check_overlap = TRUE,
              vjust = 1.5)
  }

  return(plot)
}

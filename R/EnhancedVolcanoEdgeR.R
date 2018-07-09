EnhancedVolcanoEdgeR <- function(toptable, AdjustedCutoff = 0.05, LabellingCutoff = 0.05, FCCutoff = 2.0, main = "EdgeR results", col=c("grey30", "forestgreen", "royalblue", "red2"), DrawConnectors = FALSE)
{
	if(!requireNamespace("ggplot2")) { stop( "Please install ggplot2 first.", call.=FALSE) }

	requireNamespace("ggplot2")

	toptable <- as.data.frame(toptable)

	toptable$Significance <- "NS"
	toptable$Significance[(abs(toptable$logFC) > FCCutoff)] <- "FC"
	toptable$Significance[(toptable$FDR<AdjustedCutoff)] <- "FDR"
	toptable$Significance[(toptable$FDR<AdjustedCutoff) & (abs(toptable$logFC)>FCCutoff)] <- "FC_FDR"
	table(toptable$Significance)
	toptable$Significance <- factor(toptable$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))

	plot <- ggplot2::ggplot(toptable, ggplot2::aes(x=logFC, y=-log10(FDR))) +

		#Add points:
		#	Colour based on factors set a few lines up
		#	'alpha' provides gradual shading of colour
		#	Set size of points
		ggplot2::geom_point(ggplot2::aes(color=factor(Significance)), alpha=1/2, size=0.8) +

		#Choose which colours to use; otherwise
		ggplot2::scale_color_manual(values=c(NS=col[1], FC=col[2], FDR=col[3], FC_FDR=col[4]), labels=c(NS="NS", FC=paste("LogFC>|", FCCutoff, "|", sep=""), FDR=paste("FDR Q<", AdjustedCutoff, sep=""), FC_FDR=paste("FDR Q<", AdjustedCutoff, " & LogFC>|", FCCutoff, "|", sep=""))) +

		#Set the size of the plotting window
		ggplot2::theme_bw(base_size=24) +

		#Modify various aspects of the plot text and legend
		ggplot2::theme(legend.background=ggplot2::element_rect(),
			plot.title=ggplot2::element_text(angle=0, size=12, face="bold", vjust=1),

			panel.grid.major=ggplot2::element_blank(),	#Remove gridlines
			panel.grid.minor=ggplot2::element_blank(),	#Remove gridlines

			axis.text.x=ggplot2::element_text(angle=0, size=12, vjust=1),
			axis.text.y=ggplot2::element_text(angle=0, size=12, vjust=1),
			axis.title=ggplot2::element_text(size=12),

			#Legend
			legend.position="top",			#Moves the legend to the top of the plot
			legend.key=ggplot2::element_blank(),		#removes the border
			legend.key.size=ggplot2::unit(0.5, "cm"),	#Sets overall area/size of the legend
			legend.text=ggplot2::element_text(size=8),	#Text size
			title=ggplot2::element_text(size=8),		#Title text size
			legend.title=ggplot2::element_blank()) +		#Remove the title

		#Change the size of the icons/symbols in the legend
		ggplot2::guides(colour = ggplot2::guide_legend(override.aes=list(size=2.5))) +

		#Set x- and y-axes labels
		ggplot2::xlab(bquote(~Log[2]~ "fold change")) +
		ggplot2::ylab(bquote(~-Log[10]~adjusted~italic(P))) +

		#Set the axis limits
		#xlim(-6.5, 6.5) +
		#ylim(0, 100) +

		#Set title
		ggplot2::ggtitle(main) +

		#Add a vertical line for fold change cut-offs
		ggplot2::geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +

		#Add a horizontal line for P-value cut-off
		ggplot2::geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)

		#Tidy the text labels for a subset of genes
		if (DrawConnectors) {
			plot + ggplot2::geom_text(data=subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff),
				ggplot2::aes(label=rownames(subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff))),
				size=2.25,
				segment.color="black",
				segment.size=0.01,
				check_overlap=TRUE,
				vjust=1.0)
		} else if (!DrawConnectors) {
			plot + ggplot2::geom_text(data=subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff),
				ggplot2::aes(label=rownames(subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff))),
				size=2.25,
				check_overlap=TRUE,
				vjust=1.0)
		}

	return(plot)
}

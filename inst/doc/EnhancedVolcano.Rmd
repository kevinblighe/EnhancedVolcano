# EnhancedVolcano

<h2>Vignette</h2>
Publication-ready volcano plots with enhanced colouring and labeling.

<hr>

Following tutorial of http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

```{r}

	source("https://bioconductor.org/biocLite.R")
	biocLite("airway")
	library(airway)

	library(magrittr)

	data("airway")
	airway$dex %<>% relevel("untrt")
	airway$dex

```

```{r}

	[1] untrt trt   untrt trt   untrt trt   untrt trt  
	Levels: untrt trt

```

```{r}

	library("DESeq2")
	dds <- DESeqDataSet(airway, design = ~ cell + dex)
	dds <- DESeq(dds, betaPrior=FALSE)

```

```{r}

	estimating size factors
	estimating dispersions
	gene-wise dispersion estimates
	mean-dispersion relationship
	final dispersion estimates
	fitting model and testing

```

```{r}

	res1 <- results(dds, contrast = c("dex","trt","untrt"))
	res1 <- lfcShrink(dds, contrast = c("dex","trt","untrt"), res=res1)

	res2 <- results(dds, contrast = c("cell", "N061011", "N61311"))
	res2 <- lfcShrink(dds, contrast = c("cell", "N061011", "N61311"), res=res2)

	res3 <- results(dds, contrast = c("cell", "N061011", "N052611"))
	res3 <- lfcShrink(dds, contrast = c("cell", "N061011", "N052611"), res=res3)

	res4 <- results(dds, contrast = c("cell", "N061011", "N052611"))
	res4 <- lfcShrink(dds, contrast = c("cell", "N061011", "N052611"), res=res4)

	biocLite("EnhancedVolcano")
	library(EnhancedVolcano)

	p1 <- EnhancedVolcano(res1,
		x = "log2FoldChange",
		y = "pvalue",
		pCutoff = 10e-18,
		FCcutoff = 2.0,
		title = "Treated versus untreated",
		DrawConnectors = TRUE)

	p2 <- EnhancedVolcano(res2,
		x = "log2FoldChange",
		y = "pvalue",
		title = "N061011 versus N61311",
		col = c("black", "black", "black", "red3"),
		DrawConnectors = TRUE)

	p3 <- EnhancedVolcano(res3,
		x = "log2FoldChange",
		y = "pvalue",
		title = "N061011 versus N052611",
		col = c("grey30", "limegreengreen", "darkblue", "pink"),
		DrawConnectors = FALSE)

	p4 <- EnhancedVolcano(res4,
		x = "log2FoldChange",
		y = "padj",
		ylab = bquote(~-Log[10]~adjusted~italic(P)),
		title = "N061011 versus N61311",
		col = c("black", "dodgerblue", "skyblue", "gold"),
		DrawConnectors = FALSE)

	pdf("test.pdf", width=18, height=8)
	library(gridExtra)
 	library(grid)
	grid.arrange(p1, p2, p3, p4, ncol=4, top="EnhancedVolcano")
	grid.rect(gp=gpar(fill=NA))
	dev.off()

```


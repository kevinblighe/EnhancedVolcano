Publication-ready volcano plots with enhanced colouring and labeling
================
Kevin Blighe
2018-08-17

-   [Introduction.](#introduction.)
-   [Installation.](#installation.)
    -   [1. Download the package from Bioconductor.](#download-the-package-from-bioconductor.)
    -   [2. Load the package into R session.](#load-the-package-into-r-session.)
-   [Quick start.](#quick-start.)
    -   [Plot the most basic volcano plot.](#plot-the-most-basic-volcano-plot.)
-   [Advanced features.](#advanced-features.)
    -   [Modify cut-offs for log2FC and P value; add title; adjust point and label size.](#modify-cut-offs-for-log2fc-and-p-value-add-title-adjust-point-and-label-size.)
    -   [Adjust colour and alpha for point shading.](#adjust-colour-and-alpha-for-point-shading.)
    -   [Adjust axis limits.](#adjust-axis-limits.)
    -   [Adjust cut-off lines.](#adjust-cut-off-lines.)
    -   [Adjust legend position, size, and text.](#adjust-legend-position-size-and-text.)
    -   [Plot adjusted p-values.](#plot-adjusted-p-values.)
    -   [Fit more labels by adding connectors.](#fit-more-labels-by-adding-connectors.)
    -   [Only label key transcripts.](#only-label-key-transcripts.)
    -   [Plot multiple volcanos on the same page.](#plot-multiple-volcanos-on-the-same-page.)
-   [Acknowledgments](#acknowledgments)
-   [Session info](#session-info)
    -   [References](#references)

Introduction.
=============

Volcano plots represent a useful way to visualise the results of differential expression analyses. Here, we present a highly-configurable function that produces publication-ready volcano plots (Blighe 2018). EnhancedVolcano will attempt to fit as many transcript names in the plot window as possible, thus avoiding 'clogging' up the plot with labels that could not otherwise have been read.

Installation.
=============

1. Download the package from Bioconductor.
------------------------------------------

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) 
install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
```

Note: to install development version:

``` r
devtools::install_github("kevinblighe/EnhancedVolcano")
```

2. Load the package into R session.
-----------------------------------

``` r
library(EnhancedVolcano)
```

Quick start.
============

For this example, we will follow the tutorial (from Section 3.1) of [RNA-seq workflow: gene-level exploratory analysis and differential expression](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). Specifically, we will load the 'airway' data, where different airway smooth muscle cells were treated with dexamethasone.

``` r
library(airway)

library(magrittr)

data("airway")

airway$dex %<>% relevel("untrt")
```

Conduct differential expression using DESeq2 in order to create 2 sets of results:

``` r
library("DESeq2")

dds <- DESeqDataSet(airway, design = ~cell + dex)
dds <- DESeq(dds, betaPrior = FALSE)
res1 <- results(dds, contrast = c("dex", "trt", "untrt"))
res1 <- lfcShrink(dds, contrast = c("dex", "trt", "untrt"), res = res1)
res2 <- results(dds, contrast = c("cell", "N061011", "N61311"))
res2 <- lfcShrink(dds, contrast = c("cell", "N061011", "N61311"), res = res2)
```

Plot the most basic volcano plot.
---------------------------------

For the most basic volcano plot, only a single data-frame or -matrix of test results is required, containing transcript names, log2FC, and adjusted or unadjusted P values. The default cut-off for log2FC is &gt;|2|; the default cut-off for P value is 0.05.

``` r
    EnhancedVolcano(res1,

        lab = rownames(res1),

        x = "log2FoldChange",

        y = "pvalue")
```

![Plot the most basic volcano plot.](README_files/figure-markdown_github/ex1-1.png)

Advanced features.
==================

Virtually all aspects of an EnhancedVolcano plot can be configured for the purposes of accommodating all types of statistical distributions and labelling preferences. EnhancedVolcano will only attempt to label genes that pass the thresholds that you set for statistical significance, i.e., 'pCutoff' and 'FCcutoff'. In addition, it will only label as many of these that can reasonably fit in the plot space. The user can optionally supply a vector of transcript names (as 'selectLab') that s/he wishes to label in the plot.

Modify cut-offs for log2FC and P value; add title; adjust point and label size.
-------------------------------------------------------------------------------

The default P value cut-off of 0.05 may be too relaxed for most studies, which may therefore necessitate increasing this threshold by a few orders of magnitude. Equally, the log2FC cut-offs may be too stringent, given that moderated 'shrunk' estimates of log2FC differences in differential expression analysis can now be calculated.

In this example, we also modify the point and label size, which can help to improve clarity where many transcripts went into the differential expression analysis.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-12,

        FCcutoff = 1.5,

        transcriptPointSize = 1.5,

        transcriptLabSize = 3.0,

        title = "N061011 versus N61311")
```

![Modify cut-offs for log2FC and P value; add title; adjust point and label size.](README_files/figure-markdown_github/ex2-1.png)

Adjust colour and alpha for point shading.
------------------------------------------

The default colour scheme may not be to everyone's taste. Here we make it such that only the transcripts passing both the log2FC and P value thresholds are coloured red, with everything else black. We also adjust the value for 'alpha', which controls the transparency of the plotted points: 1 = 100% opaque; 0 = 100% transparent.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-12,

        FCcutoff = 1.5,

        transcriptPointSize = 1.5,

        transcriptLabSize = 3.0,

        title = "N061011 versus N61311",

        col=c("black", "black", "black", "red3"),

        colAlpha = 1)
```

![Adjust colour and alpha for point shading.](README_files/figure-markdown_github/ex3-1.png)

Adjust axis limits.
-------------------

The x-axis limits for log2FC defaults to the max and min of the log2FC values passed to EnhancedVolcano. This can often render the plot asymmetrical; so, the user may wish to set these axis limits to the same absolute values, e.g., c(-8, 8). One can also modify the y-axis limits, but this should be a less common occurrence.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-12,

        FCcutoff = 1.5,

        transcriptPointSize = 1.5,

        transcriptLabSize = 3.0,

        title = "N061011 versus N61311",

        colAlpha = 1,

        xlim = c(-8, 8),

        ylim = c(0, -log10(10e-32)))
```

![Adjust axis limits.](README_files/figure-markdown_github/ex4-1.png)

Adjust cut-off lines.
---------------------

The lines that are drawn to indicate cut-off points are also modifiable. The parameter 'cutoffLineType' accepts the following values: "blank", "solid", "dashed", "dotted", "dotdash", "longdash", and "twodash". The colour and thickness of these can also be modified with 'cutoffLineCol' and 'cutoffLineWidth'. To disable the lines, set either cutoffLineType="blank" or cutoffLineWidth=0.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-12,

        FCcutoff = 1.5,

        transcriptPointSize = 1.5,

        transcriptLabSize = 3.0,

        title = "N061011 versus N61311",

        colAlpha = 1,

        xlim = c(-8, 8),

        ylim = c(0, -log10(10e-32)),

        cutoffLineType = "twodash",

        cutoffLineCol = "red3",

        cutoffLineWidth = 1.5)
```

![Adjust cut-off lines.](README_files/figure-markdown_github/ex5-1.png)

Adjust legend position, size, and text.
---------------------------------------

The position of the legend can also be changed to "left" or "right" (and stacked vertically), or "top" or "bottom" (stacked horizontally). The legend text, label size, and icon size can also be modified.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-12,

        FCcutoff = 1.5,

        transcriptPointSize = 1.5,

        transcriptLabSize = 3.0,

        colAlpha = 1,

        cutoffLineType = "twodash",

        cutoffLineCol = "red4",

        cutoffLineWidth = 1.0,

        legend=c("NS","Log (base 2) fold-change","P value",
            "P value & Log (base 2) fold-change"),

        legendPosition = "right",

        legendLabSize = 14,

        legendIconSize = 5.0)
```

![Adjust legend position, size, and text.](README_files/figure-markdown_github/ex6-1.png)

Note: to make the legend completely invisible, specify:

``` r
legend=c("","","",""), legendLabSize=-1, legendIconSize=-1
```

Plot adjusted p-values.
-----------------------

Volcano plots do not have to be plotted with nominal (unadjusted P values). Simply provide a column name relating to adjusted P values and you can also generate a volcano with these. In this case, the cutoff for the P value then relates to the adjusted P value. Here, we also modify the axis titles by supplying an expression via the bquote function.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "padj",

        xlab = bquote(~Log[2]~ "fold change"),

        ylab = bquote(~-Log[10]~adjusted~italic(P)),

        pCutoff = 0.0001,

        FCcutoff = 1.0,

        xlim=c(-6,6),

        transcriptLabSize = 3.0,

        colAlpha = 1,

        legend=c("NS","Log2 FC","Adjusted p-value",
            "Adjusted p-value & Log2 FC"),

        legendPosition = "bottom",

        legendLabSize = 10,

        legendIconSize = 3.0)
```

![Plot adjusted p-values.](README_files/figure-markdown_github/ex7-1.png)

Fit more labels by adding connectors.
-------------------------------------

In order to maximise free space in the plot window, one can fit more transcript labels by adding connectors from labels to points, where appropriate. The width and colour of these connectors can also be modified with widthConnectors and colConnectors, respectively.

The result may not always be desirable as it can make the plot look overcrowded.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "padj",

        xlab = bquote(~Log[2]~ "fold change"),

        ylab = bquote(~-Log[10]~adjusted~italic(P)),

        pCutoff = 0.0001,

        FCcutoff = 2.0,

        xlim = c(-6,6),

        transcriptLabSize = 3.0,

        colAlpha = 1,

        legend=c("NS","Log2 FC","Adjusted p-value",
            "Adjusted p-value & Log2 FC"),

        legendPosition = "bottom",

        legendLabSize = 10,

        legendIconSize = 3.0,

        DrawConnectors = TRUE,

        widthConnectors = 0.2,

        colConnectors = "grey30")
```

![Fit more labels by adding connectors.](README_files/figure-markdown_github/ex8-1.png)

Only label key transcripts.
---------------------------

In many situations, people may only wish to label their key transcripts / transcripts of interest. One can therefore supply a vector of these transcripts via the 'selectLab' parameter, the contents of which have to also be present in the vector passed to 'lab'. In addition, only those transcripts that pass both the cutoff for log2FC and P value will be labelled.

``` r
    EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "padj",

        selectLab = c("ENSG00000106565","ENSG00000187758"),

        xlab = bquote(~Log[2]~ "fold change"),

        ylab = bquote(~-Log[10]~adjusted~italic(P)),

        pCutoff = 0.0001,

        FCcutoff = 2.0,

        xlim = c(-6,6),

        transcriptPointSize = 1.8,

        transcriptLabSize = 5.0,

        colAlpha = 1,

        legend=c("NS","Log2 FC","Adjusted p-value",
            "Adjusted p-value & Log2 FC"),

        legendPosition = "bottom",

        legendLabSize = 10,

        legendIconSize = 3.0)
```

![Only label key transcripts.](README_files/figure-markdown_github/ex9-1.png)

Plot multiple volcanos on the same page.
----------------------------------------

One can also plot multiple volcanos on the same plot via the use of the grid and gridExtra packages.

``` r
    p1 <- EnhancedVolcano(res1,

        lab = rownames(res1),

        x = "log2FoldChange",

        y = "pvalue",

        pCutoff = 10e-24,

        FCcutoff = 2.0,

        transcriptLabSize = 2.5,

        colAlpha = 1,

        legendPosition = "bottom",

        legendLabSize = 10,

        legendIconSize = 3.0)

    p2 <- EnhancedVolcano(res2,

        lab = rownames(res2),

        x = "log2FoldChange",

        y = "padj",

        selectLab = c("ENSG00000106565","ENSG00000187758"),

        xlab = bquote(~Log[2]~ "fold change"),

        ylab = bquote(~-Log[10]~adjusted~italic(P)),

        pCutoff = 0.0001,

        FCcutoff = 2.0,

        xlim = c(-6,6),

        transcriptLabSize = 5.0,

        colAlpha = 1,

        legend=c("NS","Log2 FC","Adjusted p-value",
            "Adjusted p-value & Log2 FC"),

        legendPosition = "bottom",

        legendLabSize = 10,

        legendIconSize = 3.0)

    library(gridExtra)

    library(grid)

    grid.arrange(p1, p2, ncol=2, top="EnhancedVolcano")

    grid.rect(gp=gpar(fill=NA))
```

![Plot multiple volcanos on the same page.](README_files/figure-markdown_github/ex10-1.png)

Acknowledgments
===============

The development of *EnhancedVolcano* has benefited from contributions and suggestions from:

Sharmila Rana, [Myles Lewis](https://www.qmul.ac.uk/whri/people/academic-staff/items/lewismyles.html)

Session info
============

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3               DESeq2_1.20.0              
    ##  [3] magrittr_1.5                airway_0.114.0             
    ##  [5] SummarizedExperiment_1.10.1 DelayedArray_0.6.1         
    ##  [7] BiocParallel_1.14.1         matrixStats_0.53.1         
    ##  [9] Biobase_2.40.0              GenomicRanges_1.32.3       
    ## [11] GenomeInfoDb_1.16.0         IRanges_2.14.10            
    ## [13] S4Vectors_0.18.3            BiocGenerics_0.26.0        
    ## [15] EnhancedVolcano_0.99.8      ggrepel_0.8.0              
    ## [17] ggplot2_3.0.0               knitr_1.20                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7            splines_3.5.0          Formula_1.2-3         
    ##  [4] assertthat_0.2.0       highr_0.7              latticeExtra_0.6-28   
    ##  [7] blob_1.1.1             GenomeInfoDbData_1.1.0 yaml_2.1.19           
    ## [10] RSQLite_2.1.1          pillar_1.2.3           backports_1.1.2       
    ## [13] lattice_0.20-35        glue_1.2.0             digest_0.6.15         
    ## [16] RColorBrewer_1.1-2     XVector_0.20.0         checkmate_1.8.5       
    ## [19] colorspace_1.3-2       htmltools_0.3.6        Matrix_1.2-14         
    ## [22] plyr_1.8.4             XML_3.98-1.11          pkgconfig_2.0.1       
    ## [25] genefilter_1.62.0      zlibbioc_1.26.0        purrr_0.2.5           
    ## [28] xtable_1.8-2           scales_1.0.0           tibble_1.4.2          
    ## [31] htmlTable_1.12         annotate_1.58.0        withr_2.1.2           
    ## [34] nnet_7.3-12            lazyeval_0.2.1         survival_2.42-3       
    ## [37] memoise_1.1.0          evaluate_0.10.1        foreign_0.8-70        
    ## [40] tools_3.5.0            data.table_1.11.4      formatR_1.5           
    ## [43] stringr_1.3.1          locfit_1.5-9.1         munsell_0.5.0         
    ## [46] cluster_2.0.7-1        AnnotationDbi_1.42.1   bindrcpp_0.2.2        
    ## [49] compiler_3.5.0         rlang_0.2.1            RCurl_1.95-4.10       
    ## [52] rstudioapi_0.7         htmlwidgets_1.2        labeling_0.3          
    ## [55] bitops_1.0-6           base64enc_0.1-3        rmarkdown_1.10        
    ## [58] gtable_0.2.0           DBI_1.0.0              R6_2.2.2              
    ## [61] dplyr_0.7.5            bit_1.1-14             bindr_0.1.1           
    ## [64] Hmisc_4.1-1            rprojroot_1.3-2        stringi_1.2.3         
    ## [67] Rcpp_0.12.18           geneplotter_1.58.0     rpart_4.1-13          
    ## [70] acepack_1.4.1          tidyselect_0.2.4

References
----------

Blighe (2018)

Blighe, Kevin. 2018. “EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling.” <https://github.com/kevinblighe>.

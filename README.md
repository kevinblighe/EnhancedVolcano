# EnhancedVolcano
<h1>Publication-ready volcano plots with enhanced colouring and labeling</h1>
<img src="https://github.com/kevinblighe/EnhancedVolcano/blob/master/vignettes/Volcano.png">
<br>
<h2>Imports</h2>
<ul>
  <li>ggplot2</li>
  </ul>
<h2>Execution</h2>
<code>EnhancedVolcanoDESeq2(topTable, AdjustedCutoff, LabellingCutoff, FCCutoff, main, col)</code>
<code>EnhancedVolcanoEdgeR(topTable, AdjustedCutoff, LabellingCutoff, FCCutoff, main, col)</code>
<br>
<h2>Example</h2

<code>EnhancedVolcanoDESeq2(topTable=results, AdjustedCutoff=0.05, LabellingCutoff0.05, FCCutoff=2.0, main="DESeq2 results", col=c("grey30", "forestgreen", "royalblue", "red2"))</code>

<code>EnhancedVolcanoEdgeR(topTable=results, AdjustedCutoff=0.05, LabellingCutoff0.05, FCCutoff=2.0, main="EdgeR results", col=c("grey30", "forestgreen", "royalblue", "red2"))</code>
<br>
<h2>Parameters</h2>
<ul>
<li>toptable, data-frame of test statistics. Requires at least the following
  <ul>
    <li>gene names as rownames</li>
  <li>column named 'log2FoldChange' for DESeq2 or 'logFC' for EdgeR</li>
    <li>column named 'padj' for DESeq2 or 'FDR' for EdgeR</li>
  </ul>
<li>topTable, A data-frame of test statistics (if not a data frame, an atempt will be made to convert it to one). Requires at least the following: transcript names as rownames; a column named 'log2FoldChange' for DESeq2 or 'logFC' for EdgeR; a column named 'padj' for DESeq2 or 'FDR' for EdgeR</li>
<li>AdjustedCutoff, Adjusted p-value cut-off for statistical significance</li>
<li>LabellingCutoff, Adjusted p-value cut-off for statistical significance for labeling of transcripts</li>
<li>FCCutoff, absolute fold change cut-off for statistical significance (assumes log base 2)</li>
<li>main, Plot title</li>
<li>col, Colour shading of poins for Not Significant, >FCCutoff, <AdjustedCutoff, >FCCutoff && <AdjustedCutoff</li>
  </ul>
<br>
<h2>Credits</h2>
<ul>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Myles Lewis (Queen Mary University of London)</li>
  <li>Sharmila Rana (Queen Mary University of London)</li>
</ul>

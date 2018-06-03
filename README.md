# EnhancedVolcano
<h1>Volcano plot with enhanced colouring and labeling</h1>
<img src="https://github.com/kevinblighe/EnhancedVolcano/blob/master/Volcano.png">
<br>
<h2>Requires</h2>
<ul>
  <li>ggplot2</li>
  </ul>
<h2>Execution</h2>
<code>EnhancedVolcanoDESeq2(topTable, NominalCutoff, AdjustedCutoff, LabellingCutoff, FCCutoff, main)</code>
<code>EnhancedVolcanoEdgeR(topTable, NominalCutoff, AdjustedCutoff, LabellingCutoff, FCCutoff, main)</code>
<br>
<h2>Example</h2

<code>EnhancedVolcanoDESeq2(topTable=results, AdjustedCutoff=0.05, LabellingCutoff0.05, FCCutoff=2.0, main="DESeq2 results")</code>

<code>EnhancedVolcanoEdgeR(topTable=results, AdjustedCutoff=0.05, LabellingCutoff0.05, FCCutoff=2.0, main="EdgeR results")</code>
<br>
<h2>Parameters</h2>
<ul>
<li>toptable, data-frame of test statistics. Requires at least the following
  <ul>
    <li>gene names as rownames</li>
  <li>column named 'log2FoldChange'</li>
    <li>column named 'padj'</li>
  </ul>
<li>NominalCutoff, nominal p-value cut-off for statistical significance (obsolete)</li>
<li>AdjustedCutoff, adjusted p-value cut-off for statistical significance</li>
<li>LabellingCutoff, adjusted p-value cut-off for statistical significance for labels</li>
<li>FCCutoff, absolute log (base 2) fold change cut-off for statistical significance</li>
<li>main, title</li>
  </ul>
<br>
<h2>Credits</h2>
<ul>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Myles Lewis (Queen Mary University of London)</li>
  <li>Sharmila Rana (Queen Mary University of London)</li>
</ul>

test_that("datatypes", {
    expect_type(EnhancedVolcano(toptable),
        c("S4", "list"))
    expect_type(EnhancedVolcano(x, y, lab, selectLab, col, legend,
        legendPosition, colConnectors, cutoffLineType, cutoffLineCol),
        c("character"))
    expect_type(EnhancedVolcano(xlab, ylab, title),
        c("language", "character"))
    expect_type(EnhancedVolcano(DrawConnectors),
        c("logical"))
    expect_type(EnhancedVolcano(axisLabSize, pCutoff, pLabellingCutoff,
        FCcutoff, titleLabSize, transcriptPointSize, transcriptLabSize,
        colAlpha, legendLabSize, legendIconSize, widthConnectors,
        cutoffLineWidth),
        c("integer", "double"))
    expect_gt(EnhancedVolcano(axisLabSize, pCutoff, pLabellingCutoff,
        FCcutoff, titleLabSize, transcriptPointSize, transcriptLabSize,
        colAlpha, legendLabSize, legendIconSize, widthConnectors,
        cutoffLineWidth),
        0.0)
})

test_that("datatypes", {
    expect_type(EnhancedVolcano(toptable),
        c("S4", "list"))
    expect_type(EnhancedVolcano(x, y, lab, selectLab, col, colOverride, 
        legend, legendPosition, colConnectors, cutoffLineType,
        cutoffLineCol, border),
        c("character"))
    expect_type(EnhancedVolcano(xlab, ylab, title),
        c("language", "character"))
    expect_type(EnhancedVolcano(DrawConnectors, gridlines.major,
        gridlines.minor),
        c("logical"))
    expect_type(EnhancedVolcano(axisLabSize, pCutoff, pLabellingCutoff,
        FCcutoff, titleLabSize, transcriptPointSize, transcriptLabSize,
        colAlpha, legendLabSize, legendIconSize, widthConnectors,
        cutoffLineWidth, borderWidth),
        c("integer", "double"))
    expect_gt(EnhancedVolcano(axisLabSize, pCutoff, pLabellingCutoff,
        FCcutoff, titleLabSize, transcriptPointSize, transcriptLabSize,
        colAlpha, legendLabSize, legendIconSize, widthConnectors,
        cutoffLineWidth, borderWidth),
        0.0)
})

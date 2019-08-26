test_that("datatypes", {
  expect_type(EnhancedVolcano(toptable),
    c("S4", "list"))
  expect_type(EnhancedVolcano(lab, x, y, selectLab, shapeCustom, col,
    colCustom, legend, legendLabels, shade, shadeLabel),
    c("character"))
  expect_type(EnhancedVolcano(xlab, ylab, title, subtitle, caption,
    cutoffLineType, cutoffLineCol, transcriptLabCol, transcriptLabFace,
    legendPosition, shadeFill, typeConnectors, endsConnectors,
    colConnectors, hlineType, hlineCol, vlineType, vlineCol, border,
    borderColour),
    c("language", "character"))
  expect_type(EnhancedVolcano(boxedLabels, legendVisible, drawConnectors,
    gridlines.major, gridlines.minor),
     c("logical"))
  expect_type(EnhancedVolcano(subtitleLabSize, captionLabSize, pCutoff,
    cutoffLineWidth, pLabellingCutoff, FCcutoff, pointSize,
    labSize, labhjust, labvjust, shape,
    colAlpha, legendLabSize, legendIconSize, shadeAlpha, shadeSize,
    shadeBins, widthConnectors, lengthConnectors, hline, hlineWidth,
    vline, vlineWidth, borderWidth),
    c("integer", "double"))
  expect_gt(EnhancedVolcano(subtitleLabSize, captionLabSize, pCutoff,
    cutoffLineWidth, pLabellingCutoff, FCcutoff, transcriptPointSize,
    transcriptLabSize, labhjust, labvjust, shape,
    colAlpha, legendLabSize, legendIconSize, shadeAlpha, shadeSize,
    shadeBins, widthConnectors, lengthConnectors, hlineWidth, vlineWidth,
    borderWidth),
    0.0)
})

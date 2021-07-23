test_that('datatypes', {
  expect_type(EnhancedVolcano(toptable),
    c('S4', 'list'))
  expect_type(EnhancedVolcano(lab, x, y, pCutoffCol, selectLab, shapeCustom,
    col, colCustom, colGradient, colGradientLabels, legendLabels, encircle,
    encircleFill, encircleCol, shade), c('character'))
  expect_type(EnhancedVolcano(xlab, ylab, title, subtitle, caption,
    cutoffLineType, cutoffLineCol, legendPosition, shadeFill,
    typeConnectors, endsConnectors, colConnectors, hlineType,
    hlineCol, vlineType, vlineCol, border, borderColour, directionConnectors),
    c('language', 'character'))
  expect_type(EnhancedVolcano(boxedLabels, drawConnectors, arrowheads,
    gridlines.major, gridlines.minor, parseLabels),
     c('logical'))
  expect_type(EnhancedVolcano(subtitleLabSize, captionLabSize, pCutoff,
    cutoffLineWidth, FCcutoff, pointSize,
    labSize, labhjust, labvjust, shape,
    colAlpha, colGradientBreaks, colGradientLimits, legendLabSize,
    legendIconSize, encircleSize, encircleAlpha, shadeAlpha, shadeSize,
    shadeBins, widthConnectors, lengthConnectors, maxoverlapsConnectors, hline,
    hlineWidth, vline, vlineWidth, borderWidth, max.overlaps, min.segment.length),
    c('integer', 'double'))
  expect_gt(EnhancedVolcano(subtitleLabSize, captionLabSize, pCutoff,
    cutoffLineWidth, FCcutoff, labhjust, labvjust,
    shape, colAlpha, legendLabSize, legendIconSize, encircleSize, encircleAlpha,
    shadeAlpha, shadeSize, shadeBins, widthConnectors, lengthConnectors,
    hlineWidth, vlineWidth, borderWidth),
    0.0)
})


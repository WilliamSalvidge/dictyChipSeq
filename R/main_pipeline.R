
#' mainPipeLine
#'
#' @param x filepath to narrowPeak or broadPeak file
#' @param y filepath to 6 column un-normalised ax4 set1 readcounts file
#' @param bins number of bins
#'
#' @return
#' @export
#'
#' @examples
mainPipeLine = function(x, y, bins = 10) {

  ax4 = makeAx4TPM(y)
  ax4Bins = SplitIntoBins(ax4, bins)

  dictyPromoters = makeDictyPromotersGrangesObject()
  dictyGenes = makeDictyGenes()

  x = ChIPseeker::readPeakFile(x)

  if (grepl("narrowPeak$", x)) {

    x = add_narrowPeak_columnNames(x)

  }

  else {

    x = add_broadPeak_columnNames(x)

  }

  x = dropDictyChromosomes(x)
  x = annotatePeaksSeq(x)
  x = addExtraAnnotation(x)

  x = getTagMatricesForDiffExpLevelBins(x = ax4Bins,
                                        y = x,
                                        z = dictyPromoters)

  y = selectPromoterGRangesBasedOnTagMatrixIndex(x)

  y = addGeneIdToIndex(x = y,
                       y = dictyGenes)

  y = addExpValueAndSample(x = y,
                           y = ax4Bins,
                           z = ax4)

  y = rearrangeTagMatricesBasedOnSample(x = y,
                                        y = x)

  y = plotFunctionWillDictyChip(y)

  return (y)

}


#' mainPipeLine
#'
#' @param x filepath to narrowPeak or broadPeak file
#' @param y filepath to 6 column un-normalised ax4 set1 readcounts file
#' @param bins number of bins
#' @param subsetGenes vector of genes that you wish to perform operations on
#'
#' @return list
#' @export
mainPipeLine = function(x, y, bins = 10, subsetGenes = NULL) {

  mainPipeStartTime = Sys.time()

  ax4 = makeAx4TPM(y)
  ax4Bins = SplitIntoBins(ax4, bins)
  set1Degs = getSet1DEGS(y)

  dictyPromoters = makeDictyPromotersGrangesObject()
  dictyGenes = makeDictyGenes()

  x = ChIPseeker::readPeakFile(x)

  if (grepl("narrowPeak$", x)[1]) {

    x = add_narrowPeak_columnNames(x)

  }

  else {

    x = add_broadPeak_columnNames(x)

  }

  x = dropDictyChromosomes(x)
  x = annotatePeaksSeq(x)
  x = addExtraAnnotation(x)

  if (is.null(subsetGenes)) {

    x = getTagMatricesForDiffExpLevelBins(x = ax4Bins,
                                          y = x,
                                          z = dictyPromoters)

  }

  else if (is.vector(subsetGenes)) {

    # Convert dictyGenes GRange to a promoter GRange with geneID column
    promoterConversionGRange = promotersWithGeneId(x = dictyGenes$dictyGenesPlus,
                                               y = dictyGenes$dictyGenesNeg,
                                               promoterStart = 2500,
                                               promoterEnd = 2500)
    # subset the promoter GRange with genes of interest
    subsetPromoterGRange = promoterConversionGRange[promoterConversionGRange$geneId %in% subsetGenes, ]

    x = getTagMatricesForDiffExpLevelBins(x = ax4Bins,
                                          y = x,
                                          z = subsetPromoterGRange)

  }


  print("getTagMatricesForDiffExpLevelBins completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))


  y = selectPromoterGRangesBasedOnTagMatrixIndex(x)

  print("selectPromoterGRangesBasedOnTagMatrixIndex completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  y = addGeneIdToIndex(x = y,
                       y = dictyGenes)

  print("addGeneIdToIndex completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  y = addExpValueAndSample(x = y,
                           y = ax4Bins,
                           z = ax4)

  print("addExpValueAndSample completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  a = peakGRangesOrganisedBySample(y)

  print("peakGRangesOrganisedBySample completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  t = rearrangeTagMatricesBasedOnSample(x = y,
                                        y = x)

  print("rearrangeTagMatricesBasedOnSample completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  z = plotFunctionWillDictyChip(t)

  print("plotFunctionWillDictyChip completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  w = set1DegOverlapNumbers(a, set1_up = set1Degs[["set1_up"]], set1_down = set1Degs[["set1_down"]])

  print("set1DegOverlapNumbers completed")
  print(paste0("Current time elapsed = ", (Sys.time() - mainPipeStartTime), sep = ""))

  mainPipeEndTime = Sys.time()
  print(paste0("Total time elapsed for main pipeline = ", (mainPipeEndTime - mainPipeStartTime), sep = ""))

  return (list(plotData = z,
               CorrectTagMatrix = t,
               GRangesOrganisedBySample = a,
               Set1DegNumbers = w))

}

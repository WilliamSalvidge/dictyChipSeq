
#' selectPromoterGRangesBasedOnTagMatrixIndex
#'
#' @param x Output from getTagMatricesForDiffExpLevelBins which will likely be a list of tagMatrices and Windows / Promoters GRanges with index column
#'
#' @return A new list. Each element is based on one of the incoming tagMatrices in the input list. Each element in the new list is also a list. These elements show
#' the Promoter GRange object subsetted by the rownames of the tagMatrix. Where the rowname of the tagMatrix matches the index column in the Windows / Promoter GRange
#' object, that row is returned. This is also done in a strand specific way to generate + and - strand Promoter GRange objects.
#' @export
selectPromoterGRangesBasedOnTagMatrixIndex = function(x) {

  listOfIndices = list()

  for (i in names(x)) {

    listOfIndices[[i]] = list()
    listOfIndices[[i]][["Indices"]] = x[[i]]$Windows[x[[i]]$Windows$index %in% row.names(x[[i]]$tagMatrix), ]
    listOfIndices[[i]][["PlusIndices"]] = listOfIndices[[i]][["Indices"]][GenomicRanges::strand(listOfIndices[[i]][["Indices"]]) == "+", ]
    listOfIndices[[i]][["NegIndices"]] = listOfIndices[[i]][["Indices"]][GenomicRanges::strand(listOfIndices[[i]][["Indices"]]) == "-", ]

  }

  return (listOfIndices)

}

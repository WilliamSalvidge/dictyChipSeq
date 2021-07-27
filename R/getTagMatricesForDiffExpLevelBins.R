#' getTagMatricesForDiffExpLevelBins
#'
#' Takes in a list of dataframes. Was desgined to take in the different expression dataframes that are the output go [dictyChipSeqPackage::makeAx4TPM()].
#' It will then subset the GRanges object within the csAnno object anno slot and perform the tagMatrix generation on just that subset.
#'
#' This is actually a badly written function! When you subset based on the geneID annotation associated with high, medium or low expression
#' you do return GRanges objects that contain peaks that have been assigned to high, medium or lowly expressed genes, HOWEVER! getTagMatrix doesn't
#' care that each peak has only been assigned to 1 geneID. It just looks to see if that peak overlaps with the promoters of any genes and then returns
#' that in the tagMatrix. In the following functions we go in and take out any incorrect coverage!
#'
#' Each call to getTagMatrixWill will generate a Windows GRange with an index attached. Worthwhile checking thatt are all identical.
#'
#' It can also just accept a GRanges object in the y paramter in which case the output will contain only 1 tagMatrix and windows.
#'
#' @param x list of dataframes based on gene expression
#' @param y csAnno Object. OR if x is NULL a GRanges object.
#' @param z GRanges object of dicty promoters
#'
#' @return Output is a list with tagMatrices based on expression and Windows of promoters
#' @export
getTagMatricesForDiffExpLevelBins = function(x = NULL, y, z) {

  tagMatrixList = list()

  if (is.null(x)) {

    for (i in names(y)) {

      # y is GRangesList in this instance

      tagMatrixList[[i]] = getTagMatrixWill(y[[i]],
                                            windows=z)

      print(dim(tagMatrixList[[i]]$tagMatrix))

    }

  }

  else {

    for (i in names(x)) {

      # y is an annotated or csAnno object in this instance

      tagMatrixList[[i]] = getTagMatrixWill(y@anno[y@anno$geneId %in% row.names(x[[i]]), ],
                                            windows=z)

      print(dim(tagMatrixList[[i]]$tagMatrix))

    }

  }


  return(tagMatrixList)

}

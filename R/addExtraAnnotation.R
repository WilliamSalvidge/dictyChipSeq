
#' addExtraAnnotation
#'
#' This function looks at the annotation column which is the output of the ChIPseeker::annotatePeak functon.
#' I've wrapped the ChIPseeker::annotatePeak function in annotatePeaksSeq.
#' The annotation can be a number of things promoter, exon 1, exon 2 etc. downstream 1kb etc.
#'
#' The ChIPseeker::annotatePeak assigns a peak to a gene hierachically. First to a promoter, then exon, then downstream.
#' So if the annotation of a peak is exon 1, 2 or 3 of a gene, then the geneID present in the annotation column for said
#' peak should match that of the DDB_G0 in the geneID column.
#'
#' This function asks this question. It extracts the geneID from the annotation column using regex and the puts it in
#' a separate column. It then creates a second new column 'annotationGeneIDEqualGeneID' that tells us if the two geneIDs match.
#' If they do match the value is "Same geneID" if they don't "Different annotation".
#'
#' Other annotations such as promoter or downstream don't have a geneID in the original annotation column and so these
#' are given the value "Different annotation" in the 'annotationGeneIDEqualGeneID' column
#'
#' @param x csAnno object which is the output of annotatePeaksSeq
#'
#' @return csAnno object
#' @export
addExtraAnnotation = function(x) {

  # Create new metaData column based on extracting the geneID from the annotation (if annotation is exon)
  x@anno@elementMetadata$annotationGeneID =
    gsub("^.*/|,.*\\)$", "", x = x@anno@elementMetadata$annotation)

  # Create another new column which tells us whether the geneID for the exon is the same as assigned in the geneid column
  # by the AnnotatePeak function
  x@anno@elementMetadata$annotationGeneIDEqualGeneID = ifelse(
    grepl("^DDB_.*", x = x@anno@elementMetadata$annotationGeneID,
          perl = T),
    ifelse(x@anno@elementMetadata$annotationGeneID ==
             x@anno@elementMetadata$geneId,
           "Same geneID", "Different geneID"),
    "Different annotation"
  )

  print(data.frame(numberOfPeaks = dim(x@anno@elementMetadata)[1],
                   notExonAnnotation = length(x@anno@elementMetadata$annotationGeneIDEqualGeneID[
                     x@anno@elementMetadata$annotationGeneIDEqualGeneID == "Different annotation"]),
                   sameGeneId = length(x@anno@elementMetadata$annotationGeneIDEqualGeneID[
                     x@anno@elementMetadata$annotationGeneIDEqualGeneID == "Same geneID"]),
                   differeGeneId = length(x@anno@elementMetadata$annotationGeneIDEqualGeneID[
                     x@anno@elementMetadata$annotationGeneIDEqualGeneID == "Different geneID"])))

  return(x)

}

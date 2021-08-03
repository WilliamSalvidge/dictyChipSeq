#' getTagMatrixWill
#'
#' @param peak peak GRanges file. Will probably be in the anno GRange from the output of annotatePeaksSeq or addExtraAnnotation
#' @param weightCol ?
#' @param windows GRanges object of promoters - output from makeDictyPromotersGRangesObject()
#' @param flip_minor_strand flips tags so there all in the same direction. E.g if the output was 000000111111 and the strand is negative it will become 111111000000
#'
#' @return list of a tagMatrix, which will be used to make a plot and a GRanges object of windows / promoters.
#' This GRanges object now has the index column. The index column should map between the tagMatrix rowname and the value in this column.
#' This can be used to reverse engineer what promoter belongs to what gene.
#'
#' @export
getTagMatrixWill = function (peak, weightCol = NULL, windows, flip_minor_strand = TRUE)  {
  peak.gr <- ChIPseeker:::loadPeak(peak)
  if (!methods::is(windows, "GRanges")) {
    stop("windows should be a GRanges object...")
  }
  if (length(unique(GenomicRanges::width(windows))) != 1) {
    stop("width of windows should be equal...")
  }
  if (is.null(weightCol)) {
    peak.cov <- GenomicRanges::coverage(peak.gr)
  }
  else {
    weight <- S4Vectors::mcols(peak.gr)[[weightCol]]
    peak.cov <- GenomicRanges::coverage(peak.gr, weight = weight)
  }
  cov.len <- S4Vectors::elementNROWS(peak.cov)
  cov.width <- GenomicRanges::GRanges(seqnames = names(cov.len), IRanges::IRanges(start = rep(1,
                                                                      length(cov.len)), end = cov.len))
  # make GRanges of only promoters (windows) in the 6 chromosomes
  windows <- IRanges::subsetByOverlaps(windows, cov.width, type = "within",
                              ignore.strand = TRUE)
  # add index to windows (should be 13565)
  windows$index = seq(1, length(windows), by = 1)
  # character vector of 6 chromosome names
  chr.idx <- IRanges::intersect(names(peak.cov), unique(as.character(GenomicRanges::seqnames(windows))))

  # The subject is peak.cov which is the coverage 00000111110000 for example
  # against the promoter elements
  # Each promoter element is in there as a vector
  peakView <- IRanges::Views(peak.cov[chr.idx], methods::as(windows, "IntegerRangesList")[chr.idx])
  # Takes the vectors and puts them into a matrix of there own (one matrix per chromosome in a list)
  tagMatrixList <- lapply(peakView, function(x) t(IRanges::viewApply(x,
                                                            as.vector)))
  # bind the different matrices together (dim = 13565  5001)
  tagMatrix <- do.call("rbind", tagMatrixList)
  # length(windows) = 13565
  idx.list <- split(1:length(windows), as.factor(GenomicRanges::seqnames(windows)))
  idx <- do.call("c", idx.list)
  rownames(tagMatrix) <- idx
  tagMatrix <- tagMatrix[order(idx), ]
  if (flip_minor_strand) {
    minus.idx <- which(as.character(GenomicRanges::strand(windows)) == "-")
    tagMatrix[minus.idx, ] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
  }
  tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0, ]

  tagMatrixAndWindows = list(tagMatrix = tagMatrix,
                             Windows = windows)

  return(tagMatrixAndWindows)
}

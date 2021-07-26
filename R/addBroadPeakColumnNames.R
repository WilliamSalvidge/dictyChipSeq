

#' add broadPeak file column names to a GRanges object
#'
#' The chipSeeker readPeakFile function creates a GRanges object but the narrowPeak column names are missing.
#' This function adds them in.
#'
#' @param x GRanges object
#'
#' @return GRanges object
#' @export
add_broadPeak_columnNames = function (x) {
  # Set metadata column names to reflect narrowPeaks format
  # see http://genome.ucsc.edu/FAQ/FAQformat.html#format12
  # and https://pypi.org/project/MACS2/ for more info
  names(S4Vectors::mcols(x)) <- c("name", "score", "strand", "signalValue", "pValue", "qValue")
  return (x)

}

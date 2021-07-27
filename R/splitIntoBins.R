#' SplitIntoBins
#'
# YOU NEED TO SORT YOUR DATAFRAME BY EXPRESSION
# BEFORE RUNNING THIS FUNCTION
#'
#' @param x data frame to be split
#' @param z number - how many bins you wish to split into
#'
#' @return
#' @export
SplitIntoBins = function (x, z) {


  y = round(seq(1, dim(x)[1], by = (dim(x)[1] / z)))
  y = c(y, dim(x)[1])
  geneExpBins = list()
  for (i in 1:length(y)) {
    if (i == length(y)) {
      break()
    }
    else {
      geneExpBins[[ paste("Genes_", y[i], "_to_", (y[i + 1]), sep = "") ]] = x[y[i]:y[(i + 1)], ]
    }
  }

  return (geneExpBins)

}

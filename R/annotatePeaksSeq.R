

#' annotatePeaksSeq
#'
#' Takes in GRanges object of ChIP-seq peaks and annotates them based on a D.Discoideum TxDb object
#' coming from makeDictyGrangesfromTxDb.
#'
#' The seqnames returned by makeDictyGrangesfromTxDb are 1, 2, 3, 4, 5, 6 and so I change the input seqnames
#' to match before converting them back to what they were previously after.
#'
#' @param x GRanges object
#' @param tssRegionStart number : default is -3000bp - what is set by ChIPseeker::annotatePeak - can be changed!
#' @param tssRegionEnd number : default is 3000bp - what is set by ChIPseeker::annotatePeak - can be changed!
#'
#' @return
#' @export
annotatePeaksSeq = function(x, tssRegionStart = -3000, tssRegionEnd = 3000) {
  GenomeInfoDb::seqlevels(x) <- c("DDB0232428" = "1", "DDB0232429" = "2", "DDB0232430" = "3",
                    "DDB0232431" = "4", "DDB0232432" = "5", "DDB0232433" = "6")
  y <- ChIPseeker::annotatePeak(x, tssRegion=c(tssRegionStart, tssRegionEnd), TxDb=makeDictyGrangesfromTxDb())
  GenomeInfoDb::seqlevels(y@anno) = c("1" = "DDB0232428", "2" = "DDB0232429", "3" = "DDB0232430",
                        "4" = "DDB0232431", "5" = "DDB0232432", "6" = "DDB0232433")
  return(y)

}

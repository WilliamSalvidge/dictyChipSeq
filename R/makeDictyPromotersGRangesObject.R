
#' makeDictyPromotersGrangesObject
#'
#' Can change what is defined as a promoter by altering the Start and End parameters. I've basically wrapped the ChIPseeker::getPromoters function and then removed chromosomses 1 to 6.
#'
#' @param Start defaults to 2500 upstream
#' @param End defaults to 2500 downstream
#'
#' @return GRanges Object of D. Discoideum promoters for chromsomes 1 to 6
#' @export
makeDictyPromotersGrangesObject = function(Start = 2500, End = 2500) {

  dictyPromotersGRanges <- ChIPseeker::getPromoters(makeDictyGrangesfromTxDb(),
                                                    upstream = Start,
                                                    downstream = End)

  dictyPromotersGRanges = actualDropSeqLevels(dictyPromotersGRanges)

  return(dictyPromotersGRanges)

}







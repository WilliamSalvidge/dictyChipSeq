#' makeDictyGenesGRangesObject
#'
#' Calls makeDictyGRangesfromTxdb(). The gets the genes GRanges from that.
#'
#' @return A list of GRanges. First is a GRanges of all dicty genes.
#' Then just positive strand dicty genes and then negative strand dicty genes.
#' @export
makeDictyGenes = function() {

  dictyChromosomes1to6 = c("1", "2", "3",
                           "4", "5", "6")


  dictyGenesList = list()

  dicty_genes = GenomicFeatures::genes(makeDictyGrangesfromTxDb())

  dicty_genes = GenomeInfoDb::dropSeqlevels(dicty_genes,
                                           GenomeInfoDb::seqlevels(dicty_genes)[! ( GenomeInfoDb::seqlevels(dicty_genes) %in% dictyChromosomes1to6)],
                                           pruning.mode=c("coarse"))

  GenomeInfoDb::seqlevels(dicty_genes) <- c("1" = "DDB0232428",
                              "2" = "DDB0232429",
                              "3" = "DDB0232430",
                              "4" = "DDB0232431",
                              "5" = "DDB0232432",
                              "6" = "DDB0232433")


  dictyGenesList[["dictyGenes"]] = dicty_genes

  dictyGenesList[["dictyGenesPlus"]] = dicty_genes[GenomicRanges::strand(dicty_genes) == "+", ]

  dictyGenesList[["dictyGenesNeg"]] = dicty_genes[GenomicRanges::strand(dicty_genes) == "-", ]

  return (dictyGenesList)
}

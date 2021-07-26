#' makeDictyGenesGRangesObject
#'
#' Calls makeDictyGRangesfromTxdb(). The gets the genes GRanges from that. Then drops all but chromosome 1- 6 and returns a list of GRanges.
#'
#' @return A list of GRanges. First is a GRanges of all dicty genes.
#' Then just positive strand dicty genes and then negative strand dicty genes.
#' @export
makeDictyGenes = function() {

  dictyGenesList = list()

  dicty_genes = GenomicFeatures::genes(makeDictyGrangesfromTxDb())

  dicty_genes = actualDropSeqLevels(dicty_genes)

  dictyGenesList[["dictyGenes"]] = dicty_genes

  dictyGenesList[["dictyGenesPlus"]] = dicty_genes[GenomicRanges::strand(dicty_genes) == "+", ]

  dictyGenesList[["dictyGenesNeg"]] = dicty_genes[GenomicRanges::strand(dicty_genes) == "-", ]

  return (dictyGenesList)
}

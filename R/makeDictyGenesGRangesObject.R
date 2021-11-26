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



#' promotersWithGeneId
#'
#' @param x makeDictyGenes() $dictyGenesPlus
#' @param y makeDictyGenes() $dictyGenesNeg
#' @param promoterStart bp (in +ve numbers) upstream from TSS
#' @param promoterEnd bp (in +ve numbers) downstream from TSS
#'
#' @return GRange where makeDictyGenes is converted to a all Promoters With Gene Id GRange.
#' Use this in getTagMatricesForDiffExpLevelBins if you want to use only a subset of genes which overlap with peakFile.
#' @export
promotersWithGeneId = function (x, y, promoterStart = 2500, promoterEnd = 2500) {

  newGRangePlus = GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(x),
                                     ranges = IRanges::IRanges(start=GenomicRanges::start(x) - (promoterStart),
                                                               end= GenomicRanges::start(x) + promoterEnd,
                                                               names=x$gene_id),
                                     strand = GenomicRanges::strand(x),
                                     geneId = x$gene_id)

  newGRangeMinus = GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(y),
                                          ranges = IRanges::IRanges(start=GenomicRanges::end(y) - (promoterEnd),
                                                                    end= GenomicRanges::end(y) + promoterStart,
                                                                    names=y$gene_id),
                                          strand = GenomicRanges::strand(y),
                                          geneId = y$gene_id)

  return(c(newGRangePlus, newGRangeMinus))

}


#' addGeneIdToIndex
#'
#' Shrinks Genes GRange down to 3bp width around the TSS. The Genes GRAnge has a geneID associated with it.
#' Individual promoter GRanges are shrunk to 3bp around the TSS as well.
#' In theory the promoter and gene GRanges are the same. E.g.  DDB0232428 100-103 is the same in both dataframes.
#' We can now connect the index in the promoter GRange to a geneID!
#'
#' @param x output from selectPromoterGRangesBasedOnTagMatrixIndex()
#' @param y output from makeDictyGenesGRangesObject()
#'
#' @return List. Each List also a list containing two GRanges. The first GRange is the interesting part. Contains the Windows / Promoters GRange from
#' selectPromoterGRangesBasedOnTagMatrixIndex() but with GeneId added. You can chgeck against TestGRanges I think they should be the same.
#'
#' @export
addGeneIdToIndex = function(x, y) {

  newList = list()
  newList2 = list()

  testGRangePromsPlus = GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(y$dictyGenesPlus),
                                ranges = IRanges::IRanges(start=GenomicRanges::start(y$dictyGenesPlus),
                                               end= GenomicRanges::end(y$dictyGenesPlus) - (GenomicRanges::width(y$dictyGenesPlus) - 3),
                                               names=y$dictyGenesPlus$gene_id),
                                strand = GenomicRanges::strand(y$dictyGenesPlus),
                                geneId = y$dictyGenesPlus$gene_id)

  testGRangePromsNeg = GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(y$dictyGenesNeg),
                               ranges = IRanges::IRanges(start = GenomicRanges::end(y$dictyGenesNeg),
                                              end = GenomicRanges::end(y$dictyGenesNeg) + 3,
                                              names = y$dictyGenesNeg$gene_id),
                               strand = GenomicRanges::strand(y$dictyGenesNeg),
                               geneId = y$dictyGenesNeg$gene_id)


  testGRangeProms = c(testGRangePromsPlus, testGRangePromsNeg)

  for (i in names(x)) {


    testGRange = GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(x[[i]]$Indices),
                         ranges = IRanges::IRanges(start= GenomicRanges::start(x[[i]]$Indices) + round((GenomicRanges::width(x[[i]]$Indices)[1] -1) /2),
                                        end = GenomicRanges::end(x[[i]]$Indices) - (round((GenomicRanges::width(x[[i]]$Indices)[1] -1) /2) - 3),
                                        names = x[[i]]$Indices$index),
                         strand = GenomicRanges::strand(x[[i]]$Indices))



    testOvlp = IRanges::subsetByOverlaps(testGRangeProms, testGRange)

    testOvlp = GenomeInfoDb::sortSeqlevels(testOvlp)
    testOvlp = sort(testOvlp)

    testGRange = GenomeInfoDb::sortSeqlevels(testGRange)
    testGRange = sort(testGRange)

    testOvlp$index = S4Vectors::ROWNAMES(testGRange)

    newList[[i]] = testOvlp
    newList2[[i]] = testGRange

  }

  newList3 = list(OverlapGRanges = newList, TestGRanges = newList2)
  return (newList3)

}

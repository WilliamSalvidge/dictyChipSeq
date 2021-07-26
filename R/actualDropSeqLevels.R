actualDropSeqLevels = function(x) {

  dictyChromosomes1to6 = c("1", "2", "3",
                           "4", "5", "6")

  GRangesChr1To6 = GenomeInfoDb::dropSeqlevels(x,
                                               GenomeInfoDb::seqlevels(x)[! ( GenomeInfoDb::seqlevels(x) %in% dictyChromosomes1to6)],
                                               pruning.mode=c("coarse"))

  GenomeInfoDb::seqlevels(GRangesChr1To6) <- c("1" = "DDB0232428",
                                               "2" = "DDB0232429",
                                               "3" = "DDB0232430",
                                               "4" = "DDB0232431",
                                               "5" = "DDB0232432",
                                               "6" = "DDB0232433")


  return(GRangesChr1To6)

}

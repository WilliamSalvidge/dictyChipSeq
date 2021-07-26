
#' dropDictyChromsomes
#'
#' Function to drop any chromsomes which are not d. discoideum chromsomes 1 through 6
#' Extra chromosomes and contigs - ribosomal, mitochondrial etc. may be interesinting but this
#' is just for simplicity
#'
#' @param x Granges object
#'
#' @return Granges object
#' @export
dropDictyChromosomes = function(x) {

  dictyChromosomes1to6 = c("DDB0232428", "DDB0232429", "DDB0232430",
                           "DDB0232431", "DDB0232432", "DDB0232433")


  return (dropSeqlevels(x, seqlevels(x)[ ! ( seqlevels(x) %in% dictyChromosomes1to6)], pruning.mode=c("coarse")))

}



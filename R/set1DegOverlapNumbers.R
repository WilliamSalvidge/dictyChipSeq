#' Title
#'
#' @param x List of rearranged GRange objects. Output from The sample column should all be the same e.g. 1 / 2 ... 9. If its a mixture this isn't the right list
#' @param set1_up dataframe
#' @param set1_down dataframe
#'
#' @return
#' @export
set1DegOverlapNumbers = function (x, set1_up, set1_down) {

  set1UpRunningTotal = 0
  set1DownRunningTotal = 0

  uniqueGenesWithPeaksInBin = numeric()
  set1UpInBin = numeric()
  set1DownInBin = numeric()
  set1UpInBinGenes = character()
  set1DownInBinGenes = character()

  set1UpInBinsList = list()
  set1DownInBinsList = list()

  for (i in 1:length(x)) {

    print(paste0("# Number of unique genes in bin ", i, sep = ""))
    print(length(unique(x[[i]])))
    uniqueGenesWithPeaksInBin = c(uniqueGenesWithPeaksInBin, length(unique(x[[i]])))

    print(paste0("# Number of unique set1 up genes in bin ", i, sep = ""))
    print(length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_up), ])))
    set1UpRunningTotal = set1UpRunningTotal + length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_up), ]))
    set1UpInBin = c(set1UpInBin, length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_up), ])))
    set1UpInBinGenes = c(set1UpInBinGenes, unique(x[[i]][x[[i]]$geneId %in% row.names(set1_up), ]$geneId))
    set1UpInBinsList[[paste0("Bin ", i, sep = "")]] = unique(x[[i]][x[[i]]$geneId %in% row.names(set1_up), ]$geneId)

    print(paste0("# Number of unique set1 down genes in bin ", i, sep = ""))
    print(length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_down), ])))
    set1DownRunningTotal = set1DownRunningTotal + length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_down), ]))
    set1DownInBin = c(set1DownInBin, length(unique(x[[i]][x[[i]]$geneId %in% row.names(set1_down), ])))
    set1DownInBinGenes = c(set1DownInBinGenes, unique(x[[i]][x[[i]]$geneId %in% row.names(set1_down), ]$geneId))
    set1DownInBinsList[[paste0("Bin ", i, sep = "")]] = unique(x[[i]][x[[i]]$geneId %in% row.names(set1_down), ]$geneId)

  }

  set1DataFrame = data.frame(uniqueGenesWithPeaksInBin = uniqueGenesWithPeaksInBin,
                             set1UpInBin = set1UpInBin,
                             set1DownInBin = set1DownInBin)

  set1UpInBinPVal = numeric()
  set1DownInBinPVal = numeric()

  for (i in 1:length(set1DataFrame$uniqueGenesWithPeaksInBin)) {

    set1UpInBinPVal = c(set1UpInBinPVal,
                        stats::phyper(set1DataFrame[i, "set1UpInBin"],
                               set1DataFrame[i, "uniqueGenesWithPeaksInBin"],
                               (sum(set1DataFrame$uniqueGenesWithPeaksInBin) -  set1DataFrame[i, "uniqueGenesWithPeaksInBin"]),
                               sum(set1DataFrame$set1UpInBin), lower.tail = F))

    set1DownInBinPVal = c(set1DownInBinPVal,
                          stats::phyper(set1DataFrame[i, "set1DownInBin"],
                                 set1DataFrame[i, "uniqueGenesWithPeaksInBin"],
                                 (sum(set1DataFrame$uniqueGenesWithPeaksInBin) -  set1DataFrame[i, "uniqueGenesWithPeaksInBin"]),
                                 sum(set1DataFrame$set1DownInBin), lower.tail = F))


  }

  set1DataFrame$UpPVals = set1UpInBinPVal
  set1DataFrame$DownPVals = set1DownInBinPVal

  print("# all set1 Up genes")
  print(set1UpRunningTotal)
  print("# all set1 Down genes")
  print(set1DownRunningTotal)


  return (list(set1DataFrane = set1DataFrame,
               set1UpInBinGenes = set1UpInBinGenes,
               set1DownInBinGenes=set1DownInBinGenes,
               set1UpIndividualBins = set1UpInBinsList,
               set1DownIndividualBins = set1DownInBinsList))

}

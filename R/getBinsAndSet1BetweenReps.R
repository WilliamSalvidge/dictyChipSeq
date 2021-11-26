getBinsAndSet1BetweenReps = function(x, y) {

  coerceList = list()
  coerceListSet1Up = list()
  coerceListSet1Down = list()

  for (i in 1:length(x$GRangesOrganisedBySample)) {


    coerceList[[i]] = unique(c(x$GRangesOrganisedBySample[[i]]$geneId, y$GRangesOrganisedBySample[[i]]$geneId))

    coerceListSet1Up[[i]] = unique(c(x$Set1DegNumbers$set1UpIndividualBins[[i]],
                                     y$Set1DegNumbers$set1UpIndividualBins[[i]]))

    coerceListSet1Down[[i]] = unique(c(x$Set1DegNumbers$set1DownIndividualBins[[i]],
                                       y$Set1DegNumbers$set1DownIndividualBins[[i]]))


  }

  uniqueGenesWithPeaksInBin = numeric()
  set1UpInBin = numeric()
  set1DownInBin = numeric()

  for (i in 1:length(coerceList)) {

    uniqueGenesWithPeaksInBin = c(uniqueGenesWithPeaksInBin, length(coerceList[[i]]))
    set1UpInBin = c(set1UpInBin, length(coerceListSet1Up[[i]]))
    set1DownInBin = c(set1UpInBin, length(coerceListSet1Down[[i]]))
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

  return (set1DataFrame)

}

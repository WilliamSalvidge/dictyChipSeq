plotFunctionWillDictyChip = function(x) {

  plotList = list()

  for (i in 1:length(x)) {

    # Dig down into the plotAvgProf function
    # plotAvgProf calls ChIPseeker:::plotAvgProf.internal
    # this calls ChIPseeker:::getTagCount

    # getTagCount has the following
    # ss <- colSums(tagMatrix)
    # ss <- ss/sum(ss)
    # dd <- data.frame(pos = c(xlim[1]:xlim[2]), value = ss)

    profilePlotTemp = ChIPseeker::plotAvgProf(x[[i]], xlim=c(-2500, 2500),
                                  xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

    profilePlotTempData = profilePlotTemp$data
    profilePlotTempData$sample = i

    plotList[[i]] = profilePlotTempData

    print(paste("tagMatrix", i, "complete"))


  }

  plotTable = do.call(rbind, plotList)
  return (plotTable)

}

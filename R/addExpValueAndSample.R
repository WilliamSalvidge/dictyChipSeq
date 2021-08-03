#' addExpValueAndSample
#'
#' Adds mean expression value and sample based on expression level.
#' Internally calls lApplyFunction(), bin_search() & addSampleToNewList()
#'
#' @param x output from addGeneIdToIndex()
#' @param y output from splitIntoBins(makeAx4TPM())
#' @param z output from makeAx4TPM()
#'
#' @return Modified form of addGeneIdToIndex() output where mean and sample columns are added.
#' @export
addExpValueAndSample = function (x, y, z) {

  newList = list()
  timeStart = Sys.time()

  lApplyOutput = lapply(x$OverlapGRanges, lApplyFunction, z=z)

  newListMod = addSampleToNewList(lApplyOutput, y)

  timeEnd = Sys.time()
  print("Total time elapsed")
  print(timeEnd - timeStart)

  return (newListMod)

}

lApplyFunction = function(w, z) {

  roundStart = Sys.time()
  print("Round started")
  testEnd = w
  testEnd$mean = 0
  testEndGenes = z[z$geneId %in% testEnd$geneId, ]

  # order by geneId and then use binary search
  # should be ascending
  newData <- testEndGenes[order(testEndGenes$geneId),]


  for (j in 1:length(testEnd$geneId)) {

    tempResult = bin_search(newData$geneId, testEnd$geneId[j])
    if (tempResult != 0) {
      testEnd$mean[j] = newData$mean[tempResult]
    }


  }

  roundEnd = Sys.time()
  print(paste0("This round took ", (roundEnd - roundStart), sep = ""))
  return(testEnd)

}

bin_search = function(vectorInput, target) {

  # target would be testEnd$geneId[i]
  # vectorInput would be sorted testEndGenes$mean or newData$geneId

  lo <- 1; hi <- length(vectorInput)
  while (lo <= hi) {
    mid <- as.integer(round((lo + hi) / 2)) # always even!
    if (vectorInput[mid] == target) {
      return(mid)
    } else if (vectorInput[mid] < target) {
      lo <- mid + 1
    } else {
      hi <- mid - 1
    }
  }
  return(0)
}

addSampleToNewList = function(x, y) {

  newVector = vector()

  for (i in 1:length(y)) {

    newVector = c(newVector, max(y[[i]]$mean))

  }
  newVector = c(0, newVector)
  for (k in 1:length(x)) {
    x[[k]]$sample = 0
    for (i in 1:length(newVector)) {
      if (i < length(newVector)) {
        x[[k]][x[[k]]$mean > newVector[i] & x[[k]]$mean <= newVector[i + 1], ]$sample = i
      }
    }
  }

  return(x)

}

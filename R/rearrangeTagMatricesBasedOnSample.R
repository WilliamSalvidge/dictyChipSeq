
#' rearrangeTagMatricesBasedOnSample
#'
#' @param x output of addExpValueAndSample()
#' @param y output of getTagMatricesForDiffExpLevelBins()
#'
#' This is a list of GRanges objects based on the tagMatrices from getTagMatricesForDiffExpLevelBins()
#'
#' @return
#' @export
#'
#' @examples
rearrangeTagMatricesBasedOnSample = function(x, y) {

  # They have the index and mean as metadata columns.
  # The sample metadata column is based expression level
  # y =

  timeStart = Sys.time()
  counter = 0
  emptyList = list()
  for (i in 1:length(x)) {
    emptyList[[i]] = vector()
  }

  for (i in 1:length(x)) {
    timeCheckStart = Sys.time()


    y[[i]]$tagMatrix = cbind(y[[i]]$tagMatrix, row.names(y[[i]]$tagMatrix))
    print(paste0("Round ", i, " started", sep =""))

    for (k in 1:length(x)) {

      for (j in 1:length(x[[i]])) {
        counter = counter + 1

        if (counter %% 500 == 0) {
          print(counter)
        }

        if (x[[i]]$sample[j] == k) {

          emptyList[[k]] = rbind(emptyList[[k]], y[[i]]$tagMatrix[row.names(y[[i]]$tagMatrix) == x[[i]]$index[j], ])

        }


      }


    }

    timeCheck = Sys.time()
    print(paste0("Time since start: ", timeCheck - timeStart, sep = ""))
    print(paste0("This round took: ", timeCheck - timeCheckStart, sep = ""))

    print(dim(emptyList[[i]]))


  }

  for (i in 1:length(emptyList)) {

    row.names(emptyList[[i]]) = emptyList[[i]][,ncol(emptyList[[i]])]
    emptyList[[i]] = emptyList[[i]][ , 1:(ncol(emptyList[[i]])-1)]
    emptyList[[i]] = t(apply(emptyList[[i]], 1, as.numeric))


  }

  timeStop = Sys.time()
  print(paste0("Total time elapsed: " , (timeStop - timeStart), sep = ""))

  # Output is the correct tagMatrices!

  return(emptyList)
}

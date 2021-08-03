#' rearrangeTagMatricesBasedOnSample
#'
#' @param x output of addExpValueAndSample()
#'
#' This is a list of GRanges objects based on the tagMatrices from getTagMatricesForDiffExpLevelBins().They have the index, geneId, mean, sample as metadata columns.
#' The sample metadata column is based expression level
#'
#' @param y output of getTagMatricesForDiffExpLevelBins()
#'
#' @return output is the correct tagMatrices!
#' @export
rearrangeTagMatricesBasedOnSample = function(x, y) {

  timeStart = Sys.time()
  emptyList = list()

  # create new emptyList same length as x
  for (i in 1:length(x)) {
    emptyList[[i]] = vector()
  }

  # Search through each GRange object in list of GRange objects
  # For each GRange object cycle through every row and see what the sample number is
  # Move this row to right element of emptyList
  # E.g. if row 1000 belongs to sample 9 move this row to emptyList[[9]]
  # if row 1001 belongs to sample 2 move this row to emptyList[[2]]
  for (i in 1:length(x)) {
    timeCheckStart = Sys.time()
    counter = 0
    # make the last column of each tagMatrix the rownames. Will use these for subsetting within this for loop.
    y[[i]]$tagMatrix = cbind(y[[i]]$tagMatrix, row.names(y[[i]]$tagMatrix))
    print(paste0("Round ", i, " started", sep =""))

    for (k in 1:length(x)) {

      for (j in 1:length(x[[i]])) {
        counter = counter + 1

        if (counter %% 2000 == 0) {
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

    # make rownames whatever the last column. This is the
    row.names(emptyList[[i]]) = emptyList[[i]][,ncol(emptyList[[i]])]
    # remove the last row
    emptyList[[i]] = emptyList[[i]][ , 1:(ncol(emptyList[[i]])-1)]
    # turn every row into a numeric vector and then transform. So that indices are column names?
    emptyList[[i]] = t(apply(emptyList[[i]], 1, as.numeric))

  }

  timeStop = Sys.time()
  print(paste0("Total time elapsed: " , (timeStop - timeStart), sep = ""))

  # Output is the correct tagMatrices!

  return(emptyList)
}

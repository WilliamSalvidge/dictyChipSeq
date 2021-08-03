
peakGRangesOrganisedBySample = function(x) {

  timeStart = Sys.time()

  emptyList= list()
  t = 1
  for (w in 1:length(x)) {


    emptyList[[w]] = x[[t]][x[[t]]$sample == w, ]

  }


  for (i in 2:length(x)) {
    timeCheckStart = Sys.time()

    for (j in 1:length(x)) {

      emptyList[[j]] = c(emptyList[[j]], x[[i]][x[[i]]$sample == j, ])

    }


  }

  timeCheck = Sys.time()
  print(paste0("Time since start: ", timeCheck - timeStart, sep = ""))
  print(paste0("This round took: ", timeCheck - timeCheckStart, sep = ""))
  timeStop = Sys.time()
  print(paste0("Total time elapsed: ", (timeStop - timeStart), sep = ""))
  return(emptyList)
} # Give you the genes



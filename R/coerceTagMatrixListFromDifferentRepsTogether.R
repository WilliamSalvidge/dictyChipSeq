
coerceTagMatrixListFromDifferentRepsTogether = function(x, y) {
  coerceList = list()

   for (i in 1:length(x)) {

     coerceList[[i]] = rbind(x[[i]], y[[i]])

   }

  return (coerceList)

}

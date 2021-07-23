

add_narrowPeak_columnNames = function (x) {
  # Set metadata column names to reflect narrowPeaks format
  # see http://genome.ucsc.edu/FAQ/FAQformat.html#format12
  # and https://pypi.org/project/MACS2/ for more info
  names(mcols(x)) <- c("name", "score", "strand", "signalValue", "pValue", "qValue", "Position relative to peak start")
  return (x)

}

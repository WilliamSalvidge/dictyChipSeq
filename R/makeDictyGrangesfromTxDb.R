
#' makeDictyGrangesfromTxDb
#'
#' Not totally sure how this works so maybe come back and edit later.
#' using the biomaRt package we establish to a connection to ensembl protist database I think.
#' We then get the D. discoideum TxDb or transcript database from there.
#' What we get is a TxDb object. We can call things like the GenomicFeatures::genes function on this to get a GRanges of genes.
#'
#' @return TxDb object
#' @export
makeDictyGrangesfromTxDb = function() {

  #mart <- biomaRt::useMart(biomart="protists_mart",
                           #host="protists.ensembl.org") # make connection to db

  #biomaRt::listDatasets(mart) # ddiscoideum_eg_gene shows up in this list

  dicty_txdb <- GenomicFeatures::makeTxDbFromBiomart(biomart="protists_mart",
                                    dataset="ddiscoideum_eg_gene",
                                    host="protists.ensembl.org") # extract dicty data

  return (dicty_txdb)

}

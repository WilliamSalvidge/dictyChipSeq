% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeDictyGrangesfromTxDb.R
\name{makeDictyGrangesfromTxDb}
\alias{makeDictyGrangesfromTxDb}
\title{makeDictyGrangesfromTxDb}
\usage{
makeDictyGrangesfromTxDb()
}
\value{
TxDb object
}
\description{
Not totally sure how this works so maybe come back and edit later.
using the biomaRt package we establish to a connection to ensembl protist database I think.
We then get the D. discoideum TxDb or transcript database from there.
What we get is a TxDb object. We can call things like the GenomicFeatures::genes function on this to get a GRanges of genes.
}

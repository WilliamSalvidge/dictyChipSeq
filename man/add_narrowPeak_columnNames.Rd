% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addNarrowPeakColumnNames.R
\name{add_narrowPeak_columnNames}
\alias{add_narrowPeak_columnNames}
\title{add narrowPeak file column names to a GRanges object}
\usage{
add_narrowPeak_columnNames(x)
}
\arguments{
\item{x}{Granges object}
}
\value{
Granges object
}
\description{
The chipSeeker readPeakFile function creates a GRanges object but the narrowPeak column names are missing.
This function adds them in.
}

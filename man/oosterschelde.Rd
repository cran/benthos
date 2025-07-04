% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_oosterschelde.R
\docType{data}
\name{oosterschelde}
\alias{oosterschelde}
\title{Oosterschelde Marine Benthos Data}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 4269 rows and 8 columns.
}
\source{
Rijkswaterstaat Water, Transport and Living Environment,
     Department of Information Management, Lelystad, 
     The Netherlands (contact: servicedesk-data@rws.nl)
}
\usage{
oosterschelde
}
\description{
Oosterschelde data set. The Oosterschelde is located in the southwest of
the Netherlands.
}
\details{
The Oosterschelde data contains the following columns:
\itemize{
  \item ID sample identifier
  \item HABITAT specification of the habitat
  \item AREA sampled area
  \item DATE sampling date (YYYY-MM-DD, ISO 8601)
  \item TAXON standardized taxon code (see WoRMS-website 
     \url{https://www.marinespecies.org/})
  \item COUNT number of individuals of 'TAXON'
}
}
\note{
This is not the original data set, but a simplified version of it meant
     for didactic purposes only! For instance it only contains taxa
     identified at the species level. Other taxa have been removed.
}
\keyword{datasets}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{ambi}
\alias{ambi}
\alias{ambi_}
\alias{has_ambi}
\alias{has_ambi_}
\title{AZTI Marine Biotic Index (AMBI)}
\usage{
ambi(.data = NULL, taxon, count, group = NULL)

ambi_(.data = NULL, taxon, count, group = NULL)

has_ambi(.data = NULL, taxon, group = NULL)

has_ambi_(.data = NULL, taxon, group = NULL)
}
\arguments{
\item{.data}{data in a \code{data.frame}, \code{tibble}, 
\code{data.table}, database etc.}

\item{taxon}{species names}

\item{count}{counts of individuals (\code{numeric})}

\item{group}{sensitivity groups I, II, III, IV, or V}
}
\value{
numeric vector of length 1 containing the AMBI
}
\description{
AZTI Marine Biotic Index (AMBI) according to Borja et al. (2000)
}
\details{
The index is given by:
 \deqn{c_\mathrm{b} = \frac{3}{2} \sum_{i=2}^5 (i-1) p_i}
 where \eqn{p_i} is the proportion of species in sensitivity group \eqn{i}.
}
\section{Functions}{
\itemize{
\item \code{ambi_()}: version suitable for calling from a function
(see package \pkg{lazyeval}).

\item \code{has_ambi()}: tests if an AMBI sensitivity group is available for \code{taxon}
(returns \code{TRUE} (available) or \code{FALSE} (unavailable))

\item \code{has_ambi_()}: version suitable for calling from a function
(see package \pkg{lazyeval}).

}}
\examples{
 ambi(
     taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
     count = c(4, 6)
 )
 
     data(oosterschelde)
     has_ambi(oosterschelde, TAXON)
 
}
\references{
Borja, A., J. Franco and V. Perez, 2000. A Marine Biotic Index 
 to Establish the Ecological Quality of Soft-Bottom Benthos Within 
 European Estuarine and Coastal Environments. 
 Marine Pollution Bulletin 40:1100-1114
}

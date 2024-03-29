% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{iti}
\alias{iti}
\alias{iti_}
\alias{has_iti}
\alias{has_iti_}
\title{Infaunal Trophic Index (ITI)}
\usage{
iti(.data = NULL, taxon, count, group = NULL)

iti_(.data = NULL, taxon, count, group = NULL)

has_iti(.data = NULL, taxon, group = NULL)

has_iti_(.data = NULL, taxon, group = NULL)
}
\arguments{
\item{.data}{data in a \code{data.frame}, \code{tibble}, 
\code{data.table}, database etc.}

\item{taxon}{species names}

\item{count}{counts of individuals (\code{numeric})}

\item{group}{sensitivity groups I, II, III, or IV}
}
\value{
numeric vector of length 1 containing the ITI
}
\description{
Computes the Infaunal Trophic Index (ITI) according to 
Gittenberger & van Loon (2013).
}
\details{
The Infaunal Trophic Index (ITI) is given by
 \deqn{\mathrm{ITI} = 100 \sum_{i=1}^3 \frac{(4-i)}{3} p_i}
 where \eqn{p_i} is the proportion of species in class \eqn{i}, where
 \itemize{
     \item group   I are suspension feeders (highest quality);
     \item group  II are interface feeders
     \item group III are surface deposit feeders and
     \item group  IV are subsurface deposit feeders (lowest quality). 
 }
}
\section{Functions}{
\itemize{
\item \code{iti_()}: version suitable for calling from a function
(see package \pkg{lazyeval}).

\item \code{has_iti()}: tests if an ITI sensitivity group is available for \code{taxon}
(returns \code{TRUE} (available) or \code{FALSE} (unavailable))

\item \code{has_iti_()}: version suitable for calling from a function
(see package \pkg{lazyeval}).

}}
\examples{
     iti(taxon = c("Euspira pulchella", "Nephtys cirrosa"), count = c(4, 6))
     
 
   data(oosterschelde)
   has_iti(oosterschelde, TAXON)
 
}
\references{
Gittenberger A. and  W. van Loon, 2013. 
     Sensitivities of marine macrozoobenthos to environmental pressures 
     in the Netherlands. Nederlandse Faunistische 
     Mededelingen 41: 79-112.
}

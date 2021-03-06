% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{eqr}
\alias{eqr}
\title{Ecological Quality Ratio (EQR)}
\usage{
eqr(x, bad, ref)
}
\arguments{
\item{x}{numeric vector containing benthic indices}

\item{bad}{the value for a bad status}

\item{ref}{the value for a reference status}
}
\value{
numeric vector with EQR-values: low values indicate bad ecological
     quality and high values indicate good ecological quality.
}
\description{
The ecological quality ratio is the ratio beween a parameter value and its
reference value:  
\deqn{EQR = \frac{x-bad}{ref-bad}}{EQR = (x-bad)/(ref-bad)}
Depending on \code{bad} and \code{ref}, the EQR usually 
(but not necessarily!) varies between 0 (bad ecological quality) and 1 (
ecological quality equals the reference status).
}

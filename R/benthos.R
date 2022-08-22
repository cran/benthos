#' Marine Benthic Ecosystem Analysis
#'
#' \pkg{benthos} provides functions for facilitating the analysis of marine
#' benthos data. Examples are indicators like species abundance, 
#' species richness, Margalef's index of diversity, Shannon's Entropy, 
#' AZTI's Marine Biotic Index, and the Infaunal Trophic Index (ITI). 
#' In addition functions for data pooling,
#' genus-to-species conversion and validation and conversion of species names 
#' to those recommended by the World Register of Marine Species are provided.
#'
#' All functions are designed to work seamlessly with the
#' \pkg{dplyr}-package which implements a grammar for structured
#' data manimpulation.
#'
#' The \pkg{benthos}-package contains functions for estimating various
#' species abundance, species richness, species heterogeneity and species
#' sensitivity measures:
#' \itemize{
#'      \item total abundance (\code{\link{total_abundance}})
#'      \item abundance (\code{\link{abundance}})
#'      \item species richness (\code{\link{species_richness}})
#'      \item Margalef's index of diversity (\code{\link{margalef}})
#'      \item Rygg's index of diversity (\code{\link{rygg}})
#'      \item Hurlbert's Expected Number of Species (\code{\link{hurlbert}})
#'      \item Simpson's measure of concentration (\code{\link{simpson}})
#'      \item Hurlbert's probability of interspecific encounter (PIE) (\code{\link{hpie}})
#'      \item Shannon's index or entropy (\code{\link{shannon}})
#'      \item Hill's diversity number (\code{\link{hill}})
#'      \item AZTI Marine Biotic Index (AMBI) (\code{\link{ambi}})
#'      \item Infaunal Trophic Index (ITI) (\code{\link{iti}})
#'      \item Bray-Curtis dissimilarity (\code{\link{bray_curtis}})
#'  }
#'
#' In addition, functions are available for data preparation, e.g.:
#' \itemize{
#'      \item data pooling (\code{\link{pool}})
#'      \item genus to species conversion (\code{\link{genus_to_species}})
#' }
#'
#' For an overview of all the functions in the package click on the index link
#' at the bottom of this page.
#'
#' @author Dennis Walvoort \email{dennis.walvoort@@wur.nl}
#'
#' @seealso The \pkg{BEQI2}-package on CRAN, and the package vignettes.
"_PACKAGE"


#' Remove Redundant Spaces
#'
#' This function removes redundant spaces from character vectors
#'
#' @param x character vector
#'
#' @return  character vector without trailing or multiple spaces
strip_spaces <- 
function(x) {
    if (is.character(x)) {
    	x <- gsub(pattern = " {2,}",   replacement = " ", x = x)
    	x <- gsub(pattern = "^ +| +$", replacement = "",  x = x)
    }
	x
}



#' Test for Azoic Samples
#'
#' Case-insensitive test for taxa starting with 'azoi'
#'     
#' @details Azoic samples need special attention during data analysis. 
#'      They should be marked as 'azoic', and taken care of during
#'      analysis. Note that an azoic sample is not the same as a
#'      record where a taxon has zero counts. The latter should be
#'      removed from further analysis, whereas the former provides
#'      important information.
#'
#' @param x character vector containing taxa
#'
#' @return  logical vector, with elements \code{TRUE} for azoic samples, 
#'      and \code{FALSE} otherwise.
#'      
#' @export
is_azoic <- 
function(x) {
    grepl(pattern = "^Azoi", x = x, ignore.case = TRUE) 
}



#' Pooling
#'
#' This function randomly assigns samples to pools of approximately equal area
#'
#' @param sample_id sample identifier
#' @param area sampling area of \code{sample_id} (in the same units as 
#'     \code{target_area})
#' @param target_area vector of length 2 containing the lower and upper bound 
#'     of the pooled area (same units as \code{area})
#' @param max_try maximum number of unsuccessful pooling tries before the 
#'     algorithm gives up.
#'
#' @return vector with idenitifiers (integers) indicating the pool to which 
#'     each sample belongs (NA for samples that could not be pooled)
#' 
#' @export
pool <-
function(sample_id = 1:length(area), area, target_area, max_try = 100L) {

    # make sure that target_area is an ordered vector of length 2
    target_area <- target_area %>% range

    # construct tibble
    stopifnot(length(area) == length(sample_id))
    d_full <- tibble(sample_id, area)
    if (nrow(d_full) == 0L) {
        stop("no samples found", call. = FALSE)
    }

    # make sure all sampling units only occur once
    d_uniq <- d_full %>% distinct
    
    # test uniqueness of area for each sample_id
    if (anyDuplicated(d_uniq$sample_id)) {
        stop("Each sampling unit should have a unique area", call. = FALSE)
    }
    
    # start pooling
    d_uniq$pool_id <- .pool(sample_id = d_uniq$sample_id, area = d_uniq$area, 
                            target_area = target_area, max_try = max_try)
    
    # decompression
    d_full <- d_full %>% 
        select(sample_id) %>% 
        left_join(d_uniq, by = "sample_id")

    # return pool identifiers
    d_full$pool_id
    
}


#' @describeIn pool internal function not supposed to be called directly.
.pool <-
function(sample_id = 1:length(area), area, target_area, max_try = 100L) {

    # determine number of samples
    n_samples <- length(sample_id)
    
    # initialize pool identifier
    pool_id <- rep.int(x = NA_integer_, times = n_samples)
    
    # initialize sample index
    sample_index <- seq_len(n_samples)
    
    # check trivial cases
    big_area <- area > target_area[2]
    if (any(big_area)) {
        sample_index <- setdiff(sample_index, which(big_area))
    }
    if (all(area %>% between(target_area[1], target_area[2]))) {
        message("No pooling necessary")
        return(sample_index)
    }
    if (all(area > target_area[2])) {
        stop("No pooling possible", call. = FALSE)
        return(NULL)
    }
    s <- sum(area)
    if (s < target_area[1]) {
        stop("No pooling possible", call. = FALSE)
        return(NULL)
    }
    if (s < target_area[2]) {
        pool_id[] <- 1L
        return(pool_id)
    }
    

    # estimate the (theoretical) maximum number of samples in a pool
    # Note: the theoretical minimum is not necessary as cumsum is 
    # also computed in the repeat-loop below.
    pooled_area <- cumsum(sort(area))
    index <- which(pooled_area > target_area[2])
    max_samples <- min(index)
                  
    # initialize number of pools
    n_pools <- 0L
    
    # start pooling...
    n_try <- max_try
    repeat {

        # select samples _without_ replacement 
        # (safer than 'sample'-function: no undesired behavior when n = 1)
        index <- sample_index[
            sample.int(
                n = n_samples, 
                size = min(n_samples, max_samples), 
                replace = FALSE
            )
        ]

        # selected smallest set that is in the specified range of target areas
        pooled_area <- cumsum(area[index])
        in_range <- pooled_area %>% 
            between(target_area[1], target_area[2])
        if (any(in_range)) {

            # reset n_try
            n_try <- max_try

            pool_size <- min(which(in_range))
            index <- index[1:pool_size]
            
            # construct new pool identifier
            n_pools <- n_pools + 1L
            pool_id[index] <- n_pools
            
            # update sample indices accordingly
            sample_index <- setdiff(sample_index, index)
            
            # update remaining number of samples to pool
            n_samples <- n_samples - pool_size
            
        } else {
            # termination criterion:
            # - remaining number of samples is smaller than maximum samples in a pool
            # - remaining sample area is smaller than minimum target area.
            if (n_samples <= max_samples) {
                if (pooled_area[length(pooled_area)] < target_area[1]) {
                    break
                }
            }
            
            # termination criterion:
            # decreases and test n_try (trial-and-error)
            n_try <- n_try - 1L
            if (n_try == 0L) {
                # give up
                break
            }
        }
        
        # termination criterion:
        if (n_samples == 0L) {
            break
        }
    }
    
    # return pool identifiers
    pool_id
}



#' Ecological Quality Ratio (EQR)
#'
#' The ecological quality ratio is the ratio beween a parameter value and its
#' reference value:  
#' \deqn{EQR = \frac{x-bad}{ref-bad}}{EQR = (x-bad)/(ref-bad)}
#' Depending on \code{bad} and \code{ref}, the EQR usually 
#' (but not necessarily!) varies between 0 (bad ecological quality) and 1 (
#' ecological quality equals the reference status).
#'
#' @param x numeric vector containing benthic indices
#' @param bad the value for a bad status
#' @param ref the value for a reference status
#'
#' @return  numeric vector with EQR-values: low values indicate bad ecological
#'      quality and high values indicate good ecological quality.
#'  
#' @export
eqr <- function(x, bad, ref) {
    stopifnot((length(bad)  == 1L) | length(bad) == length(x))
    stopifnot((length(ref) == 1L) | length(ref) == length(x))
    if (all(ref > bad) | all(ref < bad)) {
        return((x - bad) / (ref - bad))
    }
    stop(
        paste0(
            "reference values for 'ref' and 'bad' status are inconsistent.\n",
            "'ref' should always be either smaller or greater than 'bad'."
        ), 
        call. = FALSE
    )
}



#' Harmonize Case
#'
#' Convert text to the most occuring case. In case of ties, the first
#' occurence in sorted order will be taken.
#'
#' @param x character vector
#'
#' @return character vector with harmonized names (i.e., same case)
#'  
#' @examples 
#'  x <- c("FOO", "Foo", "bar", "FOO", "bar", "FOO", "Bar")
#'  y <- harmonize(x)
#'  stopifnot(all.equal(y, c("FOO", "FOO", "bar", "FOO", "bar", "FOO", "bar")))
#'  
#' @export
harmonize <- function(x) {
    if (is.character(x)) {
        lx <- tolower(x)
        lut <- tapply(
            X = x, 
            INDEX = lx, 
            FUN = function(x) {
                f <- tapply(X = x, INDEX = x, FUN = length)
                names(f)[which.max(f)]
            }
        )
        x <- as.character(lut[match(x = lx, table = names(lut))])
    }
    x
}



#' Binomial Names
#'  
#' \code{is_binomial} tests for valid binomial names,
#' \code{generic_name} extracts the genus to which the species belongs,
#' \code{specific_name} extracts the species within the genus.
#'
#' @param x \code{\link{character}} vector, containing the binomial name(s) of 
#'      species (a.k.a. binomen or scientific name)
#'
#' @return character vector with either the generic name or the specific name
#'      of the species.
#'  
#' @examples 
#'  is_binomen("Venerupis corrugata") # TRUE
#'  generic_name("Venerupis corrugata") # Venerupis
#'  specific_name("Venerupis corrugata") # corrugata
#'  generic_name("venerupis corrugata") # NA (genus part should be capitalized)
#'
#' @export
is_binomen <- function(x) {
    grepl(pattern = "^([A-Z][a-z]+) +([a-z]+|sp{1,2}\\.?)$", x = x)
}

#' @describeIn is_binomen extracts the genus to which the species belongs
#' @export
generic_name <- function(x) {
    m <- regexpr(pattern = "^([A-Z][a-z]+)", text = x)
    name <- rep.int(x = NA_integer_, times = length(x))
    name[m != -1L] <- regmatches(x, m)
    name[!is_binomen(x)] <- NA_character_
    name
}

#' @describeIn is_binomen extracts the species within the genus
#' @export
specific_name <- function(x) {
    x <- strip_sp(x)
    m <- regexpr(pattern = " +([a-z]+)$", text = x)
    name <- rep.int(x = NA_integer_, times = length(x))
    name[m != -1L] <- regmatches(x, m)
    name[!is_binomen(x)] <- NA_character_
    name <- sub(pattern = "^ +", replacement = "",  x = name)
    name
}


#' @describeIn is_binomen strips postfix sp. or spp. from a binomen
#' @export
strip_sp <-
function(x) {
    ifelse(
        is_binomen(x), 
        sub(pattern = " sp{1,2}\\.?$", replacement = "", x = x),
        x
    )
}


#' Genus to Species Conversion
#'
#' This algorithm reallocates the counts of taxa, that are only identified at
#' the genus level to taxa in the same sampling unit and of the same genus 
#' but that are identified on the species level. The redistribution of counts 
#' is proportional to the number of counts at the species level.
#'
#' @param is_genus \code{logical} vector with elements \code{TRUE} if the 
#'      corresponding taxon is on the genus level, and \code{FALSE} if it is on 
#'      the  \code{species} level.
#' @param count \code{numeric} vector with elements giving the counts of each
#'      corresponding taxon.
#'
#' @note Parameters \code{is_genus} and \code{count} are of the same length and
#'      correspond to the same taxon.
#' @note The resulting counts are not necessarily integers.
#'
#' @return \code{numeric} vector with updated counts. The counts for the 
#'      taxon on the genus level have been set to zero.
#'  
#' @export
#'      
#' @examples
#'      genus_to_species(is_genus = c(TRUE, FALSE, FALSE), count = c(3, 10, 20))
#'      genus_to_species(is_genus = c(TRUE, FALSE, FALSE), count = c(1, 10, 20))
genus_to_species <- function(is_genus, count) {
    stopifnot(length(is_genus) == length(count))
    specific <- which(!is_genus)
    is_genus <- which(is_genus)
    count <- as.numeric(count)
    n_g <- count[is_genus]
    n_s <- count[specific]
    if (length(n_g) == 0L) {
        return(count)
    }
    if (length(n_s) == 0L) {
        return(count)
    }
    if (length(n_g) != 1L) {
        stop(
            "Only one taxon at the genus-level expected.\n", 
            "Use packages like 'dplyr' to group your data.",
             call. = FALSE
        )
    }
    count[specific] <- n_s + n_g * n_s / sum(n_s)
    count[is_genus] <- 0L
    count
}



#' Convert Taxon Names to Comply with WoRMS
#'  
#' @description  
#' Taxon names are standardized according to the World Register of 
#' Marine Species (WoRMS) database. The conversion is case-insensitive.
#' For this conversion, the TWN-list (Taxa Water management the Netherlands) 
#' is used, extended with species of the Southern North Sea. 
#' See references below for download locations.
#'  
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'     \code{data.table}, database etc.
#' @param taxon \code{\link{character}} vector, containing taxon names
#' @param worms an optional table usually created with \code{\link{read_twn}}.
#'  
#' @references \url{https://www.marinespecies.org/}
#' @references \url{https://taxainfo.nl/}
#'
#' @return character vector with WoRMS compliant species names
#'
#' @export
to_worms <-
function(taxon, worms = NULL) {
    .Deprecated("as_accepted")
}

#' @describeIn to_worms check if a taxon complies with WoRMS
#' @return \code{TRUE} for WoRMS compliant species names, 
#'      \code{FALSE} otherwise.
#' @export
is_worms <- 
function(.data = NULL, taxon) {
    .Deprecated("is_accepted")
}

#' @describeIn to_worms as \code{is_worms} but suitable for calling from a 
#'      function (see package \pkg{lazyeval}).  
#' @export
is_worms_ <- 
function(.data, taxon) {
    .Deprecated("is_accepted_")
}





#' Convert Taxon Names to Comply with WoRMS/TWN
#'
#' @description
#' Taxon names are standardized according to the World Register of 
#' Marine Species (WoRMS) database. The conversion is case-insensitive.
#' For this conversion, the TWN-list (Taxa Water management the Netherlands) 
#' is used, extended with species of the Southern North Sea. 
#' See references below for download locations.
#'  
#' @param taxon \code{\link{character}} vector, containing taxon names
#' @param taxa an optional table usually created with \code{\link{read_taxa}}.
#'  
#' @references \url{https://www.marinespecies.org/}
#' @references \url{https://taxainfo.nl/}
#'
#' @return character vector with WoRMS/TWN compliant species names
#'
#' @export
as_accepted <-
function(taxon, taxa = NULL) {
    taxon <- taxon %>% 
        strip_spaces %>% 
        strip_sp
    
    if (is.null(taxa)) {
        taxa <- get_taxa()
    } else {
        if (!inherits(taxa, "tbl")) {
            stop(
                "argument 'taxa' should inherit from class 'tbl'", 
                call. = FALSE
            )
        }
        if (!all(c("provided", "accepted") %in% names(taxa))) {
            stop(
                "argument 'taxa' has not the right format. See ?read_taxa().", 
                call. = FALSE
            )
        }
    }
    taxa$accepted[match(x = tolower(taxon), table = tolower(taxa$provided))]
}


#' @describeIn to_worms check if a taxon complies with WoRMS/TWN
#' @return \code{TRUE} for WoRMS/TWN compliant species names, 
#'      \code{FALSE} otherwise.
#' @export
is_accepted <- 
function(.data = NULL, taxon) {
    is_accepted_(.data, lazy(taxon))
}

#' @describeIn to_worms as \code{is_accepted} but suitable for calling from a 
#'      function (see package \pkg{lazyeval}).  
#' @export
is_accepted_ <- 
function(.data, taxon) {
    taxon <- lazy_eval(taxon, .data)
    d <- get_taxa()
    taxon %in% d$accepted
}
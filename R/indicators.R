#' Abundance
#'
#' The number of indiviuals in each taxon.
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon name of column in \code{.data} containing taxa
#' @param count name of column in \code{.data} containing counts
#'
#' @return \code{numeric} vector with abundance per taxon.
#'  
#' @note due to pooling, the abundance is not necessarily an integer
#'  
#' @examples 
#'  abundance(
#'      taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'      count = c(4, 6)
#'  )
#' @import dplyr
#' @import lazyeval  
#' 
#' @export
abundance <-
function(.data = NULL, taxon = NULL, count) {
    abundance_(.data, lazy(taxon), lazy(count))
}

#' @describeIn abundance version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
abundance_ <- 
function(.data = NULL, taxon = NULL, count) {

    # lazy evaluation of arguments
    taxon <- lazy_eval(taxon, .data)
    count <- lazy_eval(count, .data)
    
    # if 'taxon' is missing, it is assumed that each element in 'count'
    # represents an individual. 
    # 'names' are not meaningful and will be removed.
    if (is.null(taxon)) {
        return(unname(count))
    }
    
    # if '.data.' is missing, 'taxon' and 'count' should be vectors of exactly 
    # the same length
    if (is.null(.data)) {
        if (length(taxon) != length(count)) {
            stop(
                "'taxon' and 'count' should have the same length", 
                call. = FALSE
            )
        }
    }
    
    # check counts
    if (any((count < 0) | (is.na(count)))) {
        stop(
            "all elements in 'count' should be nonnegative",
            call. = FALSE
        )
    }

    # compute abundance (as a named numeric vector)
    n <- tapply(X = count, INDEX = taxon, FUN = sum)
    names_n <- names(n)
    n <- as.numeric(n)
    names(n) <- names_n
    n
}



#' Total Abundance
#'
#' The total number of individuals.
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc. 
#' @param count counts (\code{numeric})
#' @param na.rm Should missing values (including \code{NaN}) 
#'      be removed? (\code{logical})
#'      
#' @return total number of individuals (\code{integer})
#' @examples 
#'  total_abundance(count = c(4, 6))
#'      
#' @export
total_abundance <-
function(.data = NULL, count, na.rm = FALSE) {
    total_abundance_(.data, lazy(count), na.rm = na.rm)
}

#' @describeIn total_abundance version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
total_abundance_ <- 
function(.data = NULL, count, na.rm = FALSE) {
    count <- lazy_eval(count, .data)
    if (any((count < 0) | (is.na(count)))) {
        stop(
            "all elements in 'count' should be nonnegative",
            call. = FALSE
        )
    }
    sum(count, na.rm = na.rm)
}




#' @describeIn total_abundance natural log of total abundance + 1
#'      (see package \pkg{lazyeval}).
#' @export
lnn <- function(.data = NULL, count, na.rm = FALSE) {
    lnn_(.data, lazy(count), na.rm = na.rm)
}

#' @describeIn total_abundance version of lnn suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
lnn_ <- function(.data = NULL, count, na.rm = FALSE) {
    count <- lazy_eval(count, .data)
    if (any((count < 0) | (is.na(count)))) {
        stop(
            "all elements in 'count' should be nonnegative",
            call. = FALSE
        )
    }
    log(sum(count, na.rm = na.rm) + 1)
}

#' Species Richness 
#'
#' Species richness (\eqn{S}{S}) is defined as the number of taxa 
#' (lowest identification level possible) per sampling unit 
#' (data pool or box core sample).
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon taxa names (\code{character})
#' @param count number of individuals for each taxon (\code{numeric})
#'  
#' @return species richness (\code{integer} vector of length 1)
#'  
#' @examples 
#'  species_richness(
#'      taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'      count = c(4, 6)
#'  )
#'  
#' @export
species_richness <-
function(.data = NULL, taxon, count = NULL) {
    species_richness_(.data, lazy(taxon), lazy(count))
}

#' @describeIn species_richness version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
species_richness_ <- 
function(.data = NULL, taxon, count = NULL) {

    # evaluate function arguments
    taxon <- lazy_eval(taxon, .data)
    count <- lazy_eval(count, .data)
    
    if (!is.null(count)) {
        # if '.data.' is missing, 'taxon' and 'count' should be vectors
        # of exactly the same length
        if (is.null(.data)) {
            if (length(taxon) != length(count)) {
                stop(
                    "'taxon' and 'count' should have the same length", 
                    call. = FALSE
                )
            }
        }

        # check counts
        if (any((count < 0) | (is.na(count)))) {
            stop(
                "all elements in 'count' should be nonnegative",
                call. = FALSE
            )
        }
    
        # remove records without individuals
        index <- which(count == 0)
        if (length(index) > 0L) {
            taxon <- taxon[-index]
        }
    }
    
    # return species richness
    length(unique(taxon))
}



#' Margalef Index of Diversity
#'
#' Margalef Index of Diversity is given by
#' \deqn{D = \frac{S-1}{\ln(N)}}{D = (S-1)/ln(N)}
#'  
#' For \eqn{N=1}, the index is set to 0.
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon taxa names (\code{character})
#' @param count counts (\code{numeric})
#'
#' @return Margalef diversity index (\code{numeric} \code{vector} of length 1)
#'  
#' @examples 
#' margalef(
#'     taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'     count = c(4, 6)
#' )
#'  
#' @export
margalef <-
function(.data = NULL, taxon, count) {
    margalef_(.data, lazy(taxon), lazy(count))
}

#' @describeIn margalef version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
margalef_ <- 
function(.data = NULL, taxon, count) {

    # species richness
    S <- species_richness_(.data, taxon = taxon, count = count)

    # abundance
    N <- total_abundance_(.data, count = count)

    # Margalef's index of diversity
    if (N == 1L) {
        return(0)
    }
    (S-1)/log(N)
}



#' Rygg's Index of Diversity
#'
#' Rygg's index of diversity is given by
#' \deqn{SN = \frac{\ln(S)}{\ln(\ln(N))}}{SN = ln(S)/(ln(ln(N)))}
#'  
#' The adjusted version of Rygg's index which gives more consistent values
#' for smaller \code{S=2, N=2, N=3} and \code{S=3, N=3} is
#' \deqn{SN = \frac{\ln(S)}{\ln(\ln(N+1)+1)}}{SN = ln(S)/(ln(ln(N+1)+1))}
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon taxa names (\code{character})
#' @param count counts (\code{numeric})
#' @param adjusted (defaults to \code{FALSE})
#'
#' @return Rygg's index of diversity (\code{numeric} \code{vector} of length 1)
#'  
#' @note Rygg's index is not defined for \eqn{N=exp(1)}. For 
#'      \eqn{N \leq exp(1)}{N <= exp(1)}, \code{rygg} returns 
#'      \code{\link{NA_real_}}.
#'  
#' @references Rygg, B. (2006). Developing indices for quality-status 
#'      classification of marine soft-bottom fauna in Norway. 
#'      Norwegian Institute for Water Research, Oslo, Norway. 
#'      NIVA Report SNO 5208-2006.
#'  
#' @examples 
#'  rygg(
#'      taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'      count = c(4, 6)
#'  )
#'  
#' @export
rygg <-
function(.data = NULL, taxon, count, adjusted = FALSE) {
    rygg_(.data, lazy(taxon), lazy(count), adjusted = adjusted)
}

#' @describeIn rygg version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
rygg_ <- 
function(.data = NULL, taxon, count, adjusted = FALSE) {

    # species richness
    S <- species_richness_(.data, taxon = taxon, count = count)

    # abundance
    N <- total_abundance_(.data, count = count)
    
    # return Brage Rygg's index of diversity SN
    if (adjusted) {
        result <- log(S) / log1p(log1p(N))
    } else {
        # note that Rygg's index is not defined for N=e
        result <- NA_real_
        if (N > exp(1)) {
            result <- log(S)/log(log(N))
        }
    }
    result
}



#' Shannon's Index or Entropy
#'
#' Compute entropy according to Shannon (1948)
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon taxa names (\code{character})
#' @param count counts (\code{numeric})
#' @param base the base with respect to which logarithms are computed. 
#'      Defaults to 2 (unit: bits).
#'
#' @return Shannon's entropy
#'  
#' @references Shannon, C. E., 1948. A Mathematical Theory of Communication.
#'      Bell System Technical Journal 27: 379-423.
#'
#' @examples
#'  shannon(
#'      taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'      count = c(4, 6)
#'  )
#'
#' @export
shannon <- 
function(.data = NULL, taxon, count, base = 2) {
    shannon_(.data, lazy(taxon), lazy(count), base = base)
}

#' @describeIn shannon version suitable for calling from a function
#'     (see package \pkg{lazyeval}).
#' @export
shannon_ <- 
function(.data = NULL, taxon, count, base = 2) {
    if (length(base) != 1L) {
        stop(
            "'base' should be a numeric vector of length 1", 
            call. = FALSE
        )
    }
    n <- abundance_(.data, taxon = taxon, count = count)
    N <- total_abundance_(.data, count = count)
    epsilon <- .Machine$double.eps
    if (N < epsilon) {
        return(NaN)
    }
    p <- n / N
    p <- p[p > epsilon] # remove limiting case as lim p->0 (p log p = 0)
    -sum(p * log(p, base = base))
}


#' AZTI Marine Biotic Index (AMBI)
#'
#' AZTI Marine Biotic Index (AMBI) according to Borja et al. (2000)
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon species names
#' @param count counts of individuals (\code{numeric})
#' @param group sensitivity groups I, II, III, IV, or V
#'
#' @details The index is given by:
#'  \deqn{c_\mathrm{b} = \frac{3}{2} \sum_{i=2}^5 (i-1) p_i}
#'  where \eqn{p_i} is the proportion of species in sensitivity group \eqn{i}.
#'
#' @return numeric vector of length 1 containing the AMBI
#'  
#' @references Borja, A., J. Franco and V. Perez, 2000. A Marine Biotic Index 
#'  to Establish the Ecological Quality of Soft-Bottom Benthos Within 
#'  European Estuarine and Coastal Environments. 
#'  Marine Pollution Bulletin 40:1100-1114
#'  
#' @examples 
#'  ambi(
#'      taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'      count = c(4, 6)
#'  )
#'  
#' @export
ambi <- 
function(.data = NULL, taxon, count, group = NULL) {
    ambi_(.data, lazy(taxon), lazy(count), lazy(group))
}



#' @describeIn ambi version suitable for calling from a function
#'     (see package \pkg{lazyeval}).
#' @export
ambi_ <- 
function(.data = NULL, taxon, count, group = NULL) {

    # evaluate function arguments
    taxon <- lazy_eval(taxon, .data)
    count <- lazy_eval(count, .data)
    group <- lazy_eval(group, .data)
    
    # check lengths
    if (is.null(.data)) {
        if (length(taxon) != length(count)) {
            stop(
                "'taxon' and 'count' should have the same length", 
                call. = FALSE
            )
        }
        if (!is.null(group)) {
            if (length(taxon) != length(group)) {
                stop(
                    "'taxon' and 'group' should have the same length", 
                    call. = FALSE
                )
            }
        }
    }

    # check counts
    if (any((count < 0) | (is.na(count)))) {
        stop(
            "all elements in 'count' should be nonnegative",
            call. = FALSE
        )
    }
    
    # set up tibble
    d <- tibble(TAXON = taxon, COUNT = count)
    if (is.null(group)) {
        # default
        d_ambi <- .get_ambi()
        d <- d %>% left_join(d_ambi, by = "TAXON")
    } else {
        # specified group, but use default if specified group is missing
        d$GROUP <- group
        if (any(is.na(d$GROUP))) {
            d_ambi <- .get_ambi()
            d <- d %>% 
                left_join(d_ambi, by = "TAXON") %>% 
                mutate_(GROUP = ~ ifelse(is.na(GROUP.x), GROUP.y, GROUP.x))
        }
    }
    d <- d %>% 
        select_(~GROUP, ~COUNT) %>% 
        filter_(~!is.na(GROUP))
    
    # check sensitivity groups
    permissable_groups <- c("I", "II", "III", "IV", "V")
    if (!all(d$GROUP %in% permissable_groups)) {
        stop(sprintf("'group' should be one of: %s",
            paste(permissable_groups, collapse = ", ")),
            call. = FALSE
        )
    }
    
    # handle case when all COUNTs are zero
    if (sum(d$COUNT) < 1.0e-9) {
        return(NA_real_)
    }
    
    # compute AMBI
    d %>%
        group_by_(~GROUP) %>%
        summarise_(n = ~sum(COUNT)) %>%
        right_join(tibble(GROUP = permissable_groups), by = "GROUP") %>%
        arrange_(~GROUP) %>%
        mutate_(n = ~ifelse(is.na(n), 0, n), p = ~n / sum(n)) %>% 
        summarise_(AMBI = ~1.5 * sum(0:4 * p)) %>% 
        unlist(use.names = FALSE)
}
    
#' @describeIn ambi tests if an AMBI sensitivity group is available for \code{taxon}
#'      (returns \code{TRUE} (available) or \code{FALSE} (unavailable))
#' @export
has_ambi <-
function(.data = NULL, taxon, group = NULL) {
    has_ambi_(.data, lazy(taxon), lazy(group))
}

#' @describeIn ambi version suitable for calling from a function
#'  (see package \pkg{lazyeval}).
#'  
#' @examples
#'      data(oosterschelde)
#'      has_ambi(oosterschelde, TAXON)
#'  
#' @export
has_ambi_ <-
function(.data = NULL, taxon, group = NULL) {
    taxon <- lazy_eval(taxon, .data)
    group <- lazy_eval(group, .data)
    if (is.null(group)) {
        group <- rep.int(x = NA_character_, times = length(taxon))
    } else {
        if (length(taxon) != length(group)) {
            stop(
                "'taxon' and 'group' should have the same length", 
                call. = FALSE
            )
        }
    }
    (taxon %in% .get_ambi()$TAXON) | !is.na(group)
}



#' Infaunal Trophic Index (ITI)
#'  
#' @description 
#' Computes the Infaunal Trophic Index (ITI) according to 
#' Gittenberger & van Loon (2013).
#'  
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon species names
#' @param count counts of individuals (\code{numeric})
#' @param group sensitivity groups I, II, III, or IV
#'  
#' @details The Infaunal Trophic Index (ITI) is given by
#'  \deqn{\mathrm{ITI} = 100 \sum_{i=1}^3 \frac{(4-i)}{3} p_i}
#'  where \eqn{p_i} is the proportion of species in class \eqn{i}, where
#'  \itemize{
#'      \item group   I are suspension feeders (highest quality);
#'      \item group  II are interface feeders
#'      \item group III are surface deposit feeders and
#'      \item group  IV are subsurface deposit feeders (lowest quality). 
#'  }
#'
#' @return numeric vector of length 1 containing the ITI
#'  
#' @references Gittenberger A. and  W. van Loon, 2013. 
#'      Sensitivities of marine macrozoobenthos to environmental pressures 
#'      in the Netherlands. Nederlandse Faunistische 
#'      Mededelingen 41: 79-112.
#'      
#' @examples 
#'      iti(taxon = c("Euspira pulchella", "Nephtys cirrosa"), count = c(4, 6))
#'      
#'  
#' @export
iti <- 
function(.data = NULL, taxon, count, group = NULL) {
    iti_(.data, lazy(taxon), lazy(count), lazy(group))
}

#' @describeIn iti version suitable for calling from a function
#'  (see package \pkg{lazyeval}).
#' @export
iti_ <- 
function(.data = NULL, taxon, count, group = NULL) {

    # evaluate function arguments
    taxon <- lazy_eval(taxon, .data)
    count <- lazy_eval(count, .data)
    group <- lazy_eval(group, .data)
    
    # check lengths
    if (is.null(.data)) {
        if (length(taxon) != length(count)) {
            stop(
                "'taxon' and 'count' should have the same length", 
                call. = FALSE
            )
        }
        if (!is.null(group)) {
            if (length(taxon) != length(group)) {
                stop(
                    "'taxon' and 'group' should have the same length", 
                    call. = FALSE
                )
            }
        }
    }

    # check counts
    if (any((count < 0) | (is.na(count)))) {
        stop(
            "all elements in 'count' should be nonnegative",
            call. = FALSE
        )
    }

    # set up tibble
    d <- tibble(TAXON = taxon, COUNT = count)
    if (is.null(group)) {
        # default ITI
        d_iti <- get_iti()
        d <- d %>% left_join(d_iti, by = "TAXON")
    } else {
        # specified group, but default ITI if specified group is missing
        d$GROUP <- group
        if (any(is.na(d$GROUP))) {
            d_iti <- get_iti()
            d <- d %>% left_join(d_iti, by = "TAXON") %>% 
                mutate_(GROUP = ~ ifelse(is.na(GROUP.x), GROUP.y, GROUP.x))
        }
    }
    d <- d %>% 
        select_(~GROUP, ~COUNT) %>% 
        filter_(~!is.na(GROUP))
    
    # check sensitivity groups
    permissable_groups <- c("I", "II", "III", "IV")
    if (!all(d$GROUP %in% c(permissable_groups, NA_character_))) {
        stop(sprintf("'group' should be one of: %s",
            paste(permissable_groups, collapse = ", ")),
            call. = FALSE
        )
    }

    # handle case when all COUNTs are zero (or nrow(d)==0)
    if (sum(d$COUNT) < 1.0e-9) {
        return(NA_real_)
    }

    # compute ITI
    d %>%
        group_by_(~GROUP) %>%
        summarise_(n = ~sum(COUNT)) %>%
        right_join(tibble(GROUP = permissable_groups), by = "GROUP") %>%
        arrange_(~GROUP) %>%
        mutate_(n = ~ifelse(is.na(n), 0, n), p = ~n / sum(n)) %>% 
        summarise_(ITI = ~100 * sum(3:0 * p) / 3) %>% 
        unlist(use.names = FALSE)
}


#' @describeIn iti tests if an ITI sensitivity group is available for \code{taxon}
#'      (returns \code{TRUE} (available) or \code{FALSE} (unavailable))
#' @export
has_iti <-
function(.data = NULL, taxon, group = NULL) {
    has_iti_(.data, lazy(taxon), lazy(group))
}

#' @describeIn iti version suitable for calling from a function
#'  (see package \pkg{lazyeval}).
#'  
#' @examples
#'    data(oosterschelde)
#'    has_iti(oosterschelde, TAXON)
#'  
#' @export
has_iti_ <-
function(.data = NULL, taxon, group = NULL) {
    taxon <- lazy_eval(taxon, .data)
    group <- lazy_eval(group, .data)
    if (is.null(group)) {
        group <- rep.int(x = NA_character_, times = length(taxon))
    } else {
        if (length(taxon) != length(group)) {
            stop(
                "'taxon' and 'group' should have the same length", 
                call. = FALSE
            )
        }
    }
    (taxon %in% get_iti()$TAXON) | !is.na(group)
}




#'  Hurlbert's Expected Number of Species
#'
#'  The expected number of species in a sample of \code{n} individuals:
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon name of column in \code{.data} containing taxa
#' @param count name of column in \code{.data} containing counts
#' @param n number of individuals in a standard sample
#'
#' @return expected number of species in a sample of \code{n} individuals
#'  
#' @references Hurlbert, S.H., 1971. The Nonconcept of Species Diversity: 
#'      A Critique and Alternative Parameters. Ecology 52:577-586.
#'      
#' @examples 
#'      hurlbert(
#'          taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'          count = c(4, 6),
#'          n = 8
#'      )
#'      
#' @export
hurlbert <-
function(.data = NULL, taxon, count, n = 100L) {
    hurlbert_(.data, lazy(taxon), lazy(count), n = n)
}

#' @describeIn hurlbert version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#'      
#' @export
hurlbert_ <- 
function(.data = NULL, taxon, count, n = 100L) {
    if (length(n) != 1L) {
        stop(
            "'n' should be an integer vector of length 1", 
            call. = FALSE
        )
    }
    p <- abundance_(.data, taxon = taxon, count = count)
    N <- total_abundance_(.data, count = count)
    if (n > N) {
        stop(
            "'n' (=", n, ") should be less than the total abundance (=", N, ")",
            call. = FALSE
        )
    }
    
    # Let 
    #   a = ln(choose(n = N - p, k = n)
    #   b = ln(choose(n = N,     k = n))
    # and given that
    #   exp(a) / exp(b) = exp(a-b)
    # it follows that
    #   sum(1 - choose(n = N - p, k = n) / choose(n = N, k = n))
    # can be written as (less prone to overflow):
    a <- lchoose(n = N - p, k = n)
    b <- lchoose(n = N,     k = n)
    sum(1 - exp(a - b))
}



#' Hurlbert's Probability of Interspecific Encounter (PIE)
#'
#' The probability that two individuals selected at random (\emph{without} 
#' replacement) from a sample will belong to different species is given by 
#' (Hurlbert, 1971, p.579, Eq. 3):
#' \deqn{\Delta_1 = \sum_{i=1}^S (\frac{N_i}{N})(\frac{N-N_i}{N-1}) =
#' (\frac{N}{N-1})\Delta_2}{\Delta1 = \sum_{i=1}^S (N[i]/N)((N-N[i])(N-1)) =
#' (N/(N-1))\Delta2}
#' where \eqn{\Delta_2}{\Delta2} (Hurlbert, 1971, p.579, Eq. 4) is the 
#' probability that two individuals selected at random (\emph{with} 
#' replacement) from a sample will belong to different species:
#' \deqn{\Delta_2 = 1 - \sum_{i=1}^S \pi_i^2}{\Delta2 = 1 - \sum_{i=1}^S \pi[i]^2}
#' where  \eqn{N_i}{N[i]} is the number of individuals of the \emph{i}th 
#' species in the community, \eqn{N} is the total number of individuals in the
#' community, \eqn{\pi_i = N_i/N}{\pi[i] = N[i]/N}, and \eqn{S} is the number of 
#' species in the community.
#' Note that Hurlbert's PIE \code{hpie} is the complement of 
#' \code{\link{simpson}}.
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon name of column in \code{.data} containing taxa
#' @param count name of column in \code{.data} containing counts
#'
#' @return A numeric vector with the probability of interspecific encounter
#'      (PIE).
#'  
#' @references Hurlbert, S.H., 1971. The Nonconcept of Species Diversity: 
#'      A Critique and Alternative Parameters. Ecology 52:577-586.
#'      
#' @seealso \code{\link{simpson}}, \code{\link{hurlbert}}
#'  
#' @examples 
#'      hpie(
#'          taxon = c("Euspira pulchella", "Nephtys cirrosa"), 
#'          count = c(6, 12)
#'      )
#'  
#' @export
hpie <-
function(.data = NULL, taxon, count) {
    hpie_(.data, lazy(taxon), lazy(count))
}

#' @describeIn hpie suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
hpie_<- 
function(.data = NULL, taxon, count) {
    n <- abundance_(.data, taxon = taxon, count = count)
    N <- total_abundance_(.data, count = count)
    # mathematically identical to the expression in Hurlbert (1971), 
    # but less prone to round off error due (less divisions) (= 1 - simpson)
    1 - sum(n * (n-1L)) / (N * (N-1L)) # Hurlbert's Delta1
}


#' Simpson's Measure of Concentration
#'
#' The probability that two inidividuals selected at random (with replacement, 
#' Hurlbert, 1971, p.579) from a sample will belong to the same species. For
#' an infinite sample Simpson's Index is given by (Peet, 1974):
#' \deqn{\lambda = \sum_{i=1}^S p_i^2}
#' For a finite sample by:
#' \deqn{L = \sum_{i=1}^S \frac{n_i (n_i-1)}{N (N-1)}}
#' where \eqn{p_i}{p[i]} the proportion of the individuals in species \eqn{i}, 
#' \eqn{n_i}{n[i]} the number of individuals in species 
#' \eqn{i} (relative \code{\link{abundance}}), and \eqn{N} the total number 
#' of individuals (\code{\link{total_abundance}}). The finite sample case 
#' has been implemented in function \code{simpson} (and \code{simpson_}).
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon name of column in \code{.data} containing taxa
#' @param count name of column in \code{.data} containing counts
#'
#' @return The probability that two inidividuals selected at random from a 
#'      sample will belong to the same species.
#'  
#' @references Peet, R. K. 1974, The Measurement of Species Diversity. Annual 
#'      Review of Ecology and Systematics 5:285-307.
#'      
#' @seealso \code{\link{hpie}}
#'  
#' @examples 
#'      simpson(
#'          taxon = c("Euspira pulchella", "Nephtys cirrosa"),  
#'          count = c(6, 12)
#'      )
#'  
#' @export
simpson <-
function(.data = NULL, taxon, count) {
    simpson_(.data, lazy(taxon), lazy(count))
}

#' @describeIn simpson version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
simpson_ <- 
function(.data = NULL, taxon, count) {
    n <- abundance_(.data, taxon = taxon, count = count)
    N <- total_abundance_(.data, count = count)
    if (N == 1L) {
        # one individual sampled with replacement gives the same individual
        # hence: lambda = 1
        return(1.0)
    }
    sum((n * (n - 1L)) / (N * (N - 1L)))
}




#' Hill's Diversity Numbers
#'
#' According to Hill (1973): \emph{"a diversity number is figuratively a 
#' measure of how many species are present if we examine the sample down to a 
#' certain depth among its rarities. If we examine superficially (e.g.,
#' by using \eqn{N_2}{N[2]}) we shall see only the more abundant species. If we 
#' look deeply (e.g., by using \eqn{N_0}{N[0]}) we shall see all the 
#' species present."}
#'  
#' Hill's diversity numbers are given by:
#' \deqn{N_a=\sum{i=1}^S (p_i^a)^{1/(1-a)}}
#'  
#' Special cases are:
#' \describe{
#'      \item{\eqn{N_{-\infty}}{N[-Inf]}}{reciprocal of the proportional
#'          abundance of the rarest species;}
#'      \item{\eqn{N_0}{N[0]}}{total number of species present;}
#'      \item{\eqn{N_1}{N[1]}}{exp(H), where H: Shannon's index (see also 
#'          \code{\link{shannon}});}
#'      \item{\eqn{N_2}{N[2]}}{reciprocal of Simpson's index (see also 
#'          \code{\link{simpson}});}
#'      \item{\eqn{N_{\infty}}{N[Inf]}}{reciprocal of the proportional
#'          abundance of the commonest species.}
#' }
#'
#' @param .data data in a \code{data.frame}, \code{tibble}, 
#'      \code{data.table}, database etc.
#' @param taxon name of column in \code{.data} containing taxa
#' @param count name of column in \code{.data} containing counts
#' @param a exponent in Hill's diversity number (R, with special cases for 
#'      \code{a} in {0, 1, 2} (see details))
#'
#' @return numeric vector of Hill's numbers
#'  
#' @references Hill, M.O., 1973. Diversity and Evenness: 
#'      A Unifying Notation and Its Consequences. Ecology 54:427-432
#'  
#' @seealso \code{\link{species_richness}}, \code{\link{shannon}}, 
#'      \code{\link{simpson}}
#'      
#' @examples 
#'      hill(
#'          taxon = c("Euspira pulchella", "Nephtys cirrosa"),  
#'          count = c(6, 12),
#'          a = 0
#'      )
#'      hill0(
#'          taxon = c("Euspira pulchella", "Nephtys cirrosa"),  
#'          count = c(6, 12)
#'      )
#'  
#' @export
hill <-
function(.data = NULL, taxon, count, a = 0) {
    hill_(.data, lazy(taxon), lazy(count), a = a)
}

#' @describeIn hill version suitable for calling from a function
#'      (see package \pkg{lazyeval}).
#' @export
hill_ <- 
function(.data = NULL, taxon, count, a = 0) {
    if (length(a) != 1L) {
        stop(
            "'a' should be a numeric vector of length 1", 
            call. = FALSE
        )
    }
    if (a == 1) { # undefined (limit exists, see hill1)
        message("N_a(a=1) is undefined. Therefore N_a(lim a->1) will be returned")
        return(hill1_(.data, taxon = taxon, count = count))
    }
    n <- abundance_(.data, taxon = taxon, count = count)
    N <- total_abundance_(.data, count = count)
    p <- n/N
    sum(p^a)^(1/(1-a))
}


#' @describeIn hill \eqn{N_0}{N[0]}
#' @export
hill0 <-
function(.data = NULL, taxon, count) {
    hill0_(.data, lazy(taxon), lazy(count))
}

#' @describeIn hill \eqn{N_0}{N[0]}, version suitable for calling 
#'  from a function (see package \pkg{lazyeval}).
#' @export
hill0_ <- 
function(.data = NULL, taxon, count) {
    species_richness_(.data, taxon = taxon, count = count)
}

#' @describeIn hill \eqn{N_1}{N[1]}
#' @export
hill1 <-
function(.data = NULL, taxon, count) {
    hill1_(.data, lazy(taxon), lazy(count))
}

#' @describeIn hill \eqn{N_1}{N[1]}, version suitable for calling 
#' from a function (see package \pkg{lazyeval}).
#' @export
hill1_ <- 
function(.data = NULL, taxon, count) {
    H <- shannon_(.data, taxon = taxon, count = count, base = exp(1))
    exp(H)
}


#' @describeIn hill \eqn{N_2}{N[2]}
#' @export
hill2 <-
function(.data = NULL, taxon, count) {
    hill2_(.data, lazy(taxon), lazy(count))
}

#' @describeIn hill \eqn{N_2}{N[2]}, version suitable for calling 
#'  from a function (see package \pkg{lazyeval}).
#' @export
hill2_ <- 
function(.data = NULL, taxon, count) {
    hill_(.data, taxon = taxon, count = count, a = 2)
}
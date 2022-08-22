#' Read and Validate BEQI2 Input Files
#'
#' This function reads and checks BEQI2 input files. The format has been
#' specified in Van Loon (2013) and is described in the vignette of the BENMMI-package.
#'
#' @param filename name of BEQI2 input file (\code{character})
#'
#' @details The function performs the following tasks:
#' \itemize{
#'      \item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'      \item{make column names with aggregation data case-insensitive;}
#'  	\item{removes redundant spaces;}
#'      \item{checks if DATE-field adheres to ISO 8601 (YYYY-mm-dd);}
#'  	\item{constructs a unique identifier \code{ID} by concatenating 
#'          columns \code{OBJECTID} and \code{DATE};}
#'      \item{checks that each \code{ID} has a unique \code{AREA};}
#'      \item{checks azoic samples for VALUE=0;}
#'      \item{removes records with VALUE=0, not belonging to azoic samples;}
#'      \item{checks VALUE-field on missing values;}
#'      \item{checks if VALUE-field is an integer;}
#'  }
#'
#' @references Willem van Loon, 2013. BEQI2 INPUT FORMAT. See the package-vignette of the BENMMI-package.
#' 	
#' @importFrom utils read.csv
#'
#' @export
read_beqi2 <-
function(filename) {

	# check if BEQI2-file exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read BEQI2-file
	d <- read.csv(file = filename, as.is = TRUE, na.strings = c("NA", ""))

    # validation
    validate_beqi2(d)
}



#' @describeIn read_beqi2 validator for BEQI2-format
#' @param .data table in BEQI2-format
#' @import dplyr
validate_beqi2 <-
function(.data) {

    # convert to tibble
    .data <- .data %>% as_tibble

	# CHAR is not required anymore: patch for backward-compatibility
    if (!("CHAR" %in% names(.data))) {
        .data$CHAR <- NA_character_
	}
    
	# check column names (case insensitive)
	required_vars <- c("OBJECTID", "HABITAT", "SAMPLEID", "TAXON", 
                       "CHAR", "SAMPDEV", "AREA", "DATE", "VALUE")
	missing_vars <- setdiff(required_vars, toupper(names(.data)))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missing_vars)
            ),
            call. = FALSE
        )
	}

    # check required columns on missing values
	for (required_var in required_vars %>% setdiff(c("CHAR", "SAMPDEV"))) {
        is_missing <- is.na(.data[[required_var]])
        if (any(is_missing)) {
            index <- which(is_missing)
            stop(
                sprintf(
                    fmt = paste0(
                        "The value in column %s is missing for %i records.\n",
                        "The record indices are:\n%s"
                    ),
                    sQuote(required_var),
                    length(index),
                    toString(index)
                ), 
                call. = FALSE
            )
        }
	}

    # check VALUE field on non-negative values
    is_negative <- .data$VALUE < 0
    if (any(is_negative)) {
        stop(
            sprintf(
                fmt = paste0(
                    "Negative values found in column %s.\n",
                    "The record indices are: %s"
                ),
                sQuote("VALUE"),
                toString(which(is_negative))
            ), 
            call. = FALSE
        )
    }

    # check VALUE field on non-integer values
    is_integer <- (.data$VALUE - floor(.data$VALUE)) < .Machine$double.eps
    if (any(!is_integer)) {
        warning(
            sprintf(
                fmt = paste0(
                    "%s out of %s non-integer values found in column %s.\n",
                    "This is about %s%% of the records.\n",
                    "Please make sure that these records contain species abundances (and not species densities)."),
                sum(!is_integer),
                length(is_integer),
                sQuote("VALUE"),
                round(100 * sum(!is_integer) / length(is_integer))
            ), 
            call. = FALSE
        )
    }
    
	# remove redundant spaces (alternative: 'mutate_each' function did crash)
    .data[required_vars] <- .data[required_vars] %>% 
        lapply(FUN = strip_spaces) %>% 
        as_tibble

    # harmonize columns (same case)
    .data[required_vars] <- .data[required_vars] %>% 
        lapply("harmonize") %>% 
        as_tibble
    
	# coerce date to class 'Date'
	res <- try(as.Date(.data$DATE, format = "%Y-%m-%d"), silent = TRUE)
	if (inherits(res, "try-error") | any(is.na(res))) {
		stop(
			"Invalid date formats found. Please adhere to ISO 8601 (YYYY-mm-dd)",
			call. = FALSE
		)
	}
	.data$DATE <- res

    # check if areas are unique for a specific ID
    n <- .data %>%
        group_by_(~OBJECTID, ~SAMPLEID, ~DATE) %>%
        summarise_(n = ~length(unique(AREA))) %>%
        ungroup() %>%
        filter_(~n != 1L)
    if (nrow(n) != 0L) {
        stop(
            sprintf(
                fmt = paste0(
                    "AREA is not unique for %i samples. ",
                    "These samples have the following OBJECTID/SAMPLEID/DATE:\n%s"
                ),
                nrow(n),
                n %>% 
                    select_(~OBJECTID, ~SAMPLEID, ~DATE) %>% 
                    apply(MARGIN = 1, paste, collapse = "/") %>%
                    paste(collapse = "; ")
            ), 
            call. = FALSE
        )
    }

    # find duplicated records (all required columns, except for VALUE)
    d <- .data[, required_vars %>% setdiff(c("VALUE", "CHAR", "SAMPDEV"))]
    if (anyDuplicated(d)) {
        stop(
            sprintf(
                "Duplicated rows found in the benthos file.\nThe row numbers are:\n%s",
                toString(
                    sort(
                        which(
                            duplicated(d) | 
                            duplicated(d, fromLast = TRUE)
                        )
                    )
                )
            ),
            call. = FALSE
        )
    }
    
    # aggregate (by summation) the VALUE-fields of records that only differ
    # in VALUE-field value and CHAR-field. These are replicates or juvenile
    # individuals.
    # (NB: the sample areas should not be summed!)
    # .data <- .data %>% 
    #     mutate_(CHAR = ~ifelse(CHAR == "REST", "REST", NA_character_)) %>%
    #     group_by_(.dots = setdiff(names(.data), "VALUE")) %>%
    #     summarise_(VALUE = ~sum(VALUE)) %>%
    #     ungroup

    # handle azoic samples
    index <- which(is_azoic(x = .data$TAXON))
    if (length(index) != 0L) {
        if (isTRUE(any(.data$VALUE[index] != 0L))) {
            warning(
                paste(
                    "Azoic sample(s) found with non-zero VALUE-field.",
                    "These abundances will be set to zero"
                ),
                call. = FALSE
            )
        }
        .data$VALUE[index] <- 0L
    }

    # remove non-Azoic records with zero counts as these are redundant 
    # and won't affect the results
    is_zero <- (.data$VALUE == 0L) & !is_azoic(x = .data$TAXON)
    if (any(is_zero)) {
        warning(
            sprintf(
                paste(
                    "Non-azoic records (n=%s) found with zero VALUE-field.",
                    "These records are redundant and will be excluded."
                ),
                sum(is_zero)
            ),
            call. = FALSE
        )
        .data <- .data[!is_zero, ]
    }

    # final checks
    if (nrow(.data) == 0L) {
        stop("No valid records found", call. = FALSE)
    }
    .data <- .data %>% distinct_

	# return result
	.data
}



#' Read and Validate Taxa Waterbeheer Nederland (TWN) Data
#'
#' This function reads files in the Taxa Waterbeheer Nederland (TWN) format.
#'
#' @details The function adds a new column \code{taxon}. Its contents depending 
#'     on TWN-status:
#' \itemize{
#' 		\item{status = 10} {taxonname}
#'  	\item{status = 20} {prefername}
#'  	\item{status = 80} {parentname}
#' 	}
#'
#' @param filename name of TWN file (\code{character})
#'
#' @return a \code{tibble} with four columns:
#'  \itemize{
#' 		\item{GROUP} {TWN/WoRMS taxon group}
#'  	\item{LEVEL} {TWN/WoRMS taxon level}
#'  	\item{FROM} {taxon name to convert from}
#'      \item{TO} {taxon name to convert to}
#' 	}
#'
#' @references \url{https://taxainfo.nl/}
#'  
#' @importFrom utils read.csv
#' 
#' @export
read_twn <-
function(filename) {

    # raise error if filename is missing
    if (missing(filename)) {
        stop("Filename is missing", call. = FALSE)
    }
    
	# check if file exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file
	d <- try(
        read.csv(file = filename, as.is = TRUE, na.strings = c("NA", "")), 
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

    # validate TWN-file
    validate_twn(d)
}


#' @describeIn read_twn get default WoRMS list (TWN list extended with species Southern North Sea)
#' @export
get_worms <-
function() {
    .Deprecated("get_taxa")
}



#' @describeIn read_twn validator for TWN-format
#' @param .data table in TWN-format
#' @import dplyr
validate_twn <-
function(.data) {
    
    # convert to tibble
    .data <- .data %>% 
        as_tibble

	# check column names
	required_vars <- c("status", "taxonname", "taxongroup", "prefername", 
                       "parentname", "taxonlevel")
	missing_vars <- setdiff(required_vars, tolower(names(.data)))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missing_vars)
            ),
            call. = FALSE
        )
	}
	names(.data) <- tolower(names(.data))

	# select only columns of interest
	.data <- .data %>% 
	    select_(.dots = required_vars)

	# keep only status codes:
	#	10: preferred name 
	#	20: synonym
	#	80: non-taxonomic species group (groups and aggregates)
	# see also www.aquo.nl/faq/faq-twn
    .data <- .data %>% 
        filter_(~status %in% c(10L, 20L, 80L))

	# remove redundant spaces 
	# (including leading and trailing spaces)
    .data <- .data %>% 
        mutate_if(is.character, funs_(~strip_spaces))

	# construct TAXON
    .data <- bind_rows(
        .data %>% filter_(~status == 10) %>% mutate_(TAXON = ~taxonname),
        .data %>% filter_(~status == 20) %>% mutate_(TAXON = ~prefername),
        .data %>% filter_(~status == 80) %>% mutate_(TAXON = ~parentname)
    )
    
    # for status 80 also the taxon level needs to be changed 
    # to that of the parent
    id <- with(.data, 
        which(
                (status == 80) &
                (taxonlevel %in% c("Genus combi", "Subgenus", 
	                               "Species", "Species combi", "Subspecies"))
        )
    )
    .data$taxonlevel[id] <- "Genus"

	# create ordered factor of taxon levels
	.data$taxonlevel <- factor(
		x = .data$taxonlevel, 
		levels = c(
			"Regio", 
			"Regnum", 
			"Phylum", "Subphylum", 
			"Classis", "Subclassis", "Infraclassis", 
			"Ordo", "Subordo", "Infraordo", 
			"Superfamilia",	"Familia", "Subfamilia", 
			"Tribe", 
			"Genus", "Genus combi", "Subgenus", 
			"Species", "Species combi", "Subspecies", 
			"Varietas", 
			"Forma"),
		ordered = TRUE
	)

	# selection
    .data <- .data %>% 
        select_(GROUP = ~taxongroup, LEVEL = ~taxonlevel, 
                FROM = ~taxonname, TO = ~TAXON) %>%
        arrange_(~GROUP, ~FROM)

    # check if TO can be converted to FROM
    is_inconvertible <- is.na(.data$TO) | is.na(.data$FROM)
	if (any(is_inconvertible)) {
		stop(
			sprintf(
				"A total of %d taxon names in the taxon-file are inconvertible:\n%s", 
				sum(is_inconvertible),
				.data$FROM[is_inconvertible] %>% sQuote %>% toString
			), 
			call. = FALSE
		)
	}
    
    # check if all taxonomic groups are specified
    # (limits the selection of endofauna species) 
    is_missing <- is.na(.data$GROUP)
	if (any(is_missing)) {
		warning(
			sprintf(
			    paste(
				    "The taxonomic group is missing for %d records",
                    "in the taxon-file.\n",
				    "This may limit the removal of specific groups",
                    "like decapoda and insecta."
				),
				sum(is_missing)
			), 
			call. = FALSE
		)
	}
    
    # consistency check: "TO" should be part of only one "GROUP"
    d <- .data %>% 
        select_(~GROUP, ~TO) %>%
        distinct_ %>%
        group_by_(~TO) %>%
        summarise_(n = ~n()) %>%
        filter_(~n != 1L)
    is_inconsistent <- nrow(d) > 0L
    if (is_inconsistent) {
        stop(
            sprintf(
                paste(
                    "Inconsistencies found in taxa file.",
                    "Please check groups for taxa:\n %s"
                ),
                d$TO %>% sQuote %>% toString
            ),
            call. = FALSE
        )
    }

    # consistency check: "TO" should be part of only one "LEVEL"
    # not yet implemented since LEVEL refers to FROM and not yet to TO
    
    # return result
    .data
}









#' Read and Validate Taxa Data
#'
#' This function reads files in the taxa format.
#' 
#' @param filename name of taxa file
#'
#' @details Taxa files have the following format:
#'  \itemize{
#' 		\item{group} {taxonomic group}
#'  	\item{provided} {provided taxon name}
#'  	\item{accepted} {accepted taxon name}
#'      \item{level} {taxonomic level}
#' 	}
#' 	Other columns are allowed, but silently ingored.
#'
#' @importFrom readr read_csv cols_only col_character
#' 
#' @export
read_taxa <-
function(filename) {

    # raise error if filename is missing
    if (missing(filename)) {
        stop("Filename is missing", call. = FALSE)
    }
    
	# check if file exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file
	d <- try(
	    suppressMessages(read_csv(filename)),
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

    # validate taxa-file
	d %>%
        validate_taxa
}

#' @describeIn read_taxa get default taxa list (TWN list extended with species Southern North Sea)
#' @importFrom readr read_rds
#' @export
get_taxa <-
function() {
    read_rds(system.file("extdata", "taxa.rds", package = "benthos"))
}



#' @describeIn read_taxa validator for taxa-format
#' @param .data table in taxa-format
#' @import dplyr
validate_taxa <-
function(.data) {
    
    # convert to tibble
    .data <- .data %>% 
        as_tibble

	# check column names
	required_vars <- c("group", "provided", "accepted", "level")
	missing_vars <- setdiff(required_vars, tolower(names(.data)))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missing_vars)
            ),
            call. = FALSE
        )
	}
	names(.data) <- tolower(names(.data))

	# select only columns of interest
	.data <- .data %>% 
	    select_(.dots = required_vars)

	# create ordered factor of taxon levels
	.data$level <- factor(
		x = .data$level, 
		levels = c(
			"Regio", 
			"Regnum", 
			"Phylum", "Subphylum", 
			"Classis", "Subclassis", "Infraclassis", 
			"Ordo", "Subordo", "Infraordo", 
			"Superfamilia",	"Familia", "Subfamilia", 
			"Tribe", 
			"Genus", "Genus combi", "Subgenus", 
			"Species", "Species combi", "Subspecies", 
			"Varietas", 
			"Forma"),
		ordered = TRUE
	)

    # check if 'provided' can be converted to 'accepted'
    is_inconvertible <- is.na(.data$provided) | is.na(.data$accepted)
	if (any(is_inconvertible)) {
		stop(
			sprintf(
				"A total of %d taxon names in the taxon-file are inconvertible:\n%s", 
				sum(is_inconvertible),
				.data$provided[is_inconvertible] %>% sQuote %>% toString
			), 
			call. = FALSE
		)
	}
    
    # check if all taxonomic groups are specified
    # (limits the selection of endofauna species) 
    is_missing <- is.na(.data$group)
	if (any(is_missing)) {
		warning(
			sprintf(
			    paste(
				    "The taxonomic group is missing for %d records",
                    "in the taxon-file.\n",
				    "This may limit the removal of specific groups",
                    "like decapoda and insecta."
				),
				sum(is_missing)
			), 
			call. = FALSE
		)
	}
    
    # consistency check: "accepted" should be part of only one "group"
    d <- .data %>% 
        select_(~group, ~accepted) %>%
        distinct_ %>%
        group_by_(~accepted) %>%
        summarise_(n = ~n()) %>%
        filter_(~n != 1L)
    is_inconsistent <- nrow(d) > 0L
    if (is_inconsistent) {
        stop(
            sprintf(
                paste(
                    "Inconsistencies found in taxa file.",
                    "Please check groups for taxa:\n %s"
                ),
                d$accepted %>% 
                    sQuote %>% 
                    toString
            ),
            call. = FALSE
        )
    }

    # consistency check: 'accepted' should be part of only one 'level'
    # not yet implemented since 'level' refers to 'provided' and not yet to 'level'
    # (will be checked by RWS)
    
    # return result
    .data
}




#' Read and Validate AMBI Sensitivity Data
#'
#' This function reads and checks files with AMBI sensitivity data. 
#' The data should be stored in 'comma separated values' format (csv) 
#' consisting of two columns:
#' \itemize{
#' 		\item{TAXON} {species name;}
#'  	\item{GROUP} {Roman numeral (I, II, III, IV, V) giving the sensitivity
#'      group}
#' }
#'
#' @param filename name of the AMBI sensitivity file (character)
#'
#' @details The function performs the following tasks:
#' \itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces;}
#'  	\item{removes duplicated records.}
#'  }
#'
#' @references Borja, A., J. Franco and V. Perez, 2000. A Marine Biotic Index 
#'  to Establish the Ecological Quality of Soft-Bottom Benthos Within 
#'  European Estuarine and Coastal Environments. 
#'  Marine Pollution Bulletin 40:1100-1114
#'  
#' @importFrom utils read.csv
#' 
#' @export
read_ambi <-
function(filename) {

    # check filename
    if (missing(filename)) {
        stop("Filename is missing", call. = FALSE)
    }

    # read AMBI-file
    if (!(file.exists(filename))) {
    	stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}
	d <- try(
        read.csv(file = filename, as.is = TRUE, na.strings = c("NA", "")), 
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

    # validate AMBI-file
    validate_ambi(d)
}


#' Get Supplementary AMBI Sensitivity Groups
#'
#' This function gets sensitivity groups that are supplementary to the AMBI of
#' Borja et al., (2000)
#'  
#' @param which which AMBI supplement? Currently only the Dutch supplement is
#'     available (\code{which} = "NL")
#'  
#' @return a data frame with columns \code{TAXON} containing taxa and 
#'     \code{GROUP} containing Dutch AMBI-groups
#'  
#' @references Borja, A., J. Franco and V. Perez, 2000. A Marine Biotic Index 
#' to Establish the Ecological Quality of Soft-Bottom Benthos Within 
#' European Estuarine and Coastal Environments. 
#' Marine Pollution Bulletin 40:1100-1114
#'
#' @importFrom utils read.csv
#'
#' @export
get_ambi <-
function(which = "NL") {
    switch(which,
        NL = readRDS(system.file("extdata", "ambi_nl.rds", package = "benthos")),
        stop("This AMBI-supplement is not available", call. = FALSE)
    )
}

.get_ambi <-
function() {
    frame_number <- sys.nframe()
    if (frame_number == 1L) {
        stop("Internal function, not to be called directly", call. = FALSE)
    }
    readRDS(system.file("extdata", "azti.rds", package = "benthos"))
}


#' @describeIn read_ambi validator for AMBI-format
#' @param .data table in AMBI-format
validate_ambi <-
function(.data) {
    .validate_groups(
        .data = .data, 
        permissable_groups = c("I", "II", "III", "IV", "V")
    )
}

#' @import dplyr
#' @import lazyeval
.validate_groups <-
function(.data, permissable_groups) {

    # convert to tibble
    .data <- .data %>% as_tibble

    # check column names (case insensitive)
    names(.data) <- toupper(names(.data))
	required_vars <- c("TAXON", "GROUP")
	missing_vars <- setdiff(required_vars, names(.data))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missing_vars)
            ),
            call. = FALSE
        )
	}
    
    # remove sp. from binomen with unknown species
    .data$TAXON <- .data$TAXON %>% strip_sp

    # remove redundant spaces
    .data <- .data %>% 
        mutate_if(is.character, funs_(~strip_spaces))

    # remove duplicated records
    n <- nrow(.data)
    .data <- .data %>% distinct
    n <- n - nrow(.data)
    if (n != 0L) {
		message(sprintf(fmt = "Number of duplicated records: %i", n))
		message("These will be removed")
	}

    # check on duplicated taxa
    if (anyDuplicated(.data$TAXON)) {
    	index <- which(duplicated(.data$TAXON))
		stop(
            sprintf(
                fmt = "Duplicated taxonnames found: %s", 
                toString(.data$TAXON[index])
            ),
            call. = FALSE
		)
    }

    # check if all groups are allowed
    if (!all(.data$GROUP %in% permissable_groups)) {
		stop(
            sprintf(
                fmt = "sensitivity groups should be in {%s}",
                toString(permissable_groups)
            ),
            call. = FALSE
		)
    }

    # return result
	.data
}



#' Get Infaunal Trophic Index
#'
#' This function gets the sensitivity groups to estimate the infaunal 
#' trophic index of Gittenberger et al., (2011)
#'  
#' @return a data frame with columns \code{TAXON} containing taxa and 
#'     \code{GROUP} containing the ITI-groups of Gittenberger & Van Loon (2013).
#'  
#' @references Gittenberger A. and  W. van Loon, 2013. 
#'     Sensitivities of marine macrozoobenthos to environmental pressures 
#'     in the Netherlands. Nederlandse Faunistische 
#'     Mededelingen 41: 79-112.
#'
#' @export
get_iti <-
function() {
    readRDS(system.file("extdata", "iti.rds", package = "benthos"))
}


#' Read and Validate Infaunal Trophic Index Files
#'
#' This function reads and checks files containing Infaunal Trophic Index 
#' (ITI) data (Gittenberger & Van Loon, 2013)
#'
#' @param filename name of the ITI file (character).
#'
#' @details The function performs the following tasks:
#' \itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive), 
#'          i.e., TAXON and GROUP;}
#'  	\item{removes redundant spaces;}
#'  	\item{removes duplicated records.}
#'      \item{checks if all ITI classes are I, II, III, or IV}
#' }
#' The column 'GROUP' contains the Roman numerals I, II, III, and IV, with
#' the following meaning:
#' \itemize{
#'     	\item{  I: } {suspension feeders;}
#'  	\item{ II: } {interface feeders;}
#'  	\item{III: } {surface deposit feeders;}
#'      \item{ IV: } {subsurface deposit feeders.}
#' }
#'  
#' @return A data frame with columns \code{TAXON} containing taxa and 
#'      \code{GROUP} containing user-defined ITI-groups
#'       (see Gittenberger & Van Loon, 2013).
#'  
#' @references Gittenberger A. and  W. van Loon, 2013. 
#'      Sensitivities of marine macrozoobenthos to environmental pressures 
#'      in the Netherlands. Nederlandse Faunistische 
#'      Mededelingen 41: 79-112.
#'
#' @importFrom utils read.csv
#'
#' @export
read_iti <-
function(filename) {

    # check if file exists
    if (missing(filename)) {
        stop("Filename is missing", call. = FALSE)
    }
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file
	d <- try(
        read.csv(file = filename, as.is = TRUE, na.strings = c("NA", "")), 
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}
    
    # validate ITI-file
    validate_iti(d)
}



#' @describeIn read_iti validator for ITI-format
#' @param .data table in ITI-format
validate_iti <- 
function(.data) {
    .validate_groups(
        .data = .data, 
        permissable_groups = c("I", "II", "III", "IV")
    )
}



#' Read and Validate Habitat References Files
#'
#' This function reads and checks files with reference values
#'
#' @param filename name of the habitat reference file (\code{character})
#' @param indicators indicators to be processed (\code{character}, 
#'      see details)
#'
#' @details The function performs the following tasks:
#' \itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces}
#'  	\item{removes duplicated records}
#'  }
#'  
#' Argument \code{indicators} is a \code{character} vector of additional benthic 
#' indicators to be checked for. For example, if \code{indicators = "ITI"}, then
#' the habitat reference file should also contain columns ITIREF and ITIBAD.
#' Implemented indicators are N, LNN, S, D, SN, SNA, H, L, AMBI, ITI, PIE, N2 (see package vignette).
#'  
#' The format of the habitat reference file is documented in the 
#' BEQI2-package vignette.
#'
#' @references Van Loon, W, 2013. Loon2013-BEQI2-Specs-Ecotopes-27nov.doc
#'
#' @importFrom utils read.csv
#'      
#' @export
read_ref <-
function(filename, indicators = c("S", "H", "AMBI")) {

	# check if 'filename' exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file
	d <- try(
        read.csv(file = filename, as.is = TRUE, na.strings = c("NA", "")),  
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}
    
    # validate REF-file
    validate_ref(d, indicators = indicators)
}



#' @describeIn read_ref validator for REF-format
#' @param .data table in REF-format
validate_ref <- 
function(.data, indicators = c("S", "H", "AMBI")) {
    
    # convert to tibble
    .data <- .data %>% as_tibble
    
    # check indicators
    valid_indicators <- c("N", "LNN", "S", "D", "SN", "SNA", "H", "L", 
                          "AMBI", "ITI", "PIE", "N2")
    indicators <- toupper(indicators)
    if (length(indicators) == 0L) {
    	stop(
            sprintf(
                fmt = "No indicators specified. Select a subset from:\n%s", 
                toString(sQuote(valid_indicators))
            ),
            call. = FALSE
        )
    }
    invalid_indicators <- indicators[!(indicators %in% valid_indicators)]
    if (length(invalid_indicators) > 0L) {
        stop(
            sprintf(
                fmt = "Invalid indicators found:\n%s", 
                toString(sQuote(invalid_indicators))
            ),
            call. = FALSE
        )
    }

	# check column names (case insensitive)
    names(.data) <- toupper(names(.data))
    required_vars <- c(c("OBJECTID", "RELAREA", "HABITAT"),
        paste0(rep(indicators, each = 2), c("REF", "BAD"))
    )
	missing_vars <- setdiff(required_vars, names(.data))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(sQuote(missing_vars))
            ),
            call. = FALSE
        )
	}
    
    # keep only required variables
    .data <- .data[, required_vars]

    # remove redundant spaces
    .data <- .data %>% 
        mutate_if(is.character, funs_(~strip_spaces))

    # remove duplicated records
    n <- nrow(.data)
    .data <- .data %>% distinct
    n <- n - nrow(.data)
    if (n != 0L) {
		message(sprintf(fmt = "Number of duplicated records: %i", n))
		message("These will be removed")
	}

    # harmonize columns to make case-insensitive matching possible
    .data <- .data %>% 
        mutate_if(is.character, funs_(~harmonize))

    # return result
	.data
}
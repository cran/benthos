#' Bray-Curtis Dissimilarity 
#'
#' @param n1 abundances of species at site 1
#' @param n2 abundances of species at site 2
#' 
#' @note species in n1 and n2 need to be aligned
#'
#' @return Bray-Curtis dissimilarity (0..1, 0 = equal, 1 = different)
#' 
#' @examples
#'      n1 <- c(11,  0, 7,  8, 0)
#'      n2 <- c(24, 37, 5, 18, 1)
#'      bray_curtis(n1, n2)
#' @export
bray_curtis <- 
function(n1, n2) {
    sum(abs(n1-n2)) / sum(n1+n2)
}

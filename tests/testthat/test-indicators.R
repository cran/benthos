context("Biodiversity Indicators")



test_that("'abundance' has been correctly implemented", {

    # reference data: unique taxa
    taxon <- c("a", "b", "c", "d", "e")
    count <- c( 3,   7,  11,  17,  13)
    expected_unnamed <- count
    d <- data.frame(TAXON = taxon, COUNT = count)
    names(count) <- taxon
    expected_named <- count

    # unit tests: unique taxa
    expect_identical(abundance(count = count), expected_unnamed)
    expect_identical(abundance(taxon = taxon, count = count), expected_named)
    expect_identical(abundance(d, count = COUNT), expected_unnamed)
    expect_identical(abundance(d, taxon = TAXON, count = COUNT), expected_named)
    expect_identical(abundance(d, TAXON, COUNT), expected_named)

    # reference data: duplicated taxa
    taxon <- c("a", "b", "c", "a", "b")
    count <- c( 3,   7,  11,  17,  13)
    expected_unnamed <- count
    d <- data.frame(TAXON = taxon, COUNT = count)
    names(count) <- taxon
    tmp <- tapply(X = count, INDEX = taxon, FUN = sum)
    expected_named <- as.numeric(tmp)
    names(expected_named) <- names(tmp)
    
    # unit tests: duplicated taxa
    expect_identical(abundance(count = count), expected_unnamed)
    expect_identical(abundance(taxon = taxon, count = count), expected_named)
    expect_identical(abundance(d, count = COUNT), expected_unnamed)
    expect_identical(abundance(d, taxon = TAXON, count = COUNT), expected_named)
    expect_identical(abundance(d, TAXON, COUNT), expected_named)
    

    # exception handling
    expect_error(
        abundance(taxon = c("a", "b", "c"), count = c(7, 11)),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        abundance(taxon = c("a", "b", "c"), count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        abundance(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'total abundance' has been correctly implemented", {
    expect_identical(total_abundance(count = c(3, 7, 11, 17, 13)), 51)
    expect_error(
        total_abundance(count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        total_abundance(count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'species richness' has been correctly implemented", {
    expect_identical(
        species_richness(taxon = c("a", "b", "c", "a"), count = c(3, 7, 1, 9)),
        3L
    )
    expect_identical(
        species_richness(taxon = c("a", "b", "c", "a")),
        3L
    )
    expect_identical(
        species_richness(taxon = c("a", "b", "c"), count = c(3, 0, 1)),
        2L
    )
    expect_error(
        species_richness(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9)),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        species_richness(taxon = c("a", "b", "c"), count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        species_richness(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'Margalef index of diversity' has been correctly implemented", {
    expect_equal(
        margalef(taxon = c("a", "b", "c", "a"), count = c(3, 7, 1, 2)),
        (3 - 1)/log(13)
    )
    expect_error(
        margalef(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9)),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        margalef(taxon = c("a", "b", "c"), count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        margalef(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'Margalef index of diversity' has been correctly implemented", {
    expect_equal(
        rygg(taxon = c("a", "b", "c", "a"), count = c(3, 7, 1, 2)),
        log(3)/log(log(13))
    )
    expect_identical(
        rygg(taxon = c("a"), count = 0),
        NA_real_
    )
    expect_identical(
        rygg(taxon = c("a"), count = 1),
        NA_real_
    )
    expect_identical(
        rygg(taxon = c("a"), count = 2),
        NA_real_
    )
    expect_identical(
        rygg(taxon = c("a"), count = 3),
        0
    )
    expect_identical(
        rygg(taxon = c("a"), count = 3, adjusted = TRUE),
        0
    )
    expect_equal(
        rygg(taxon = c("a", "b", "c", "a"), count = c(3, 7, 1, 2), adjusted = TRUE),
        log(3)/log1p(log1p(13))
    )
    expect_error(
        rygg(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9)),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        rygg(taxon = c("a", "b", "c"), count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        rygg(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'Shannon index' has been correctly implemented", {
    expect_equal(
        shannon(taxon = c("a", "b", "c"), count = c( 3,   7,  11)),
        -sum(
             3/21 * log2( 3/21) + 
             7/21 * log2( 7/21) + 
            11/21 * log2(11/21)
        )
    )
    expect_error(
        shannon(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9)),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        shannon(taxon = c("a", "b", "c"), count = c(7, 11, -1)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        shannon(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'AMBI-index' has been correctly implemented", {
    expect_equal(
        ambi(
            taxon = letters[1:5],
            count = c( 3,    7,    11,   17,   13),
            group = c("I", "II", "III", "IV", "III")
        ),
        1.5 * sum((0:4) * c(I = 3, II = 7, III = 24, IV = 17, V = 0) / 51)
    )
    expect_error(
        ambi(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9), 
             group = c("I", "II", "III", "IV", "III")),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        ambi(taxon = c("a", "b", "c"), count = c(3, 7, 1), 
             group = c("X", "Y", "III")),
       "'group' should be one of: I, II, III, IV, V"
    )
    expect_error(
        ambi(taxon = c("a", "b", "c"), count = c(7, 11, -1), 
             group = c("I", "II", "III")),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        ambi(taxon = c("a", "b", "c"), count = c(7, NA, 3), 
             group = c("I", "II", "III")),
        "all elements in 'count' should be nonnegative"
    )
    expect_true(has_ambi(taxon = "a", group = "I"))
    expect_false(has_ambi(taxon = "a"))
})



test_that("'ITI-index' has been correctly implemented", {
    expect_equal(
        iti(
            taxon = letters[1:5],
            count = c( 3,    7,    11,   17,   13),
            group = c("I", "II", "III", "IV", "III")
        ),
        100/3 * sum((3:0) * c(I = 3, II = 7, III = 24, IV = 17) / 51)
    )
    expect_error(
        iti(taxon = c("a", "b", "c"), count = c(3, 7, 1, 9), 
             group = c("I", "II", "III", "IV", "III")),
       "'taxon' and 'count' should have the same length"
    )
    expect_error(
        iti(taxon = c("a", "b", "c"), count = c(3, 7, 1), 
             group = c("X", "Y", "III")),
       "'group' should be one of: I, II, III, IV"
    )
    expect_error(
        iti(taxon = c("a", "b", "c"), count = c(7, 11, -1), 
             group = c("I", "II", "III")),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        iti(taxon = c("a", "b", "c"), count = c(7, NA, 3), 
             group = c("I", "II", "III")),
        "all elements in 'count' should be nonnegative"
    )
    expect_true(has_iti(taxon = "a", group = "I"))
    expect_false(has_iti(taxon = "a"))
})



test_that("'Hurlbert index' has been correctly implemented", {
    expect_equal(
        hurlbert(taxon = c("a", "b", "c", "d", "e"), count = c(96, 1, 1, 1, 1)),
        5
    )
    expect_equal(
        hurlbert(taxon = c("a", "b", "c"), count = c(100, 100, 100), n = 1),
        1
    )
    expect_equal(
        hurlbert(taxon = c("a", "b", "c"), count = c(100, 100, 100), n = 2),
        1.66889632107
    )
    expect_equal(
        hurlbert(taxon = c("a", "b", "c"), count = c(100, 100, 100), n = 5),
        2.61155015984
    )
    expect_equal(
        hurlbert(taxon = c("a", "b", "c"), count = c(100, 0, 100), n = 100),
        2
    )
    expect_error(
        hurlbert(taxon = c("a", "b", "c"), count = c(7, -1, 3)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        hurlbert(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'Simpson's index' has been correctly implemented", {
    expect_equal(
        simpson(taxon = c("a", "b", "c"), count = c(1000, 1000, 1000)),
        3 * (1/3)^2, tolerance = 0.001
    )
    expect_gt(
        simpson(taxon = c("a", "b", "c"), count = c(1000, 1, 1)),
        0.99
    )
    expect_error(
        simpson(taxon = c("a", "b", "c"), count = c(7, -1, 3)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        simpson(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'Hill's diversity number' has been correctly implemented", {
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = -Inf),
        1
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(1, 10, 100), a = -Inf),
        1
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = Inf),
        1
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(1, 10, 100), a = Inf),
        1
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = 0),
        3
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = 0),
        hill0(taxon = c("a", "b", "c"), count = c(100, 100, 100))
    )
    expect_equal(
        hill0(taxon = c("a", "b", "c"), count = c(100, 100, 100)),
        species_richness(taxon = c("a", "b", "c"), count = c(100, 100, 100))
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(1, 10, 100), a = 0),
        3
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = 1),
        hill1(taxon = c("a", "b", "c"), count = c(100, 100, 100))
    )
    expect_message(
        hill(taxon = c("a", "b", "c"), count = c(100, 100, 100), a = 1),
        "N_a(a=1) is undefined. Therefore N_a(lim a->1) will be returned",
        fixed = TRUE
    )
    expect_equal(
        hill1(taxon = c("a", "b", "c"), count = c(100, 100, 100)),
        exp(shannon(taxon = c("a", "b", "c"), count = c(100, 100, 100), base = exp(1)))
    )
    expect_equal(
        hill1(taxon = c("a", "b", "c"), count = c(1, 10, 100)),
        exp(shannon(taxon = c("a", "b", "c"), count = c(1, 10, 100), base = exp(1)))
    )
    expect_equal(
               hill(taxon = c("a", "b", "c"), count = c(1000, 1000, 1000), a = 2),
        1 / simpson(taxon = c("a", "b", "c"), count = c(1000, 1000, 1000)),
        tolerance = 0.001
    )
    expect_equal(
               hill(taxon = c("a", "b", "c"), count = c(1, 10, 100), a = 2),
        1 / simpson(taxon = c("a", "b", "c"), count = c(1, 10, 100)),
        tolerance = 0.01
    )
    expect_equal(
        hill(taxon = c("a", "b", "c"), count = c(1000, 1000, 1000), a = 2),
        hill2(taxon = c("a", "b", "c"), count = c(1000, 1000, 1000)),
        tolerance = 0.001
    )
    expect_equal(
        hill( taxon = c("a", "b", "c"), count = c(10, 100, 1000), a = 2),
        hill2(taxon = c("a", "b", "c"), count = c(10, 100, 1000)),
        tolerance = 0.001
    )
    expect_error(
        hill(taxon = c("a", "b", "c"), count = c(7, -1, 3)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        hill(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'hpie' has been correctly implemented", {
    expect_equal(
        hpie(taxon = c("a", "b"), count = c(100, 100)),
        sum((c(100, 100) / (100 + 100)) * 
                ((100 + 100) -  c(100, 100)) / 
                ((100 + 100) - 1)) # Hurlbert's Delta1
    )
    expect_equal(
        hpie(taxon = c("a", "b"), count = c(13, 71)),
        sum((c(13, 71) / (13 + 71)) * 
                ((13 + 71) - c(13, 71)) / 
                ((13 + 71) - 1)) # Hurlbert's Delta1
    )
    expect_equal(
        hpie(taxon = c("a", "b"), count = c(10, 100)),
        1 - simpson(taxon = c("a", "b"), count = c(10, 100))
    )
    expect_error(
        hpie(taxon = c("a", "b", "c"), count = c(7, -1, 3)),
        "all elements in 'count' should be nonnegative"
    )
    expect_error(
        hpie(taxon = c("a", "b", "c"), count = c(7, NA, 3)),
        "all elements in 'count' should be nonnegative"
    )
})



test_that("'hpie' has been correctly implemented (test by simulation)", {
    skip_on_cran()
    set.seed(314)
    N <- 1000000L
    d <- rep.int(x = 1:2, times = c(100, 100))
    r <- replicate(n = N, {sample(x = d, size = 2L)})
    expect_equal(
        hpie(taxon = c("a", "b"), count = c(100, 100)),
        mean(r[1, ] !=  r[2, ]),
        tolerance = 1.0e-3
    )
})


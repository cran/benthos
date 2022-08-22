context("Utilities")


test_that("binomial names are correctly tested and split", {

    expect_true(is_binomen("Venerupis corrugata"))
    expect_false(is_binomen("venerupis corrugata"))
    expect_false(is_binomen("Venerupis"))
    expect_that(
        is_binomen(c("Venerupis corrugata", "Urothoe poseidonis")),
        is_identical_to(c(TRUE, TRUE))
    )
    expect_that(
        is_binomen(c("Venerupis", "Urothoe poseidonis")),
        is_identical_to(c(FALSE, TRUE))
    )
    expect_that(
        is_binomen(c("Venerupis", "Urothoe")),
        is_identical_to(c(FALSE, FALSE))
    )
    expect_true(
        is_binomen("Venerupis sp."),
        info = "genus known, species unsure"
    )
    expect_true(
        is_binomen("Venerupis spp."),
        info = "Species pluralis, multiple species"
    )
    expect_false(is_binomen("Venerupis xxx."))
    expect_false(is_binomen("Venerupis sppp."))
    expect_true(
        is_binomen("Venerupis sp"),
        info = "genus known, species unsure, period is missing"
    )
    expect_true(
        is_binomen("Venerupis spp"),
        info = "Species pluralis, multiple species, period is missing"
    )
    expect_false(is_binomen("Venerupis sppp."))
    
    expect_that(
        generic_name("Venerupis"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        generic_name("Venerupis corrugata"), 
        is_identical_to("Venerupis")
    )
    expect_that(
        generic_name("venerupis corrugata"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        generic_name(c("Venerupis corrugata", "Urothoe poseidonis")), 
        is_identical_to(c("Venerupis", "Urothoe"))
    )
    expect_that(
        generic_name(c("OLIGOCHAETA", "Pygospio elegans")), 
        is_identical_to(c(NA_character_, "Pygospio"))
    )
    
    expect_that(
        specific_name("Venerupis"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name("Venerupis sp."), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name("Venerupis spp."), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name("Venerupis sp"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name("Venerupis spp"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name("Venerupis spx"), 
        is_identical_to("spx")
    )
    expect_that(
        specific_name("Venerupis corrugata"), 
        is_identical_to("corrugata")
    )
    expect_that(
        specific_name("venerupis corrugata"), 
        is_identical_to(NA_character_)
    )
    expect_that(
        specific_name(c("Venerupis corrugata", "Urothoe poseidonis")), 
        is_identical_to(c("corrugata", "poseidonis"))
    )
    expect_that(
        specific_name(c("OLIGOCHAETA", "Pygospio elegans")), 
        is_identical_to(c(NA_character_, "elegans"))
    )

    expect_that(
        strip_sp("Venerupis sp"), 
        is_identical_to("Venerupis")
    )
    expect_that(
        strip_sp("Venerupis sp."), 
        is_identical_to("Venerupis")
    )
    expect_that(
        strip_sp("Venerupis spp."), 
        is_identical_to("Venerupis")
    )
    expect_that(
        strip_sp("Venerupis spp.x"), 
        is_identical_to("Venerupis spp.x")
    )
    
})




test_that("stripping spaces and harmonization", {

    expect_that(
        strip_spaces(c(" Hello  World  ", NA_character_)), 
        is_identical_to(c("Hello World", NA_character_))
    )

    expect_that(
        harmonize(c("FOO", "Foo", "bar", "FOO", "bar", "FOO", "Bar")), 
        is_identical_to(c("FOO", "FOO", "bar", "FOO", "bar", "FOO", "bar"))
    )
})




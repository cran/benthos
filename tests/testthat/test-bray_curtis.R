context("Bray Curtis")



test_that("'bray_curtis dissimilarity' has been correctly implemented", {
    
    # test example https://www.statisticshowto.datasciencecentral.com/bray-curtis-dissimilarity/
    # (but different implementation) 
    expect_equal(bray_curtis(c(6, 7, 4), c(10, 0, 6)), 0.39, tolerance = 0.01)
    
    # test extremes
    expect_identical(bray_curtis(c(0, 0, 0), c(5, 5, 5)), 1)
    expect_identical(bray_curtis(c(5, 5, 5), c(5, 5, 5)), 0)
})

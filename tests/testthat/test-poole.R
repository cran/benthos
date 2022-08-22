context("Pooling")



test_that("no pooling occurs if all sample areas are within the target area.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 10L
        area <- rep.int(x = 0.10, times = n_samples)
        sample_id <- 1:n_samples
    
    	expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            shows_message("No pooling necessary.")
        )
        expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            is_identical_to(sample_id)
        )
})



test_that("no pooling occurs if all sample areas are greater than the target area.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 10L
        area <- rep.int(x = 0.20, times = n_samples)
        sample_id <- 1:n_samples
    
        expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            throws_error("No pooling possible")
        )
})



test_that("pooling occurs for trivial cases.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 10L
        area <- rep.int(x = 0.01, times = n_samples)
        sample_id <- 1:n_samples
    
        expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            is_identical_to(rep.int(x = 1L, times = n_samples))
        )
})



test_that("no pooling occurs if the sum of all sample areas is below the target area.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 10L
        area <- rep.int(x = 0.001, times = n_samples)
        sample_id <- 1:n_samples
    
        expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            throws_error("No pooling possible")
        )
})



test_that("no pooling occurs if all sample areas are slightly smaller than the target area.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 10L
        area <- rep.int(x = 0.08, times = n_samples)
        sample_id <- 1:n_samples
    
        expect_that(
            pool(sample_id = sample_id, area = area, target_area = target_area),
            is_identical_to(rep.int(x = NA_integer_, times = n_samples))
        )
})



test_that("pooling processes all samples if possible.", {

        target_area <- c(0.09, 0.11)
        n_samples <- 8L
        area <- rep.int(x = 0.025, times = n_samples)
        sample_id <- 1:n_samples
    
        for (i in 1:10) {
            expect_false(
                any(is.na(pool(sample_id = sample_id, area = area, target_area = target_area)))
            )
        }
})



test_that("areas of pools are in target interval.", {
        target_area <- c(0.09, 0.11)
        for (i in 1:10) {
            n_samples <- sample(x = 10:250, size = 1)
            sample_id <- 1:n_samples
            area <- runif(n = n_samples, min = 0.01, max = 0.04)
            index <- pool(sample_id = sample_id, area = area, target_area = target_area)
            expect_true(
                all(tapply(X = area, INDEX = index, FUN = sum) %>% 
                        between(target_area[1], target_area[2]))
            )
        }
})
test_that("validate_ena_sra_run is present", {
  expect_true(is.function(validate_ena_sra_run))
})

test_that("validate_gsa_run is present", {
  expect_true(is.function(validate_gsa_run))
})

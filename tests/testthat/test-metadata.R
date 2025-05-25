test_that("validate_and_query_ena recognizes accession types", {
  expect_true(grepl("study_accession", validate_and_query_ena("SRP123456")))
  expect_true(grepl("sample_accession", validate_and_query_ena("SRS123456")))
  expect_error(validate_and_query_ena("INVALID"))
})

test_that("get_ena_sra_metadata is present", {
  expect_true(is.function(get_ena_sra_metadata))
})

test_that("get_gsa_metadata is present", {
  expect_true(is.function(get_gsa_metadata))
})

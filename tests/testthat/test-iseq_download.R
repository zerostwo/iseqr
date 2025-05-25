test_that("iseq_download is present", {
  expect_true("iseq_download" %in% ls("package:iseqr"))
})

test_that("validate_and_init_params handles input and output_dir", {
  tmpdir <- tempdir()
  input <- c("SRR123456", "SRR123457")
  res <- validate_and_init_params(input, tmpdir, merge = "none", database = "ena", threads = 2, parallel = 0, retries = 1)
  expect_true(is.list(res))
  expect_true("accessions" %in% names(res))
  expect_true(dir.exists(res$output_dir))
})

test_that("init_software_paths returns a list", {
  res <- init_software_paths(gzip = FALSE, fastq = FALSE, merge = "none", parallel = 0, aspera = FALSE)
  expect_true(is.list(res))
})

test_that("process_metadata is present", {
  expect_true(is.function(process_metadata))
})

test_that("process_run is present", {
  expect_true(is.function(process_run))
})

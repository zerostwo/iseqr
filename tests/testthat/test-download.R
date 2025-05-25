test_that("execute_download is present", {
  expect_true(is.function(execute_download))
})

test_that("execute_aspera is present", {
  expect_true(is.function(execute_aspera))
})

test_that("get_aspera_key_path is present", {
  expect_true(is.function(get_aspera_key_path))
})

test_that("download_ena_sra_run is present", {
  expect_true(is.function(download_ena_sra_run))
})

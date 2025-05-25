# get_software_path

test_that("get_software_path returns path for existing software", {
  expect_true(nzchar(get_software_path("echo", required = TRUE)))
})

test_that("get_software_path returns NULL or error for non-existent software", {
  expect_error(get_software_path("nonexistent_software", required = TRUE))
  expect_null(get_software_path("nonexistent_software", required = FALSE))
})

# run_command

test_that("run_command runs a simple command and returns correct output", {
  res <- run_command("echo", args = c("hello world"), error_on_status = TRUE)
  expect_equal(res$status, 0)
  expect_match(res$stdout, "hello world")
})

test_that("run_command returns error on bad command", {
  expect_error(run_command("nonexistent_command", args = c("foo"), error_on_status = TRUE))
})

# check_file_valid

test_that("check_file_valid works as expected", {
  # Case 1: File does not exist
  nonexistent_file <- tempfile()
  expect_warning(
    expect_false(check_file_valid(nonexistent_file, "fake_accession")),
    regexp = "File not found"
  )

  # Case 2: File exists but is empty
  empty_file <- tempfile()
  file.create(empty_file)
  expect_warning(
    expect_false(check_file_valid(empty_file, "test_accession")),
    regexp = "File is empty"
  )
  unlink(empty_file)

  # Case 3: File exists and has content
  valid_file <- tempfile()
  writeLines("dummy content", valid_file)
  expect_true(check_file_valid(valid_file, "test_accession"))
  unlink(valid_file)
})

test_that("check_file_valid returns TRUE for non-empty file", {
  tf <- tempfile()
  writeLines("abc", tf)
  expect_true(check_file_valid(tf, "test_accession"))
  unlink(tf)
})

# log_status & is_already_successful

test_that("log_status and is_already_successful work as expected", {
  tmpdir <- tempdir()
  run_id <- paste0("TEST", as.integer(Sys.time()))
  log_status(run_id, status = "success", log_dir = tmpdir)
  expect_true(is_already_successful(run_id, log_dir = tmpdir))
  # Negative case
  expect_false(is_already_successful("NOT_EXIST", log_dir = tmpdir))
})

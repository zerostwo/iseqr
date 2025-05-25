test_that("count_lines works for a small file", {
  tf <- tempfile(fileext = ".fastq")
  writeLines(rep("A", 8), tf)
  expect_equal(count_lines(tf), 8)
  unlink(tf)
})

test_that("count_lines works for empty file", {
  tf <- tempfile(fileext = ".fastq")
  file.create(tf)
  expect_equal(count_lines(tf), 0)
  unlink(tf)
})

test_that("check_fastq_integrity detects valid files", {
  tf1 <- tempfile(fileext = ".fastq")
  tf2 <- tempfile(fileext = ".fastq")
  writeLines(rep("A", 4), tf1)
  writeLines(rep("A", 4), tf2)
  res <- check_fastq_integrity(c(tf1, tf2))
  expect_true(res$valid)
  expect_equal(res$lines, c(4, 4))
  unlink(c(tf1, tf2))
})

test_that("check_fastq_integrity detects invalid files (not multiple of 4)", {
  tf1 <- tempfile(fileext = ".fastq")
  writeLines(rep("A", 5), tf1)
  res <- check_fastq_integrity(tf1)
  expect_false(res$valid)
  unlink(tf1)
})

test_that("check_fastq_integrity detects invalid files (unequal lines)", {
  tf1 <- tempfile(fileext = ".fastq")
  tf2 <- tempfile(fileext = ".fastq")
  writeLines(rep("A", 4), tf1)
  writeLines(rep("A", 8), tf2)
  res <- check_fastq_integrity(c(tf1, tf2))
  expect_false(res$valid)
  unlink(c(tf1, tf2))
})

# convert_sra_to_fastq, compress_fastq, merge_fastq_files 需要依赖外部软件和大文件，建议mock或跳过

test_that("convert_sra_to_fastq is present", {
  expect_true(is.function(convert_sra_to_fastq))
})

test_that("compress_fastq is present", {
  expect_true(is.function(compress_fastq))
})

test_that("merge_fastq_files is present", {
  expect_true(is.function(merge_fastq_files))
})

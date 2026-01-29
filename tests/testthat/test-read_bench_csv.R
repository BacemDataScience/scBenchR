test_that("read_bench_csv reads .csv files and adds provenance", {
  d <- file.path(tempdir(), "scBenchR_test")
  dir.create(d, showWarnings = FALSE)

  readr::write_csv(
    tibble::tibble(
      dataset_id  = "pbmc3k",
      method      = "seurat_v5",
      run_id      = "run_001",
      step        = "step7",
      elapsed_sec = 1.0
    ),
    file.path(d, "bench.csv")
  )

  out <- read_bench_csv(d)
  expect_equal(nrow(out), 1)
  expect_true(".source_file" %in% names(out))
})

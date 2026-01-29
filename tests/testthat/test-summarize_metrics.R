test_that("summarize_metrics summarizes elapsed_sec by dataset/method/step", {
  x <- tibble::tibble(
    dataset_id  = c("pbmc3k","pbmc3k"),
    method      = c("seurat_v5","seurat_v5"),
    run_id      = c("run_001","run_002"),
    step        = c("step7","step7"),
    elapsed_sec = c(10, 14)
  )

  s <- summarize_metrics(x)
  expect_equal(nrow(s), 1)
  expect_equal(s$n, 2)
  expect_equal(s$elapsed_mean, 12)
  expect_equal(s$elapsed_min, 10)
  expect_equal(s$elapsed_max, 14)
})

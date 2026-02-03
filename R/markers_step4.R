# R/markers_step4.R
# Step 4 — Marker detection + benchmarking (Seurat v5 SAFE)
# Design:
#   SMALL  → RNA, layer = "data", serial vs parallel
#   LARGE  → RNA, layer = "data", serial vs parallel
#   MERGED → RNA, layer = "counts", SERIAL ONLY (memory-safe)
#
# NOTE:
# - Requires helper: scbenchr_counts_layers() from R/utils_layers.R
# - Requires timing helpers: scbenchr_time_seconds() from R/utils_timing.R
# - Requires paths helper: scbenchr_init_outdirs() from R/paths.R

#' Find markers per cluster (Seurat v5 layer-aware)
#'
#' @param seu Seurat object with Idents set.
#' @param parallel Use multisession parallelism.
#' @param workers Worker count (only used if parallel=TRUE).
#' @param assay Assay name.
#' @param layer Layer name ("data" for normalized; "counts" for raw).
#' @return data.table with markers for all clusters.
#' @export
scbenchr_find_markers_clusterwise <- function(
    seu,
    parallel = FALSE,
    workers = 4,
    assay = "RNA",
    layer = "data"
) {
  Seurat::DefaultAssay(seu) <- assay
  clusters <- levels(Seurat::Idents(seu))
  stopifnot(length(clusters) > 1)

  one_cluster <- function(cl) {
    mk <- Seurat::FindMarkers(
      object = seu,
      ident.1 = cl,
      assay = assay,
      layer = layer,
      test.use = "wilcox",
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    mk <- data.table::as.data.table(mk, keep.rownames = "gene")
    mk$cluster <- cl
    mk
  }

  if (!parallel) {
    res <- lapply(clusters, one_cluster)
    return(data.table::rbindlist(res, use.names = TRUE, fill = TRUE))
  }

  future::plan(future::multisession, workers = workers)
  res <- future.apply::future_lapply(clusters, one_cluster, future.seed = TRUE)
  future::plan(future::sequential)

  data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
}

#' Top markers per cluster helper
#'
#' @param dt data.table with cluster and avg_log2FC/avg_logFC column.
#' @param n number of genes per cluster.
#' @return data.table of top n markers per cluster.
#' @export
scbenchr_top_markers_per_cluster <- function(dt, n = 10) {
  dt <- data.table::as.data.table(dt)
  fc_col <- grep("avg_log2FC|avg_logFC|log2FC|avg_logFC", names(dt), value = TRUE)[1]
  if (is.na(fc_col)) {
    stop("No fold-change column found. Columns are: ", paste(names(dt), collapse = ", "))
  }
  dt <- dt[order(cluster, -dt[[fc_col]])]
  dt[, utils::head(.SD, n), by = cluster]
}

#' Merge multi-layer counts into a single sparse matrix (for safe DE)
#'
#' @param seu Integrated object with Assay5 RNA having counts.* layers.
#' @param assay Assay name.
#' @return Sparse dgCMatrix of merged counts.
#' @export
scbenchr_merge_counts_layers <- function(seu, assay = "RNA") {
  cnt_layers <- scbenchr_counts_layers(seu, assay = assay)
  stopifnot(length(cnt_layers) >= 2)

  mats <- lapply(cnt_layers, function(ly) {
    Seurat::GetAssayData(seu, assay = assay, layer = ly)
  })

  all_genes <- Reduce(union, lapply(mats, rownames))

  pad_to_union <- function(mat, genes_union) {
    mat <- mat[intersect(rownames(mat), genes_union), , drop = FALSE]
    missing <- setdiff(genes_union, rownames(mat))
    if (length(missing) > 0) {
      zeros <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(mat), sparse = TRUE)
      rownames(zeros) <- missing
      colnames(zeros) <- colnames(mat)
      mat <- rbind(mat, zeros)
    }
    mat[genes_union, , drop = FALSE]
  }

  mats2 <- lapply(mats, pad_to_union, genes_union = all_genes)
  counts_merged <- do.call(cbind, mats2)

  colnames(counts_merged) <- make.unique(colnames(counts_merged))
  counts_merged
}

#' Run Step 4 marker benchmarking for SMALL/LARGE (serial vs parallel) and MERGED (serial counts only)
#'
#' @param seu_small Step3 small object (Idents set).
#' @param seu_large Step3 large object (Idents set).
#' @param seu_merged_integrated Integrated object (Idents set).
#' @param out_dir Output directory.
#' @param workers_use Parallel worker count (if NULL, auto-select).
#' @return Invisible list with worker count and timing table.
#' @export
scbenchr_run_step4 <- function(
    seu_small,
    seu_large,
    seu_merged_integrated,
    out_dir,
    workers_use = NULL
) {
  scbenchr_init_outdirs(out_dir)

  # ---- Windows-safe worker cap ----
  avail <- parallel::detectCores()
  if (is.null(workers_use)) workers_use <- max(2, min(8, avail - 1))
  options(future.globals.maxSize = 6 * 1024^3)

  cat("\nSTEP 4 workers available:", avail, "\n")
  cat("STEP 4 workers used:", workers_use, "\n")

  # ============================================================
  # STEP 4A — SMALL (serial vs parallel, RNA data)
  # ============================================================
  Seurat::DefaultAssay(seu_small) <- "RNA"
  stopifnot("data" %in% SeuratObject::Layers(seu_small[["RNA"]]))


  tmp <- scbenchr_time_seconds(
    scbenchr_find_markers_clusterwise(seu_small, parallel = FALSE, layer = "data")
  )
  markers_small_serial <- tmp$value
  t_small_serial <- tmp$sec

  tmp <- scbenchr_time_seconds(
    scbenchr_find_markers_clusterwise(seu_small, parallel = TRUE, workers = workers_use, layer = "data")
  )
  markers_small_parallel <- tmp$value
  t_small_parallel <- tmp$sec

  data.table::fwrite(markers_small_serial,
                     file.path(out_dir, "markers", "markers_SMALL_serial.csv"))
  data.table::fwrite(markers_small_parallel,
                     file.path(out_dir, "markers", paste0("markers_SMALL_parallel_w", workers_use, ".csv")))
  data.table::fwrite(scbenchr_top_markers_per_cluster(data.table::copy(markers_small_serial), 10),
                     file.path(out_dir, "markers", "top10_markers_SMALL_serial.csv"))

  tim_small <- data.table::data.table(
    dataset = "SMALL",
    mode = c("serial", paste0("parallel_w", workers_use)),
    workers = c(1, workers_use),
    elapsed_sec = c(t_small_serial, t_small_parallel)
  )
  serial_small <- tim_small$elapsed_sec[tim_small$mode == "serial"][1]
  tim_small$speedup_vs_serial <- serial_small / tim_small$elapsed_sec

  # ============================================================
  # STEP 4B — LARGE (serial vs parallel, RNA data)
  # ============================================================
  Seurat::DefaultAssay(seu_large) <- "RNA"
  stopifnot("data" %in% SeuratObject::Layers(seu_large[["RNA"]]))

  tmp <- scbenchr_time_seconds(
    scbenchr_find_markers_clusterwise(seu_large, parallel = FALSE, layer = "data")
  )
  markers_large_serial <- tmp$value
  t_large_serial <- tmp$sec

  tmp <- scbenchr_time_seconds(
    scbenchr_find_markers_clusterwise(seu_large, parallel = TRUE, workers = workers_use, layer = "data")
  )
  markers_large_parallel <- tmp$value
  t_large_parallel <- tmp$sec

  data.table::fwrite(markers_large_serial,
                     file.path(out_dir, "markers", "markers_LARGE_serial.csv"))
  data.table::fwrite(markers_large_parallel,
                     file.path(out_dir, "markers", paste0("markers_LARGE_parallel_w", workers_use, ".csv")))
  data.table::fwrite(scbenchr_top_markers_per_cluster(data.table::copy(markers_large_serial), 10),
                     file.path(out_dir, "markers", "top10_markers_LARGE_serial.csv"))

  tim_large <- data.table::data.table(
    dataset = "LARGE",
    mode = c("serial", paste0("parallel_w", workers_use)),
    workers = c(1, workers_use),
    elapsed_sec = c(t_large_serial, t_large_parallel)
  )
  serial_large <- tim_large$elapsed_sec[tim_large$mode == "serial"][1]
  tim_large$speedup_vs_serial <- serial_large / tim_large$elapsed_sec

  # ============================================================
  # STEP 4C — MERGED integrated (SERIAL ONLY, RNA counts)
  # Robust fix for Seurat v5 Assay5 multi-layer counts
  # ============================================================
  seu_merged <- seu_merged_integrated
  Seurat::DefaultAssay(seu_merged) <- "RNA"

  counts_merged <- scbenchr_merge_counts_layers(seu_merged, assay = "RNA")

  cat("Merged counts matrix:", nrow(counts_merged), "genes x", ncol(counts_merged), "cells\n")

  # Build a NEW single-layer assay for DE only (keeps your cluster labels)
  seu_merged[["RNA_DE"]] <- Seurat::CreateAssayObject(counts = counts_merged)
  Seurat::DefaultAssay(seu_merged) <- "RNA_DE"

  future::plan(future::sequential)

  tmp <- scbenchr_time_seconds(
    scbenchr_find_markers_clusterwise(
      seu_merged,
      parallel = FALSE,
      assay = "RNA_DE",
      layer = "counts"
    )
  )

  markers_merged_serial <- tmp$value
  t_merged_serial <- tmp$sec

  data.table::fwrite(markers_merged_serial,
                     file.path(out_dir, "markers", "markers_MERGED_INTEGRATED_serial.csv"))
  data.table::fwrite(scbenchr_top_markers_per_cluster(data.table::copy(markers_merged_serial), 10),
                     file.path(out_dir, "markers", "top10_markers_MERGED_INTEGRATED_serial.csv"))

  tim_merged <- data.table::data.table(
    dataset = "MERGED_INTEGRATED",
    mode = "serial",
    workers = 1,
    elapsed_sec = t_merged_serial,
    speedup_vs_serial = 1
  )

  cat("MERGED integrated markers completed in", round(t_merged_serial, 2), "seconds\n")

  # ============================================================
  # STEP 4D — Combine timing (publication table)
  # ============================================================
  tim_all <- data.table::rbindlist(list(tim_small, tim_large, tim_merged), use.names = TRUE, fill = TRUE)
  data.table::fwrite(tim_all,
                     file.path(out_dir, "markers", "timing_step4_markers_summary.csv"))

  cat("\nSTEP 4 COMPLETE\n")
  cat(" -", file.path(out_dir, "markers", "timing_step4_markers_summary.csv"), "\n")
  cat(" -", file.path(out_dir, "markers", "markers_*.csv"), "\n")
  cat(" -", file.path(out_dir, "markers", "top10_markers_*.csv"), "\n")

  invisible(list(
    workers_use = workers_use,
    timing = tim_all
  ))
}

#' Worker-grid benchmark (algorithmic insight)
#' Creates: out_dir/bench/worker_grid_markers.csv
#' SMALL + LARGE only (MERGED excluded by design)
#'
#' @param seu_small Step3 small object.
#' @param seu_large Step3 large object.
#' @param out_dir Output directory.
#' @param worker_grid Vector of worker counts to benchmark.
#' @return data.table of benchmark results (invisibly).
#' @export
scbenchr_worker_grid_benchmark <- function(
    seu_small,
    seu_large,
    out_dir,
    worker_grid = c(1, 2, 4, 6, 8)
) {
  scbenchr_init_outdirs(out_dir)

  grid_bench_one <- function(seu, tag, assay = "RNA", layer = "data") {
    Seurat::DefaultAssay(seu) <- assay

    res <- lapply(worker_grid, function(w) {
      if (w == 1) {
        future::plan(future::sequential)
        tmp <- scbenchr_time_seconds(
          scbenchr_find_markers_clusterwise(
            seu,
            parallel = FALSE,
            assay = assay,
            layer = layer
          )
        )
        data.table::data.table(dataset = tag, workers = 1, elapsed_sec = tmp$sec)
      } else {
        future::plan(future::multisession, workers = w)
        tmp <- scbenchr_time_seconds(
          scbenchr_find_markers_clusterwise(
            seu,
            parallel = TRUE,
            workers = w,
            assay = assay,
            layer = layer
          )
        )
        future::plan(future::sequential)
        data.table::data.table(dataset = tag, workers = w, elapsed_sec = tmp$sec)
      }
    })

    dt <- data.table::rbindlist(res)
    serial_sec <- dt$elapsed_sec[dt$workers == 1][1]
    dt$speedup <- serial_sec / dt$elapsed_sec
    dt
  }

  bench_small <- grid_bench_one(seu_small, "SMALL", layer = "data")
  bench_large <- grid_bench_one(seu_large, "LARGE", layer = "data")

  bench_all <- data.table::rbindlist(list(bench_small, bench_large), use.names = TRUE, fill = TRUE)

  data.table::fwrite(bench_all, file.path(out_dir, "bench", "worker_grid_markers.csv"))

  cat(" -", file.path(out_dir, "bench", "worker_grid_markers.csv"), "\n")

  invisible(bench_all)
}

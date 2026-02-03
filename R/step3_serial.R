#' Run serial Seurat pipeline (Normalize -> HVG -> Scale -> PCA -> Neighbors -> Clusters -> UMAP)
#'
#' @param seu Seurat object.
#' @param prefix String prefix for timing step names.
#' @param npcs Number of PCs.
#' @param resolution Clustering resolution.
#' @return list(seu=SeuratObject, timing=data.frame)
#' @export
scbenchr_run_serial_step3 <- function(seu, prefix = "DATA", npcs = 30, resolution = 0.6) {
  future::plan(future::sequential)
  options(future.globals.maxSize = NULL)

  timings <- list()

  tmp <- scbenchr_time_it(paste0(prefix, "_NormalizeData"), Seurat::NormalizeData(seu, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_FindVariableFeatures"),
                          Seurat::FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_ScaleData"),
                          Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu), verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_RunPCA"),
                          Seurat::RunPCA(seu, npcs = npcs, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_FindNeighbors"),
                          Seurat::FindNeighbors(seu, dims = 1:npcs, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_FindClusters"),
                          Seurat::FindClusters(seu, resolution = resolution, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  tmp <- scbenchr_time_it(paste0(prefix, "_RunUMAP"),
                          Seurat::RunUMAP(seu, dims = 1:npcs, verbose = FALSE))
  seu <- tmp$value; timings[[length(timings)+1]] <- tmp$timing

  list(seu = seu, timing = do.call(rbind, timings))
}

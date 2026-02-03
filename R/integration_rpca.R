#' RPCA integration of two Seurat objects (serial)
#'
#' @param seu_small Filtered small object (counts state).
#' @param seu_large Filtered large object (counts state).
#' @param npcs Number of PCs / dims.
#' @param resolution Clustering resolution.
#' @return Integrated Seurat object (DefaultAssay integrated).
#' @export
scbenchr_integrate_rpca <- function(seu_small, seu_large, npcs = 30, resolution = 0.6) {
  future::plan(future::sequential)

  seu_small$batch <- "SMALL"
  seu_large$batch <- "LARGE"
  obj_list <- list(SMALL = seu_small, LARGE = seu_large)

  obj_list <- lapply(obj_list, function(x) {
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, nfeatures = 3000, verbose = FALSE)
    x
  })

  features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

  obj_list <- lapply(obj_list, function(x) {
    x <- Seurat::ScaleData(x, features = features, verbose = FALSE)
    x <- Seurat::RunPCA(x, features = features, npcs = npcs, verbose = FALSE)
    x
  })

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = features,
    reduction = "rpca",
    dims = 1:npcs
  )

  seu_int <- Seurat::IntegrateData(anchorset = anchors, dims = 1:npcs)
  Seurat::DefaultAssay(seu_int) <- "integrated"

  seu_int <- Seurat::ScaleData(seu_int, verbose = FALSE)
  seu_int <- Seurat::RunPCA(seu_int, npcs = npcs, verbose = FALSE)
  seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:npcs, verbose = FALSE)
  seu_int <- Seurat::FindClusters(seu_int, resolution = resolution, verbose = FALSE)
  seu_int <- Seurat::RunUMAP(seu_int, dims = 1:npcs, verbose = FALSE)

  seu_int
}

#' Export merged metadata (cell-level) for paper step
#' @param seu Integrated Seurat object.
#' @param out_dir Output directory.
#' @param filename Output filename in integrated/ subfolder.
#' @export
scbenchr_export_merged_metadata <- function(seu, out_dir, filename = "merged_metadata.csv") {
  scbenchr_init_outdirs(out_dir)
  meta_out <- data.table::as.data.table(seu@meta.data, keep.rownames = "cell")
  meta_out$seurat_clusters <- as.character(Seurat::Idents(seu))

  if (!("batch" %in% names(meta_out))) {
    if ("orig.ident" %in% names(meta_out)) {
      meta_out$batch <- as.character(meta_out$orig.ident)
    } else if ("dataset" %in% names(meta_out)) {
      meta_out$batch <- as.character(meta_out$dataset)
    }
  }
  p <- file.path(out_dir, "integrated", filename)
  data.table::fwrite(meta_out, p)
  invisible(p)
}

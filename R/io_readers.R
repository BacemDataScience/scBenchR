#' Read expression matrix safely (fix duplicated genes/cells)
#'
#' @param path Path to a gzipped or plain text matrix with first column genes and remaining columns cells.
#' @param cell_prefix Prefix to guarantee cell name uniqueness across inputs.
#' @return A numeric matrix with unique rownames and colnames.
#' @export
scbenchr_read_expMat_safe <- function(path, cell_prefix) {
  stopifnot(file.exists(path))
  dt <- data.table::fread(path, data.table = FALSE)

  genes <- make.unique(as.character(dt[[1]]))
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes

  if (is.null(colnames(mat)) || any(colnames(mat) == "")) {
    colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  }
  colnames(mat) <- make.unique(paste0(cell_prefix, "_", colnames(mat)))
  mat
}

#' Load PDAC datasets (GSE111672 small, GSE154778 large) into Seurat objects
#'
#' @param base_dir Base directory containing subfolders GSE111672 and GSE154778.
#' @param min_features_small Minimum features for CreateSeuratObject (small).
#' @param min_features_large Minimum features for CreateSeuratObject (large).
#' @return A list with seu_small and seu_large.
#' @export
scbenchr_load_pdac <- function(
    base_dir = scbenchr_default_base_dir(),
    min_features_small = 200,
    min_features_large = 300
) {
  # ---- SMALL ----
  small_dir <- file.path(base_dir, "GSE111672")
  path_A <- file.path(small_dir, "GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz")
  path_B <- file.path(small_dir, "GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz")

  mat_A <- scbenchr_read_expMat_safe(path_A, "GSE111672_PDACA")
  mat_B <- scbenchr_read_expMat_safe(path_B, "GSE111672_PDACB")
  colnames(mat_B) <- gsub(" ", "", colnames(mat_B))

  all_genes <- union(rownames(mat_A), rownames(mat_B))

  mat_A2 <- matrix(0, nrow = length(all_genes), ncol = ncol(mat_A),
                   dimnames = list(all_genes, colnames(mat_A)))
  mat_B2 <- matrix(0, nrow = length(all_genes), ncol = ncol(mat_B),
                   dimnames = list(all_genes, colnames(mat_B)))

  mat_A2[rownames(mat_A), ] <- mat_A
  mat_B2[rownames(mat_B), ] <- mat_B
  counts_small <- cbind(mat_A2, mat_B2)

  stopifnot(!anyDuplicated(rownames(counts_small)))
  stopifnot(!anyDuplicated(colnames(counts_small)))

  seu_small <- Seurat::CreateSeuratObject(
    counts = counts_small,
    project = "PDAC_GSE111672",
    min.cells = 3,
    min.features = min_features_small
  )
  seu_small$dataset <- "GSE111672"

  # ---- LARGE ----
  large_dir <- file.path(base_dir, "GSE154778")
  large_path <- file.path(large_dir, "GSE154778_dgeMtx.csv.gz")
  stopifnot(file.exists(large_path))

  dt_large <- data.table::fread(large_path, data.table = FALSE)
  genes_large <- make.unique(as.character(dt_large[[1]]))
  mat_large <- as.matrix(dt_large[, -1, drop = FALSE])
  storage.mode(mat_large) <- "numeric"
  rownames(mat_large) <- genes_large

  if (is.null(colnames(mat_large)) || any(colnames(mat_large) == "")) {
    colnames(mat_large) <- paste0("cell", seq_len(ncol(mat_large)))
  }
  colnames(mat_large) <- make.unique(paste0("GSE154778_", colnames(mat_large)))

  seu_large <- Seurat::CreateSeuratObject(
    counts = mat_large,
    project = "PDAC_GSE154778",
    min.cells = 3,
    min.features = min_features_large
  )
  seu_large$dataset <- "GSE154778"

  list(seu_small = seu_small, seu_large = seu_large)
}

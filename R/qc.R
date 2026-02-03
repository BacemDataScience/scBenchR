#' Add percent.mt to a Seurat object
#' @param seu Seurat object.
#' @param pattern Regex to identify mitochondrial genes.
#' @export
scbenchr_add_percent_mt <- function(seu, pattern = "^MT-") {
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = pattern)
  seu
}

#' Filter Seurat object using QC thresholds
#' @param seu Seurat object.
#' @param min_features Minimum nFeature_RNA.
#' @param max_features Maximum nFeature_RNA.
#' @param min_counts Minimum nCount_RNA.
#' @param max_mt Maximum percent.mt.
#' @export
scbenchr_qc_filter <- function(seu, min_features, max_features, min_counts, max_mt) {
  subset(
    seu,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      nCount_RNA   >= min_counts &
      percent.mt   <= max_mt
  )
}


#' QC summary (printed)
#' @param seu Seurat object.
#' @param name Label printed in output.
#' @export
scbenchr_qc_summary <- function(seu, name = "object") {
  md <- seu@meta.data
  cat("\n====================\n")
  cat("QC summary:", name, "\n")
  cat("====================\n")
  cat("Cells:", ncol(seu), " Genes:", nrow(seu), "\n\n")
  cat("nFeature_RNA:\n"); print(summary(md$nFeature_RNA))
  cat("\nnCount_RNA:\n");  print(summary(md$nCount_RNA))
  cat("\npercent.mt:\n");  print(summary(md$percent.mt))
  invisible(TRUE)
}

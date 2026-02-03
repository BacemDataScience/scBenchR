#' Get counts layers for an Assay5 RNA assay
#' @keywords internal
scbenchr_counts_layers <- function(seu, assay = "RNA") {
  grep("^counts\\.", SeuratObject::Layers(seu[[assay]]), value = TRUE)
}

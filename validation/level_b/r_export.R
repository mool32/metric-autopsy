# Level B, step S5 (plan_s5.md, section 4): a deposit held only in an R object (a Seurat or SingleCellExperiment
# object), described or exported by R. Nothing of the engine runs; nothing here computes a claim's metric.
#
#   Rscript validation/level_b/r_export.R describe <object.rds> <out dir>
#       its class, slots, assays and dimensions, and its cells' metadata: each column's type, number of distinct
#       values and its 30 commonest values with their counts (describe.json, meta_summary.tsv)
#   Rscript validation/level_b/r_export.R export <object.rds> <out dir>
#       the counts (the RNA assay's counts layer or slot; genes by cells) as counts.mtx.gz, the genes as
#       genes.tsv.gz, the cells' metadata as cells.csv.gz in the counts' order (data_run.py's r_export loader)
#
# A gzip of an RDS file (GEO's .rds.gz) is read after one gunzip, which the caller does.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3, args[1] %in% c("describe", "export"))
mode <- args[1]; path <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(Matrix))
obj <- readRDS(path)
cls <- class(obj)[1]

meta_of <- function(o) {
  if (inherits(o, "Seurat")) return(o@meta.data)
  if (inherits(o, "SingleCellExperiment")) return(as.data.frame(SummarizedExperiment::colData(o)))
  stop("not a Seurat or SingleCellExperiment object: ", class(o)[1])
}
counts_of <- function(o) {
  if (inherits(o, "Seurat")) {
    suppressPackageStartupMessages(library(SeuratObject))
    a <- o[["RNA"]]
    if (inherits(a, "Assay5")) return(SeuratObject::LayerData(a, layer = "counts"))
    return(SeuratObject::GetAssayData(a, slot = "counts"))
  }
  return(SummarizedExperiment::assay(o, "counts"))
}

meta <- meta_of(obj)
if (mode == "describe") {
  rows <- list()
  for (col in names(meta)) {
    v <- meta[[col]]
    tab <- sort(table(as.character(v), useNA = "ifany"), decreasing = TRUE)
    top <- head(tab, 30)
    rows[[length(rows) + 1]] <- data.frame(column = col, type = class(v)[1], distinct = length(tab),
                                           value = names(top), cells = as.integer(top), stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, rows)
  write.table(summary, file.path(out, "meta_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  info <- list(class = cls, cells = nrow(meta), meta_columns = names(meta),
               slots = methods::slotNames(obj),
               assays = if (inherits(obj, "Seurat")) names(obj@assays) else SummarizedExperiment::assayNames(obj))
  m <- tryCatch(counts_of(obj), error = function(e) NULL)
  info$counts_dim <- if (is.null(m)) "none" else dim(m)
  info$counts_integer <- if (is.null(m)) NA else all(m@x == round(m@x)) && all(m@x >= 0)
  writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE), file.path(out, "describe.json"))
  cat(readLines(file.path(out, "describe.json")), sep = "\n")
} else {
  m <- as(counts_of(obj), "CsparseMatrix")
  stopifnot(identical(colnames(m), rownames(meta)))
  tmp <- file.path(out, "counts.mtx")
  Matrix::writeMM(m, tmp)
  system2("gzip", c("-n", tmp))
  gz <- gzfile(file.path(out, "genes.tsv.gz"), "w"); writeLines(rownames(m), gz); close(gz)
  cells <- cbind(cell = rownames(meta), meta)
  gz <- gzfile(file.path(out, "cells.csv.gz"), "w"); write.csv(cells, gz, row.names = FALSE); close(gz)
  cat("exported", nrow(m), "genes by", ncol(m), "cells\n")
}

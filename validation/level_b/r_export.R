# Level B, step S5 (plan_s5.md, section 4): a deposit held only in an R object (a Seurat or SingleCellExperiment
# object), described or exported by R. Nothing of the engine runs; nothing here computes a claim's metric.
#
#   Rscript validation/level_b/r_export.R describe <object.rds> <out dir>
#       its class, slots, assays and dimensions, and its cells' metadata: each column's type, number of distinct
#       values and its 30 commonest values with their counts (describe.json, meta_summary.tsv)
#   Rscript validation/level_b/r_export.R export <object.rds> <out dir>
#       the counts (the RNA assay's counts layer or slot, or the "counts" assay; genes by cells; a dense matrix made
#       sparse 5,000 cells at a time) as counts.mtx.gz, the genes as genes.tsv.gz (and a SingleCellExperiment's
#       rowData as rowdata.csv.gz), the cells' metadata as cells.csv.gz in the counts' order (data_run.py's r_export
#       loader)
#
# A gzip of an RDS file (GEO's .rds.gz) is read after one gunzip, which the caller does.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3, args[1] %in% c("describe", "export"))
mode <- args[1]; path <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(Matrix))
obj <- readRDS(path)
cls <- class(obj)[1]

# A Seurat object is read by its slots as attributes (attr), by its class's name: an old object's class belongs to
# the package Seurat (not installed), whose definition inherits(), `@` and slotNames() would need (2026-11-11, J25).
is_seurat <- function(o) class(o)[1] == "Seurat"
meta_of <- function(o) {
  if (is_seurat(o)) return(attr(o, "meta.data", exact = TRUE))
  if (inherits(o, "SingleCellExperiment")) return(as.data.frame(SummarizedExperiment::colData(o)))
  stop("not a Seurat or SingleCellExperiment object: ", class(o)[1])
}
counts_of <- function(o) {
  if (is_seurat(o)) {
    a <- attr(o, "assays", exact = TRUE)[["RNA"]]
    if (is.null(a)) stop("no RNA assay")
    if (class(a)[1] == "Assay5") {  # Seurat 5: the layer, named by the assay's features and cells
      m <- attr(a, "layers", exact = TRUE)[["counts"]]
      if (is.null(rownames(m))) rownames(m) <- rownames(attr(a, "features", exact = TRUE))
      if (is.null(colnames(m))) colnames(m) <- rownames(attr(a, "cells", exact = TRUE))
      return(m)
    }
    return(attr(a, "counts", exact = TRUE))
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
               slots = if (is_seurat(obj)) setdiff(names(attributes(obj)), "class") else methods::slotNames(obj),
               assays = if (is_seurat(obj)) names(attr(obj, "assays", exact = TRUE)) else SummarizedExperiment::assayNames(obj))
  m <- tryCatch(counts_of(obj), error = function(e) NULL)
  info$counts_dim <- if (is.null(m)) "none" else dim(m)
  info$counts_class <- if (is.null(m)) "none" else class(m)[1]
  # non-negative integers: every stored value of a sparse matrix; a dense one's first 1,000 cells (said so)
  info$counts_integer <- if (is.null(m)) NA else if (methods::is(m, "sparseMatrix")) {
    all(m@x == round(m@x)) && all(m@x >= 0)
  } else {
    s <- m[, seq_len(min(ncol(m), 1000)), drop = FALSE]; all(s == round(s)) && all(s >= 0)
  }
  info$counts_integer_of <- if (!is.null(m) && !methods::is(m, "sparseMatrix")) "the first 1,000 cells" else "every value"
  writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE), file.path(out, "describe.json"))
  cat(readLines(file.path(out, "describe.json")), sep = "\n")
} else {
  m <- counts_of(obj)
  cn <- if (is_seurat(obj)) rownames(meta) else colnames(obj)
  rn <- if (is_seurat(obj)) rownames(m) else rownames(obj)
  rd <- if (!is_seurat(obj) && inherits(obj, "SingleCellExperiment")) as.data.frame(SummarizedExperiment::rowData(obj)) else NULL
  rm(obj); invisible(gc())
  cat("counts:", class(m)[1], nrow(m), "genes by", ncol(m), "cells\n")
  if (methods::is(m, "sparseMatrix")) {
    m <- as(m, "CsparseMatrix")
  } else {  # a dense matrix, made sparse 5,000 cells at a time (memory)
    blocks <- split(seq_len(ncol(m)), ceiling(seq_len(ncol(m)) / 5000))
    m <- do.call(cbind, lapply(blocks, function(j) as(m[, j, drop = FALSE], "CsparseMatrix")))
    invisible(gc())
  }
  if (is.null(colnames(m))) colnames(m) <- cn
  if (is.null(rownames(m))) rownames(m) <- rn
  stopifnot(identical(colnames(m), rownames(meta)))
  tmp <- file.path(out, "counts.mtx")
  Matrix::writeMM(m, tmp)
  system2("gzip", c("-n", tmp))
  gz <- gzfile(file.path(out, "genes.tsv.gz"), "w"); writeLines(rownames(m), gz); close(gz)
  cat("genes, the first three:", head(rownames(m), 3), "\n")
  if (!is.null(rd) && ncol(rd) > 0) {  # the genes' own table (a symbol beside an id, say), in the counts' order
    cat("rowData's columns:", names(rd), "\n")
    gz <- gzfile(file.path(out, "rowdata.csv.gz"), "w"); write.csv(cbind(gene = rownames(m), rd), gz, row.names = FALSE); close(gz)
  }
  cells <- cbind(cell = rownames(meta), meta)
  gz <- gzfile(file.path(out, "cells.csv.gz"), "w"); write.csv(cells, gz, row.names = FALSE); close(gz)
  cat("exported", nrow(m), "genes by", ncol(m), "cells\n")
}

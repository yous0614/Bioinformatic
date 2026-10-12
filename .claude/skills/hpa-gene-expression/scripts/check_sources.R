#!/usr/bin/env Rscript
# Check: --source=local (computed from bulk TSVs) reproduces --source=web (HPA page chart data) for single_cell plots.
# Runs hpa_plot.R both ways and compares the written tables. Expected: bars exact or within 0.1 (TSV cluster nCPM is
# rounded to 1 decimal); dot colors may differ only where the value is ~1 nCPM.
# Usage (any cwd): Rscript check_sources.R [GENE ...] [--tissues=brain,liver,...]

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(file.path(script_dir, "..", "..", "..", ".."))
args <- commandArgs(trailingOnly = TRUE)
tissues <- strsplit(sub("^--tissues=", "", c(grep("^--tissues=", args, value = TRUE), "--tissues=brain,eye,liver,lung,blood"))[1], ",")[[1]]
genes <- grep("^--", args, value = TRUE, invert = TRUE)
if (!length(genes)) genes <- c("MS4A4A", "CX3CR1", "PTPRC", "APOE", "GFAP", "CD3E", "ALB", "INS", "MYH7", "SNAP25")

run <- function(...) {
  a <- c(file.path(script_dir, "hpa_plot.R"), ...)
  o <- system2("Rscript", a, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(o, "status"))) stop("hpa_plot.R failed: ", paste(a[-1], collapse = " "), "\n", paste(o, collapse = "\n"))
  sub("^Saved: (.*)\\.\\{pdf,png,tsv\\} *$", "\\1", grep("^Saved:", o, value = TRUE))
}
cmp <- function(x, y, by) {
  m <- merge(x, y, by = by, all = TRUE, suffixes = c(".local", ".web"))
  v <- grep("^(value|ncpm)\\.", names(m), value = TRUE)
  diff <- abs(m[[v[1]]] - m[[v[2]]])
  c(n = nrow(m), missing = sum(is.na(diff)), exact = sum(diff < 0.05, na.rm = TRUE), within_0.1 = sum(diff < 0.15, na.rm = TRUE),
    max_diff = max(c(0, diff), na.rm = TRUE),
    color_diff = if ("colored.local" %in% names(m)) sum(m$colored.local != m$colored.web, na.rm = TRUE) else NA)
}
read <- function(f) read.delim(f, check.names = FALSE, encoding = "UTF-8")

res <- do.call(rbind, lapply(genes, function(g) {
  cat(g, "")
  l <- run(paste0("--gene=", g), "--dataset=single_cell"); w <- run(paste0("--gene=", g), "--dataset=single_cell", "--source=web")
  rows <- list(c(gene = g, plot = "overview bars", cmp(read(paste0(l, ".tsv")), read(paste0(w, ".tsv")), "label")),
               c(gene = g, plot = "overview dots", cmp(read(paste0(l, "_dots.tsv")), read(paste0(w, "_dots.tsv")), c("tissue", "group"))))
  for (t in tissues) {
    l <- run(paste0("--gene=", g), "--dataset=single_cell", paste0("--tissue=", t))
    w <- run(paste0("--gene=", g), "--dataset=single_cell", paste0("--tissue=", t), "--source=web")
    rows[[length(rows) + 1]] <- c(gene = g, plot = t, cmp(read(paste0(l, ".tsv")), read(paste0(w, ".tsv")), "label"))
  }
  do.call(rbind, rows)
}))
res <- as.data.frame(res); options(width = 200)
cat("\n\n"); print(res, row.names = FALSE)
num <- function(x) as.numeric(res[[x]])
cat(sprintf("\nTOTAL %d values: missing %d, exact %d (%.1f%%), within 0.1 %d (%.1f%%), max diff %.1f, dot color diffs %d\n",
            sum(num("n")), sum(num("missing")), sum(num("exact")), 100 * sum(num("exact")) / sum(num("n")),
            sum(num("within_0.1")), 100 * sum(num("within_0.1")) / sum(num("n")), max(num("max_diff")), sum(num("color_diff"), na.rm = TRUE)))
ok <- sum(num("missing")) == 0 && max(num("max_diff")) <= 0.1 + 1e-9
cat(if (ok) "OK\n" else "FAIL\n")
quit(status = if (ok) 0 else 1)

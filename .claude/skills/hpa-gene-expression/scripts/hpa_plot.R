#!/usr/bin/env Rscript
# Plot one gene's Human Protein Atlas RNA expression in the HPA website's layout (vertical bars, HPA order and colors).
# Usage (any cwd):
#   Rscript hpa_plot.R --gene=MS4A4A --dataset=brain_sn|tissue|single_cell [--tissue=liver] [--source=local|web]
# --source=local (default): compute from HPA bulk downloads in data/. --source=web: read the chart data embedded in
# the gene's HPA page (single_cell only; fallback, and the reference for check_sources.R).
# Writes figures/<GENE>_<dataset>[_<tissue>][_web].{pdf,png,tsv} (+ _dots.tsv for the single_cell overview).

suppressPackageStartupMessages({ library(ggplot2); library(jsonlite); library(patchwork) })

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(file.path(script_dir, "..", "..", "..", ".."))  # project root: data/ and figures/ live there
skill_dir <- file.path(script_dir, "..")

HPA_VERSION <- "v25.1"
CAPTION <- paste("Data: Human Protein Atlas", HPA_VERSION, "(proteinatlas.org)")

opt <- list(gene = NA, dataset = NA, tissue = NA, source = "local")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (!kv[1] %in% names(opt) || length(kv) != 2) stop("Bad argument: ", a)
  opt[[kv[1]]] <- kv[2]
}
if (is.na(opt$gene) || !grepl("^[A-Za-z0-9._-]+$", opt$gene)) stop("--gene must be a gene symbol or Ensembl ID")
if (!opt$dataset %in% c("brain_sn", "tissue", "single_cell")) stop("--dataset must be one of: brain_sn, tissue, single_cell")
if (!opt$source %in% c("local", "web")) stop("--source must be local or web")
if (opt$source == "web" && opt$dataset != "single_cell") stop("--source=web only applies to --dataset=single_cell")
if (!is.na(opt$tissue) && opt$dataset != "single_cell") stop("--tissue only applies to --dataset=single_cell")
if (!is.na(opt$tissue) && !grepl("^[A-Za-z &-]+$", opt$tissue)) stop("--tissue must be an HPA tissue name, e.g. liver")

# ---------- data access
# HPA bulk downloads, fetched once into data/.
data_file <- function(name, url = sprintf("https://www.proteinatlas.org/download/tsv/%s.tsv.zip", name)) {
  f <- file.path("data", paste0(name, ".tsv"))
  if (!file.exists(f)) {
    message("Downloading ", name, " (first use only)...")
    dir.create("data", showWarnings = FALSE); z <- tempfile(fileext = ".zip")
    download.file(url, z, mode = "wb", quiet = TRUE)
    unzip(z, exdir = "data")
  }
  f
}
# Rows of a big TSV for one Ensembl ID, via grep (the cluster file is 1.4 GB). col: 1-based column holding the ID.
grep_rows <- function(f, ensg, col = 1) {
  pat <- if (col == 1) sprintf("'^%s\t'", ensg) else sprintf("-F '\t%s\t'", ensg)
  read.delim(pipe(sprintf("(head -1 '%s'; grep %s '%s')", f, pat, f)), check.names = FALSE, quote = "\"", encoding = "UTF-8")
}
table_file <- function(name) read.delim(file.path(skill_dir, name), encoding = "UTF-8")
# Gene pages embed chart data as JSON: $('#id').barChart([...]) / $("#id").bubble_plot([...]).
# ponytail: HTML scraping, breaks if HPA redesigns; only --source=web uses it.
fetch <- function(u) paste(readLines(url(u), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
chart <- function(html, id) {
  m <- regmatches(html, regexec(sprintf("(?s)\\$\\(['\"]#%s['\"]\\)\\.\\w+\\((\\[.*?\\])\\s*[,)]", id), html, perl = TRUE))[[1]]
  if (length(m) < 2) return(NULL)
  fromJSON(m[2])
}

# ---------- plots
# HPA-style vertical barplot. d: label, value, group, color, already in display order.
bars <- function(d, unit, title, subtitle = NULL, show_x = TRUE) {
  d$label <- factor(d$label, levels = unique(d$label))
  d$group <- factor(d$group, levels = unique(d$group))
  ggplot(d, aes(label, value, fill = group)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = setNames(d$color, d$group)[!duplicated(d$group)], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = unit, title = title, subtitle = subtitle) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 1),
          legend.position = "right", legend.key.size = unit(0.35, "cm"),
          plot.title = element_text(face = "bold"), plot.title.position = "plot")
}
# Single cell overview: group bars on top, tissue x group dots below.
# dots: tissue (factor, bottom-up HPA order), group, ncpm, colored (>= 1 nCPM), color.
overview <- function(b, dots, title, subtitle) {
  dots$x <- factor(dots$group, levels = b$label)
  top <- bars(b, "nCPM", title, subtitle, show_x = FALSE)
  bottom <- ggplot(dots, aes(x, tissue)) +
    geom_point(data = dots[!dots$colored, ], colour = "#C8C8C8", size = 1.6) +
    geom_point(data = dots[dots$colored, ], aes(colour = color, size = ncpm)) +
    scale_colour_identity() + scale_size(range = c(1.8, 7), guide = "none") +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +  # two layers: keep HPA order, not layer-union order
    labs(x = NULL, y = NULL, caption = paste("Dot: colored = detected (>= 1 nCPM), size ~ nCPM; grey = present, < 1 nCPM; * single nuclei.", CAPTION)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.3))
  top / bottom + plot_layout(heights = c(1, 1.7))
}

# ---------- gene
genes <- read.delim(data_file("rna_tissue_consensus"), check.names = FALSE, quote = "", encoding = "UTF-8")
g <- genes[toupper(genes$`Gene name`) == toupper(opt$gene) | genes$Gene == opt$gene, ]
if (nrow(g) == 0) stop("Gene '", opt$gene, "' not found in HPA")
gene <- g$`Gene name`[1]; ensg <- g$Gene[1]
dots <- NULL

if (opt$dataset %in% c("brain_sn", "tissue")) {
  src <- list(brain_sn = c("rna_single_nuclei_cluster_type", "Cluster type", "nCPM", "Brain single nuclei RNA-seq"),
              tissue   = c("rna_tissue_consensus", "Tissue", "nTPM", "Tissue RNA (consensus)"))[[opt$dataset]]
  d <- if (opt$dataset == "tissue") g else grep_rows(data_file(src[1]), ensg)
  pal <- table_file("hpa_palette.tsv"); pal <- pal[pal$dataset == opt$dataset, ]
  m <- match(tolower(d[[src[2]]]), pal$key)
  if (anyNA(m)) stop("Not in hpa_palette.tsv (rerun build_tables.R?): ", paste(d[[src[2]]][is.na(m)], collapse = ", "))
  d <- data.frame(label = tools::toTitleCase(d[[src[2]]]), value = d[[src[3]]], group = pal$group[m], color = pal$color[m], pos = pal$pos[m])
  d <- d[order(d$pos), ]
  p <- bars(d, src[3], paste0(gene, ": ", src[4])) + labs(caption = CAPTION)
  w <- max(6, 3 + 0.2 * nrow(d)); h <- 5

} else if (opt$source == "local") {
  grp <- table_file("sc_groups.tsv")
  cl <- grep_rows(data_file("rna_single_cell_cluster"), ensg)
  names(cl) <- c("gene", "name", "tissue", "cluster", "cell_type", "reads", "ncpm")
  cls <- table_file("sc_clusters.tsv")
  cl <- merge(cl, cls[cls$included, c("web_tissue", "tissue", "cluster", "cells")], by = c("tissue", "cluster"))
  # HPA's per-tissue cell type value = cell-count-weighted mean of nCPM over that cell type's clusters.
  wmean <- function(by) { a <- aggregate(as.formula(paste("cbind(w = ncpm * cells, cells) ~", by)), cl, sum); a$v <- a$w / a$cells; a }

  if (is.na(opt$tissue)) {
    sct <- grep_rows(data_file("rna_single_cell_type"), ensg)
    v <- sct$nCPM[match(grp$fine, sct$`Cell type`)]
    gv <- tapply(v, grp$group, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))  # genes can lack some cell types
    d <- unique(grp[c("group", "pos", "legend", "color")])
    d <- data.frame(label = d$group, value = as.numeric(gv[d$group]), group = d$legend, color = d$color, pos = d$pos)
    d <- d[!is.na(d$value), ]; d <- d[order(d$pos), ]
    # Dots: eye and retina stay separate here (as on HPA's dot tooltips); group value = max over its cell types.
    f <- wmean("tissue + web_tissue + cell_type")
    f$group <- grp$group[match(f$cell_type, grp$fine)]
    dots <- aggregate(v ~ web_tissue + group, f, max)
    tis <- table_file("sc_tissues.tsv")
    dots <- data.frame(tissue = factor(tis$label[match(dots$web_tissue, tis$web_tissue)], levels = tis$label),
                       group = dots$group, ncpm = round(dots$v, 1), colored = round(dots$v, 1) >= 1,
                       color = d$color[match(dots$group, d$label)])
    pa <- grep_rows(data_file("proteinatlas", "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"), ensg, col = 3)
    sp <- pa$`RNA single cell type group specificity`[1]
    top <- strsplit(strsplit(pa$`RNA single cell type group specific nCPM`[1], ";", fixed = TRUE)[[1]], ": ", fixed = TRUE)
    top <- if (length(top)) sapply(top, `[`, 1)[order(-as.numeric(sapply(top, `[`, 2)))] else character()
    top <- d$label[match(tolower(top), tolower(d$label))]
    spec <- if (is.na(sp) || sp == "") NA else paste0(sp, if (length(top)) paste0(" (", paste(top, collapse = ", "), ")"))
  } else {
    web <- tolower(trimws(opt$tissue))
    f <- wmean("web_tissue + cell_type")  # tissue page merges eye + retina per cell type
    f <- f[f$web_tissue == web, ]
    if (nrow(f) == 0) stop("No single cell data for tissue '", opt$tissue, "'. HPA tissues: ", paste(sort(unique(cls$web_tissue)), collapse = ", "))
    m <- match(f$cell_type, grp$fine)
    d <- data.frame(label = grp$label[m], value = round(f$v, 1), group = grp$legend[m], color = grp$color[m], pos = grp$pos[m] * 100 + grp$fine_pos[m])
    d <- d[order(d$pos), ]
  }

} else {  # --source=web
  base <- sprintf("https://www.proteinatlas.org/%s/single+cell", ensg)
  if (is.na(opt$tissue)) {
    html <- fetch(base)
    b <- chart(html, "ctChart_110"); dw <- chart(html, "cell_tissue_plot_110")
    if (is.null(b) || is.null(dw)) stop("HPA single cell charts not found for ", gene, "; page layout may have changed")
    spec <- trimws(regmatches(html, regexec("RNA single cell type group specificity:</b>&nbsp;([^<]+)<", html))[[1]][2])
    d <- data.frame(label = b$label, value = as.numeric(b$value), group = b$legend, color = b$color, pos = b$pos)
    d <- d[order(d$pos), ]
    # value: -1 = absent, 0 = present but < 1 nCPM (grey), > 0 = detected (proportional to nCPM). Tooltips list nCPM per cell type.
    lv <- unique(dw$cell_type)
    dw <- dw[dw$value >= 0, ]
    ncpm <- sapply(regmatches(dw$tooltip, gregexpr("[0-9.]+(?= nCPM)", dw$tooltip, perl = TRUE)), function(x) max(as.numeric(x), 0))
    dots <- data.frame(tissue = factor(dw$cell_type, levels = lv),
                       group = d$label[match(tolower(dw$tissue), tolower(d$label))],  # case differs: NK-cells vs Nk-cells
                       ncpm = ncpm, colored = dw$value > 0)
    dots$color <- d$color[match(dots$group, d$label)]
  } else {
    slug <- gsub(" ", "+", tolower(trimws(opt$tissue)))
    x <- chart(tryCatch(fetch(paste0(base, "/", slug)), error = function(e) ""), "sc_detail_summary")
    if (is.null(x)) stop("No HPA single cell chart for tissue '", opt$tissue, "'")
    d <- data.frame(label = x$label, value = as.numeric(x$value), group = x$legend, color = x$color, pos = x$pos)
    d <- d[order(d$pos), ]
  }
}

if (opt$dataset == "single_cell" && is.na(opt$tissue)) {
  p <- overview(d, dots, paste0(gene, ": Single cell type groups"), if (!is.na(spec)) paste("RNA single cell type group specificity:", spec))
  w <- 15; h <- 10
} else if (opt$dataset == "single_cell") {
  p <- bars(d, "nCPM", paste0(gene, ": Single cell types in ", tools::toTitleCase(tolower(trimws(opt$tissue))))) + labs(caption = CAPTION)
  w <- max(6, 3 + 0.3 * nrow(d)); h <- 5
}

dir.create("figures", showWarnings = FALSE)
out <- file.path("figures", paste(c(gene, opt$dataset, if (!is.na(opt$tissue)) gsub("\\W+", "_", tolower(opt$tissue)),
                                    if (opt$source == "web") "web"), collapse = "_"))
ggsave(paste0(out, ".pdf"), p, width = w, height = h)
ggsave(paste0(out, ".png"), p, width = w, height = h, dpi = 150)
write.table(d[c("label", "value", "group")], paste0(out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(dots)) write.table(dots[c("tissue", "group", "ncpm", "colored")], paste0(out, "_dots.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(head(d[order(-d$value), c("label", "value", "group")], 5), stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:", paste0(out, ".{pdf,png,tsv}"), "\n")

#!/usr/bin/env Rscript
# One-off: scrape the gene-independent display metadata that HPA's bulk TSVs lack. Rerun only after a new HPA release.
#   hpa_palette.tsv : brain_sn cluster types and tissues -> HPA order, group, color
#   sc_groups.tsv   : fine single cell type -> 53 display groups, HPA order, colors, display label
#   sc_clusters.tsv : per tissue cluster cell counts + "included in category" (for cell-weighted means)
#   sc_tissues.tsv  : single cell tissues in HPA dot-plot order (* = single nuclei)
# ponytail: HTML scraping of one reference gene's pages; gene expression itself always comes from the TSVs.
# Usage (any cwd): Rscript .claude/skills/hpa-gene-expression/scripts/build_tables.R

suppressPackageStartupMessages(library(jsonlite))

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(file.path(script_dir, "..", "..", "..", ".."))
out <- function(df, name) {
  write.table(df, file.path(script_dir, "..", name), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Wrote", name, ":", nrow(df), "rows\n")
}

ref <- "https://www.proteinatlas.org/ENSG00000110079-MS4A4A"  # any gene works; MS4A4A is the reference
fetch <- function(u) paste(readLines(url(u), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
chart <- function(html, id, required = TRUE) {
  m <- regmatches(html, regexec(sprintf("(?s)\\$\\(['\"]#%s['\"]\\)\\.\\w+\\((\\[.*?\\])\\s*[,)]", id), html, perl = TRUE))[[1]]
  if (length(m) < 2) { if (required) stop("Chart '", id, "' not found; HPA page layout may have changed"); return(NULL) }
  fromJSON(m[2])
}
get_tsv <- function(name) {
  f <- file.path("data", paste0(name, ".tsv"))
  if (!file.exists(f)) {
    dir.create("data", showWarnings = FALSE); z <- tempfile(fileext = ".zip")
    download.file(sprintf("https://www.proteinatlas.org/download/tsv/%s.tsv.zip", name), z, mode = "wb", quiet = TRUE)
    unzip(z, exdir = "data")
  }
  f
}

# --- brain_sn + tissue palette
main <- fetch(paste0(ref, "/single+cell"))
sn <- chart(main, "ctChart_103")
tis <- chart(fetch(paste0(ref, "/tissue")), "RNAChart_52")
out(rbind(
  data.frame(dataset = "brain_sn", key = tolower(sn$label), pos = sn$pos, group = sn$legend, color = sn$color),
  data.frame(dataset = "tissue", key = tolower(tis$label), pos = tis$pos,
             group = sub("(?s).*Organ: ([^<]+).*", "\\1", tis$tooltip, perl = TRUE), color = tis$color)),
  "hpa_palette.tsv")

# --- 53 groups: order/colors from the bar chart, members (in display order) from the hidden per-group popups.
bars <- chart(main, "ctChart_110")
js <- gsub("\\/", "/", gsub('\\"', '"', main, fixed = TRUE), fixed = TRUE)
pops <- strsplit(js, 'class="ct_popup rounded"', fixed = TRUE)[[1]][-1]
pops <- pops[grepl('^\\s*id="ct_[^"]*_110"', pops)]
groups <- do.call(rbind, lapply(pops, function(p) {
  hdr <- strsplit(p, "</tr>", fixed = TRUE)[[1]][1]
  fine <- regmatches(hdr, gregexpr("(?<=<td >)[^<]*?(?=&nbsp</td>)", hdr, perl = TRUE))[[1]]
  data.frame(fine = tolower(fine), label = fine, group = regmatches(p, regexec("<h4>(.*?)&nbsp;-&nbsp;", p))[[1]][2],
             fine_pos = seq_along(fine))
}))
m <- match(tolower(groups$group), tolower(bars$label))
if (anyNA(m)) stop("Popup groups missing from bar chart: ", paste(unique(groups$group[is.na(m)]), collapse = ", "))
groups <- cbind(groups[c("fine", "label", "fine_pos")], group = bars$label[m], pos = bars$pos[m], legend = bars$legend[m], color = bars$color[m])

# --- dot-plot tissue order (HPA lists rows bottom-up)
dots <- chart(main, "cell_tissue_plot_110")
tl <- unique(dots$cell_type)
out(data.frame(web_tissue = tolower(sub("*", "", tl, fixed = TRUE)), label = tl, pos = seq_along(tl)), "sc_tissues.tsv")

# --- cluster cell counts from each tissue subpage; also collect display labels (e.g. "NK-cells").
clu <- read.delim(pipe(sprintf("grep -F '\tMS4A4A\t' %s", get_tsv("rna_single_cell_cluster"))), header = FALSE, quote = "",
                  col.names = c("gene", "name", "tissue", "cluster", "cell_type", "reads", "ncpm"))
slugs <- unique(regmatches(main, gregexpr("(?<=-MS4A4A/single\\+cell/)[a-z+-]+(?=\")", main, perl = TRUE))[[1]])
tsv_tissues <- function(web) switch(web, brain = c("brain neurons", "brain non-neurons"), eye = c("eye", "retina"), web)
labels <- character()
cl <- do.call(rbind, lapply(slugs, function(s) {
  web <- gsub("+", " ", s, fixed = TRUE)
  cand <- clu[clu$tissue %in% tsv_tissues(web), ]
  if (nrow(cand) == 0) { message("skip ", web, ": not in cluster TSV"); return(NULL) }
  h <- fetch(paste0(ref, "/single+cell/", s))
  lab <- chart(h, "sc_detail_summary", required = FALSE)$label
  labels <<- c(labels, setNames(lab, tolower(lab)))
  r <- regmatches(h, gregexpr('<tr class="group\\d+ (c-\\d+)" title="[^"]*"[^>]*data-cell_type="[^"]*"', h))[[1]]
  x <- unique(data.frame(
    cluster   = sub('.*class="group\\d+ (c-\\d+)".*', "\\1", r),
    cell_type = tolower(sub('.*data-cell_type="([^"]*)".*', "\\1", r)),
    cells     = as.integer(sub(".*Cell count: (\\d+).*", "\\1", r)),
    reads     = sub(".*Read count: (\\d+).*", "\\1", r),
    included  = grepl("Included in category: Yes", r, fixed = TRUE)))
  # Page rows -> TSV tissue via (cluster, cell type, reference-gene read count); read count separates eye vs retina c-5.
  # Unmatched page rows belong to other datasets on the page (e.g. Tabula Sapiens) and are dropped.
  key <- paste(cand$cluster, tolower(cand$cell_type), cand$reads)
  x$tissue <- vapply(paste(x$cluster, x$cell_type, x$reads), function(k) {
    t <- unique(cand$tissue[key == k]); if (length(t) == 1) t else NA_character_ }, "")
  miss <- setdiff(key, paste(x$cluster, x$cell_type, x$reads))
  if (length(miss)) stop(web, ": TSV clusters not found on page: ", paste(head(miss, 5), collapse = "; "))
  cbind(web_tissue = web, x[!is.na(x$tissue), c("tissue", "cluster", "cell_type", "cells", "included")])
}))
if (nrow(cl) != nrow(clu)) stop("Matched ", nrow(cl), " clusters but TSV has ", nrow(clu))
out(cl, "sc_clusters.tsv")

l <- labels[groups$fine]; groups$label <- ifelse(is.na(l), groups$label, l)
out(groups, "sc_groups.tsv")

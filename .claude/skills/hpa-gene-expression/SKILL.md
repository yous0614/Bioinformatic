---
name: hpa-gene-expression
description: Query a gene's RNA expression from the Human Protein Atlas (brain single nuclei, single cell types overall or within one tissue, or tissues) and draw it in the HPA website's chart layout with R/ggplot2. Use when the user runs /hpa-gene-expression or asks which brain cell / cell type / tissue expresses a gene most according to HPA.
argument-hint: <GENE> [dataset or tissue]
---

# HPA gene expression plots

Input: `/hpa-gene-expression <GENE> [more words]`, e.g. `MS4A4A`, `MS4A4A brain single nuclei`, `MS4A4A liver cells`.

## 1. Walk the option tree

No gene given: ask for it first (plain text). Then walk the tree below **one level per AskUserQuestion call**, in the user's language, until a leaf (= one script call) is reached.
Skip any level the user's words already settle (e.g. `MS4A4A liver cells` → straight to the liver leaf; `MS4A4A tissues` → tissue leaf). If the words are ambiguous ("brain" fits both A→Brain and B), ask that level; never guess.
For a level with more than 4 choices, offer the 4 listed and name the rest in the question text so the user can type one via "Other".

```
L1 Dataset
├─ A Single cell scRNA-seq (cell types across organs)
│   L2 Scope
│   ├─ Whole-body overview ................. --dataset=single_cell
│   └─ One tissue → L3 Organ system
│       ├─ Nervous / sensory / endocrine → L4: Brain | Eye | Pituitary gland | Adrenal gland
│       ├─ Digestive                     → L4: Liver | Pancreas | Stomach | Colon
│       │                                     (Other: Tongue, Salivary gland, Esophagus, Small intestine, Rectum)
│       ├─ Immune / blood                → L4: Blood | Bone marrow | Lymph node | Spleen   (Other: Thymus)
│       └─ Other organs → L4 Subsystem
│           ├─ Respiratory / urinary     → L5: Lung | Kidney | Urinary bladder
│           ├─ Reproductive              → L5: Breast | Prostate | Testis | Placenta
│           │                                  (Other: Epididymis, Ovary, Fallopian tube, Endometrium)
│           ├─ Muscle / vascular         → L5: Heart muscle | Skeletal muscle | Vasculature
│           └─ Skin / adipose            → L5: Skin | Adipose tissue
│       leaf: --dataset=single_cell --tissue=<tissue, lowercase>
├─ B Brain single nuclei snRNA-seq (34 brain cluster types) ...... --dataset=brain_sn
└─ C Tissue bulk RNA (51 tissues, colored by organ) ............. --dataset=tissue
```

Option descriptions should say what the figure shows: A-overview = 53 cell type groups + tissue dot matrix; A-tissue = that tissue's cell types; B = brain cluster types; C = tissue bars.

Not supported (do not offer, do not substitute): IHC protein, tissue MS proteomics, blood/Olink, cell lines, brain subregions, spatial. If asked, say so and stop.

## 2. Run

```bash
Rscript .claude/skills/hpa-gene-expression/scripts/hpa_plot.R --gene=<GENE> --dataset=<D> [--tissue=<T>] [--source=web]
```
- Gene symbol (case-insensitive) or Ensembl ID. Works from any cwd; writes to the project's `figures/`.
- Default `--source=local`: everything is computed from HPA bulk downloads in `data/` (fetched on first use; the single cell cluster file is 1.4 GB, first run takes a few minutes).
- `--source=web` (single_cell only) reads the chart data embedded in the gene's HPA page. Use it **only** if the user asks for it, or local fails (e.g. a missing file can't be downloaded); say which source was used.
- Output: `figures/<GENE>_<dataset>[_<tissue>][_web].{pdf,png,tsv}` (+ `_dots.tsv` for the overview); stdout prints the top 5 by value.
- Bars follow HPA's display order (not sorted by value), colors are HPA's.

## 3. Report

Read the PNG to check it rendered, then answer in 2–4 lines: top cell type/tissue with value and unit, how far ahead of #2, the file path. For the single_cell overview also quote the specificity line in the subtitle and name tissues with the largest dots in the top group.

## Reading the single_cell overview

- Top bars = max nCPM over the fine cell types in each of 53 groups (e.g. Mononuclear phagocytes = max of macrophages, Kupffer, Hofbauer, monocytes).
- Dots: colored = max nCPM of that group in that tissue >= 1, size ~ nCPM; grey = present but < 1 nCPM; none = absent. `*` = single nuclei. HPA help calls colored dots "enriched", but the data is a plain >= 1 nCPM threshold: never describe them as enrichment.
- Per-tissue values = cell-count-weighted mean of cluster nCPM. Tissue pages merge eye + retina per cell type; dot tooltips keep them separate (the script does both).
- Local vs website: bars of the overview are exact; per-tissue values can differ by 0.1 (the cluster TSV is rounded to 1 decimal); a dot right at ~1 nCPM can differ in color.

## Naming gotchas

- brain_sn has no "microglia" label: microglia = `Central Nervous System Macrophage` (HPA groups it under Glial cells).
- In single_cell, microglia (Glial cells) is separate from macrophages/monocytes (Blood and immune cells).

## Maintenance (new HPA release)

1. Delete `data/` (re-downloads), update `HPA_VERSION` in hpa_plot.R.
2. `Rscript .claude/skills/hpa-gene-expression/scripts/build_tables.R` rebuilds the gene-independent tables (`hpa_palette.tsv`, `sc_groups.tsv`, `sc_clusters.tsv`, `sc_tissues.tsv`) by scraping one reference gene's pages.
3. `Rscript .claude/skills/hpa-gene-expression/scripts/check_sources.R` (~10 min) must print OK: local vs web for 10 genes × overview + 5 tissues.

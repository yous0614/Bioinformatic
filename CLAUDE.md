# bioinformatic-tool

Toolkit of Claude Code skills/agents/harness config for fetching public omics data and plotting it (mostly R/ggplot2). See README.md for layout and how to add a tool.

## Conventions (all tools)

- R 4.6 (Homebrew). Base R + ggplot2/jsonlite/patchwork; don't add tidyverse or other packages for small things.
- Scripts run from any cwd and write to repo-root `data/` (downloads, gitignored, re-creatable) and `figures/` (outputs, gitignored).
- Never commit large or re-downloadable data. Commit only small reference tables together with the script that rebuilds them.
- Expression values come from official bulk downloads; HTML scraping only for display metadata or explicit fallbacks.
- Every figure caption cites its data source.
- Each tool keeps one runnable check, listed below.

## Tools

### hpa-gene-expression (`.claude/skills/hpa-gene-expression/`)
- Human Protein Atlas v25.1; usage and data rules in its SKILL.md, method/validation in its README.md.
- Bulk files: `https://www.proteinatlas.org/download/tsv/<name>.tsv.zip` (not listed on the download page but live); `proteinatlas.tsv.zip` is at `/download/` (no `tsv/`).
- SEA-AD disease-progression questions go to the `seaad` MCP, not HPA.
- Quick check: `Rscript .claude/skills/hpa-gene-expression/scripts/hpa_plot.R --gene=MS4A4A --dataset=brain_sn | sed -n 2p | grep -q "^Central Nervous System Macrophage" && echo OK`
- Full check (~10 min, single cell local vs website): `Rscript .claude/skills/hpa-gene-expression/scripts/check_sources.R`

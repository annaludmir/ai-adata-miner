# Shared context for generated analyses (step 3)

Prepend this to any prompt that asks a model to write a downstream analysis.

## Your inputs

You work **only** on CSVs under `csv_exports/`. You must not open the `.h5ad`
files — they are 300k and 1.67M cells and the CSV layer exists precisely so you
do not have to.

Discover what is available by reading `csv_exports/<dataset>/_manifest.csv`.
It has one row per CSV with `csv`, `script`, `n_rows`, `n_cols`, `columns` and a
`description`. **Read the manifest before assuming a filename exists.** If the
table you need is not in the manifest, say so rather than inventing it.

Datasets: `cortex` (cortex/forebrain, 297,927 cells) and `human_dev`
(whole brain, 1,665,937 cells). `csv_exports/_cross_dataset/` holds the keys
that make comparing them valid.

## Non-negotiable statistical rules

1. **Donors are the replicate unit, not cells.** There are 15 donors in `cortex`
   and 26 in `human_dev`. Never compute a p-value over ~10^6 cells. Aggregate to
   donor or sample first (`sample__*` and `donor__*` pseudobulk tables exist for
   exactly this).
2. **Age is nested in chemistry.** `05_confounds/confound_warnings.csv` quantifies
   it. Any developmental claim must state this limitation or control for it.
3. **Composition data are compositional.** Use `log2oe_*` tables for enrichment,
   not raw `fractions_*`, which cannot rise independently.
4. **Pick the right expression table.**
   - `*__pseudobulk_counts.csv` → counts, for edgeR/DESeq2 only
   - `*__mean_lognorm.csv` → comparing groups of different depth
   - `*__detection_fraction.csv` → how *broadly* a gene is on
   - `*__cpm.csv` → depth-independent magnitude
   Using counts where lognorm is meant produces a depth artefact.
   Before comparing groups, check `14_normalization/normalization_summary.csv`:
   a wide TMM spread or a flagged group means CPM alone is unsafe there, and the
   TMM factor must be applied. Never compute RPKM/FPKM — this is 10x 3' UMI data
   with no transcript-length bias to correct.
5. **Structure before significance.** `15_relationships/` already holds the
   correlation matrices, dendrograms and PCA. If groups cluster by donor or
   chemistry rather than by biology, say so before reporting any DE result.
6. **Absent ≠ zero.** Check `11_panels/panel_coverage.csv` before concluding a
   gene is not expressed.
7. **Small groups are flagged, not dropped.** Respect `below_min_cells` and
   `n_cells` columns; a group of 7 cells should not drive a conclusion.
8. **Gene panels are seed lists**, not authoritative releases (`panels/README.md`).

## Output contract

- Write scripts to `downstream_analyses/`, results to `downstream_analyses/results/`.
- Use only numpy / pandas / scipy / matplotlib unless told otherwise.
- Every script states its inputs at the top and fails loudly if they are missing.
- Report effect sizes with uncertainty, and always report the n *of donors*.
- Describe what the result does **not** support as readily as what it does.

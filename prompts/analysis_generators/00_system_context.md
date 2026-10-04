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

Before writing anything, read `downstream_analyses/REPORT.md` and
`downstream_analyses/results/01_data_audit/SUMMARY.md`: they record what the
data can support. Reuse `downstream_analyses/_common.py` (loading, TMM log
CPM, exact permutation tests, signed Stouffer replication, SUMMARY.md writer)
rather than re-implementing it.

## Non-negotiable statistical rules

1. **Donors are the replicate unit, not cells.** There are 15 donors in `cortex`
   and 26 in `human_dev`. Never compute a p-value over ~10^6 cells. Aggregate to
   donor or sample first (`sample__*` and `donor__*` pseudobulk tables exist for
   exactly this).
2. **Chemistry is a confound with age, and the CSVs are stratified for it.**
   `csv_exports/<dataset>__v2/` and `__v3/` hold the per-chemistry results;
   `csv_exports/<dataset>/` holds only the pooled-by-design steps (00, 07, 13, 18).
   Work *within* one chemistry folder. Never compare a group from `__v2/`
   against one from `__v3/` and report the difference as biological — v2 and v3
   cover near-disjoint age ranges. If a contrast must span chemistries, restrict
   to the ages that `18_chemistry/age_chemistry_overlap.csv` marks `comparable`,
   and say how many cells that leaves. Donors are nested in chemistry in both
   files, so v2 and v3 are **independent donor sets: use them to replicate each
   other** (same direction in both, then combine), never pool them.
   Within a chemistry, each donor has exactly one age and most ages have one
   donor, so an age trend is a trend across 5-15 people.
3. **cortex is not an independent cohort.** All 15 cortex donors are human_dev
   donors (IDs differ: cortex `XHU:1966:307` = human_dev `XHU:307`; use
   `_common.normalise_donor`). Agreement between the files is reproducibility,
   not replication. human_dev cell classes pool brain regions whose dissection
   changes with age, so its within-class age trends are region-confounded;
   cortex is the cleaner test.
4. **Donor sex is not annotated** (`obs['sex']` is "unknown") but is inferred in
   `results/01_data_audit/donor_sex_inferred.csv`, and is unevenly spread over
   age. Flag chrY genes and XIST/TSIX; check any age result against sex.
5. **Composition data are compositional.** Use `log2oe_*` tables for enrichment,
   not raw `fractions_*`, which cannot rise independently.
6. **Pick the right expression table.**
   - `*__pseudobulk_counts.csv` → the input for any comparison *between groups*:
     TMM-normalise to log CPM (`_common.tmm_log_cpm`) or use edgeR/DESeq2
   - `*__mean_lognorm.csv` → describing a group (mean log1p CP10K per cell).
     **Not for comparing groups that differ in depth**: it drifts with UMIs per
     cell. In cortex, UMIs per cell fall with age, and mean_lognorm made the
     median gene "decrease with age" (rho ≈ -0.5); TMM log CPM removed it.
     Always check the median gene-level statistic is near 0.
   - `*__detection_fraction.csv` → how *broadly* a gene is on
   - `*__cpm.csv` → plain CPM; fails when a few genes dominate (erythrocytes)
   Exports made before the empty-group fix contain all-zero columns for levels
   that exist only in the other chemistry: drop groups by `n_cells` from
   `*__group_summary.csv` (`_common.group_matrix` does this).
   Before comparing groups, check `14_normalization/normalization_summary.csv`:
   a wide TMM spread or a flagged group means CPM alone is unsafe there, and the
   TMM factor must be applied. Never compute RPKM/FPKM — this is 10x 3' UMI data
   with no transcript-length bias to correct.
7. **Structure before significance.** `15_relationships/` already holds the
   correlation matrices, dendrograms and PCA. If groups cluster by donor or
   chemistry rather than by biology, say so before reporting any DE result.
8. **Match the null on expression level** for any gene-set test. NDD genes are
   long, neuronal and highly expressed; highly expressed genes are also measured
   more precisely, so they pass significance tests more often. An unmatched
   null calls both of those things "enrichment".
9. **Absent ≠ zero.** Check `11_panels/panel_coverage.csv` before concluding a
   gene is not expressed.
10. **Small groups are flagged, not dropped.** Respect `below_min_cells` and
   `n_cells` columns; a group of 7 cells should not drive a conclusion.
11. **Gene panels are seed lists**, not authoritative releases (`panels/README.md`).

## Output contract

- Write scripts to `downstream_analyses/<nn>_<slug>.py`, results to
  `downstream_analyses/results/<nn>_<slug>/` (`_common.Output` does this), then
  run `downstream_analyses/build_report.py` so REPORT.md includes the new one.
- Use only numpy / pandas / scipy / matplotlib unless told otherwise.
- Every script states its inputs at the top and fails loudly if they are missing.
- Report effect sizes with uncertainty, and always report the n *of donors*.
- Describe what the result does **not** support as readily as what it does.

# Step 3: downstream analyses

Analyses that read only `csv_exports/`, never an `.h5ad`, so they run on a
laptop. Each one writes a folder under `results/` holding a `SUMMARY.md`,
tables and figures. `build_report.py` collects the summaries into
[`REPORT.md`](REPORT.md).

```bash
./downstream_analyses/run_all.sh                 # everything, about a minute
PYTHON=~/envs/sc/bin/python ./downstream_analyses/run_all.sh
python3 downstream_analyses/03_age_trends_within_cell_class.py   # one analysis
```

Needs numpy, pandas and scipy. matplotlib is optional; without it the figures
are skipped. To read exports from somewhere else, set `AI_ADATA_OUT_ROOT`, and
to write results elsewhere, set `AIM_DOWNSTREAM_OUT`.

## The analyses

| # | Question | Replicate unit | Main output |
|---|---|---|---|
| 01 | What can the data support? Replicate structure, donor matching across files, inferred donor sex, marker sanity, export artefacts | donor | `replicate_structure.csv`, `donor_sex_inferred.csv` |
| 02 | Which cell classes expand or shrink with age, in both donor sets? | donor | `composition_trends_replicated.csv` + per-scope figures |
| 03 | Which genes rise or fall with age *within* a cell class, in both donor sets? | age point (≈ donor) | `age_trends_combined.csv`, `ndd_genes_with_age_trends.csv` |
| 04 | Where are NDD panels expressed, against expression-matched random genes, and are they enriched among age trends? | age point (≈ donor) | `panel_class_preference_combined.csv`, `panel_age_trend_enrichment.csv` |
| 05 | Does each cell class mean the same thing in v2 and v3, and in cortex and human_dev? | class profile | `identity_reproducibility.csv` |

Order matters: 02 reads the sex calls from 01, and 04 reads the trends from 03.
`run_all.sh` runs them in number order.

## Rules every analysis here follows

These come from the audit (01), not from habit:

- **Strata stay separate.** v2 and v3 have disjoint donors and mostly disjoint
  ages. Each is analysed on its own, and the two serve as each other's replicate.
- **Tiers:**
  - *replicated*: same direction in both chemistries, each nominally
    significant on its own (one-sided p < 0.05), and combined
    (signed Stouffer) BH q < 0.05;
  - *supported*: passes the combined test, but one chemistry is weak.
- **Exact permutation p-values.** With 5–9 points the t-approximation behind
  `scipy.stats.spearmanr` is unreliable, so `_common.spearman_perm_p` enumerates
  every ordering (Monte Carlo above 9! orderings).
- **Comparisons between groups use TMM log CPM from pseudobulk counts**, not
  `mean_lognorm`. The latter drifts with UMIs per cell, and in cortex those fall
  with age. 03 records the median gene-level rho per class as a drift check.
- **Gene-set nulls are matched on expression level.** NDD genes are highly
  expressed, and highly expressed genes pass tests more often.
- **cortex is primary for within-class questions.** human_dev classes pool
  regions whose sampling changes with age.

## Adding an analysis

1. Copy the shape of an existing script: a module docstring stating the
   question, method and inputs, then `out = C.Output("<nn>_<slug>")`, then
   `out.used(...)` for every input, `out.write(...)` for every table, and
   `out.summary(...)` at the end.
2. Load tables through `_common` (`group_matrix`, `tmm_log_cpm`, `donor_ages`,
   `panels`), which fail loudly when an input is missing.
3. Compute findings from the results instead of writing them by hand, so
   `SUMMARY.md` stays true when the data change.
4. Run it, then `build_report.py`.

The generator prompt in `prompts/analysis_generators/` encodes the same rules,
so an AI-generated analysis follows them too.

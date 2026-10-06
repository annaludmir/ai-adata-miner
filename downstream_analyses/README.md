# Step 3: downstream analyses

Analyses that read only `csv_exports/`, never an `.h5ad`, so they run on a
laptop. Each one writes a folder under `results/` holding a `SUMMARY.md`,
tables and figures. `build_report.py` collects the summaries into
[`REPORT.md`](REPORT.md).

```bash
./downstream_analyses/run_all.sh                 # everything, a few minutes
PYTHON=~/envs/sc/bin/python ./downstream_analyses/run_all.sh
python3 downstream_analyses/03_age_trends_within_cell_class.py   # one analysis
```

Needs numpy, pandas and scipy. matplotlib is optional; without it the figures
are skipped. To read exports from somewhere else, set `AI_ADATA_OUT_ROOT`, and
to write results elsewhere, set `AIM_DOWNSTREAM_OUT`.

Every analysis first checks that numpy computes correctly, and stops with
instructions if it does not. numpy 2.2.6 on macOS 26 (the Mac's system Python)
returned rank correlations above 1 on large arrays; numpy 2.5 and the cluster's
numpy are fine. Correlation helpers also spot-check a few values against a
pure-Python recomputation on every call. If the check stops a run, use a fresh
environment:

```bash
python3 -m venv ~/aim-env && ~/aim-env/bin/pip install -U numpy scipy pandas matplotlib
PYTHON=~/aim-env/bin/python ./downstream_analyses/run_all.sh
```

## The analyses

| # | Question | Replicate unit | Main output |
|---|---|---|---|
| 01 | What can the data support? Replicate structure, donor matching across files, inferred donor sex, marker sanity, export artefacts | donor | `replicate_structure.csv`, `donor_sex_inferred.csv` |
| 02 | Which cell classes expand or shrink with age, in both donor sets? | donor | `composition_trends_replicated.csv` + per-scope figures |
| 03 | Which genes rise or fall with age *within* a cell class, in both donor sets? | age point (≈ donor) | `age_trends_combined.csv`, `ndd_genes_with_age_trends.csv` |
| 04 | Where are NDD panels expressed, against expression-matched random genes, and are they enriched among age trends? | age point (≈ donor) | `panel_class_preference_combined.csv`, `panel_age_trend_enrichment.csv` |
| 05 | Does each cell class mean the same thing in v2 and v3, and in cortex and human_dev? | class profile | `identity_reproducibility.csv` |
| 06 | Gene lists: what the data see of each list, how the lists overlap, where each is expressed, and whether a list moves with age as a group | age point (≈ donor) | `list_coverage.csv`, `list_class_preference_combined.csv`, `list_age_coordination_combined.csv` |
| 07 | Do a list's genes co-express (across clusters, and within cell classes)? Which members carry that, how does the list split into sub-modules, which outside genes track it? | cluster profile; sub-modules cross-validated between donor sets | `coherence_combined.csv`, `submodules.csv`, `candidate_members.csv` |
| 08 | Data-driven co-expression modules: which gene groups exist without any list, where they peak, how they change with age, which lists they hold; each module labelled from the seed panels it concentrates (G2/M phase, S phase, cell-class markers, patterning) | cluster profile; modules rediscovered in each chemistry | `modules.csv` (`annotation` column), `module_annotation.csv`, `module_gene_set_enrichment.csv` |
| 09 | Cell-cycle programs: proliferation and G1/S/G2M shares over development within progenitor types, raw and compared at matched sequencing depth; outer vs ventricular radial glia (script 19): oRG share over age, each sub-type's phase shares, and the radial-glia trend at a fixed sub-type mix; which genes, lists and modules follow proliferation and lean to S or G2/M; whether list age trends in progenitors survive once proliferation-linked genes are set aside | age point (≈ donor); cluster profile | `proliferation_trajectories_combined.csv`, `gene_phase_map.csv`, `set_phase_profile.csv`, `age_trends_cycle_independent.csv` |
| 10 | Splicing dynamics (cortex): in which cell classes genes, lists and modules are being switched on or off (nascent vs mature RNA), and whether their unspliced share rises or falls with age alongside expression | age point (≈ donor); cell-class pseudobulk | `set_induction_by_cell_class.csv`, `unspliced_age_trends.csv`, `set_unspliced_age_trends.csv` |
| 11 | Which fine clusters (sub-types) each list concentrates in, replicated per cluster; lists along neuron sub-type axes (excitatory vs inhibitory, deep vs upper layer) | cluster profile; clusters paired across chemistries | `cluster_scores_combined.csv`, `neuron_axis_tests_combined.csv` |
| 12 | Which sub-types (clusters) expand or shrink with age within their cell class | donor | `cluster_abundance_combined.csv` |
| 13 | Regional differences of lists and modules within a cell type, with each region's age offset (human_dev) | class x region pseudobulk | `region_contrast_combined.csv` |
| 14 | Sex differences of lists within cell types, age-adjusted; Y genes / XIST as positive control (exploratory: few male donors) | age point (≈ donor) | `set_sex_differences.csv`, `gene_sex_differences.csv` |
| 15 | Between-donor variability of lists after the age trend, against expression-matched genes (dosage control) | age point (≈ donor) | `set_variability_combined.csv`, `consistently_tight_genes.csv` |
| 16 | Age clock: predict a cell type's age from expression, trained on one donor set and tested on the other; clock genes; lists that predict age better than random genes | class x age pseudobulk | `clock_accuracy.csv`, `clock_gene_weights.csv`, `list_clocks.csv` |
| 17 | human_dev's 50 latent factors: what each marks, which lists sit among its top genes, which change with age beyond cell-type mix | donor | `factor_identity.csv`, `lists_in_factors.csv`, `factor_age_trends_combined.csv` |
| 18 | Co-expression rewiring: list coherence in early vs late clusters, same class mix | cluster profile | `rewiring_combined.csv` |
| 19 | Robustness: gene-length-matched nulls for 06/07 tests; modules tracking QC metrics; dissociation-stress score vs age | as the tests checked | `length_check_summary.csv`, `module_qc.csv`, `stress_vs_age.csv` |
| 20 | Agreement of CellClass with the files' other annotations (cortex 'classes', human_dev Cell Ontology terms); needs stage 1 from this version | cells | `annotation_agreement.csv`, `annotation_mapping.csv` |

Order matters: 02 and 14 read the sex calls from 01; 04, 06, 09, 10 and 19
read the trends from 03; and 09, 10, 13, 15, 18 and 19 also use 08's modules.

06 and 07 need gene lists. They read the folder named by `config.GENE_LISTS_DIR`:
`AIM_GENE_LISTS` if set, otherwise the cluster folder when it exists, otherwise
the git-ignored `gene_lists/` in the repo. With no lists, they write a short
summary saying so. 08 runs either way and also tests the seed NDD panels.
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
- **GWAS lists count each locus once.** Lists whose name matches
  `AIM_COLLAPSE_LISTS` (default `GWAS`) keep one gene per locus, because
  neighbouring genes at one association signal are often co-regulated and
  would otherwise inflate coherence and enrichment.
- **Co-expression is measured across fine clusters, and nulls are matched on
  level and spread.** Clusters need at least 100 cells and must not be
  dominated by one donor. Random comparison sets draw each gene from the same
  bin of mean level × spread as the real gene.
- **A group that was found in the data is never tested on the data that found
  it.** Sub-modules (07) and modules (08) are discovered in one chemistry and
  tested in the other.
- **"The same cluster" in v2 and v3 must be checked, not assumed.** cortex's
  fine clusterings (`Clusters`, `ClustersModularity`, `ClustersSurprise`) were
  computed per chemistry and reuse labels for different cells; `leiden_scVI`,
  `louvain` and human_dev's `cluster_id` were computed once.
  `_common.cluster_label_identity` tests this (a label's v2 profile must best
  match the same label in v3); 11 pairs cortex clusters by expression instead,
  and 12 uses `leiden_scVI`.
- **A null must vary what the question varies.** For sex differences (14),
  random genes are the wrong comparison: donors differ in ways that are not
  sex, so sex labels are shuffled among donors too, and a set must beat both.
- **Sub-modules are compared with what the same procedure finds in random
  lists.** During development, a random list split into sub-modules that also
  held up in the other donor set, because the transcriptome has strong,
  reproducible structure. So 07 asks whether a list's sub-modules are tighter
  than the sub-modules carved out of matched random lists.

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

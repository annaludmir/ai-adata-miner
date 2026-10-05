# Data audit: what the exports can support

_Generated 2026-10-05 10:57 UTC by `downstream_analyses/01_data_audit.py` from `csv_exports/`._

## Question

Before any biology: how many independent units stand behind each comparison, which covariates are entangled, and do the labels and tables behave?

## Inputs

- `_cross_dataset/gene_id_map.csv`
- `_cross_dataset/label_overlap_donor.csv`
- `cortex/18_chemistry/age_chemistry_overlap.csv`
- `cortex/18_chemistry/chemistry_summary.csv`
- `cortex__v2/02_composition/counts_cell_class_by_donor.csv`
- `cortex__v2/05_confounds/crosstab_age_x_donor.csv`
- `cortex__v2/09_pseudobulk/donor__mean_lognorm.csv`
- `cortex__v2/10_markers/top_markers_cell_class.csv`
- `cortex__v2/11_panels/marker_label_check.csv`
- `cortex__v3/02_composition/counts_cell_class_by_donor.csv`
- `cortex__v3/05_confounds/crosstab_age_x_donor.csv`
- `cortex__v3/09_pseudobulk/donor__mean_lognorm.csv`
- `cortex__v3/10_markers/top_markers_cell_class.csv`
- `cortex__v3/11_panels/marker_label_check.csv`
- `human_dev/18_chemistry/age_chemistry_overlap.csv`
- `human_dev/18_chemistry/chemistry_summary.csv`
- `human_dev__v2/02_composition/counts_cell_class_by_donor.csv`
- `human_dev__v2/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v2/09_pseudobulk/donor__mean_lognorm.csv`
- `human_dev__v2/10_markers/top_markers_cell_class.csv`
- `human_dev__v2/11_panels/marker_label_check.csv`
- `human_dev__v3/02_composition/counts_cell_class_by_donor.csv`
- `human_dev__v3/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v3/09_pseudobulk/donor__mean_lognorm.csv`
- `human_dev__v3/10_markers/top_markers_cell_class.csv`
- `human_dev__v3/11_panels/marker_label_check.csv`

## Method

- Replicate structure from the age x donor crosstabs of each chemistry stratum.
- Cross-dataset donor matching after normalising ID formats; checked on chemistry, age and cell count (cortex cells <= human_dev cells).
- Donor sex inferred from donor-level pseudobulk: male if mean Y-gene log1p(CP10K) > 0.1 and XIST <= 0.1, female if the reverse.
- Sex-age association: Spearman between male indicator and age, exact permutation p.
- Marker sanity from 11_panels/marker_label_check.csv; sex-linked markers from the top 50 per class.

## Key findings

- **Donor is the replicate unit, and age is nested in donor.** Every donor has exactly one age (max ages per donor = 1), and 85% of age points rest on a single donor. Donors per stratum: cortex__v2 7, cortex__v3 8, human_dev__v2 15, human_dev__v3 9. An age effect is therefore indistinguishable from a donor effect at that age.
- **cortex: chemistries cover different ages.** v2 only: 7.5|8|8.5|9.2|9.5|10; v3 only: 5|5.5|11.5|12|13.25|14; both: 6.9 (10% of cells). Donors nested in chemistry: True. v2 and v3 are therefore independent donor sets that can replicate each other, but not be pooled as if equivalent.
- **human_dev: chemistries cover different ages.** v2 only: 6.6|6.7|7.5|8|8.1|9.2|9.5|10; v3 only: 7|11.5|12|13|14; both: 6|6.9|8.5 (42% of cells). Donors nested in chemistry: True. v2 and v3 are therefore independent donor sets that can replicate each other, but not be pooled as if equivalent.
- **cortex is not an independent cohort.** 13/15 cortex donors are human_dev donors (5 under a differently written ID, e.g. `XHU:1966:307` = `XHU:307`; `13_cross_dataset_keys` joins them). The other 2 (XDD:348, XDD:400) fall under a human_dev exclusion rule, so they have no human_dev counterpart in these exports. Agreement between the two files is reproducibility of processing, not replication. Age annotations disagree for XDD:359 (13.25 vs 13 pcw).
- **Donor sex is recoverable from expression** (obs['sex'] is 'unknown'): 39/39 donor entries call cleanly.
- **Sex is unevenly spread over age** in cortex__v2 (male ages 6.9|8, female ages 7.5|8.5|9.2|9.5|10; rho = -0.63, exact p = 0.19); cortex__v3 (male ages 6.9|12|14, female ages 5|5.5|6.9|11.5|13.25; rho = 0.45, exact p = 0.29). With this few donors the association is not significant, but it does not need to be to matter: a sex-differential gene can look like an age trend. Later analyses flag sex-linked genes and check trends against sex.
- **Sex-linked genes rank among cell-class markers**, which happens when a class is drawn unevenly from male and female donors: RPS4Y1 in Radial glia (cortex__v2, rank 10); RPS4Y1 in Radial glia (human_dev__v2, rank 34); EIF1AY in Erythrocyte (human_dev__v3, rank 43). See cell_class_sex_balance.csv for the imbalance behind each.
- **Cell-class labels recover known markers**: 20/30 seed marker panels score highest in the class they name. Misses: radial_glia->Glioblast (cortex__v2), neuroblast->Neuron (cortex__v3), radial_glia->Glioblast (cortex__v3), glioblast_opc->Oligo (human_dev__v2), neuron->Placodes (human_dev__v2), oligo->Neural crest (human_dev__v2), radial_glia->Glioblast (human_dev__v2), glioblast_opc->Oligo (human_dev__v3), oligo->Neural crest (human_dev__v3), radial_glia->Glioblast (human_dev__v3). Panels naming a class absent from a stratum are not scored. The misses are neighbouring lineages: glioblasts carry radial-glia genes, placode-derived sensory neurons carry pan-neuronal genes, and the short seed lists cannot separate them.
- **Exclusions (exclusions.csv)**: human_dev age=5.0 (The 'without_week_5' file still holds cells at 5.0 pcw (donor XDD:348; 62,786 cells)); human_dev age=5.5 (The 'without_week_5' file still holds cells at 5.5 pcw (donor XDD:400; 59,667 cells)). These exports already omit them.
- **XIST reads ~6.1x higher in human_dev than in cortex for the same female donors**, so the two files probably count reads differently (XIST is largely nuclear and intronic). Compare genes across the files by rank or within-file contrast, not by absolute level.

## Limitations

- Sex calls rest on expression and need confirming against sample records.
- Per-donor composition reflects which tissue was dissected from that donor, so a class imbalance between sexes may be dissection, not sex biology.
- Seed marker panels are short, hand-picked lists (panels/README.md).

## What would strengthen this

- Fill obs['sex'] from records and re-export, so sex can enter models as a covariate.
- Fix the donor-ID join in 13_cross_dataset_keys (strip the middle field).
- Re-run stage 2 so the pseudobulk exports carry no empty groups.

## Output files

- `replicate_structure.csv` -- Per chemistry stratum: donors, ages, donors per age
- `chemistry_age_overlap.csv` -- Which ages each chemistry covers
- `donor_map_cortex_to_human_dev.csv` -- Each cortex donor matched to its human_dev ID (IDs are written differently)
- `donor_sex_inferred.csv` -- Donor sex from Y-gene vs XIST expression (obs['sex'] is 'unknown')
- `sex_age_entanglement.csv` -- Is donor sex associated with age within a stratum? (exact permutation p)
- `cell_class_sex_balance.csv` -- Fraction of each cell class contributed by male donors vs the stratum overall
- `sex_linked_genes_in_top_markers.csv` -- Y-linked genes and XIST among the top 50 markers of a cell class
- `marker_panel_recovery.csv` -- Does each seed marker panel score highest in the cell class it names?
- `pseudobulk_group_sizes.csv` -- Groups per pseudobulk table; empty = level that exists only in the other chemistry
- `exclusions_status.csv` -- Rules in exclusions.csv: cells covered and whether exports already omit them
- `xist_cortex_vs_human_dev.csv` -- XIST level for the same female donors in each file

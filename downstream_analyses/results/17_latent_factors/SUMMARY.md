# Latent factors (human_dev): identity, gene lists, and age trends within cell types

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/17_latent_factors.py` from `csv_exports/`._

## Question

Which cell types and regions does each latent factor mark, which gene lists concentrate among its top genes, and which factors change with age beyond the shifting cell-type mix?

## Inputs

- `_cross_dataset/gene_id_map.csv`
- `human_dev/07_modules/module_genes_long.csv`
- `human_dev__v2/02_composition/counts_cell_class_by_donor.csv`
- `human_dev__v2/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v2/08_factor_activity/factor_activity_Factors_by_cell_class.csv`
- `human_dev__v2/08_factor_activity/factor_activity_Factors_by_donor.csv`
- `human_dev__v2/08_factor_activity/factor_activity_Factors_by_region.csv`
- `human_dev__v2/08_factor_activity/factor_qc_correlation_Factors.csv`
- `human_dev__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/02_composition/counts_cell_class_by_donor.csv`
- `human_dev__v3/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v3/08_factor_activity/factor_activity_Factors_by_cell_class.csv`
- `human_dev__v3/08_factor_activity/factor_activity_Factors_by_donor.csv`
- `human_dev__v3/08_factor_activity/factor_activity_Factors_by_region.csv`
- `human_dev__v3/08_factor_activity/factor_qc_correlation_Factors.csv`
- `human_dev__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`

## Method

- A: class whose pseudobulk (log2 TMM-CPM, Z across classes) expresses a factor's top positive genes most; top major class (>= 2% of cells) by mean activity; v2/v3 agreement; 08's QC-correlation flags.
- B: each factor's 100 top positive and negative genes (07_modules); overlap with each set vs an expectation where each set gene has its total-UMI decile's rate of pole membership; exact Poisson-binomial p; BH over all tests.
- C: donor mean activity vs donor age (Spearman, exact permutation), raw and after subtracting each donor's expected activity from its class mix; v2 x v3 signed Stouffer (weights sqrt(donors)), BH, tiered.

## Key findings

- **Factor identity**: 35 of 50 factors' top genes peak in the same cell class in both chemistries (median activity class-profile r 0.94); 1 flagged as likely technical and 1 as tracking the cell cycle by 08. Class where each factor's top genes peak (v2): Erythrocyte 1, Fibroblast 6, Glioblast 2, Immune 2, Neuroblast 2, Neuron 29, Oligo 4, Radial glia 1, Vascular 3.
- **No gene set concentrates among any factor's top genes** beyond expression-matched chance.
- **Factors changing with age (raw)**: F5 (Glioblast genes) rises with age (rho +0.74/+0.85; replicated); F10 (Glioblast genes) rises with age (rho +0.62/+0.91; replicated); F23 (Neuron genes) falls with age (rho -0.75/-0.40; supported); F12 (Radial glia genes) falls with age (rho -0.63/-0.67; replicated); F6 (Neuron genes) rises with age (rho +0.59/+0.63; replicated); F28 (Neuron genes) falls with age (rho -0.51/-0.74; replicated); F17 (Oligo genes) rises with age (rho +0.21/+0.95; supported); F4 (Neuron genes) [cell cycle] rises with age (rho +0.63/+0.43; supported).
- **Factors changing with age (composition-adjusted)**: F5 (Glioblast genes) rises with age (rho +0.75/+0.83; replicated); F4 (Neuron genes) [cell cycle] rises with age (rho +0.88/+0.33; supported); F28 (Neuron genes) falls with age (rho -0.63/-0.74; replicated); F12 (Radial glia genes) falls with age (rho -0.56/-0.66; replicated); F17 (Oligo genes) rises with age (rho +0.20/+0.95; supported); F6 (Neuron genes) rises with age (rho +0.57/+0.63; replicated); F7 (Vascular genes) falls with age (rho -0.54/-0.63; replicated); F10 (Glioblast genes) rises with age (rho +0.55/+0.61; replicated) (+2 more).
- **Explained by cell-type mix** (age trend in raw activity that vanishes after composition adjustment): none.

## Limitations

- Factor activity per donor pools all of the donor's cells; the composition adjustment uses class means from all donors, so within-class differences between donors that are not age remain in the residual.
- Top-gene sets are fixed at 100 genes per pole by the export; a factor whose signal is spread over many genes is under-represented.
- human_dev donors were dissected differently; region mix is not adjusted (only class mix).

## What would strengthen this

- Export class x donor factor activity to test age within each class directly.

## Output files

- `factor_identity.csv` -- Per factor: class whose pseudobulk expresses its top positive genes most; top major class by mean activity; top region; v2/v3 agreement; QC flags from 08
- `lists_in_factors.csv` -- Per set x factor pole: overlap with the factor's 100 top genes vs the expression-matched expectation; exact Poisson-binomial p; BH over all tests
- `factor_age_trends_per_stratum.csv` -- Per stratum x factor: donor mean activity vs age, raw and composition-adjusted (exact permutation p)
- `factor_age_trends_combined.csv` -- v2 x v3 combined factor age trends; tier

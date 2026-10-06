# Sex differences in gene lists within cell types (age-adjusted, exploratory)

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/14_sex_differences.py` from `csv_exports/`._

## Question

Within a cell type and at the same age, are NDD gene lists expressed differently in male and female donors?

## Inputs

- `_cross_dataset/gene_id_map.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `results/01_data_audit/donor_sex_inferred.csv`

## Method

- Donor sex from 01 (Y genes, XIST); age points kept when all their donors share one sex.
- Per dataset x class, both chemistries pooled (>= 3 age points per sex): genes >= 5 CPM regressed on age and chemistry; residuals scaled by their SD; set score = mean scaled residual; male - female difference vs 2,000 random sets matched on expression decile; sex-linked genes excluded. BH per dataset; second null shuffles sex labels among age points within chemistry (all relabellings or 5,000); p = the larger of the two; 'consistent' = q < 0.05 and same sign within each chemistry.
- Gene level: Welch t on residuals; sex-linked genes as positive control.

## Key findings

- **Design** (age points male/female, pooled over chemistries): cortex: Glioblast 3/7, Neuroblast 4/8, Neuron 4/7, Neuronal IPC 4/8, Radial glia 4/9; human_dev: Erythrocyte 5/9, Glioblast 5/10, Immune 6/7, Neuroblast 6/10, Neuron 6/10, Neuronal IPC 6/10, Oligo 3/5, Radial glia 6/10, Vascular 6/9. Sex vs age: cortex rho +0.10, human_dev rho -0.11 (median over classes; far from 0 means sex and age are entangled). Possible sex relabellings per class: cortex Glioblast 36, cortex Neuroblast 210, cortex Neuron 126, cortex Neuronal IPC 210, cortex Radial glia 315, human_dev Erythrocyte 560, human_dev Glioblast 720, human_dev Immune 315, human_dev Neuroblast 1260, human_dev Neuron 1260, human_dev Neuronal IPC 1260, human_dev Oligo 10, human_dev Radial glia 1260, human_dev Vascular 980 -- with fewer than 20, the label-shuffle p cannot reach 0.05.
- **Positive control** (best rank of a Y gene or XIST among all genes, per class): cortex Glioblast #1; cortex Neuroblast #1; cortex Neuron #1; cortex Neuronal IPC #1; cortex Radial glia #1; human_dev Glioblast #1; human_dev Immune #1; human_dev Neuroblast #1; human_dev Neuron #1; human_dev Neuronal IPC #1; human_dev Oligo #1; human_dev Radial glia #1; human_dev Vascular #1. Rank 1-3 means the design detects real sex differences.
- **Autosomal / X genes differing by sex at q < 0.05**: cortex Neuroblast: KDM5C; cortex Neuronal IPC: PABPC1, LINC01551; cortex Radial glia: AGAP6; human_dev Glioblast: JPX; human_dev Neuron: KDM5C, JPX; human_dev Neuronal IPC: JPX.
- **Gene sets differing by sex** (passing both nulls; effect vs random genes in null SDs; label-shuffle p; within-chemistry differences v2/v3): none at q < 0.05. 24 set x class tests beat random genes (p < 0.01) but not shuffled sex labels (p >= 0.05): there the donors differ, but not by sex.

## Limitations

- Very few male donors (2-4 age points per chemistry): one male donor's peculiarities can pass for sex. The random-gene null removes donor-wide shifts, not donor-specific programs.
- Pooling chemistries means v2 and v3 are not independent replicates here; the within-chemistry signs are a weak consistency check.
- cortex donors are a subset of human_dev donors.

## What would strengthen this

- A larger donor panel, or per-donor pseudobulks from a dataset with balanced sexes.

## Output files

- `design.csv` -- Per dataset x class: age points by inferred sex and chemistry; Spearman of sex (male = 1) with age (entanglement)
- `set_sex_differences.csv` -- Per dataset x class x set: male - female mean scaled residual (age, chemistry adjusted) vs expression-matched random sets (perm_p) and vs shuffled sex labels (label_perm_p); within-chemistry differences; tier on the larger p
- `gene_sex_differences.csv` -- Per gene x class: Welch t of age/chemistry residuals, male vs female (positive control: sex-linked genes should rank first)

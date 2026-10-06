# Between-donor variability of NDD genes, against expression-matched genes

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/15_donor_variability.py` from `csv_exports/`._

## Question

After removing the age trend, do NDD gene lists vary less (or more) between donors than genes of the same expression level, within a cell class and in both donor sets?

## Inputs

- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- Per stratum x class with >= 5 age points: log2 TMM-CPM per class x age; genes >= 5 CPM; residual SD after a linear age fit; log SD as a percentile within 20 expression bins.
- Set statistic = mean percentile vs 5,000 random sets matched bin for bin; v2 x v3 signed Stouffer (weights sqrt(age points)), BH per dataset, tiered.

## Key findings

- **cortex: sets whose between-donor variability differs from matched genes** (0.5 = like matched genes): more variable: module:CX03 in Radial glia (mean percentile 0.64/0.69; replicated); module:CXw03 in Radial glia (mean percentile 0.65/0.71; replicated); module:CXw08 in Radial glia (mean percentile 0.78/0.81; replicated); module:CX02 in Neuroblast (mean percentile 0.65/0.81; replicated); module:CX03 in Neuroblast (mean percentile 0.67/0.65; replicated); module:CX05 in Neuroblast (mean percentile 0.71/0.75; replicated); module:CX07 in Neuroblast (mean percentile 0.70/0.74; replicated); module:CX09 in Neuroblast (mean percentile 0.74/0.82; replicated) (+41 more).
- **human_dev: sets whose between-donor variability differs from matched genes** (0.5 = like matched genes): less variable: seed:chromatin_transcription_regulators in Neuroblast (mean percentile 0.40/0.40; replicated); seed:chromatin_transcription_regulators in Radial glia (mean percentile 0.46/0.35; supported); seed:chromatin_transcription_regulators in Vascular (mean percentile 0.44/0.38; supported); seed:chromatin_transcription_regulators in Neuronal IPC (mean percentile 0.43/0.40; supported) | more variable: module:HD01 in Vascular (mean percentile 0.69/0.64; replicated); module:HD02 in Vascular (mean percentile 0.61/0.59; replicated); module:HDw01 in Vascular (mean percentile 0.74/0.69; replicated); module:HD01 in Immune (mean percentile 0.64/0.66; replicated); module:HD02 in Immune (mean percentile 0.66/0.61; replicated); module:HDw01 in Immune (mean percentile 0.66/0.69; replicated); module:HD02 in Neuroblast (mean percentile 0.73/0.71; replicated); module:HD03 in Neuroblast (mean percentile 0.70/0.65; replicated) (+65 more).
- **Consistently tight genes** (tightest fifth of their expression bin in both chemistries): cortex Neuroblast 490; cortex Neuron 476; cortex Neuronal IPC 451; cortex Radial glia 450; human_dev Erythrocyte 96; human_dev Glioblast 732; human_dev Immune 508; human_dev Neuroblast 898; human_dev Neuron 875; human_dev Neuronal IPC 809; human_dev Radial glia 821; human_dev Vascular 588. See consistently_tight_genes.csv.

## Limitations

- Residual SD mixes biological between-donor variation with dissection and sampling differences, and with any non-linear age trend.
- Co-regulated programmes (cell cycle, neuronal genes, most modules) move together when a donor's pseudobulk holds more cycling or more mature cells, so they read as 'more variable'; that is sub-type mix varying between donors, not loose control of each gene. Genes are treated as independent in the null, which makes such sets reach significance easily.
- An age point is about one donor, so residuals are between donors only approximately (two donors of one age are pooled).
- Dosage sensitivity is one reason for low variability; housekeeping-like stable expression is another.

## What would strengthen this

- C3: compare with gnomAD LOEUF -- do the tightest genes concentrate among loss-of-function intolerant genes?

## Output files

- `set_variability_per_stratum.csv` -- Per stratum x class x set: mean variability percentile of its genes vs matched random sets
- `set_variability_combined.csv` -- v2 x v3 combined; tier; direction vs matched genes
- `gene_variability.csv` -- Per gene x class x stratum: residual SD after age, percentile within expression bin (0 = least variable)
- `consistently_tight_genes.csv` -- Genes in the tightest 20% of their expression bin in both chemistries, per class
- `donor_variability_cortex.png` -- Between-donor variability of gene sets, cortex
- `donor_variability_human_dev.png` -- Between-donor variability of gene sets, human_dev

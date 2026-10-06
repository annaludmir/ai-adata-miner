# Transcriptomic age clock within cell types, validated across donor sets

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/16_age_clock.py` from `csv_exports/`._

## Question

Can a cell type's developmental age be predicted from expression by a model trained on other donors; which genes carry the signal; and do gene lists predict age better than random genes of the same expression?

## Inputs

- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`

## Method

- Class x age pseudobulks (log2 TMM-CPM), classes with >= 4 ages in both chemistries, genes >= 5 CPM; genes and ages centred within chemistry x class; genes scaled by training SD.
- Ridge regression (dual), penalty by leave-one-age-out CV within the training chemistry; train v2 -> test v3 and back; Spearman of predicted vs true (class-centred) age. Null: ages shuffled within class in training (200 times).
- List clocks: the list's genes only, fixed relative penalty 1.0; mean accuracy of both directions vs 200 random sets matched on expression decile; one-sided p, BH per dataset.

## Key findings

- **cortex: age is predictable across donor sets** -- train v2 -> test v3 rho +0.93 (shuffle p = 0.005), v3 -> v2 +0.78 (p = 0.025); per class (v2->v3 / v3->v2): Glioblast +0.80/+0.77, Neuroblast +0.94/+0.79, Neuron +0.90/+0.79, Neuronal IPC +0.94/+0.82, Radial glia +0.93/+0.79.
- **cortex: consensus clock genes** (5 in the top 200 of both models, same sign) -- rising with age: RDX, SCYL3, C2orf40; falling: NEFM, FAM207A.
- **cortex: gene sets that predict age better than matched random genes** (mean accuracy vs random): none at q < 0.05. Weakest relative to random: seed:id_dd_dominant +0.72 vs +0.74, seed:asd_high_confidence +0.77 vs +0.73, seed:s_phase +0.81 vs +0.78.
- **human_dev: age is predictable across donor sets** -- train v2 -> test v3 rho +0.93 (shuffle p = 0.005), v3 -> v2 +0.91 (p = 0.005); per class (v2->v3 / v3->v2): Erythrocyte +0.93/+0.92, Fibroblast +0.80/+1.00, Glioblast +0.96/+0.90, Immune +0.86/+0.95, Neuroblast +1.00/+0.89, Neuron +0.98/+0.91, Neuronal IPC +0.98/+0.89, Oligo +0.70/+0.80, Radial glia +0.98/+0.89, Vascular +0.90/+0.95.
- **human_dev: consensus clock genes** (21 in the top 200 of both models, same sign) -- rising with age: IFITM2, RARRES2, LGALS1, HBB, AC104389.6, ITM2A, C9orf24, APOE, DPP7, AC131571.1, EIF2AK2, APPL2; falling: HBE1, LIN28B, NR2F2, AC090204.1.
- **human_dev: gene sets that predict age better than matched random genes** (mean accuracy vs random): none at q < 0.05. Weakest relative to random: seed:g2m_phase +0.52 vs +0.72, seed:s_phase +0.64 vs +0.68, seed:synaptic_and_channels +0.63 vs +0.59.

## Limitations

- Few training samples (age points x classes, ~20-40 per chemistry): the clock is a ranking device, not a calibrated age estimate.
- The two chemistries cover different age ranges; within-class centring makes accuracy a rank agreement inside each test chemistry.
- human_dev classes pool regions whose sampling changes with age, so its clock can partly read region.

## What would strengthen this

- B5 pseudotime would separate developmental age from differentiation state within a class.

## Output files

- `clock_accuracy.csv` -- Cross-chemistry age prediction: Spearman of predicted vs true age (pooled over classes with class-centred ages, and per class); shuffle-null p
- `clock_gene_weights.csv` -- Gene weights (coefficient x training SD) per training chemistry, and consensus genes (top by |weight| in both with the same sign)
- `list_clocks.csv` -- Per gene set: mean cross-chemistry accuracy of a clock built on its genes only vs random sets matched on expression decile
- `age_clock_accuracy.png` -- Cross-chemistry age-clock accuracy per class

# Regional differences in gene lists within cell types (human_dev)

_Generated 2026-10-06 08:12 UTC by `downstream_analyses/13_regional_differences.py` from `csv_exports/`._

## Question

Within a cell type, do gene lists and modules differ between brain regions, in both donor sets, and could the difference be age?

## Inputs

- `human_dev__v2/05_confounds/crosstab_age_x_region.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_region__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/05_confounds/crosstab_age_x_region.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_region__pseudobulk_counts.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- Pseudobulk per class x region; log2 TMM-CPM within class; genes >= 5 CPM; Z across the class's regions; set score per region minus mean of the class's other regions; 2,000 random sets matched on expression decile; classes with >= 3 regions.
- v2 x v3 signed Stouffer, BH, tiered. Region age offsets from 05_confounds (cell-weighted mean age of the region in the stratum).

## Key findings

- **Tested**: 10 cell classes with >= 3 regions, 15 gene sets, 912 class x region x set tests in both chemistries; 507 tiered (356 replicated), 383 of them in neural classes (Glioblast, Neuroblast, Neuron, Neuronal IPC, Oligo, Radial glia).
- **Gene lists: no region difference within a class replicates.**
- **Seed NDD panels that differ between regions within a neural class** (strongest three per set; effect in null SDs v2/v3; region age offset vs the class's other regions): seed:chromatin_transcription_regulators (30 tiered): higher in Glioblast of Telencephalon (+3.9/+3.9; age +1.2/+0.9 wk), lower in Neuroblast of Medulla (-3.9/-3.6; age -0.8/-0.5 wk), lower in Neuroblast of Pons (-5.1/-3.6; age -0.1/-0.5 wk) | seed:id_dd_dominant (27 tiered): lower in Glioblast of Medulla (-4.8/-3.9; age -0.8/-1.0 wk), lower in Neuroblast of Medulla (-4.3/-4.0; age -0.8/-0.5 wk), higher in Neuroblast of Telencephalon (+4.6/+3.9; age +1.2/+1.3 wk) | seed:asd_high_confidence (27 tiered): higher in Neuroblast of Forebrain (+5.6/+5.1; age +1.0/+2.6 wk), lower in Neuroblast of Medulla (-3.7/-5.3; age -0.8/-0.5 wk), lower in Neuroblast of Pons (-3.7/-4.1; age -0.1/-0.5 wk) | seed:synaptic_and_channels (19 tiered): higher in Neuroblast of Forebrain (+4.0/+5.2; age +1.0/+2.6 wk), lower in Radial glia of Cerebellum (-2.8/-3.8; age -0.7/+0.5 wk), higher in Radial glia of Forebrain (+2.2/+3.9; age +1.0/+2.6 wk) | seed:epilepsy_dee (16 tiered): higher in Neuroblast of Forebrain (+2.8/+3.5; age +1.0/+2.6 wk), lower in Oligo of Medulla (-2.5/-4.0; age -1.3/-1.0 wk), higher in Oligo of Forebrain (+2.5/+2.7; age +0.7/+2.2 wk).
- **Modules that differ between regions within a neural class** (strongest three per set; effect in null SDs v2/v3; region age offset vs the class's other regions): module:HD01 (27 tiered): lower in Glioblast of Diencephalon (-13.3/-7.4; age -0.3/-1.2 wk), higher in Neuronal IPC of Forebrain (+11.3/+9.2; age +1.0/+2.6 wk), higher in Radial glia of Telencephalon (+10.0/+17.7; age +1.2/+1.3 wk) | module:HD02 (26 tiered): higher in Glioblast of Pons (+5.6/+6.9; age -0.1/-0.9 wk), higher in Neuronal IPC of Cerebellum (+6.8/+8.6; age -0.7/+0.5 wk), higher in Radial glia of Pons (+3.7/+5.2; age -0.1/-0.5 wk) | module:HD03 (25 tiered): lower in Glioblast of Telencephalon (-5.2/-6.0; age +1.2/+0.9 wk), lower in Neuroblast of Cerebellum (-14.2/-7.8; age -0.7/+0.5 wk), higher in Neuroblast of Diencephalon (+7.8/+10.3; age -0.3/-0.7 wk) | module:HD05 (24 tiered): higher in Glioblast of Diencephalon (+6.8/+8.0; age -0.3/-1.2 wk), higher in Glioblast of Pons (+4.8/+5.3; age -0.1/-0.9 wk), lower in Glioblast of Telencephalon (-8.7/-9.8; age +1.2/+0.9 wk) | module:HD11 (24 tiered): lower in Glioblast of Telencephalon (-4.6/-5.0; age +1.2/+0.9 wk), higher in Neuroblast of Diencephalon (+4.6/+4.6; age -0.3/-0.7 wk), lower in Neuroblast of Telencephalon (-7.2/-9.2; age +1.2/+1.3 wk) | module:HDw01 (28 tiered): lower in Glioblast of Diencephalon (-14.0/-7.9; age -0.3/-1.2 wk), higher in Neuronal IPC of Forebrain (+12.0/+9.4; age +1.0/+2.6 wk), higher in Radial glia of Telencephalon (+8.2/+14.9; age +1.2/+1.3 wk) | module:HDw02 (28 tiered): higher in Glioblast of Diencephalon (+6.9/+7.9; age -0.3/-1.2 wk), higher in Glioblast of Midbrain (+6.5/+4.5; age -0.4/-0.9 wk), higher in Glioblast of Pons (+5.2/+5.7; age -0.1/-0.9 wk) | module:HDw03 (26 tiered): lower in Glioblast of Forebrain (-4.3/-4.0; age +1.0/+2.2 wk), lower in Glioblast of Telencephalon (-4.7/-5.1; age +1.2/+0.9 wk), higher in Neuroblast of Pons (+6.2/+5.3; age -0.1/-0.5 wk) | module:HDw04 (38 tiered): higher in Glioblast of Cerebellum (+3.8/+4.6; age -0.7/+0.1 wk), lower in Neuron of Telencephalon (-5.6/-7.6; age +1.2/+1.3 wk), lower in Radial glia of Telencephalon (-7.7/-10.2; age +1.2/+1.3 wk) | module:HD04 (18 tiered): lower in Radial glia of Telencephalon (-4.3/-3.8; age +1.2/+1.3 wk), lower in Neuron of Telencephalon (-3.7/-4.1; age +1.2/+1.3 wk), higher in Glioblast of Cerebellum (+3.4/+3.6; age -0.7/+0.1 wk).
- **Non-neural classes, likely ambient RNA** (not interpreted): Erythrocyte 5, Fibroblast 14, Immune 14, Vascular 16 tiered list / panel results (47% 'higher'); neural gene sets in these cells mostly measure ambient neural RNA, so they are listed in region_contrast_combined.csv (neural_class = False) but not interpreted.
- **Age check**: no set's region effects track region age (|rho| >= 0.5 in both chemistries), so the differences above are not simply age.

## Limitations

- Non-neural cells carry ambient neural RNA from dissociation, so their region differences for neural gene sets track how much neural tissue the region had; they are reported separately and not interpreted.
- Each class x region pseudobulk pools several donors of varying ages; regions were sampled at different ages, so a region difference can be partly age (see the offsets and the age check).
- Region labels are dissection labels; 'Forebrain' overlaps Telencephalon and Diencephalon.
- Pseudobulks pool cells, not donors, so the per-chemistry p-values treat genes as the random unit; replication across donor sets is the safeguard.

## What would strengthen this

- B1: class x region x age pseudobulks would compare regions at matched ages directly.

## Output files

- `region_contrast_per_stratum.csv` -- Per stratum x class x region x set: set score minus other regions of the class, vs expression-matched random sets; region age offset
- `region_contrast_combined.csv` -- Per class x region x set: v2 x v3 combined; tier; region age offset per chemistry; neural_class = False marks classes where ambient neural RNA is likely
- `region_effect_vs_age.csv` -- Per set x stratum: Spearman between region effects and region age offsets (a strong value means region differences may be age differences)
- `regional_differences.png` -- Gene sets by region within cell class (combined Z)

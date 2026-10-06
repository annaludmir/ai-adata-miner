# Robustness checks: gene length, quality metrics, dissociation stress

_Generated 2026-10-06 04:44 UTC by `downstream_analyses/19_robustness_checks.py` from `csv_exports/`._

## Question

Are the list-level results explained by gene length, by quality metrics, or by dissociation stress?

## Inputs

- `_cross_dataset/gene_id_map.csv`
- `cortex__v2/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v3/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `human_dev__v2/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `results/03_age_trends_within_cell_class/age_trends_per_stratum.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- 1. Length = genomic span; each test run with random genes matched on expression only and on expression x length tertile, in the same code: class preference (06), age coordination (06), across-cluster coherence (07); v2 x v3 combined and tiered for both.
- 2. Module score per cluster (mean Z) vs the cluster's median QC values, raw and within class; flag |within-class rho| >= 0.5 in both chemistries.
- 3. Stress score = mean Z of FOS, FOSB, JUN, JUNB, EGR1, IER2, IER3, ATF3, DUSP1, ZFP36, HSPA1A, HSPA1B, HSPA8, HSP90AA1 per class x age; Spearman with age; list age trends plain and partial on the stress score.

## Key findings

- **Gene length** (median length percentile among expressed genes; 0.5 = typical): cortex: seed:synaptic_and_channels 0.84, seed:epilepsy_dee 0.79, seed:asd_high_confidence 0.79, seed:chromatin_transcription_regulators 0.78, seed:id_dd_dominant 0.76; human_dev: seed:synaptic_and_channels 0.80, seed:epilepsy_dee 0.74, seed:asd_high_confidence 0.74, seed:chromatin_transcription_regulators 0.72, seed:id_dd_dominant 0.71.
- **Length check, class preference**: 52 tiered results with expression-matched nulls; 7 lost when length is matched too (cortex seed:epilepsy_dee in Neuroblast; cortex seed:id_dd_dominant in Neuroblast; cortex seed:id_dd_dominant in Radial glia; human_dev seed:chromatin_transcription_regulators in Neural crest; human_dev seed:chromatin_transcription_regulators in Vascular; human_dev seed:id_dd_dominant in Neural crest; human_dev seed:id_dd_dominant in Neuronal IPC), 1 gained.
- **Length check, age coordination**: 31 tiered results with expression-matched nulls; 4 lost when length is matched too (human_dev seed:asd_high_confidence in Glioblast; human_dev seed:epilepsy_dee in Neuroblast; human_dev seed:synaptic_and_channels in Neuronal IPC; human_dev seed:asd_high_confidence in Vascular), 0 gained.
- **Length check, coherence**: 10 tiered results with expression-matched nulls; 0 lost when length is matched too, 0 gained.
- **Modules tracking quality metrics within cell classes** (|rho| >= 0.5 in both chemistries): none.
- **Dissociation-stress score vs age** (rho v2/v3): human_dev Neuron -0.88/-0.45 (falls with age, supported); human_dev Neuroblast -0.79/-0.40 (falls with age, supported); human_dev Neuronal IPC -0.73/-0.19; cortex Neuronal IPC -0.61/-0.71; cortex Neuroblast -0.71/-0.37; human_dev Erythrocyte -0.52/+0.36; human_dev Glioblast -0.04/+0.46; human_dev Radial glia -0.23/+0.02; cortex Radial glia -0.04/-0.29; human_dev Vascular -0.10/+0.33; human_dev Immune -0.02/-0.14; cortex Neuron -0.57/+0.90.
- **List age trends that shrink by half or more when stress is controlled** (classes whose stress score changes with age; mean rho plain -> partial): none.

## Limitations

- Genomic span is a proxy for what matters in 3' data (intron content, transcript length); tertiles are coarse.
- Cluster QC medians summarise many cells; a module can still track QC inside clusters.
- The stress gene set also responds to real activity (immediate-early genes in neurons), so a stress trend is not necessarily technical.

## What would strengthen this

- B8: re-export with stricter per-cell QC and rerun the key analyses.

## Output files

- `list_gene_length.csv` -- Per set: median gene length vs expressed genes (percentile 0.5 = typical)
- `length_check_class_preference.csv` -- Class preference (as 06) with expression-only and expression x length matched nulls
- `length_check_age_coordination.csv` -- Age coordination (as 06) with both nulls
- `length_check_coherence.csv` -- Across-cluster coherence (as 07) with both nulls
- `length_check_summary.csv` -- Per test: tier with expression-only vs expression x length nulls; lost / gained with length
- `module_qc_per_stratum.csv` -- Per module x QC metric x stratum: Spearman across clusters, raw and within class
- `module_qc.csv` -- Per module x QC metric: v2/v3 rhos; flag = |within-class rho| >= 0.5 in both, same sign
- `stress_vs_age.csv` -- Dissociation-stress score vs age per class; v2 x v3 combined
- `list_age_trends_given_stress.csv` -- Per list x class x stratum: mean rho of its genes with age, plain and partial on the stress score

# Co-expression modules across fine clusters

_Generated 2026-10-05 06:37 UTC by `downstream_analyses/08_coexpression_modules.py` from `csv_exports/`._

## Question

Without starting from a list: which groups of genes co-vary across cell clusters robustly enough to be found again in independent donors, where do they peak, how do they change with age, and which gene lists concentrate in them?

## Inputs

- `cortex__v2/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v3/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v3/11_panels/panel_coverage.csv`
- `human_dev__v2/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v3/11_panels/panel_coverage.csv`

## Method

- Clusters as in 07 (cluster_ClustersSurprise, cluster_cluster_id; >= 100 cells; donor-dominated dropped; log2 TMM-CPM; Spearman).
- Universe: the 3,000 most variable genes expressed in both chemistries, per context (across clusters; within cell classes after removing class means).
- Modules: average linkage on 1 - r of the v2/v3 consensus (Fisher z), cut at the 99% quantile of the consensus correlations, bounded to 0.3-0.5 (within-class residuals correlate far less than profiles do); modules > 250 genes re-split at +0.1 (to 0.8); < 20 genes unassigned. Named CX/HD + number (w = within-class context).
- Robust = the best-matching module found in v2 alone is preserved in v3, and vice versa (held-out coherence vs 1,000 matched random sets, BH q < 0.05; Jaccard >= 0.3).
- Age trends: module score (mean z of members over a class's age points, TMM log CPM) vs age, exact permutation p, v2/v3 signed Stouffer, tiered as elsewhere.
- Enrichment: overlap vs random sets matched on level x spread within the universe (empirical p; beyond the null's resolution, a normal approximation from its mean and SD); BH per dataset x context.
- Labels in [brackets]: the seed reference panels (cell-cycle phase, cell-class markers, patterning) a module concentrates (q < 0.05, >= 2x matched expectation; top two).

## Key findings

- **cortex / across clusters: 11 modules, 7 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - CX01 [G2/M phase] (225 genes; peak Neuronal IPC; hubs HMGB2, NUSAP1, SMC4, PIMREG, KIFC1)
  - CX02 (221 genes; peak Neuron; hubs MAPT, XPR1, SHTN1, SCG5, SPTAN1; holds seed:epilepsy_dee (6, 7.0x), seed:synaptic_and_channels (7, 3.4x))
  - CX03 [radial glia markers] (165 genes; peak Glioblast; hubs GULP1, SOX9, GLI3, GATM, CLU)
  - CX04 [S phase] (156 genes; peak Glioblast/Neuronal IPC; hubs GMNN, PCNA, CENPK, TYMS, HELLS)
  - CX05 (113 genes; peak Neuron; hubs FGF12, GPR22, GAS7, ZBTB38, FXYD7)
  - CX07 (79 genes; peak Neuroblast; hubs GPC2, MLLT11, TAGLN3, DISP3, RASGEF1B)
  - CX09 [neuronal ipc markers] (30 genes; peak Neuronal IPC; hubs NHLH1, ELAVL2, EOMES, INSM1, NEUROD4; falls with age in Neuronal IPC)
- **cortex / within class: 9 modules, 6 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - CXw01 (241 genes; peak Neuron; hubs MAPT, RAB3A, RTN1, SEZ6L2, TTC9B; holds seed:synaptic_and_channels (10, 3.8x), seed:asd_high_confidence (5, 3.9x))
  - CXw02 [G2/M phase] (199 genes; peak Neuronal IPC; hubs NUSAP1, HMGB2, CCNA2, GTSE1, TPX2)
  - CXw03 [radial glia markers] (169 genes; peak Glioblast; hubs SOX2, CLU, GULP1, SOX9, PON2)
  - CXw04 [S phase] (104 genes; peak Neuronal IPC; hubs TYMS, GMNN, CENPK, PCNA, ORC6)
  - CXw05 (81 genes; peak Neuron; hubs CNTN1, MEF2C, DACT1, SCN2A, VSTM2L)
  - CXw08 [neuronal ipc markers] (26 genes; peak Neuronal IPC; hubs NHLH1, ELAVL2, NEUROD4, EOMES, HES6; falls with age in Neuronal IPC, rises with age in Radial glia)
- **human_dev / across clusters: 15 modules, 6 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - HD01 [G2/M phase; S phase] (272 genes; peak Neuronal IPC, Telencephalon/Forebrain; hubs SMC4, CENPK, MAD2L1, KIF11, CKS1B)
  - HD02 (261 genes; peak Fibroblast/Vascular, Medulla/Hindbrain; hubs ANXA5, SUCLG2, SERPINH1, STK3, PLIN3; rises with age in Vascular)
  - HD03 (225 genes; peak Neuron, Pons; hubs TMEM35A, SV2A, SCN3B, PTPN5, JPH4; holds seed:synaptic_and_channels (10, 4.1x))
  - HD04 (72 genes; peak Neuroblast, Cerebellum/Pons; hubs TUBB3, CD24, GPC2, MLLT11, SMPD3; falls with age in Neuroblast, falls with age in Neuron)
  - HD05 (43 genes; peak Glioblast, Diencephalon/Hindbrain; hubs NME5, C9orf116, SPATA17, CFAP54, ENKUR; falls with age in Neuronal IPC)
  - HD11 (25 genes; peak Neuron, Diencephalon/Pons; hubs MTUS2, TENM2, TAC1, AC011369.1, PLEKHA6; falls with age in Neuroblast, falls with age in Neuron)
- **human_dev / within class: 6 modules, 4 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.44)
  - HDw01 [G2/M phase; S phase] (197 genes; peak Neuronal IPC, Telencephalon/Forebrain; hubs TOP2A, CDK1, KIFC1, KIF11, GTSE1; rises with age in Neuronal IPC)
  - HDw02 (48 genes; peak Glioblast, Diencephalon/Hindbrain; hubs NME5, SPATA17, KIF9, MOK, C5orf49; falls with age in Neuronal IPC)
  - HDw03 (29 genes; peak Neuron, Pons/Hindbrain; hubs MAP1LC3A, LY6H, PCSK1N, DDX25, FUCA1)
  - HDw04 (28 genes; peak Glioblast/Radial glia, Pons/Midbrain; hubs IRX2, IRX3, IRX1, C5orf38, IRX5; falls with age in Neuronal IPC, falls with age in Radial glia, falls with age in Neuroblast, falls with age in Neuron)
- **Gene sets concentrated in robust modules** (q < 0.05, >= 2x matched expectation): seed:asd_high_confidence: cortex CXw01 (5 genes, 3.9x) | seed:epilepsy_dee: cortex CX02 (6 genes, 7.0x) | seed:synaptic_and_channels: human_dev HD03 (10 genes, 4.1x); cortex CXw01 (10 genes, 3.8x); cortex CX02 (7 genes, 3.4x).
- **cortex and human_dev share 2 of 20 cortex modules** (Jaccard >= 0.3 with a human_dev module; same donors, so this is reproducibility of processing and regional pooling, not replication).

## Limitations

- Modules depend on the cut and the universe: they are a summary of correlation structure, not discrete biological units. The tree path in modules.csv shows which modules came from re-splitting one larger module.
- Cluster pseudobulks are not independent samples; preservation is judged against matched random gene sets in held-out donors, not by absolute correlation.
- Gene-set nulls ignore correlation among set members; replication and fold >= 2 are the guards.
- human_dev modules pool brain regions; cortex and human_dev share donors.

## What would strengthen this

- Annotate robust modules with GO / TF-target enrichment (needs annotation files).
- Compare module eigengenes between regions once a region x cluster export exists.

## Output files

- `module_annotation.csv` -- Overlap of each module with seed reference panels (cell-cycle phase, cell-class markers, patterning) vs matched random sets; used to label modules
- `modules.csv` -- Per module: size, robustness (independent discovery + held-out preservation), peak class/region, top clusters, hub genes, members
- `module_membership.csv` -- Gene -> module, with kME per chemistry and hub flag
- `module_age_trends_per_stratum.csv` -- Module score vs age within each cell class, per chemistry (exact permutation p)
- `module_age_trends_combined.csv` -- v2 x v3 combined module age trends per class; tier replicated / supported
- `module_gene_set_enrichment.csv` -- Overlap of each user list / seed NDD panel with each module vs matched random sets
- `cortex_vs_human_dev_module_match.csv` -- Best-matching human_dev module for each cortex module (Jaccard; same donors)
- `modules_cortex_across_clusters.png` -- Module eigengene by cell class, cortex / across_clusters
- `modules_cortex_within_class.png` -- Module eigengene by cell class, cortex / within_class
- `modules_human_dev_across_clusters.png` -- Module eigengene by cell class, human_dev / across_clusters
- `modules_human_dev_within_class.png` -- Module eigengene by cell class, human_dev / within_class
- `module_enrichment_cortex_across_clusters.png` -- Gene-set enrichment in modules, cortex / across_clusters
- `module_enrichment_cortex_within_class.png` -- Gene-set enrichment in modules, cortex / within_class
- `module_enrichment_human_dev_across_clusters.png` -- Gene-set enrichment in modules, human_dev / across_clusters
- `module_enrichment_human_dev_within_class.png` -- Gene-set enrichment in modules, human_dev / within_class

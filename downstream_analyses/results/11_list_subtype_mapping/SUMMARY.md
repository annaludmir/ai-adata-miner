# Gene lists at sub-type resolution: fine clusters and neuron sub-type axes

_Generated 2026-10-06 04:42 UTC by `downstream_analyses/11_list_subtype_mapping.py` from `csv_exports/`._

## Question

Within broad cell classes, which fine clusters (sub-types) does each gene list concentrate in, consistently in both donor sets; and among neurons, does the list follow the excitatory-inhibitory and deep-upper layer axes?

## Inputs

- `cortex__v2/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v2/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v2/10_markers/top_markers_cluster_ClustersSurprise.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v3/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v3/10_markers/top_markers_cluster_ClustersSurprise.csv`
- `human_dev__v2/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/10_markers/top_markers_cluster_cluster_id.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v3/10_markers/top_markers_cluster_cluster_id.csv`

## Method

- A: cluster pseudobulks (07/08 rules), gene Z-scores across clusters, set score = mean Z; 1,000 random sets matched member by member on mean level x spread; per cluster effect in null SDs and empirical p; the same cluster in v2 and v3 (shared label, or mutual best expression match with r >= 0.5 where labels differ between chemistries) combined by signed Stouffer, BH per dataset over set x cluster tests, tiered. Classes / regions of a set's replicated clusters compared with all tested clusters.
- B: neuron clusters (purity >= 0.6); axis = mean Z of one marker group minus the other (excitatory SLC17A6, SLC17A7, NEUROD2, NEUROD6 vs inhibitory GAD1, GAD2, SLC32A1, DLX5, human_dev; deep BCL11B, TBR1, FEZF2, SOX5 vs upper SATB2, CUX2, POU3F2, POU3F3, cortex and human_dev telencephalic excitatory clusters). Spearman of set score with the axis vs matched random sets, axis markers removed; v2 x v3 combined.
- Seed NDD panels and user lists (GWAS lists one gene per locus).

## Key findings

- **cortex: clusters where each set concentrates** (28 clusters in both chemistries, paired by expression match; effect in null SDs v2/v3): seed:asd_high_confidence: 4 replicated clusters, top cluster 127~2 (Neuron; DLX2, FAM19A2; +3.7/+3.8); cluster 207~151 (Neuron; +2.8/+5.0); cluster 48~112 (Neuron; EBF1, GREM2; +2.7/+3.8) -- over-represented: Neuron (3.5x) | seed:chromatin_transcription_regulators: 3 replicated clusters, top cluster 29~30 (Neuroblast; +3.4/+4.0); cluster 48~112 (Neuron; EBF1, GREM2; +2.0/+2.7); cluster 137~67 (Neuron; +1.7/+2.2) -- over-represented: Neuron (2.3x) | seed:synaptic_and_channels: 9 replicated clusters, top cluster 48~112 (Neuron; EBF1, GREM2; +4.1/+4.6); cluster 127~2 (Neuron; DLX2, FAM19A2; +4.3/+3.4); cluster 207~151 (Neuron; +5.0/+5.2) -- over-represented: Neuron (2.7x) | seed:epilepsy_dee: 5 replicated clusters, top cluster 207~151 (Neuron; +2.7/+3.5); cluster 89~16 (Neuron; +2.3/+2.8); cluster 18~43 (Neuron; +1.9/+2.9) -- over-represented: Neuron (3.5x) | seed:id_dd_dominant: 3 replicated clusters, top cluster 29~30 (Neuroblast; +3.0/+3.6); cluster 105~41 (Radial glia; +2.4/+1.7); cluster 137~67 (Neuron; +1.7/+1.7).
- **human_dev: clusters where each set concentrates** (447 clusters in both chemistries, paired by shared label; effect in null SDs v2/v3): seed:asd_high_confidence: 115 replicated clusters, top cluster 301 (Radial glia, Telencephalon; RAB11FIP3, ZNF785; +5.1/+3.9); cluster 476 (Neuron, Forebrain; +5.5/+4.1); cluster 376 (Neuron, Midbrain; MIR124-1HG, ADAMTS20; +5.8/+5.0) -- over-represented: Neuron (2.1x) | seed:chromatin_transcription_regulators: 112 replicated clusters, top cluster 188 (Neuron, Forebrain; AC092691.1, ZSWIM5; +5.9/+5.4); cluster 350 (Neuron, Forebrain; AC092422.1, RARB; +4.1/+4.0); cluster 352 (Neuron, Forebrain; BCL11B, KLF3-AS1; +5.3/+4.9) -- over-represented: Neuronal IPC (1.8x), Neuron (1.7x) | seed:epilepsy_dee: 120 replicated clusters, top cluster 314 (Neuron, Forebrain; SYN2, GRIA1; +3.5/+3.9); cluster 350 (Neuron, Forebrain; AC092422.1, RARB; +3.3/+4.3); cluster 363 (Neuron, Telencephalon; +3.6/+3.5) -- over-represented: Neuron (2.5x) | seed:id_dd_dominant: 89 replicated clusters, top cluster 177 (Radial glia, Forebrain; ASPM, SGO2; +4.1/+4.0); cluster 359 (Neuron, Forebrain; +5.0/+4.5); cluster 311 (Neuron, Forebrain; SLC24A2, EPHA5; +3.8/+3.4) -- over-represented: Neuronal IPC (2.4x), Immune (1.7x) | seed:synaptic_and_channels: 182 replicated clusters, top cluster 306 (Neuroblast, Telencephalon; OCA2, SH3RF3; +4.5/+3.7); cluster 494 (Neuron, Telencephalon; COL25A1, GABRB2; +3.8/+4.5); cluster 492 (Neuron, Diencephalon; +4.6/+4.7) -- over-represented: Neuron (2.2x).
- **Cluster profiles reproduce across donor sets** (Spearman of set scores v2 vs v3 over paired clusters): cortex median 0.83 (range 0.62-0.94); human_dev median 0.89 (range 0.78-0.97).
- **cortex, deep vs upper layer** (38/31 neuron clusters v2/v3; marker groups correlate rho [-0.35, -0.36]): seed:epilepsy_dee toward deep layer (rho +0.53/+0.71; replicated); seed:synaptic_and_channels toward deep layer (rho +0.56/+0.60; replicated).
- **human_dev, excitatory vs inhibitory** (212/168 neuron clusters v2/v3; marker groups correlate rho [-0.61, -0.67]): no set leans to either pole in both donor sets.
- **human_dev, deep vs upper layer** (15/15 neuron clusters v2/v3; marker groups correlate rho [-0.42, 0.36]): no set leans to either pole in both donor sets.

## Limitations

- Clusters are not donors: a cluster's v2 and v3 cells come from different people, which is what makes the per-cluster replication meaningful, but within a chemistry a cluster can be dominated by few donors (one-donor clusters are excluded).
- Cluster scores are relative to the other clusters of the same file; a list 'enriched' in a cluster is high there compared with the rest of the dataset, not necessarily specific.
- Neuron axes are confounded with maturation: deep-layer neurons are born first, so a list leaning deep may simply be higher in older neurons.

## What would strengthen this

- Pair with B5 (pseudotime) to separate sub-type identity from maturation.
- Name clusters with a reference atlas instead of their top markers.

## Output files

- `cluster_pairs.csv` -- How clusters were matched between v2 and v3: shared label, or mutual best expression match (r = Pearson of gene-centred profiles)
- `cluster_scores_per_stratum.csv` -- Per set x cluster x stratum: mean Z of set genes, effect vs matched random sets, p
- `cluster_scores_combined.csv` -- Per set x cluster: v2/v3 combined enrichment, tier, cluster annotation and top markers
- `list_cluster_concentration.csv` -- Per set: share of its replicated enriched clusters in each class / region vs share of tested clusters (ratio > 1 = over-represented)
- `cluster_profile_agreement.csv` -- Per set: Spearman of its cluster scores between v2 and v3 over paired clusters
- `neuron_axes.csv` -- Neuron sub-type axes: clusters and markers used per stratum
- `neuron_axis_tests_per_stratum.csv` -- Per set x axis x stratum: Spearman of set score with the axis across neuron clusters, vs matched random sets
- `neuron_axis_tests_combined.csv` -- Per set x axis: v2/v3 combined association with the axis; direction = the pole the set leans to
- `cluster_concentration_cortex.png` -- Replicated cluster enrichment by cell class, cortex
- `cluster_concentration_human_dev.png` -- Replicated cluster enrichment by cell class, human_dev
- `neuron_axes.png` -- Gene sets along neuron sub-type axes (combined Z)

# Co-expression rewiring: coherence of gene lists in early vs late clusters

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/18_coexpression_rewiring.py` from `csv_exports/`._

## Question

Do gene lists and modules co-express more or less tightly in late than in early clusters, within the same mix of cell classes, in both donor sets?

## Inputs

- `cortex__v2/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v2/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/04_clusters/cluster_profile_ClustersSurprise.csv`
- `cortex__v3/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `human_dev__v2/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- Cluster pseudobulks (07/08 rules). Per class with >= 10 clusters, clusters split at the class's median cluster age; halves pooled over classes.
- Coherence = mean pairwise Spearman over one half's clusters; contexts across clusters and within class (class means removed inside each half). Late - early vs 500 random sets matched on mean level x spread; v2 x v3 signed Stouffer, BH per dataset x context, tiered.

## Key findings

- **Age split** (median cluster age early -> late, per class): cortex v2: Glioblast 9.2->10, Neuroblast 9.2->10, Neuron 9.2->10, Neuronal IPC 9.2->9.75, Radial glia 9.2->9.5; cortex v3: Glioblast 11.5->12, Neuroblast 11.5->12, Neuron 11.5->11.75, Neuronal IPC 11.5->12, Radial glia 11.5->12; human_dev v2: Glioblast 8->9.2, Neuroblast 6.7->8.05, Neuron 6.9->8, Neuronal IPC 7.2->9.2, Radial glia 6.7->8, Vascular 6.9->9.2; human_dev v3: Glioblast 12->13, Neuroblast 6.9->11.5, Neuron 6.9->8.5, Neuronal IPC 6.9->12, Radial glia 6.9->6.9, Vascular 11.75->12.5.
- **cortex, across clusters: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:chromatin_transcription_regulators more coherent late (0.15->0.18 / 0.16->0.21; replicated); seed:asd_high_confidence less coherent late (0.14->0.12 / 0.16->0.15; replicated).
- **cortex, within class: sets whose coherence changes with age** (coherence early -> late, v2 / v3): module:CX07 more coherent late (0.43->0.46 / 0.38->0.51; replicated); seed:chromatin_transcription_regulators more coherent late (0.10->0.12 / 0.16->0.21; replicated).
- **human_dev, across clusters: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:id_dd_dominant more coherent late (0.14->0.21 / 0.11->0.23; replicated); module:HD02 less coherent late (0.70->0.67 / 0.72->0.68; replicated); module:HD03 less coherent late (0.77->0.73 / 0.77->0.71; replicated); module:HD04 less coherent late (0.68->0.63 / 0.68->0.65; replicated); module:HDw02 more coherent late (0.47->0.51 / 0.52->0.59; replicated); seed:synaptic_and_channels less coherent late (0.48->0.41 / 0.48->0.43; replicated); module:HDw03 less coherent late (0.71->0.68 / 0.66->0.63; replicated); module:HDw04 more coherent late (0.46->0.50 / 0.48->0.57; replicated) (+3 more).
- **human_dev, within class: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:g2m_phase more coherent late (0.59->0.67 / 0.64->0.68; replicated); seed:chromatin_transcription_regulators more coherent late (0.26->0.32 / 0.21->0.30; replicated); seed:id_dd_dominant more coherent late (0.17->0.23 / 0.12->0.24; replicated); module:HD03 less coherent late (0.29->0.26 / 0.31->0.27; replicated); module:HD11 more coherent late (0.21->0.46 / 0.22->0.47; replicated); module:HDw02 more coherent late (0.42->0.50 / 0.50->0.57; replicated); module:HDw04 more coherent late (0.45->0.54 / 0.50->0.61; replicated); module:HD02 less coherent late (0.18->0.17 / 0.23->0.20; replicated) (+4 more).

## Limitations

- A cluster's age is the median age of its cells; clusters pooling several donors blur the split, which biases changes towards zero.
- Early and late halves can differ in sub-type mix within a class, which changes co-expression without any rewiring inside cells.
- Coherence over ~half the clusters is noisier than 07's estimate over all of them.

## What would strengthen this

- B4 (within-cell co-expression per age) would test rewiring inside one cell type directly.

## Output files

- `age_split.csv` -- Per stratum x class: clusters per half and median cluster age of each half
- `rewiring_per_stratum.csv` -- Per stratum x context x set: coherence in early and late clusters, late - early vs matched random sets
- `rewiring_combined.csv` -- v2 x v3 combined change in coherence; BH per dataset within context; tier
- `coexpression_rewiring.png` -- Change in gene-set coherence, late vs early clusters (combined Z)

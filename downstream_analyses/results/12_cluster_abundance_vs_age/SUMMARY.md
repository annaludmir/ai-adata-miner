# Sub-type (cluster) abundance over development, within cell classes

_Generated 2026-10-06 04:43 UTC by `downstream_analyses/12_cluster_abundance_vs_age.py` from `csv_exports/`._

## Question

Within each cell class, which sub-types (clusters) expand or shrink with age, consistently in both donor sets?

## Inputs

- `cortex__v2/04_clusters/cluster_composition_leiden_scVI_by_donor.csv`
- `cortex__v2/04_clusters/cluster_profile_leiden_scVI.csv`
- `cortex__v2/05_confounds/crosstab_age_x_donor.csv`
- `cortex__v2/09_pseudobulk/cluster_leiden_scVI__pseudobulk_counts.csv`
- `cortex__v2/10_markers/top_markers_cluster_leiden_scVI.csv`
- `cortex__v3/04_clusters/cluster_composition_leiden_scVI_by_donor.csv`
- `cortex__v3/04_clusters/cluster_profile_leiden_scVI.csv`
- `cortex__v3/05_confounds/crosstab_age_x_donor.csv`
- `cortex__v3/09_pseudobulk/cluster_leiden_scVI__pseudobulk_counts.csv`
- `cortex__v3/10_markers/top_markers_cluster_leiden_scVI.csv`
- `human_dev__v2/04_clusters/cluster_composition_cluster_id_by_donor.csv`
- `human_dev__v2/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v2/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/10_markers/top_markers_cluster_cluster_id.csv`
- `human_dev__v3/04_clusters/cluster_composition_cluster_id_by_donor.csv`
- `human_dev__v3/04_clusters/cluster_profile_cluster_id.csv`
- `human_dev__v3/05_confounds/crosstab_age_x_donor.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v3/10_markers/top_markers_cluster_cluster_id.csv`

## Method

- Clusterings with labels shared across chemistries: cortex leiden_scVI, human_dev cluster_id (share of labels whose v2 expression profile best matches the same label in v3: cortex 100%, human_dev 90%).
- Cells per cluster x donor from 04's cluster sizes and donor fractions; parent = dominant class (cortex) or class x dominant region (human_dev); clusters with class purity >= 0.5; donors with >= 50 cells in the parent, >= 5 donors.
- Per donor, CLR of counts across the parent's clusters; Spearman with donor age, exact permutation p; v2 x v3 signed Stouffer (weights sqrt(donors)), BH per dataset, tiered.

## Key findings

- **cortex: sub-types whose share within their class changes with age** (2 of 9 clusters tiered; share youngest -> oldest donor, v2 / v3): cluster 2 in Neuroblast expands [PPP1R17, EPHB6, SNCB, NEUROD2, PRKX] (25%->74% / 1%->81%; replicated); cluster 9 in Neuroblast shrinks [LHX1, MAB21L1, RELN, TP73, PGF] (22%->0% / 20%->0%; supported).
- **human_dev (parent = class x region): sub-types whose share within their class changes with age** (131 of 293 clusters tiered; share youngest -> oldest donor, v2 / v3): cluster 227 in Radial glia | Telencephalon expands [STK17A, ESR2, CBFA2T2, SYNE2, ZC3H12C] (0%->12% / 0%->26%; replicated); cluster 488 in Neuron | Telencephalon shrinks [LEMD1, SRP14] (12%->0% / 5%->0%; replicated); cluster 347 in Neuron | Telencephalon shrinks [PLEKHA8, RFX7, STX16-NPEPL1] (5%->0% / 4%->0%; replicated); cluster 469 in Neuroblast | Telencephalon shrinks [KCNC2, ISLR2, FAM155A, RBM19, CASTOR3] (78%->0% / 40%->0%; replicated); cluster 357 in Neuron | Telencephalon shrinks [AC004943.2, RAPGEF2] (6%->0% / 4%->0%; replicated); cluster 573 in Neuroblast | Midbrain expands [EBF2, LHX2] (0%->45% / 0%->81%; replicated); cluster 225 in Neuronal IPC | Telencephalon shrinks [EMX1, RNASEH2A, ORC6, DNA2, IVNS1ABP] (33%->1% / 30%->0%; replicated); cluster 297 in Neuroblast | Telencephalon expands [DDAH2] (0%->41% / 0%->7%; replicated); cluster 228 in Neuronal IPC | Telencephalon expands [AIM2, FRYL, ZNF620] (0%->23% / 0%->32%; replicated); cluster 54 in Radial glia | Telencephalon shrinks [AC092958.1, SMS, AP002026.1, QTRT2] (12%->0% / 9%->0%; replicated); cluster 116 in Neuroblast | Cerebellum shrinks [KIRREL2, NPHS1, PRMT8, CT75, SPSB4] (10%->0% / 23%->0%; replicated); cluster 181 in Radial glia | Telencephalon expands (0%->15% / 0%->8%; replicated); cluster 566 in Neuroblast | Cerebellum shrinks [AC007130.1, IGFBPL1, C15orf41, LNPK, MAP1A] (14%->0% / 42%->0%; replicated); cluster 139 in Radial glia | Telencephalon shrinks [TMTC4] (9%->0% / 2%->0%; replicated); cluster 376 in Neuron | Midbrain expands [MIR124-1HG, ADAMTS20] (0%->4% / 0%->21%; replicated); cluster 35 in Radial glia | Midbrain expands [XPO4, HEY1, AC004470.2, SLIT2, PLOD2] (0%->2% / 0%->20%; replicated) (+115 more).

## Limitations

- Clusters were defined on all ages together; a sub-type that exists only late appears as an expanding cluster, which is the intended reading, but cluster boundaries are not independent of age.
- human_dev regions were dissected differently per donor; class x region parents reduce, but do not remove, that confound.
- 5-15 donors per chemistry: a trend is across that many people.

## What would strengthen this

- B6: Milo neighbourhoods (cortex nhoods) test abundance at finer resolution without fixed cluster boundaries.

## Output files

- `cluster_abundance_per_stratum.csv` -- Per stratum x parent x cluster: Spearman of the cluster's CLR share (within parent) with donor age, exact permutation p
- `cluster_shares_by_donor.csv` -- Per stratum x parent: each cluster's share per donor
- `cluster_abundance_combined.csv` -- Per parent x cluster: v2 x v3 combined age trend of its share within the parent; tier
- `cluster_abundance_cortex.png` -- Sub-types whose within-class share changes with age, cortex
- `cluster_abundance_human_dev.png` -- Sub-types whose within-class share changes with age, human_dev

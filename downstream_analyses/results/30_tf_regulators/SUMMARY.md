# Transcription factors: content, target enrichment and regulon activity in modules and NDD lists

_Generated 2026-10-07 11:40 UTC by `downstream_analyses/30_tf_regulators.py` from `csv_exports/`._

## Question

Which TFs sit in the modules and NDD lists, whose targets are over-represented in them, and is that TF's expression tied to its targets' in these data?

## Inputs

- `annotations/tf/TF_names_v_1.01.txt`
- `annotations/tf/collectri.tsv`
- `cortex__v2/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cluster_ClustersSurprise__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cluster_cluster_id__pseudobulk_counts.csv`
- `results/03_age_trends_within_cell_class/age_trends_combined.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- TFs from Lambert et al. 2018; targets from CollecTRI (OmniPath), signed where known; symbols mapped to each dataset; universe = expressed genes (mean CPM >= 1).
- A: hypergeometric TF content. B: hypergeometric target enrichment for regulons with >= 10 expressed targets, BH within each set.
- C: Spearman of TF expression with its targets' signed mean Z across clusters vs 200 random target sets matched on mean level x spread; v2 x v3 signed Stouffer, BH, tiered.

## Key findings

- **Gene sets rich in TFs** (TFs / expected; q < 0.05): human_dev module:HDw04 11 / 2.3 (DMBX1, EN1, EN2, IRX1, IRX2, IRX3); cortex seed:chromatin_transcription_regulators 8 / 2.4 (ADNP, ASH1L, BCL11A, BCL11B, KDM5B, KMT2A); cortex seed:id_dd_dominant 8 / 2.4 (FOXG1, GATAD2B, KMT2A, PURA, SATB2, SON); cortex module:CXw03 26 / 14.0 (ARX, CREB5, DACH1, GLI3, HES1, NFATC4); human_dev seed:chromatin_transcription_regulators 8 / 2.4 (ADNP, ASH1L, BCL11A, BCL11B, KDM5B, KMT2A); human_dev seed:id_dd_dominant 8 / 2.4 (FOXG1, GATAD2B, KMT2A, PURA, SATB2, SON); cortex module:CX09 8 / 2.5 (EOMES, HES6, INSM1, LHX9, NEUROD4, NEUROG2); cortex module:CXw08 7 / 2.2 (EOMES, HES6, INSM1, LHX9, NEUROD4, NEUROG2); cortex module:CX03 23 / 13.7 (ARX, CREB5, DACH1, GLI3, HES1, HES4); human_dev seed:asd_high_confidence 8 / 3.1 (ADNP, ASH1L, DEAF1, FOXP1, KDM5B, MECP2).
- **cortex: regulons whose TF tracks its targets across clusters** (32 of 221 tested; strongest): NEUROG2 (+0.67/+0.55), NFKB1 (+0.29/+0.39), NFKB2 (+0.34/+0.50), REST (+0.53/+0.52), SP1 (+0.36/+0.65), E2F1 (+0.63/+0.75), JUND (+0.41/+0.46), JUN (+0.37/+0.53), TFDP1 (+0.65/+0.66), FOS (+0.51/+0.51), E2F4 (+0.43/+0.43), E2F3 (+0.78/+0.83).
- **human_dev: regulons whose TF tracks its targets across clusters** (50 of 272 tested; strongest): EGR1 (+0.57/+0.67), HIF1A (+0.56/+0.38), JUN (+0.79/+0.79), NEUROG2 (+0.47/+0.53), SOX11 (+0.79/+0.81), XBP1 (+0.78/+0.69), FOXO3 (+0.27/+0.53), DLX2 (+0.37/+0.51), PAX2 (+0.55/+0.60), RARA (+0.33/+0.38), SRSF2 (+0.88/+0.87), RELA (+0.36/+0.44).
- **Candidate drivers** (TF targets enriched in the set, q < 0.05, and the TF tracks its targets in both donor sets; overlap, fold): cortex seed:chromatin_transcription_regulators: SP1 (8, 4.5x), MYB (3, 13.1x) | cortex seed:epilepsy_dee: REST (5, 62.9x) | cortex seed:synaptic_and_channels: REST (6, 56.6x) | cortex module:CX01: E2F4 (30, 10.7x), E2F1 (24, 5.3x), E2F3 (10, 9.7x), FOXM1 (8, 9.5x), E2F2 (7, 9.3x) | cortex module:CX02: REST (8, 8.2x) | cortex module:CX03: SOX2 (9, 14.6x), SP1 (29, 2.9x), JUN (18, 3.0x), SP3 (12, 3.9x), NFKB2 (17, 2.9x) | cortex module:CX04: E2F4 (29, 14.9x), E2F1 (34, 10.8x), E2F2 (13, 25.0x), E2F3 (12, 16.8x), FOXM1 (8, 13.7x) | cortex module:CX09: NEUROG2 (4, 72.8x) | cortex module:CXw01: REST (12, 11.3x) | cortex module:CXw02: E2F4 (29, 11.7x), E2F1 (25, 6.2x), E2F3 (10, 11.0x), FOXM1 (8, 10.7x), E2F2 (7, 10.6x) | cortex module:CXw03: SP1 (36, 3.5x), SOX2 (9, 14.2x), NFKB2 (22, 3.7x), NFKB1 (22, 3.2x), JUN (20, 3.2x) | cortex module:CXw04: E2F4 (29, 22.3x), E2F1 (31, 14.7x), E2F2 (13, 37.5x), E2F3 (12, 25.2x), TFDP1 (6, 28.9x) | cortex module:CXw08: NEUROG2 (4, 84.0x) | human_dev seed:asd_high_confidence: REST (3, 20.2x) | human_dev seed:chromatin_transcription_regulators: MBD2 (3, 27.6x), SP1 (8, 4.6x) | human_dev seed:epilepsy_dee: REST (6, 73.1x), JUND (5, 9.6x).

## Limitations

- CollecTRI targets come mostly from other tissues and cell lines; a regulon enriched here is a hypothesis about regulation in developing brain, not evidence of it.
- A TF tracking its targets across clusters can reflect shared cell identity (both high in one cell type) rather than regulation; within-class tests would be stricter.
- TF mRNA is a weak proxy for TF activity.

## What would strengthen this

- Within-class regulon activity (as 07's within_class context); motif or ATAC evidence for direct binding.

## Output files

- `tf_content.csv` -- Per gene set: TFs among its expressed genes vs the expressed genome
- `regulon_enrichment.csv` -- Per gene set x TF: CollecTRI targets in the set (q < 0.25, >= 3); with the TF's regulon-activity tier
- `regulon_activity_per_stratum.csv` -- Per stratum x TF: Spearman of TF expression with its targets' signed mean Z across clusters, vs matched random target sets
- `regulon_activity_combined.csv` -- Per TF: v2 x v3 combined; tier; the TF's own replicated age trends (03)

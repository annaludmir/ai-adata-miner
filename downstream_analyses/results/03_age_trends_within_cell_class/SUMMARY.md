# Within-cell-class gene expression across age

_Generated 2026-10-05 06:05 UTC by `downstream_analyses/03_age_trends_within_cell_class.py` from `csv_exports/`._

## Question

Holding cell identity fixed, which genes rise or fall across development, consistently in two independent donor sets?

## Inputs

- `_cross_dataset/gene_id_map.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__detection_fraction.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__group_summary.csv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__detection_fraction.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__group_summary.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v3/11_panels/panel_coverage.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__detection_fraction.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__group_summary.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__detection_fraction.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__group_summary.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v3/11_panels/panel_coverage.csv`

## Method

- Input: per-(cell class, age) pseudobulk counts from 09_pseudobulk, TMM-normalised within each class to log2 CPM (mean log1p(CP10K) drifts with UMIs per cell, which falls with age in cortex). Age points with < 50 cells dropped; classes need >= 5 age points per chemistry.
- Genes tested if detected in >= 10% of cells at >= 2 age points and averaging >= 5 CPM.
- Replicate unit: the age point (one donor, occasionally 2-3 at the same age). Spearman rho vs age; two-sided p exact over all age orderings up to 9 points, Monte Carlo (400,000) above.
- v2 and v3 combined by signed Stouffer (weights sqrt(n ages)), BH per dataset x class. Replicated = same direction, q < 0.05, each chemistry nominal (one-sided p < 0.05).
- Sex-linked genes (chrY, XIST/TSIX) flagged and excluded from gene lists.

## Key findings

- **cortex (primary)**: replicated within-class age trends per class -- Neuroblast 61 up / 84 down of 3937 tested; Neuron 0 up / 2 down of 4558 tested; Neuronal IPC 50 up / 47 down of 4530 tested; Radial glia 179 up / 204 down of 4371 tested.
- **human_dev (region-confounded, see limitations)**: replicated within-class age trends per class -- Erythrocyte 45 up / 31 down of 938 tested; Glioblast 76 up / 68 down of 8239 tested; Immune 72 up / 56 down of 6226 tested; Neuroblast 896 up / 674 down of 7040 tested; Neuron 891 up / 704 down of 6974 tested; Neuronal IPC 745 up / 494 down of 7705 tested; Radial glia 582 up / 532 down of 8032 tested; Vascular 183 up / 84 down of 7952 tested.
- **cortex Neuroblast** -- strongest replicated, rising: NFIC, EIF1AX, MLLT3, PTPRD, EPHA5, USP46, RBBP4, TRIM2, EIF4B, CALCOCO1, SETD2, LINC01560, TNRC6B, CEP170, CCND2; falling: CHST8, REM2, FBXO7, NNAT, EYA2, TMEM163, PPA2, CIAPIN1, PUF60, CREG2, EBF1, CDKN1C, KMT5C, PLPPR2, PEPD.
- **cortex Neuron** -- strongest replicated, rising: none; falling: MTCH1, BAD.
- **cortex Neuronal IPC** -- strongest replicated, rising: SSX2IP, NFIC, DCAF16, FAM171B, SET, AFF3, CCP110, KIF18A, MEGF9, CEP350, CDC42BPA, SCRN1, MMS22L, HELLS, SUZ12; falling: NKAIN4, CD200, TAGLN3, DEDD2, PCBP4, PSMD4, PHLDB2, EIF4A1, BLOC1S4, NOVA2, ZBTB17, CCND1, CCNJL, DDIT4, CHN1.
- **cortex Radial glia** -- strongest replicated, rising: NFIC, SPAG9, EIF2AK2, ITGA2, TSPO, KIF15, DESI2, UBXN2A, PPM1K, FER, SSX2IP, CKAP2, LYRM7, LGALS1, JAK1; falling: HMGA2, CDC23, FBL, CNTNAP2, HMGA1, NCALD, CDK4, OGFOD1, CHAF1B, NDUFS3, MED22, SNU13, BBS4, CENPM, PRMT1.
- **Residual global drift after TMM** (median gene rho beyond +/-0.1) in cortex Radial glia (v2 -0.14, v3 +0.07); human_dev Erythrocyte (v2 -0.54, v3 +0.07). Trends there are partly a whole-transcriptome shift (in erythrocytes, haemoglobin taking over the transcriptome as they mature) -- prefer genes whose |rho| clearly exceeds that offset.
- **NDD-panel genes with replicated trends in cortex**: CNTNAP2 down in Radial glia (asd_high_confidence|synaptic_and_channels); ASH1L up in Radial glia (asd_high_confidence|chromatin_transcription_regulators); CREBBP up in Neuronal IPC (asd_high_confidence|chromatin_transcription_regulators|id_dd_dominant); SLC2A1 down in Radial glia (epilepsy_dee); SLC2A1 down in Neuroblast (epilepsy_dee); CREBBP up in Neuroblast (asd_high_confidence|chromatin_transcription_regulators|id_dd_dominant); KAT6B up in Radial glia (id_dd_dominant). Full list: ndd_genes_with_age_trends.csv.
- **cortex and human_dev agree on within-class trends** (Spearman of combined Z, shared classes): Neuroblast 0.60; Neuron 0.49; Neuronal IPC 0.62; Radial glia 0.61. Same donors, so this measures how much whole-brain pooling and processing change the answer -- not independent replication.

## Limitations

- Each age point is effectively one donor: a trend is a trend across 5-9 people per chemistry. A gene specific to one unusual donor can still pass if that donor sits at an extreme age.
- v2 covers ~6.9-10 pcw and v3 ~5-14 pcw; replicated means monotonic across both windows, so transient (rise-then-fall) programmes are missed by design.
- Within-class means still mix sub-types: a 'Neuron' trend can be a shift in which neuron types are present (e.g. more deep-layer vs upper-layer neurons with age).
- human_dev classes pool regions whose dissection changes with age (see 01); treat its trends as hypotheses, and cortex as the cleaner test.
- Spearman captures monotonic trends; with 5-9 points, small effects are underpowered.

## What would strengthen this

- Repeat at cluster resolution within a class, to separate sub-type shifts from within-cell-type change.
- Add donors at ages each chemistry lacks; with >1 donor per age, model donor as random.
- Check top genes against an external developmental atlas (e.g. BrainSpan).

## Output files

- `age_trends_per_stratum.csv` -- Every expressed gene x class x stratum: Spearman vs age, exact permutation p
- `age_trends_combined.csv` -- v2 and v3 combined per gene x class; tier replicated / supported
- `trend_counts.csv` -- Genes per class by tier and direction; median_gene_rho_* near 0 = no global drift
- `ndd_genes_with_age_trends.csv` -- NDD-panel genes whose within-class age trend is replicated or supported
- `cortex_vs_human_dev_concordance.csv` -- Correlation of combined Z between files, shared classes (same donors)
- `trends_cortex_Neuroblast.png` -- Top replicated age trends in cortex Neuroblast
- `trends_cortex_Neuron.png` -- Top replicated age trends in cortex Neuron
- `trends_cortex_Neuronal_IPC.png` -- Top replicated age trends in cortex Neuronal IPC
- `trends_cortex_Radial_glia.png` -- Top replicated age trends in cortex Radial glia

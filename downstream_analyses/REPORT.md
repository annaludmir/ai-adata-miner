# Step 3 report: downstream analyses

_Built 2026-10-05 06:05 UTC from 8 analyses over `csv_exports/`._

Read this first. Three facts about the data constrain every result below
(details in 01):

1. **Donors are the replicate unit, and each donor has one age.** Within a
   chemistry there are 7-15 donors, mostly one per age. A trend over age is a
   trend across those few people.
2. **v2 and v3 are disjoint donor sets covering different ages.** They are
   analysed separately and used to replicate each other; nothing here pools
   them. "Replicated" means the same direction and nominal significance in
   both, plus a combined FDR < 0.05.
3. **cortex is a subset of the human_dev donors.** Agreement between the two
   files is not independent evidence. cortex is the primary dataset for
   within-cell-type questions; human_dev classes pool brain regions whose
   sampling changes with age.

## Contents

- [01. Data audit: what the exports can support](#01-data-audit-what-the-exports-can-support)
- [02. Cell-class composition across age](#02-cell-class-composition-across-age)
- [03. Within-cell-class gene expression across age](#03-within-cell-class-gene-expression-across-age)
- [04. NDD gene panels across cell classes and age](#04-ndd-gene-panels-across-cell-classes-and-age)
- [05. Cell-class identity across chemistries and files](#05-cell-class-identity-across-chemistries-and-files)
- [06. Gene lists: coverage, overlap, cell-class preference and age coordination](#06-gene-lists-coverage-overlap-cell-class-preference-and-age-coordination)
- [07. Gene-list co-expression: do list genes work as a group?](#07-gene-list-co-expression-do-list-genes-work-as-a-group)
- [08. Co-expression modules across fine clusters](#08-co-expression-modules-across-fine-clusters)

## 01. Data audit: what the exports can support

**Question.** Before any biology: how many independent units stand behind each comparison, which covariates are entangled, and do the labels and tables behave?

- **Donor is the replicate unit, and age is nested in donor.** Every donor has exactly one age (max ages per donor = 1), and 85% of age points rest on a single donor. Donors per stratum: cortex__v2 7, cortex__v3 8, human_dev__v2 15, human_dev__v3 9. An age effect is therefore indistinguishable from a donor effect at that age.
- **cortex: chemistries cover different ages.** v2 only: 7.5|8|8.5|9.2|9.5|10; v3 only: 5|5.5|11.5|12|13.25|14; both: 6.9 (10% of cells). Donors nested in chemistry: True. v2 and v3 are therefore independent donor sets that can replicate each other, but not be pooled as if equivalent.
- **human_dev: chemistries cover different ages.** v2 only: 6.6|6.7|7.5|8|8.1|9.2|9.5|10; v3 only: 7|11.5|12|13|14; both: 6|6.9|8.5 (42% of cells). Donors nested in chemistry: True. v2 and v3 are therefore independent donor sets that can replicate each other, but not be pooled as if equivalent.
- **cortex is not an independent cohort.** 13/15 cortex donors are human_dev donors (5 under a differently written ID, e.g. `XHU:1966:307` = `XHU:307`; `13_cross_dataset_keys` joins them). The other 2 (XDD:348, XDD:400) fall under a human_dev exclusion rule, so they have no human_dev counterpart in these exports. Agreement between the two files is reproducibility of processing, not replication. Age annotations disagree for XDD:359 (13.25 vs 13 pcw).
- **Donor sex is recoverable from expression** (obs['sex'] is 'unknown'): 39/39 donor entries call cleanly.
- **Sex is unevenly spread over age** in cortex__v2 (male ages 6.9|8, female ages 7.5|8.5|9.2|9.5|10; rho = -0.63, exact p = 0.19); cortex__v3 (male ages 6.9|12|14, female ages 5|5.5|6.9|11.5|13.25; rho = 0.45, exact p = 0.29). With this few donors the association is not significant, but it does not need to be to matter: a sex-differential gene can look like an age trend. Later analyses flag sex-linked genes and check trends against sex.
- **Sex-linked genes rank among cell-class markers**, which happens when a class is drawn unevenly from male and female donors: RPS4Y1 in Radial glia (cortex__v2, rank 10); RPS4Y1 in Radial glia (human_dev__v2, rank 34); EIF1AY in Erythrocyte (human_dev__v3, rank 43). See cell_class_sex_balance.csv for the imbalance behind each.
- **Cell-class labels recover known markers**: 20/30 seed marker panels score highest in the class they name. Misses: radial_glia->Glioblast (cortex__v2), neuroblast->Neuron (cortex__v3), radial_glia->Glioblast (cortex__v3), glioblast_opc->Oligo (human_dev__v2), neuron->Placodes (human_dev__v2), oligo->Neural crest (human_dev__v2), radial_glia->Glioblast (human_dev__v2), glioblast_opc->Oligo (human_dev__v3), oligo->Neural crest (human_dev__v3), radial_glia->Glioblast (human_dev__v3). Panels naming a class absent from a stratum are not scored. The misses are neighbouring lineages: glioblasts carry radial-glia genes, placode-derived sensory neurons carry pan-neuronal genes, and the short seed lists cannot separate them.
- **Exclusions (exclusions.csv)**: human_dev age=5.0 (File is human_dev_without_week_5 but still holds 62,786 cells at 5.0 pcw (donor XDD:348)); human_dev age=5.5 (File is human_dev_without_week_5 but still holds 59,667 cells at 5.5 pcw (donor XDD:400)). These exports already omit them.
- **XIST reads ~6.1x higher in human_dev than in cortex for the same female donors**, so the two files probably count reads differently (XIST is largely nuclear and intronic). Compare genes across the files by rank or within-file contrast, not by absolute level.

Method, limitations and output files: [results/01_data_audit/SUMMARY.md](results/01_data_audit/SUMMARY.md)

## 02. Cell-class composition across age

**Question.** Which cell classes expand or shrink with developmental age, consistently in two independent donor sets (v2 and v3 chemistry)?

- **cortex / all cells** (7 + 7 donors), replicated: Radial glia down (rho v2 -0.79, v3 -0.95, q = 0.0025); Neuron up (rho v2 +0.79, v3 +0.88, q = 0.0046). Supported by the combined test but weak in one chemistry: Neuroblast down (rho v2 -0.71, v3 -0.59).
- **human_dev / Cerebellum** (6 + 6 donors): no class trend replicates in both chemistries. Strongest: Glioblast (rho v2 +0.93, v3 +0.64, combined p = 0.0093).
- **human_dev / Midbrain** (8 + 7 donors), replicated: Glioblast up (rho v2 +0.95, v3 +0.99, q = 4.3e-05); Radial glia down (rho v2 -0.82, v3 -0.95, q = 9.9e-04); Neuroblast down (rho v2 -0.69, v3 -0.88, q = 0.0069). Supported by the combined test but weak in one chemistry: Other up (rho v2 +0.53, v3 +0.99).
- **human_dev / Telencephalon** (7 + 6 donors), replicated: Radial glia down (rho v2 -0.81, v3 -0.93, q = 0.012).
- **Sex does not explain the replicated trends.** Where a class tracks donor sex, it does so in one chemistry only, or in opposite directions (e.g. cortex/all cells Radial glia: rho with male +0.79 in v2, -0.29 in v3; cortex/all cells Neuron: rho with male -0.79 in v2, +0.29 in v3; human_dev/Midbrain Neuroblast: rho with male +0.62 in v2, -0.14 in v3), while the age trend keeps its direction in both.
- **Trends that flip between chemistries** (significant in one, opposite sign in the other) -- not replicated, possibly donor- or window-specific: human_dev/Cerebellum Other; human_dev/Cerebellum Neuron; human_dev/Midbrain Neuron.

Method, limitations and output files: [results/02_composition_vs_age/SUMMARY.md](results/02_composition_vs_age/SUMMARY.md)

## 03. Within-cell-class gene expression across age

**Question.** Holding cell identity fixed, which genes rise or fall across development, consistently in two independent donor sets?

- **cortex (primary)**: replicated within-class age trends per class -- Neuroblast 61 up / 84 down of 3937 tested; Neuron 0 up / 2 down of 4558 tested; Neuronal IPC 50 up / 47 down of 4530 tested; Radial glia 179 up / 204 down of 4371 tested.
- **human_dev (region-confounded, see limitations)**: replicated within-class age trends per class -- Erythrocyte 45 up / 31 down of 938 tested; Glioblast 76 up / 68 down of 8239 tested; Immune 72 up / 56 down of 6226 tested; Neuroblast 896 up / 674 down of 7040 tested; Neuron 891 up / 704 down of 6974 tested; Neuronal IPC 745 up / 494 down of 7705 tested; Radial glia 582 up / 532 down of 8032 tested; Vascular 183 up / 84 down of 7952 tested.
- **cortex Neuroblast** -- strongest replicated, rising: NFIC, EIF1AX, MLLT3, PTPRD, EPHA5, USP46, RBBP4, TRIM2, EIF4B, CALCOCO1, SETD2, LINC01560, TNRC6B, CEP170, CCND2; falling: CHST8, REM2, FBXO7, NNAT, EYA2, TMEM163, PPA2, CIAPIN1, PUF60, CREG2, EBF1, CDKN1C, KMT5C, PLPPR2, PEPD.
- **cortex Neuron** -- strongest replicated, rising: none; falling: MTCH1, BAD.
- **cortex Neuronal IPC** -- strongest replicated, rising: SSX2IP, NFIC, DCAF16, FAM171B, SET, AFF3, CCP110, KIF18A, MEGF9, CEP350, CDC42BPA, SCRN1, MMS22L, HELLS, SUZ12; falling: NKAIN4, CD200, TAGLN3, DEDD2, PCBP4, PSMD4, PHLDB2, EIF4A1, BLOC1S4, NOVA2, ZBTB17, CCND1, CCNJL, DDIT4, CHN1.
- **cortex Radial glia** -- strongest replicated, rising: NFIC, SPAG9, EIF2AK2, ITGA2, TSPO, KIF15, DESI2, UBXN2A, PPM1K, FER, SSX2IP, CKAP2, LYRM7, LGALS1, JAK1; falling: HMGA2, CDC23, FBL, CNTNAP2, HMGA1, NCALD, CDK4, OGFOD1, CHAF1B, NDUFS3, MED22, SNU13, BBS4, CENPM, PRMT1.
- **Residual global drift after TMM** (median gene rho beyond +/-0.1) in cortex Radial glia (v2 -0.14, v3 +0.07); human_dev Erythrocyte (v2 -0.54, v3 +0.07). Trends there are partly a whole-transcriptome shift (in erythrocytes, haemoglobin taking over the transcriptome as they mature) -- prefer genes whose |rho| clearly exceeds that offset.
- **NDD-panel genes with replicated trends in cortex**: CNTNAP2 down in Radial glia (asd_high_confidence|synaptic_and_channels); ASH1L up in Radial glia (asd_high_confidence|chromatin_transcription_regulators); CREBBP up in Neuronal IPC (asd_high_confidence|chromatin_transcription_regulators|id_dd_dominant); SLC2A1 down in Radial glia (epilepsy_dee); SLC2A1 down in Neuroblast (epilepsy_dee); CREBBP up in Neuroblast (asd_high_confidence|chromatin_transcription_regulators|id_dd_dominant); KAT6B up in Radial glia (id_dd_dominant). Full list: ndd_genes_with_age_trends.csv.
- **cortex and human_dev agree on within-class trends** (Spearman of combined Z, shared classes): Neuroblast 0.60; Neuron 0.49; Neuronal IPC 0.62; Radial glia 0.61. Same donors, so this measures how much whole-brain pooling and processing change the answer -- not independent replication.

Method, limitations and output files: [results/03_age_trends_within_cell_class/SUMMARY.md](results/03_age_trends_within_cell_class/SUMMARY.md)

## 04. NDD gene panels across cell classes and age

**Question.** Which cell classes preferentially express each neurodevelopmental-disorder panel, consistently across donors and chemistries -- and are NDD genes enriched among genes that change with age?

- **cortex: NDD panels enriched in a cell class, replicated in both donor sets** (effect v2 / v3 vs expression-matched random genes): synaptic_and_channels in Neuroblast (+4.1 / +4.3 null SDs; holds at 100% of age points or more); synaptic_and_channels in Neuron (+6.1 / +5.8 null SDs; holds at 100% of age points or more); asd_high_confidence in Neuron (+3.1 / +4.0 null SDs; holds at 86% of age points or more); chromatin_transcription_regulators in Neuroblast (+3.1 / +2.9 null SDs; holds at 100% of age points or more); asd_high_confidence in Neuroblast (+3.0 / +2.9 null SDs; holds at 100% of age points or more); epilepsy_dee in Neuron (+1.7 / +2.3 null SDs; holds at 86% of age points or more).
- **cortex: NDD panels depleted in a cell class, replicated in both donor sets** (effect v2 / v3 vs expression-matched random genes): synaptic_and_channels in Radial glia (-5.0 / -4.7 null SDs; holds at 100% of age points or more); chromatin_transcription_regulators in Glioblast (-4.7 / -3.7 null SDs; holds at 100% of age points or more); synaptic_and_channels in Neuronal IPC (-3.7 / -3.9 null SDs; holds at 100% of age points or more); synaptic_and_channels in Glioblast (-3.5 / -2.6 null SDs; holds at 100% of age points or more); asd_high_confidence in Glioblast (-3.4 / -2.5 null SDs; holds at 100% of age points or more); id_dd_dominant in Glioblast (-3.2 / -3.2 null SDs; holds at 100% of age points or more) (+4 more).
- **cortex, supported by the combined test only**: chromatin_transcription_regulators enriched in Neuron; asd_high_confidence depleted in Neuronal IPC; id_dd_dominant enriched in Neuroblast; epilepsy_dee enriched in Neuroblast; id_dd_dominant depleted in Radial glia.
- **human_dev: NDD panels enriched in a cell class, replicated in both donor sets** (effect v2 / v3 vs expression-matched random genes): synaptic_and_channels in Oligo (+4.6 / +5.2 null SDs; holds at 100% of age points or more); asd_high_confidence in Neuron (+4.3 / +4.5 null SDs; holds at 100% of age points or more); synaptic_and_channels in Neuroblast (+6.1 / +5.6 null SDs; holds at 100% of age points or more); synaptic_and_channels in Neuron (+7.2 / +7.1 null SDs; holds at 100% of age points or more); chromatin_transcription_regulators in Neuroblast (+4.2 / +4.0 null SDs; holds at 100% of age points or more); chromatin_transcription_regulators in Neuron (+3.4 / +4.2 null SDs; holds at 100% of age points or more) (+6 more).
- **human_dev: NDD panels depleted in a cell class, replicated in both donor sets** (effect v2 / v3 vs expression-matched random genes): synaptic_and_channels in Immune (-5.8 / -6.1 null SDs; holds at 100% of age points or more); synaptic_and_channels in Vascular (-4.1 / -4.1 null SDs; holds at 100% of age points or more); chromatin_transcription_regulators in Erythrocyte (-3.3 / -2.1 null SDs; holds at 100% of age points or more); chromatin_transcription_regulators in Oligo (-3.0 / -2.1 null SDs; holds at 100% of age points or more); epilepsy_dee in Immune (-2.1 / -2.8 null SDs; holds at 100% of age points or more); asd_high_confidence in Oligo (-2.4 / -1.8 null SDs; holds at 100% of age points or more) (+5 more).
- **human_dev, supported by the combined test only**: id_dd_dominant depleted in Erythrocyte; chromatin_transcription_regulators depleted in Fibroblast; epilepsy_dee depleted in Fibroblast; chromatin_transcription_regulators depleted in Neural crest; id_dd_dominant depleted in Neural crest; chromatin_transcription_regulators depleted in Vascular (+2 more).
- **cortex: NDD panels are not over-represented among age-trending genes** once expression level is matched (no panel x class x direction at q < 0.05; strongest epilepsy_dee down in Neuroblast: 1 vs 0.1 expected, p = 0.12).
- **human_dev (region-confounded trends, see 03): NDD panels over-represented among replicated age trends**: id_dd_dominant up in Neuron (12 vs 3.6 expected from expression-matched genes, q = 0.014; ASXL1, FOXG1, GATAD2B, KAT6A, KAT6B, KMT2A, MED13L, PACS1, SATB2, SMARCA2, TCF4, ZEB2).
- **GSEA (17) NDD enrichments significant in both chemistries**: cortex synaptic_and_channels in Neuron (cell_class); cortex epilepsy_dee in Neuron (cell_class); cortex id_dd_dominant in Neuron (cell_class); human_dev synaptic_and_channels in Neuron (cell_class). Across all panels: 12 of 35 enrichments significant in either chemistry hold in both.
- **Per-gene cell-class specificity is chemistry-robust for NDD genes**: the most specific class agrees between v2 and v3 for 82% in cortex, 66% in human_dev (ndd_gene_top_class.csv).

Method, limitations and output files: [results/04_ndd_panel_landscape/SUMMARY.md](results/04_ndd_panel_landscape/SUMMARY.md)

## 05. Cell-class identity across chemistries and files

**Question.** Does each cell class carry the same expression identity in v2 and v3, and in cortex and human_dev -- i.e. can class-level results be compared across strata?

- **cortex__v2 vs cortex__v3** (different chemistry, independent donors): same-class r 0.96-0.97 (median 0.96) vs best other class 0.80 at most; top-100 marker Jaccard median 0.39. Closest neighbours (margin < 0.2): Glioblast (r 0.96 vs Radial glia 0.77); Radial glia (r 0.96 vs Glioblast 0.80).
- **human_dev__v2 vs human_dev__v3** (different chemistry, independent donors): same-class r 0.82-0.98 (median 0.95) vs best other class 0.88 at most; top-100 marker Jaccard median 0.47. Closest neighbours (margin < 0.2): Neuroblast (r 0.98 vs Neuron 0.88); Neuron (r 0.98 vs Neuroblast 0.84).
- **cortex__v2 vs human_dev__v2** (same donors, different file): same-class r 0.90-0.93 (median 0.91) vs best other class 0.79 at most; top-100 marker Jaccard median 0.14. Closest neighbours (margin < 0.2): Glioblast (r 0.90 vs Radial glia 0.79); Neuroblast (r 0.93 vs Neuron 0.79); Neuron (r 0.91 vs Neuroblast 0.72); Radial glia (r 0.93 vs Glioblast 0.77).
- **cortex__v3 vs human_dev__v3** (same donors, different file): same-class r 0.87-0.91 (median 0.89) vs best other class 0.79 at most; top-100 marker Jaccard median 0.09. Closest neighbours (margin < 0.2): Neuroblast (r 0.91 vs Neuron 0.79); Neuron (r 0.89 vs Neuroblast 0.74).
- **Class identity reproduces across chemistries**: a class matches itself first in 16/16 chemistry comparisons, so v2 and v3 results about the same class describe the same cell type. The smallest margins fall between adjacent lineages (radial glia / glioblast, neuroblast / neuron), as expected along a differentiation continuum. Marker lists overlap less than profiles agree, because v2's lower sensitivity reorders the tail of each list.
- **Across files, profiles agree but marker lists do not** (median Jaccard 0.14): human_dev classes span the whole brain, so their top markers include regional genes that cortex cells never express. Compare the files by profile, not by marker list.

Method, limitations and output files: [results/05_identity_reproducibility/SUMMARY.md](results/05_identity_reproducibility/SUMMARY.md)

## 06. Gene lists: coverage, overlap, cell-class preference and age coordination

**Question.** Skipped: no gene lists.

- **No gene lists found in `gene_lists`.** Copy them there or set AIM_GENE_LISTS.

Method, limitations and output files: [results/06_gene_list_landscape/SUMMARY.md](results/06_gene_list_landscape/SUMMARY.md)

## 07. Gene-list co-expression: do list genes work as a group?

**Question.** Skipped: no gene lists.

- **No gene lists found in `gene_lists`.** Copy them there or set AIM_GENE_LISTS.

Method, limitations and output files: [results/07_gene_list_coherence/SUMMARY.md](results/07_gene_list_coherence/SUMMARY.md)

## 08. Co-expression modules across fine clusters

**Question.** Without starting from a list: which groups of genes co-vary across cell clusters robustly enough to be found again in independent donors, where do they peak, how do they change with age, and which gene lists concentrate in them?

- **cortex / across clusters: 11 modules, 7 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - CX01 (225 genes; peak Neuronal IPC; hubs HMGB2, NUSAP1, SMC4, PIMREG, KIFC1)
  - CX02 (221 genes; peak Neuron; hubs MAPT, XPR1, SHTN1, SCG5, SPTAN1; holds seed:epilepsy_dee (6, 7.0x), seed:synaptic_and_channels (7, 3.4x))
  - CX03 (165 genes; peak Glioblast; hubs GULP1, SOX9, GLI3, GATM, CLU)
  - CX04 (156 genes; peak Glioblast/Neuronal IPC; hubs GMNN, PCNA, CENPK, TYMS, HELLS)
  - CX05 (113 genes; peak Neuron; hubs FGF12, GPR22, GAS7, ZBTB38, FXYD7)
  - CX07 (79 genes; peak Neuroblast; hubs GPC2, MLLT11, TAGLN3, DISP3, RASGEF1B)
  - CX09 (30 genes; peak Neuronal IPC; hubs NHLH1, ELAVL2, EOMES, INSM1, NEUROD4; falls with age in Neuronal IPC)
- **cortex / within class: 9 modules, 6 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - CXw01 (241 genes; peak Neuron; hubs MAPT, RAB3A, RTN1, SEZ6L2, TTC9B; holds seed:synaptic_and_channels (10, 3.8x), seed:asd_high_confidence (5, 3.9x))
  - CXw02 (199 genes; peak Neuronal IPC; hubs NUSAP1, HMGB2, CCNA2, GTSE1, TPX2)
  - CXw03 (169 genes; peak Glioblast; hubs SOX2, CLU, GULP1, SOX9, PON2)
  - CXw04 (104 genes; peak Neuronal IPC; hubs TYMS, GMNN, CENPK, PCNA, ORC6)
  - CXw05 (81 genes; peak Neuron; hubs CNTN1, MEF2C, DACT1, SCN2A, VSTM2L)
  - CXw08 (26 genes; peak Neuronal IPC; hubs NHLH1, ELAVL2, NEUROD4, EOMES, HES6; falls with age in Neuronal IPC, rises with age in Radial glia)
- **human_dev / across clusters: 15 modules, 6 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.50)
  - HD01 (272 genes; peak Neuronal IPC, Telencephalon/Forebrain; hubs SMC4, CENPK, MAD2L1, KIF11, CKS1B)
  - HD02 (261 genes; peak Fibroblast/Vascular, Medulla/Hindbrain; hubs ANXA5, SUCLG2, SERPINH1, STK3, PLIN3; rises with age in Vascular)
  - HD03 (225 genes; peak Neuron, Pons; hubs TMEM35A, SV2A, SCN3B, PTPN5, JPH4)
  - HD04 (72 genes; peak Neuroblast, Cerebellum/Pons; hubs TUBB3, CD24, GPC2, MLLT11, SMPD3; falls with age in Neuroblast, falls with age in Neuron)
  - HD05 (43 genes; peak Glioblast, Diencephalon/Hindbrain; hubs NME5, C9orf116, SPATA17, CFAP54, ENKUR; falls with age in Neuronal IPC)
  - HD11 (25 genes; peak Neuron, Diencephalon/Pons; hubs MTUS2, TENM2, TAC1, AC011369.1, PLEKHA6; falls with age in Neuroblast, falls with age in Neuron)
- **human_dev / within class: 6 modules, 4 robust** (re-discovered in each chemistry alone and preserved in the other; cut r >= 0.44)
  - HDw01 (197 genes; peak Neuronal IPC, Telencephalon/Forebrain; hubs TOP2A, CDK1, KIFC1, KIF11, GTSE1; rises with age in Neuronal IPC)
  - HDw02 (48 genes; peak Glioblast, Diencephalon/Hindbrain; hubs NME5, SPATA17, KIF9, MOK, C5orf49; falls with age in Neuronal IPC)
  - HDw03 (29 genes; peak Neuron, Pons/Hindbrain; hubs MAP1LC3A, LY6H, PCSK1N, DDX25, FUCA1)
  - HDw04 (28 genes; peak Glioblast/Radial glia, Pons/Midbrain; hubs IRX2, IRX3, IRX1, C5orf38, IRX5; falls with age in Neuronal IPC, falls with age in Radial glia, falls with age in Neuroblast, falls with age in Neuron)
- **Gene sets concentrated in robust modules** (q < 0.05, >= 2x matched expectation): seed:asd_high_confidence: cortex CXw01 (5 genes, 3.9x) | seed:epilepsy_dee: cortex CX02 (6 genes, 7.0x) | seed:synaptic_and_channels: cortex CXw01 (10 genes, 3.8x); cortex CX02 (7 genes, 3.4x).
- **cortex and human_dev share 2 of 20 cortex modules** (Jaccard >= 0.3 with a human_dev module; same donors, so this is reproducibility of processing and regional pooling, not replication).

Method, limitations and output files: [results/08_coexpression_modules/SUMMARY.md](results/08_coexpression_modules/SUMMARY.md)


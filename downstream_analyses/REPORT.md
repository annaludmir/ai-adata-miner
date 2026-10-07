# Step 3 report: downstream analyses

_Built 2026-10-07 11:40 UTC from 34 analyses over `csv_exports/`._

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
- [09. Cell-cycle programmes: proliferation over development and a phase map of genes](#09-cell-cycle-programmes-proliferation-over-development-and-a-phase-map-of-genes)
- [10. Splicing dynamics: where and when genes are being switched on or off](#10-splicing-dynamics-where-and-when-genes-are-being-switched-on-or-off)
- [11. Gene lists at sub-type resolution: fine clusters and neuron sub-type axes](#11-gene-lists-at-sub-type-resolution-fine-clusters-and-neuron-sub-type-axes)
- [12. Sub-type (cluster) abundance over development, within cell classes](#12-sub-type-cluster-abundance-over-development-within-cell-classes)
- [13. Regional differences in gene lists within cell types (human_dev)](#13-regional-differences-in-gene-lists-within-cell-types-humandev)
- [14. Sex differences in gene lists within cell types (age-adjusted, exploratory)](#14-sex-differences-in-gene-lists-within-cell-types-age-adjusted-exploratory)
- [15. Between-donor variability of NDD genes, against expression-matched genes](#15-between-donor-variability-of-ndd-genes-against-expression-matched-genes)
- [16. Transcriptomic age clock within cell types, validated across donor sets](#16-transcriptomic-age-clock-within-cell-types-validated-across-donor-sets)
- [17. Latent factors (human_dev): identity, gene lists, and age trends within cell types](#17-latent-factors-humandev-identity-gene-lists-and-age-trends-within-cell-types)
- [18. Co-expression rewiring: coherence of gene lists in early vs late clusters](#18-co-expression-rewiring-coherence-of-gene-lists-in-early-vs-late-clusters)
- [19. Robustness checks: gene length, quality metrics, dissociation stress](#19-robustness-checks-gene-length-quality-metrics-dissociation-stress)
- [20. Agreement between CellClass and the files' other cell-type annotations](#20-agreement-between-cellclass-and-the-files-other-cell-type-annotations)
- [21. Age trends within one brain region (human_dev): the region confound removed](#21-age-trends-within-one-brain-region-humandev-the-region-confound-removed)
- [22. Expression measured in each cell-cycle phase: S vs G2/M, cycling vs not, and age within a phase](#22-expression-measured-in-each-cell-cycle-phase-s-vs-g2m-cycling-vs-not-and-age-within-a-phase)
- [23. Gene lists at single-cell level: broad or subset activity, by age and phase](#23-gene-lists-at-single-cell-level-broad-or-subset-activity-by-age-and-phase)
- [24. Within-cell co-expression of gene lists inside cell types](#24-within-cell-co-expression-of-gene-lists-inside-cell-types)
- [25. Gene lists along differentiation pseudotime, and age at matched differentiation](#25-gene-lists-along-differentiation-pseudotime-and-age-at-matched-differentiation)
- [26. Neighbourhood (Milo) abundance over age within cell classes (cortex)](#26-neighbourhood-milo-abundance-over-age-within-cell-classes-cortex)
- [27. Cell states over development: truncated radial glia, neuron sub-types, glial precursors](#27-cell-states-over-development-truncated-radial-glia-neuron-sub-types-glial-precursors)
- [28. Robustness to stricter cell QC: standard vs strict-QC results](#28-robustness-to-stricter-cell-qc-standard-vs-strict-qc-results)
- [29. Functional enrichment of modules, sub-modules, age-trend genes and phase classes](#29-functional-enrichment-of-modules-sub-modules-age-trend-genes-and-phase-classes)
- [30. Transcription factors: content, target enrichment and regulon activity in modules and NDD lists](#30-transcription-factors-content-target-enrichment-and-regulon-activity-in-modules-and-ndd-lists)
- [31. Mutation intolerance (gnomAD LOEUF) and Mendelian disease genes in lists and modules](#31-mutation-intolerance-gnomad-loeuf-and-mendelian-disease-genes-in-lists-and-modules)
- [32. Ligand-receptor signalling between cell types over development](#32-ligand-receptor-signalling-between-cell-types-over-development)
- [33. Gene-list symbols recovered through HGNC, and entries still missing](#33-gene-list-symbols-recovered-through-hgnc-and-entries-still-missing)
- [34. External validation: cortex age trends vs BrainSpan neocortex](#34-external-validation-cortex-age-trends-vs-brainspan-neocortex)

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
- **Exclusions (exclusions.csv)**: human_dev age=5.0 (The 'without_week_5' file still holds cells at 5.0 pcw (donor XDD:348; 62,786 cells)); human_dev age=5.5 (The 'without_week_5' file still holds cells at 5.5 pcw (donor XDD:400; 59,667 cells)). These exports already omit them.
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

Method, limitations and output files: [results/08_coexpression_modules/SUMMARY.md](results/08_coexpression_modules/SUMMARY.md)

## 09. Cell-cycle programmes: proliferation over development and a phase map of genes

**Question.** How does proliferation change over development within progenitor types; which genes and gene groups follow proliferation, and which lean to S or G2/M; and do list-level age trends in progenitors survive once proliferation-linked genes are set aside?

- **cortex: proliferation over development in progenitors** -- fraction cycling falls with age in Radial glia (0.88->0.763 in v2, 0.746->0.6 in v3; replicated).
- **human_dev: no replicated proliferation trend in progenitors.**
- **G1 lengthening not testable yet**: phase fractions per cell class x age are exported by script 03 from this version on; re-run stage 1 to add them.
- **Radial-glia sub-types not available yet**: script 19 (stage 2) exports oRG vs vRG calls; re-run stage 2 to add them.
- **cortex gene phase map** (10,439 genes): 2323 proliferative (251 lean S, 397 lean G2/M), 2342 anti-proliferative, 2292 cycle-independent.
- **human_dev gene phase map** (12,018 genes): 3129 proliferative, 2903 anti-proliferative, 2655 cycle-independent.
- **Sanity check against the seed phase panels** (partly circular: phase calls come from such genes): cortex g2m_phase: lean -16.3/-18.2 null SDs (+ = S), 0 S / 47 G2/M genes; cortex s_phase: lean +12.1/+12.0 null SDs (+ = S), 33 S / 0 G2/M genes.
- **cortex: gene sets tied to proliferation, replicated** -- seed:g2m_phase proliferative, leans G2/M (rho +0.58/+0.65 vs +0.07 matched; 0 S-lean, 47 G2/M-lean genes); seed:s_phase proliferative, leans S (rho +0.75/+0.74 vs +0.04 matched; 33 S-lean, 0 G2/M-lean genes); seed:asd_high_confidence anti-proliferative (rho -0.30/-0.28 vs +0.01 matched); seed:chromatin_transcription_regulators anti-proliferative (rho -0.25/-0.22 vs +0.00 matched; 0 S-lean, 1 G2/M-lean genes); seed:epilepsy_dee anti-proliferative (rho -0.30/-0.24 vs +0.03 matched); seed:id_dd_dominant anti-proliferative (rho -0.14/-0.15 vs +0.00 matched); seed:synaptic_and_channels anti-proliferative (rho -0.53/-0.44 vs +0.04 matched).
- **cortex: proliferative modules with a replicated phase lean** -- CX01 leans G2/M; CX04 leans S; CXw02 leans G2/M; CXw04 leans S.
- **human_dev: gene sets tied to proliferation, replicated** -- seed:g2m_phase proliferative (rho +0.81/+0.78 vs +0.04 matched); seed:s_phase proliferative (rho +0.75/+0.71 vs +0.06 matched); seed:asd_high_confidence anti-proliferative (rho -0.26/-0.26 vs -0.03 matched); seed:epilepsy_dee anti-proliferative (rho -0.24/-0.24 vs +0.03 matched); seed:synaptic_and_channels anti-proliferative (rho -0.48/-0.45 vs -0.04 matched).

Method, limitations and output files: [results/09_cell_cycle_programs/SUMMARY.md](results/09_cell_cycle_programs/SUMMARY.md)

## 10. Splicing dynamics: where and when genes are being switched on or off

**Question.** In which cell classes are genes and gene sets being switched on or off (nascent vs mature RNA), and inside a class, is their unspliced share rising or falling with age alongside their expression?

- **Induction scores reproduce across donor sets** (Spearman v2 vs v3 per class): Glioblast 0.59, Neuroblast 0.57, Neuron 0.71, Neuronal IPC 0.53, Radial glia 0.54. Low values would mean the unspliced signal is noise at this depth.
- **Where gene sets are being switched on or off along the lineage** (replicated; + = more nascent RNA than their mRNA explains) -- seed:s_phase being switched off in Radial glia (-3.2/-5.4 null SDs); seed:s_phase being switched on in Neuronal IPC (+4.5/+2.3 null SDs); seed:g2m_phase being switched off in Radial glia (-4.8/-2.0 null SDs); seed:g2m_phase being switched on in Glioblast (+4.0/+2.0 null SDs); seed:epilepsy_dee being switched on in Neuronal IPC (+1.9/+1.9 null SDs); seed:asd_high_confidence being switched on in Radial glia (+1.7/+1.9 null SDs). Modules: CXw04 being switched on in Neuronal IPC; CXw02 being switched off in Radial glia; CXw01 being switched on in Radial glia; CXw02 being switched on in Glioblast; CXw01 being switched on in Glioblast; CX02 being switched on in Neuronal IPC; CX02 being switched on in Radial glia; CX01 being switched on in Glioblast.
- **Genes by dynamics within classes** (unspliced-share trend and expression trend both replicated, or expression flat): Radial glia: 2 induction ahead of expression, 3 shut-down ahead of expression.
- **Gene sets whose unspliced share changes with age inside a class** (replicated) -- seed:s_phase unspliced share rises in Radial glia (mean rho +0.33/+0.28); seed:g2m_phase unspliced share rises in Neuronal IPC (mean rho +0.24/+0.38); seed:asd_high_confidence unspliced share falls in Neuroblast (mean rho -0.13/-0.16); seed:chromatin_transcription_regulators unspliced share falls in Neuroblast (mean rho -0.11/-0.16). Read with 06/09: a list rising with age whose unspliced share also rises is being actively induced; one whose unspliced share is flat is accumulating mRNA.
- **Level coupling removed**: across groups, a gene's log U/S falls as its expression rises (pooled within-gene slope per log2 CPM: v2 cell_class -0.125, v3 cell_class -0.196, v2 cell_class_x_age +0.010, v3 cell_class_x_age -0.035); scores above are net of it, and genes are compared only where clearly on.
- **Overall unspliced level vs age** (removed before every test above; Spearman): Glioblast +0.54/-0.80, Neuroblast +0.36/+0.71, Neuron +0.39/+1.00, Neuronal IPC +0.36/+0.64, Radial glia +0.43/+0.64.

Method, limitations and output files: [results/10_splicing_dynamics/SUMMARY.md](results/10_splicing_dynamics/SUMMARY.md)

## 11. Gene lists at sub-type resolution: fine clusters and neuron sub-type axes

**Question.** Within broad cell classes, which fine clusters (sub-types) does each gene list concentrate in, consistently in both donor sets; and among neurons, does the list follow the excitatory-inhibitory and deep-upper layer axes?

- **cortex: clusters where each set concentrates** (28 clusters in both chemistries, paired by expression match; effect in null SDs v2/v3): seed:asd_high_confidence: 4 replicated clusters, top cluster 127~2 (Neuron; DLX2, FAM19A2; +3.7/+3.8); cluster 207~151 (Neuron; +2.8/+5.0); cluster 48~112 (Neuron; EBF1, GREM2; +2.7/+3.8) -- over-represented: Neuron (3.5x) | seed:chromatin_transcription_regulators: 3 replicated clusters, top cluster 29~30 (Neuroblast; +3.4/+4.0); cluster 48~112 (Neuron; EBF1, GREM2; +2.0/+2.7); cluster 137~67 (Neuron; +1.7/+2.2) -- over-represented: Neuron (2.3x) | seed:synaptic_and_channels: 9 replicated clusters, top cluster 48~112 (Neuron; EBF1, GREM2; +4.1/+4.6); cluster 127~2 (Neuron; DLX2, FAM19A2; +4.3/+3.4); cluster 207~151 (Neuron; +5.0/+5.2) -- over-represented: Neuron (2.7x) | seed:epilepsy_dee: 5 replicated clusters, top cluster 207~151 (Neuron; +2.7/+3.5); cluster 89~16 (Neuron; +2.3/+2.8); cluster 18~43 (Neuron; +1.9/+2.9) -- over-represented: Neuron (3.5x) | seed:id_dd_dominant: 3 replicated clusters, top cluster 29~30 (Neuroblast; +3.0/+3.6); cluster 105~41 (Radial glia; +2.4/+1.7); cluster 137~67 (Neuron; +1.7/+1.7).
- **human_dev: clusters where each set concentrates** (447 clusters in both chemistries, paired by shared label; effect in null SDs v2/v3): seed:asd_high_confidence: 115 replicated clusters, top cluster 301 (Radial glia, Telencephalon; RAB11FIP3, ZNF785; +5.1/+3.9); cluster 476 (Neuron, Forebrain; +5.5/+4.1); cluster 376 (Neuron, Midbrain; MIR124-1HG, ADAMTS20; +5.8/+5.0) -- over-represented: Neuron (2.1x) | seed:chromatin_transcription_regulators: 112 replicated clusters, top cluster 188 (Neuron, Forebrain; AC092691.1, ZSWIM5; +5.9/+5.4); cluster 350 (Neuron, Forebrain; AC092422.1, RARB; +4.1/+4.0); cluster 352 (Neuron, Forebrain; BCL11B, KLF3-AS1; +5.3/+4.9) -- over-represented: Neuronal IPC (1.8x), Neuron (1.7x) | seed:epilepsy_dee: 120 replicated clusters, top cluster 314 (Neuron, Forebrain; SYN2, GRIA1; +3.5/+3.9); cluster 350 (Neuron, Forebrain; AC092422.1, RARB; +3.3/+4.3); cluster 363 (Neuron, Telencephalon; +3.6/+3.5) -- over-represented: Neuron (2.5x) | seed:id_dd_dominant: 89 replicated clusters, top cluster 177 (Radial glia, Forebrain; ASPM, SGO2; +4.1/+4.0); cluster 359 (Neuron, Forebrain; +5.0/+4.5); cluster 311 (Neuron, Forebrain; SLC24A2, EPHA5; +3.8/+3.4) -- over-represented: Neuronal IPC (2.4x), Immune (1.7x) | seed:synaptic_and_channels: 182 replicated clusters, top cluster 306 (Neuroblast, Telencephalon; OCA2, SH3RF3; +4.5/+3.7); cluster 494 (Neuron, Telencephalon; COL25A1, GABRB2; +3.8/+4.5); cluster 492 (Neuron, Diencephalon; +4.6/+4.7) -- over-represented: Neuron (2.2x).
- **Cluster profiles reproduce across donor sets** (Spearman of set scores v2 vs v3 over paired clusters): cortex median 0.83 (range 0.62-0.94); human_dev median 0.89 (range 0.78-0.97).
- **cortex, deep vs upper layer** (38/31 neuron clusters v2/v3; marker groups correlate rho [-0.35, -0.36]): seed:epilepsy_dee toward deep layer (rho +0.53/+0.71; replicated); seed:synaptic_and_channels toward deep layer (rho +0.56/+0.60; replicated).
- **human_dev, excitatory vs inhibitory** (212/168 neuron clusters v2/v3; marker groups correlate rho [-0.61, -0.67]): no set leans to either pole in both donor sets.
- **human_dev, deep vs upper layer** (15/15 neuron clusters v2/v3; marker groups correlate rho [-0.42, 0.36]): no set leans to either pole in both donor sets.

Method, limitations and output files: [results/11_list_subtype_mapping/SUMMARY.md](results/11_list_subtype_mapping/SUMMARY.md)

## 12. Sub-type (cluster) abundance over development, within cell classes

**Question.** Within each cell class, which sub-types (clusters) expand or shrink with age, consistently in both donor sets?

- **cortex: sub-types whose share within their class changes with age** (2 of 9 clusters tiered; share youngest -> oldest donor, v2 / v3): cluster 2 in Neuroblast expands [PPP1R17, EPHB6, SNCB, NEUROD2, PRKX] (25%->74% / 1%->81%; replicated); cluster 9 in Neuroblast shrinks [LHX1, MAB21L1, RELN, TP73, PGF] (22%->0% / 20%->0%; supported).
- **human_dev (parent = class x region): sub-types whose share within their class changes with age** (131 of 293 clusters tiered; share youngest -> oldest donor, v2 / v3): cluster 227 in Radial glia | Telencephalon expands [STK17A, ESR2, CBFA2T2, SYNE2, ZC3H12C] (0%->12% / 0%->26%; replicated); cluster 488 in Neuron | Telencephalon shrinks [LEMD1, SRP14] (12%->0% / 5%->0%; replicated); cluster 347 in Neuron | Telencephalon shrinks [PLEKHA8, RFX7, STX16-NPEPL1] (5%->0% / 4%->0%; replicated); cluster 469 in Neuroblast | Telencephalon shrinks [KCNC2, ISLR2, FAM155A, RBM19, CASTOR3] (78%->0% / 40%->0%; replicated); cluster 357 in Neuron | Telencephalon shrinks [AC004943.2, RAPGEF2] (6%->0% / 4%->0%; replicated); cluster 573 in Neuroblast | Midbrain expands [EBF2, LHX2] (0%->45% / 0%->81%; replicated); cluster 225 in Neuronal IPC | Telencephalon shrinks [EMX1, RNASEH2A, ORC6, DNA2, IVNS1ABP] (33%->1% / 30%->0%; replicated); cluster 297 in Neuroblast | Telencephalon expands [DDAH2] (0%->41% / 0%->7%; replicated); cluster 228 in Neuronal IPC | Telencephalon expands [AIM2, FRYL, ZNF620] (0%->23% / 0%->32%; replicated); cluster 54 in Radial glia | Telencephalon shrinks [AC092958.1, SMS, AP002026.1, QTRT2] (12%->0% / 9%->0%; replicated); cluster 116 in Neuroblast | Cerebellum shrinks [KIRREL2, NPHS1, PRMT8, CT75, SPSB4] (10%->0% / 23%->0%; replicated); cluster 181 in Radial glia | Telencephalon expands (0%->15% / 0%->8%; replicated); cluster 566 in Neuroblast | Cerebellum shrinks [AC007130.1, IGFBPL1, C15orf41, LNPK, MAP1A] (14%->0% / 42%->0%; replicated); cluster 139 in Radial glia | Telencephalon shrinks [TMTC4] (9%->0% / 2%->0%; replicated); cluster 376 in Neuron | Midbrain expands [MIR124-1HG, ADAMTS20] (0%->4% / 0%->21%; replicated); cluster 35 in Radial glia | Midbrain expands [XPO4, HEY1, AC004470.2, SLIT2, PLOD2] (0%->2% / 0%->20%; replicated) (+115 more).

Method, limitations and output files: [results/12_cluster_abundance_vs_age/SUMMARY.md](results/12_cluster_abundance_vs_age/SUMMARY.md)

## 13. Regional differences in gene lists within cell types (human_dev)

**Question.** Within a cell type, do gene lists and modules differ between brain regions, in both donor sets, and could the difference be age?

- **Tested**: 10 cell classes with >= 3 regions, 15 gene sets, 912 class x region x set tests in both chemistries; 507 tiered (356 replicated), 383 of them in neural classes (Glioblast, Neuroblast, Neuron, Neuronal IPC, Oligo, Radial glia).
- **Gene lists: no region difference within a class replicates.**
- **Seed NDD panels that differ between regions within a neural class** (strongest three per set; effect in null SDs v2/v3; region age offset vs the class's other regions): seed:chromatin_transcription_regulators (30 tiered): higher in Glioblast of Telencephalon (+3.9/+3.9; age +1.2/+0.9 wk), lower in Neuroblast of Medulla (-3.9/-3.6; age -0.8/-0.5 wk), lower in Neuroblast of Pons (-5.1/-3.6; age -0.1/-0.5 wk) | seed:id_dd_dominant (27 tiered): lower in Glioblast of Medulla (-4.8/-3.9; age -0.8/-1.0 wk), lower in Neuroblast of Medulla (-4.3/-4.0; age -0.8/-0.5 wk), higher in Neuroblast of Telencephalon (+4.6/+3.9; age +1.2/+1.3 wk) | seed:asd_high_confidence (27 tiered): higher in Neuroblast of Forebrain (+5.6/+5.1; age +1.0/+2.6 wk), lower in Neuroblast of Medulla (-3.7/-5.3; age -0.8/-0.5 wk), lower in Neuroblast of Pons (-3.7/-4.1; age -0.1/-0.5 wk) | seed:synaptic_and_channels (19 tiered): higher in Neuroblast of Forebrain (+4.0/+5.2; age +1.0/+2.6 wk), lower in Radial glia of Cerebellum (-2.8/-3.8; age -0.7/+0.5 wk), higher in Radial glia of Forebrain (+2.2/+3.9; age +1.0/+2.6 wk) | seed:epilepsy_dee (16 tiered): higher in Neuroblast of Forebrain (+2.8/+3.5; age +1.0/+2.6 wk), lower in Oligo of Medulla (-2.5/-4.0; age -1.3/-1.0 wk), higher in Oligo of Forebrain (+2.5/+2.7; age +0.7/+2.2 wk).
- **Modules that differ between regions within a neural class** (strongest three per set; effect in null SDs v2/v3; region age offset vs the class's other regions): module:HD01 (27 tiered): lower in Glioblast of Diencephalon (-13.3/-7.4; age -0.3/-1.2 wk), higher in Neuronal IPC of Forebrain (+11.3/+9.2; age +1.0/+2.6 wk), higher in Radial glia of Telencephalon (+10.0/+17.7; age +1.2/+1.3 wk) | module:HD02 (26 tiered): higher in Glioblast of Pons (+5.6/+6.9; age -0.1/-0.9 wk), higher in Neuronal IPC of Cerebellum (+6.8/+8.6; age -0.7/+0.5 wk), higher in Radial glia of Pons (+3.7/+5.2; age -0.1/-0.5 wk) | module:HD03 (25 tiered): lower in Glioblast of Telencephalon (-5.2/-6.0; age +1.2/+0.9 wk), lower in Neuroblast of Cerebellum (-14.2/-7.8; age -0.7/+0.5 wk), higher in Neuroblast of Diencephalon (+7.8/+10.3; age -0.3/-0.7 wk) | module:HD05 (24 tiered): higher in Glioblast of Diencephalon (+6.8/+8.0; age -0.3/-1.2 wk), higher in Glioblast of Pons (+4.8/+5.3; age -0.1/-0.9 wk), lower in Glioblast of Telencephalon (-8.7/-9.8; age +1.2/+0.9 wk) | module:HD11 (24 tiered): lower in Glioblast of Telencephalon (-4.6/-5.0; age +1.2/+0.9 wk), higher in Neuroblast of Diencephalon (+4.6/+4.6; age -0.3/-0.7 wk), lower in Neuroblast of Telencephalon (-7.2/-9.2; age +1.2/+1.3 wk) | module:HDw01 (28 tiered): lower in Glioblast of Diencephalon (-14.0/-7.9; age -0.3/-1.2 wk), higher in Neuronal IPC of Forebrain (+12.0/+9.4; age +1.0/+2.6 wk), higher in Radial glia of Telencephalon (+8.2/+14.9; age +1.2/+1.3 wk) | module:HDw02 (28 tiered): higher in Glioblast of Diencephalon (+6.9/+7.9; age -0.3/-1.2 wk), higher in Glioblast of Midbrain (+6.5/+4.5; age -0.4/-0.9 wk), higher in Glioblast of Pons (+5.2/+5.7; age -0.1/-0.9 wk) | module:HDw03 (26 tiered): lower in Glioblast of Forebrain (-4.3/-4.0; age +1.0/+2.2 wk), lower in Glioblast of Telencephalon (-4.7/-5.1; age +1.2/+0.9 wk), higher in Neuroblast of Pons (+6.2/+5.3; age -0.1/-0.5 wk) | module:HDw04 (38 tiered): higher in Glioblast of Cerebellum (+3.8/+4.6; age -0.7/+0.1 wk), lower in Neuron of Telencephalon (-5.6/-7.6; age +1.2/+1.3 wk), lower in Radial glia of Telencephalon (-7.7/-10.2; age +1.2/+1.3 wk) | module:HD04 (18 tiered): lower in Radial glia of Telencephalon (-4.3/-3.8; age +1.2/+1.3 wk), lower in Neuron of Telencephalon (-3.7/-4.1; age +1.2/+1.3 wk), higher in Glioblast of Cerebellum (+3.4/+3.6; age -0.7/+0.1 wk).
- **Non-neural classes, likely ambient RNA** (not interpreted): Erythrocyte 5, Fibroblast 14, Immune 14, Vascular 16 tiered list / panel results (47% 'higher'); neural gene sets in these cells mostly measure ambient neural RNA, so they are listed in region_contrast_combined.csv (neural_class = False) but not interpreted.
- **Age check**: no set's region effects track region age (|rho| >= 0.5 in both chemistries), so the differences above are not simply age.

Method, limitations and output files: [results/13_regional_differences/SUMMARY.md](results/13_regional_differences/SUMMARY.md)

## 14. Sex differences in gene lists within cell types (age-adjusted, exploratory)

**Question.** Within a cell type and at the same age, are NDD gene lists expressed differently in male and female donors?

- **Design** (age points male/female, pooled over chemistries): cortex: Glioblast 3/7, Neuroblast 4/8, Neuron 4/7, Neuronal IPC 4/8, Radial glia 4/9; human_dev: Erythrocyte 5/9, Glioblast 5/10, Immune 6/7, Neuroblast 6/10, Neuron 6/10, Neuronal IPC 6/10, Oligo 3/5, Radial glia 6/10, Vascular 6/9. Sex vs age: cortex rho +0.10, human_dev rho -0.11 (median over classes; far from 0 means sex and age are entangled). Possible sex relabellings per class: cortex Glioblast 36, cortex Neuroblast 210, cortex Neuron 126, cortex Neuronal IPC 210, cortex Radial glia 315, human_dev Erythrocyte 560, human_dev Glioblast 720, human_dev Immune 315, human_dev Neuroblast 1260, human_dev Neuron 1260, human_dev Neuronal IPC 1260, human_dev Oligo 10, human_dev Radial glia 1260, human_dev Vascular 980 -- with fewer than 20, the label-shuffle p cannot reach 0.05.
- **Positive control** (best rank of a Y gene or XIST among all genes, per class): cortex Glioblast #1; cortex Neuroblast #1; cortex Neuron #1; cortex Neuronal IPC #1; cortex Radial glia #1; human_dev Glioblast #1; human_dev Immune #1; human_dev Neuroblast #1; human_dev Neuron #1; human_dev Neuronal IPC #1; human_dev Oligo #1; human_dev Radial glia #1; human_dev Vascular #1. Rank 1-3 means the design detects real sex differences.
- **Autosomal / X genes differing by sex at q < 0.05**: cortex Neuroblast: KDM5C; cortex Neuronal IPC: PABPC1, LINC01551; cortex Radial glia: AGAP6; human_dev Glioblast: JPX; human_dev Neuron: KDM5C, JPX; human_dev Neuronal IPC: JPX.
- **Gene sets differing by sex** (passing both nulls; effect vs random genes in null SDs; label-shuffle p; within-chemistry differences v2/v3): none at q < 0.05. 24 set x class tests beat random genes (p < 0.01) but not shuffled sex labels (p >= 0.05): there the donors differ, but not by sex.

Method, limitations and output files: [results/14_sex_differences/SUMMARY.md](results/14_sex_differences/SUMMARY.md)

## 15. Between-donor variability of NDD genes, against expression-matched genes

**Question.** After removing the age trend, do NDD gene lists vary less (or more) between donors than genes of the same expression level, within a cell class and in both donor sets?

- **cortex: sets whose between-donor variability differs from matched genes** (0.5 = like matched genes): more variable: module:CX03 in Radial glia (mean percentile 0.64/0.69; replicated); module:CXw03 in Radial glia (mean percentile 0.65/0.71; replicated); module:CXw08 in Radial glia (mean percentile 0.78/0.81; replicated); module:CX02 in Neuroblast (mean percentile 0.65/0.81; replicated); module:CX03 in Neuroblast (mean percentile 0.67/0.65; replicated); module:CX05 in Neuroblast (mean percentile 0.71/0.75; replicated); module:CX07 in Neuroblast (mean percentile 0.70/0.74; replicated); module:CX09 in Neuroblast (mean percentile 0.74/0.82; replicated) (+41 more).
- **human_dev: sets whose between-donor variability differs from matched genes** (0.5 = like matched genes): less variable: seed:chromatin_transcription_regulators in Neuroblast (mean percentile 0.40/0.40; replicated); seed:chromatin_transcription_regulators in Radial glia (mean percentile 0.46/0.35; supported); seed:chromatin_transcription_regulators in Vascular (mean percentile 0.44/0.38; supported); seed:chromatin_transcription_regulators in Neuronal IPC (mean percentile 0.43/0.40; supported) | more variable: module:HD01 in Vascular (mean percentile 0.69/0.64; replicated); module:HD02 in Vascular (mean percentile 0.61/0.59; replicated); module:HDw01 in Vascular (mean percentile 0.74/0.69; replicated); module:HD01 in Immune (mean percentile 0.64/0.66; replicated); module:HD02 in Immune (mean percentile 0.66/0.61; replicated); module:HDw01 in Immune (mean percentile 0.66/0.69; replicated); module:HD02 in Neuroblast (mean percentile 0.73/0.71; replicated); module:HD03 in Neuroblast (mean percentile 0.70/0.65; replicated) (+65 more).
- **Consistently tight genes** (tightest fifth of their expression bin in both chemistries): cortex Neuroblast 490; cortex Neuron 476; cortex Neuronal IPC 451; cortex Radial glia 450; human_dev Erythrocyte 96; human_dev Glioblast 732; human_dev Immune 508; human_dev Neuroblast 898; human_dev Neuron 875; human_dev Neuronal IPC 809; human_dev Radial glia 821; human_dev Vascular 588. See consistently_tight_genes.csv.

Method, limitations and output files: [results/15_donor_variability/SUMMARY.md](results/15_donor_variability/SUMMARY.md)

## 16. Transcriptomic age clock within cell types, validated across donor sets

**Question.** Can a cell type's developmental age be predicted from expression by a model trained on other donors; which genes carry the signal; and do gene lists predict age better than random genes of the same expression?

- **cortex: age is predictable across donor sets** -- train v2 -> test v3 rho +0.93 (shuffle p = 0.005), v3 -> v2 +0.78 (p = 0.025); per class (v2->v3 / v3->v2): Glioblast +0.80/+0.77, Neuroblast +0.94/+0.79, Neuron +0.90/+0.79, Neuronal IPC +0.94/+0.82, Radial glia +0.93/+0.79.
- **cortex: consensus clock genes** (5 in the top 200 of both models, same sign) -- rising with age: RDX, SCYL3, C2orf40; falling: NEFM, FAM207A.
- **cortex: gene sets that predict age better than matched random genes** (mean accuracy vs random): none at q < 0.05. Weakest relative to random: seed:id_dd_dominant +0.72 vs +0.74, seed:asd_high_confidence +0.77 vs +0.73, seed:s_phase +0.81 vs +0.78.
- **human_dev: age is predictable across donor sets** -- train v2 -> test v3 rho +0.93 (shuffle p = 0.005), v3 -> v2 +0.91 (p = 0.005); per class (v2->v3 / v3->v2): Erythrocyte +0.93/+0.92, Fibroblast +0.80/+1.00, Glioblast +0.96/+0.90, Immune +0.86/+0.95, Neuroblast +1.00/+0.89, Neuron +0.98/+0.91, Neuronal IPC +0.98/+0.89, Oligo +0.70/+0.80, Radial glia +0.98/+0.89, Vascular +0.90/+0.95.
- **human_dev: consensus clock genes** (21 in the top 200 of both models, same sign) -- rising with age: IFITM2, RARRES2, LGALS1, HBB, AC104389.6, ITM2A, C9orf24, APOE, DPP7, AC131571.1, EIF2AK2, APPL2; falling: HBE1, LIN28B, NR2F2, AC090204.1.
- **human_dev: gene sets that predict age better than matched random genes** (mean accuracy vs random): none at q < 0.05. Weakest relative to random: seed:g2m_phase +0.52 vs +0.72, seed:s_phase +0.64 vs +0.68, seed:synaptic_and_channels +0.63 vs +0.59.

Method, limitations and output files: [results/16_age_clock/SUMMARY.md](results/16_age_clock/SUMMARY.md)

## 17. Latent factors (human_dev): identity, gene lists, and age trends within cell types

**Question.** Which cell types and regions does each latent factor mark, which gene lists concentrate among its top genes, and which factors change with age beyond the shifting cell-type mix?

- **Factor identity**: 35 of 50 factors' top genes peak in the same cell class in both chemistries (median activity class-profile r 0.94); 1 flagged as likely technical and 1 as tracking the cell cycle by 08. Class where each factor's top genes peak (v2): Erythrocyte 1, Fibroblast 6, Glioblast 2, Immune 2, Neuroblast 2, Neuron 29, Oligo 4, Radial glia 1, Vascular 3.
- **No gene set concentrates among any factor's top genes** beyond expression-matched chance.
- **Factors changing with age (raw)**: F5 (Glioblast genes) rises with age (rho +0.74/+0.85; replicated); F10 (Glioblast genes) rises with age (rho +0.62/+0.91; replicated); F23 (Neuron genes) falls with age (rho -0.75/-0.40; supported); F12 (Radial glia genes) falls with age (rho -0.63/-0.67; replicated); F6 (Neuron genes) rises with age (rho +0.59/+0.63; replicated); F28 (Neuron genes) falls with age (rho -0.51/-0.74; replicated); F17 (Oligo genes) rises with age (rho +0.21/+0.95; supported); F4 (Neuron genes) [cell cycle] rises with age (rho +0.63/+0.43; supported).
- **Factors changing with age (composition-adjusted)**: F5 (Glioblast genes) rises with age (rho +0.75/+0.83; replicated); F4 (Neuron genes) [cell cycle] rises with age (rho +0.88/+0.33; supported); F28 (Neuron genes) falls with age (rho -0.63/-0.74; replicated); F12 (Radial glia genes) falls with age (rho -0.56/-0.66; replicated); F17 (Oligo genes) rises with age (rho +0.20/+0.95; supported); F6 (Neuron genes) rises with age (rho +0.57/+0.63; replicated); F7 (Vascular genes) falls with age (rho -0.54/-0.63; replicated); F10 (Glioblast genes) rises with age (rho +0.55/+0.61; replicated) (+2 more).
- **Explained by cell-type mix** (age trend in raw activity that vanishes after composition adjustment): none.

Method, limitations and output files: [results/17_latent_factors/SUMMARY.md](results/17_latent_factors/SUMMARY.md)

## 18. Co-expression rewiring: coherence of gene lists in early vs late clusters

**Question.** Do gene lists and modules co-express more or less tightly in late than in early clusters, within the same mix of cell classes, in both donor sets?

- **Age split** (median cluster age early -> late, per class): cortex v2: Glioblast 9.2->10, Neuroblast 9.2->10, Neuron 9.2->10, Neuronal IPC 9.2->9.75, Radial glia 9.2->9.5; cortex v3: Glioblast 11.5->12, Neuroblast 11.5->12, Neuron 11.5->11.75, Neuronal IPC 11.5->12, Radial glia 11.5->12; human_dev v2: Glioblast 8->9.2, Neuroblast 6.7->8.05, Neuron 6.9->8, Neuronal IPC 7.2->9.2, Radial glia 6.7->8, Vascular 6.9->9.2; human_dev v3: Glioblast 12->13, Neuroblast 6.9->11.5, Neuron 6.9->8.5, Neuronal IPC 6.9->12, Radial glia 6.9->6.9, Vascular 11.75->12.5.
- **cortex, across clusters: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:chromatin_transcription_regulators more coherent late (0.15->0.18 / 0.16->0.21; replicated); seed:asd_high_confidence less coherent late (0.14->0.12 / 0.16->0.15; replicated).
- **cortex, within class: sets whose coherence changes with age** (coherence early -> late, v2 / v3): module:CX07 more coherent late (0.43->0.46 / 0.38->0.51; replicated); seed:chromatin_transcription_regulators more coherent late (0.10->0.12 / 0.16->0.21; replicated).
- **human_dev, across clusters: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:id_dd_dominant more coherent late (0.14->0.21 / 0.11->0.23; replicated); module:HD02 less coherent late (0.70->0.67 / 0.72->0.68; replicated); module:HD03 less coherent late (0.77->0.73 / 0.77->0.71; replicated); module:HD04 less coherent late (0.68->0.63 / 0.68->0.65; replicated); module:HDw02 more coherent late (0.47->0.51 / 0.52->0.59; replicated); seed:synaptic_and_channels less coherent late (0.48->0.41 / 0.48->0.43; replicated); module:HDw03 less coherent late (0.71->0.68 / 0.66->0.63; replicated); module:HDw04 more coherent late (0.46->0.50 / 0.48->0.57; replicated) (+3 more).
- **human_dev, within class: sets whose coherence changes with age** (coherence early -> late, v2 / v3): seed:g2m_phase more coherent late (0.59->0.67 / 0.64->0.68; replicated); seed:chromatin_transcription_regulators more coherent late (0.26->0.32 / 0.21->0.30; replicated); seed:id_dd_dominant more coherent late (0.17->0.23 / 0.12->0.24; replicated); module:HD03 less coherent late (0.29->0.26 / 0.31->0.27; replicated); module:HD11 more coherent late (0.21->0.46 / 0.22->0.47; replicated); module:HDw02 more coherent late (0.42->0.50 / 0.50->0.57; replicated); module:HDw04 more coherent late (0.45->0.54 / 0.50->0.61; replicated); module:HD02 less coherent late (0.18->0.17 / 0.23->0.20; replicated) (+4 more).

Method, limitations and output files: [results/18_coexpression_rewiring/SUMMARY.md](results/18_coexpression_rewiring/SUMMARY.md)

## 19. Robustness checks: gene length, quality metrics, dissociation stress

**Question.** Are the list-level results explained by gene length, by quality metrics, or by dissociation stress?

- **Gene length** (median length percentile among expressed genes; 0.5 = typical): cortex: seed:synaptic_and_channels 0.84, seed:epilepsy_dee 0.79, seed:asd_high_confidence 0.79, seed:chromatin_transcription_regulators 0.78, seed:id_dd_dominant 0.76; human_dev: seed:synaptic_and_channels 0.80, seed:epilepsy_dee 0.74, seed:asd_high_confidence 0.74, seed:chromatin_transcription_regulators 0.72, seed:id_dd_dominant 0.71.
- **Length check, class preference**: 52 tiered results with expression-matched nulls; 7 lost when length is matched too (cortex seed:epilepsy_dee in Neuroblast; cortex seed:id_dd_dominant in Neuroblast; cortex seed:id_dd_dominant in Radial glia; human_dev seed:chromatin_transcription_regulators in Neural crest; human_dev seed:chromatin_transcription_regulators in Vascular; human_dev seed:id_dd_dominant in Neural crest; human_dev seed:id_dd_dominant in Neuronal IPC), 1 gained.
- **Length check, age coordination**: 31 tiered results with expression-matched nulls; 4 lost when length is matched too (human_dev seed:asd_high_confidence in Glioblast; human_dev seed:epilepsy_dee in Neuroblast; human_dev seed:synaptic_and_channels in Neuronal IPC; human_dev seed:asd_high_confidence in Vascular), 0 gained.
- **Length check, coherence**: 10 tiered results with expression-matched nulls; 0 lost when length is matched too, 0 gained.
- **Modules tracking quality metrics within cell classes** (|rho| >= 0.5 in both chemistries): none.
- **Dissociation-stress score vs age** (rho v2/v3): human_dev Neuron -0.88/-0.45 (falls with age, supported); human_dev Neuroblast -0.79/-0.40 (falls with age, supported); human_dev Neuronal IPC -0.73/-0.19; cortex Neuronal IPC -0.61/-0.71; cortex Neuroblast -0.71/-0.37; human_dev Erythrocyte -0.52/+0.36; human_dev Glioblast -0.04/+0.46; human_dev Radial glia -0.23/+0.02; cortex Radial glia -0.04/-0.29; human_dev Vascular -0.10/+0.33; human_dev Immune -0.02/-0.14; cortex Neuron -0.57/+0.90.
- **List age trends that shrink by half or more when stress is controlled** (classes whose stress score changes with age; mean rho plain -> partial): none.

Method, limitations and output files: [results/19_robustness_checks/SUMMARY.md](results/19_robustness_checks/SUMMARY.md)

## 20. Agreement between CellClass and the files' other cell-type annotations

**Question.** Skipped: the cross-tabulations are not exported yet.

- **Not run**: these tables are written by stage 1 from this version of the repo on (config role alt_cell_class, and cell_type_id added to 02_composition's groupings). Re-run stage 1 (slurm_01_metadata.sh), then this analysis.

Method, limitations and output files: [results/20_annotation_agreement/SUMMARY.md](results/20_annotation_agreement/SUMMARY.md)

## 21. Age trends within one brain region (human_dev): the region confound removed

**Question.** Skipped: no class x region x age pseudobulk yet.

- **Not run**: script 09 writes this grouping from this version of the repo on; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis.

Method, limitations and output files: [results/21_region_age_trends/SUMMARY.md](results/21_region_age_trends/SUMMARY.md)

## 22. Expression measured in each cell-cycle phase: S vs G2/M, cycling vs not, and age within a phase

**Question.** Skipped: no class x phase pseudobulk yet.

- **Not run**: script 09 writes the phase groupings from this version of the repo on; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis.

Method, limitations and output files: [results/22_phase_resolved_expression/SUMMARY.md](results/22_phase_resolved_expression/SUMMARY.md)

## 23. Gene lists at single-cell level: broad or subset activity, by age and phase

**Question.** Skipped: no per-cell programme scores yet.

- **Not run**: script 22 is new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis.

Method, limitations and output files: [results/23_cell_level_activity/SUMMARY.md](results/23_cell_level_activity/SUMMARY.md)

## 24. Within-cell co-expression of gene lists inside cell types

**Question.** Skipped: no within-cell co-expression yet.

- **Not run**: script 22 is new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis.

Method, limitations and output files: [results/24_within_cell_coexpression/SUMMARY.md](results/24_within_cell_coexpression/SUMMARY.md)

## 25. Gene lists along differentiation pseudotime, and age at matched differentiation

**Question.** Skipped: no pseudotime pseudobulks yet.

- **Not run**: scripts 20 and 09's pseudotime groupings are new; re-run stages 1 and 2, then this analysis.

Method, limitations and output files: [results/25_pseudotime_programs/SUMMARY.md](results/25_pseudotime_programs/SUMMARY.md)

## 26. Neighbourhood (Milo) abundance over age within cell classes (cortex)

**Question.** Skipped: no Milo neighbourhood counts.

- **Not run**: script 21 is new; re-run stage 1 (slurm_01_metadata.sh), then this analysis.

Method, limitations and output files: [results/26_milo_abundance/SUMMARY.md](results/26_milo_abundance/SUMMARY.md)

## 27. Cell states over development: truncated radial glia, neuron sub-types, glial precursors

**Question.** Skipped: no tRG flags or state scores yet.

- **Not run**: the tRG flag (script 19) and state scores (script 22) are new; re-run stage 2 (slurm_02_pseudobulk.sh), then this analysis.

Method, limitations and output files: [results/27_cell_states/SUMMARY.md](results/27_cell_states/SUMMARY.md)

## 28. Robustness to stricter cell QC: standard vs strict-QC results

**Question.** Skipped: no strict-QC results at /nonexistent.

- **Not run**: no strict-QC results at /nonexistent. To produce them: `AIM_STRICT=true ./submit_all.sh` (exclusions_strict.csv; exports to csv_exports_strict/, results to aim_downstream_strict/), then rerun stage 4 here.

Method, limitations and output files: [results/28_strict_qc_comparison/SUMMARY.md](results/28_strict_qc_comparison/SUMMARY.md)

## 29. Functional enrichment of modules, sub-modules, age-trend genes and phase classes

**Question.** Which biological processes, components, functions and pathways do the modules, list sub-modules, age-trend genes, proliferation classes and user lists carry?

- **Libraries**: GO_Biological_Process_2023 (5,407 terms), GO_Cellular_Component_2023 (474 terms), GO_Molecular_Function_2023 (1,147 terms), KEGG_2021_Human (320 terms), Reactome_2022 (1,818 terms); query groups tested: 54.
- **Modules -- top terms** (overlap, fold, q): cortex CX01 (G2/M phase): Cell Cycle, Mitotic R-HSA-69278 (80, 9.1x, q 1.7e-52); M Phase R-HSA-68886 (63, 9.9x, q 2.2e-42); Mitotic Prometaphase R-HSA-68877 (46, 13.7x, q 8.4e-37) | cortex CX02: Neuron Projection (35, 5.1x, q 5.7e-12); Chemical Synaptic Transmission (24, 8.0x, q 5.7e-12); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.7x, q 8.7e-12) | cortex CX03 (radial glia markers): Regulation Of Epithelial Cell Proliferation (9, 13.6x, q 7.2e-05); Collagen-Containing Extracellular Matrix (13, 6.3x, q 3.6e-04); Extracellular Matrix Organization R-HSA-1474244 (12, 6.2x, q 8.4e-04) | cortex CX04 (S phase): DNA Metabolic Process (39, 11.8x, q 1.4e-27); Cell Cycle, Mitotic R-HSA-69278 (46, 7.6x, q 7.1e-25); DNA-templated DNA Replication (23, 26.8x, q 2.7e-24) | cortex CX05: Cardiac Conduction R-HSA-5576891 (7, 10.6x, q 0.016); Muscle Contraction R-HSA-397014 (8, 8.0x, q 0.016) | cortex CXw01: Neuronal System R-HSA-112316 (31, 6.1x, q 1.8e-12); Neuron Projection (36, 4.8x, q 9.5e-12); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.1x, q 6.2e-11) | cortex CXw02 (G2/M phase): Cell Cycle, Mitotic R-HSA-69278 (75, 9.7x, q 4.1e-51); M Phase R-HSA-68886 (59, 10.5x, q 4.9e-41); Mitotic Prometaphase R-HSA-68877 (44, 14.8x, q 1.3e-36) | cortex CXw03 (radial glia markers): Collagen-Containing Extracellular Matrix (16, 7.5x, q 1.7e-06); Regulation Of Epithelial Cell Proliferation (9, 13.3x, q 4.4e-05); Negative Regulation Of Cell Differentiation (13, 6.6x, q 1.3e-04) | cortex CXw04 (S phase): DNA Metabolic Process (35, 15.9x, q 1.8e-29); Cell Cycle, Mitotic R-HSA-69278 (41, 10.1x, q 1.2e-27); Cell Cycle Checkpoints R-HSA-69620 (32, 15.0x, q 3.4e-26) | cortex CXw05: Neuron Projection (14, 5.5x, q 9.0e-04); Central Nervous System Development (9, 6.4x, q 0.026); Inorganic Cation Import Across Plasma Membrane (5, 15.4x, q 0.027) | human_dev HD01 (G2/M phase; S phase): Cell Cycle, Mitotic R-HSA-69278 (102, 9.8x, q 3.4e-72); Cell Cycle Checkpoints R-HSA-69620 (67, 12.6x, q 7.3e-53); Mitotic Prometaphase R-HSA-68877 (48, 12.0x, q 1.1e-35) | human_dev HD02: Collagen-Containing Extracellular Matrix (17, 4.7x, q 5.0e-04); Hemostasis R-HSA-109582 (24, 3.3x, q 5.0e-04); Vesicle (16, 4.7x, q 5.0e-04) | human_dev HD03: Neuron Projection (41, 6.0x, q 3.5e-17); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.6x, q 1.9e-11); Neuronal System R-HSA-112316 (27, 5.9x, q 1.3e-10) | human_dev HD05: Cilium Movement (9, 96.8x, q 5.6e-13); Cilium (12, 19.6x, q 1.2e-09); Axoneme Assembly (5, 66.6x, q 1.4e-05) | human_dev HDw01 (G2/M phase; S phase): Cell Cycle, Mitotic R-HSA-69278 (79, 10.5x, q 2.4e-57); Cell Cycle Checkpoints R-HSA-69620 (56, 14.5x, q 3.9e-47); M Phase R-HSA-68886 (54, 10.0x, q 3.6e-36) | human_dev HDw02: Cilium (14, 20.5x, q 1.2e-11); Cilium Movement (8, 77.1x, q 1.2e-10); Axoneme Assembly (6, 71.6x, q 2.4e-07) | (+1 more groups).
- **Age trends -- top terms** (overlap, fold, q): cortex Neuroblast: up: Ventricular Septum Development (4, 30.3x, q 0.035); Cohesin Loading Onto Chromatin R-HSA-2470946 (3, 59.0x, q 0.035); Outflow Tract Septum Morphogenesis (3, 49.2x, q 0.042) | cortex Neuronal IPC: down: Transcriptional Regulation By RUNX3 R-HSA-8878159 (5, 16.0x, q 0.036); Ubiquitin-dependent Degradation Of Cyclin D R-HSA-75815 (4, 23.2x, q 0.036); Hh Mutants Abrogate Ligand Secretion R-HSA-5387390 (4, 22.2x, q 0.036) | cortex Radial glia: down: Glycolysis / Gluconeogenesis (8, 12.7x, q 3.2e-04); S Phase R-HSA-69242 (14, 5.7x, q 3.2e-04); Synthesis Of DNA R-HSA-69239 (12, 6.7x, q 3.2e-04) | cortex Radial glia: up: Microtubule Binding (21, 7.3x, q 4.5e-09); Tubulin Binding (22, 5.8x, q 7.3e-08); Microtubule Cytoskeleton (22, 5.4x, q 1.7e-07) | human_dev Erythrocyte: up: RNA Polymerase III Transcription Termination R-HSA-73980 (4, 53.4x, q 0.0038); RNA Polymerase III Abortive And Retractive Initiation R-HSA-749476 (4, 28.1x, q 0.028) | human_dev Immune: up: Regulation Of Actin Filament Polymerization (5, 20.4x, q 0.02) | human_dev Neuroblast: up: Neuronal System R-HSA-112316 (44, 2.4x, q 1.5e-04); Axon Guidance (25, 2.8x, q 0.0042); Axonogenesis (29, 2.5x, q 0.0052) | human_dev Neuron: up: Selenocysteine Synthesis R-HSA-2408557 (28, 4.3x, q 5.2e-08); Peptide Chain Elongation R-HSA-156902 (26, 4.2x, q 2.2e-07); Selenoamino Acid Metabolism R-HSA-2408522 (29, 3.8x, q 2.2e-07) | human_dev Neuronal IPC: up: Neuron Projection Development (27, 3.0x, q 0.0012); Regulation Of Dendrite Development (11, 5.7x, q 0.0031); Cytoskeleton (51, 1.9x, q 0.0097) | human_dev Radial glia: up: Tubulin Binding (30, 2.4x, q 0.031) | human_dev Vascular: down: Purine-Containing Compound Biosynthetic Process (4, 35.8x, q 0.018); Cholesterol Biosynthesis R-HSA-191273 (4, 26.0x, q 0.034); Sterol Biosynthetic Process (4, 23.9x, q 0.034).
- **Proliferations -- top terms** (overlap, fold, q): cortex anti-proliferative: Neuron Projection (169, 2.3x, q 1.8e-26); Neuronal System R-HSA-112316 (116, 2.4x, q 1.2e-18); Transmission Across Chemical Synapses R-HSA-112315 (85, 2.6x, q 2.6e-16) | cortex proliferative: Cell Cycle, Mitotic R-HSA-69278 (226, 2.5x, q 1.2e-43); M Phase R-HSA-68886 (147, 2.2x, q 3.2e-21); Mitotic Prometaphase R-HSA-68877 (94, 2.7x, q 2.3e-20) | cortex proliferative, G2/M-leaning: Cell Cycle, Mitotic R-HSA-69278 (97, 6.3x, q 8.3e-48); M Phase R-HSA-68886 (79, 7.0x, q 4.1e-42); Mitotic Prometaphase R-HSA-68877 (54, 9.1x, q 3.3e-34) | cortex proliferative, S-leaning: DNA Metabolic Process (49, 9.2x, q 4.4e-30); DNA-templated DNA Replication (29, 21.0x, q 2.1e-28); DNA Strand Elongation R-HSA-69190 (19, 33.7x, q 1.9e-23) | human_dev anti-proliferative: Neuron Projection (188, 2.1x, q 9.1e-27); Transmission Across Chemical Synapses R-HSA-112315 (105, 2.6x, q 8.0e-23); Neuronal System R-HSA-112316 (135, 2.3x, q 8.0e-23) | human_dev proliferative: Cell Cycle, Mitotic R-HSA-69278 (247, 2.1x, q 1.7e-34); DNA Metabolic Process (143, 2.2x, q 4.6e-23); Processing Of Capped Intron-Containing Pre-mRNA R-HSA-72203 (136, 2.2x, q 1.0e-22).

Method, limitations and output files: [results/29_functional_enrichment/SUMMARY.md](results/29_functional_enrichment/SUMMARY.md)

## 30. Transcription factors: content, target enrichment and regulon activity in modules and NDD lists

**Question.** Which TFs sit in the modules and NDD lists, whose targets are over-represented in them, and is that TF's expression tied to its targets' in these data?

- **Gene sets rich in TFs** (TFs / expected; q < 0.05): human_dev module:HDw04 11 / 2.3 (DMBX1, EN1, EN2, IRX1, IRX2, IRX3); cortex seed:chromatin_transcription_regulators 8 / 2.4 (ADNP, ASH1L, BCL11A, BCL11B, KDM5B, KMT2A); cortex seed:id_dd_dominant 8 / 2.4 (FOXG1, GATAD2B, KMT2A, PURA, SATB2, SON); cortex module:CXw03 26 / 14.0 (ARX, CREB5, DACH1, GLI3, HES1, NFATC4); human_dev seed:chromatin_transcription_regulators 8 / 2.4 (ADNP, ASH1L, BCL11A, BCL11B, KDM5B, KMT2A); human_dev seed:id_dd_dominant 8 / 2.4 (FOXG1, GATAD2B, KMT2A, PURA, SATB2, SON); cortex module:CX09 8 / 2.5 (EOMES, HES6, INSM1, LHX9, NEUROD4, NEUROG2); cortex module:CXw08 7 / 2.2 (EOMES, HES6, INSM1, LHX9, NEUROD4, NEUROG2); cortex module:CX03 23 / 13.7 (ARX, CREB5, DACH1, GLI3, HES1, HES4); human_dev seed:asd_high_confidence 8 / 3.1 (ADNP, ASH1L, DEAF1, FOXP1, KDM5B, MECP2).
- **cortex: regulons whose TF tracks its targets across clusters** (32 of 221 tested; strongest): NEUROG2 (+0.67/+0.55), NFKB1 (+0.29/+0.39), NFKB2 (+0.34/+0.50), REST (+0.53/+0.52), SP1 (+0.36/+0.65), E2F1 (+0.63/+0.75), JUND (+0.41/+0.46), JUN (+0.37/+0.53), TFDP1 (+0.65/+0.66), FOS (+0.51/+0.51), E2F4 (+0.43/+0.43), E2F3 (+0.78/+0.83).
- **human_dev: regulons whose TF tracks its targets across clusters** (50 of 272 tested; strongest): EGR1 (+0.57/+0.67), HIF1A (+0.56/+0.38), JUN (+0.79/+0.79), NEUROG2 (+0.47/+0.53), SOX11 (+0.79/+0.81), XBP1 (+0.78/+0.69), FOXO3 (+0.27/+0.53), DLX2 (+0.37/+0.51), PAX2 (+0.55/+0.60), RARA (+0.33/+0.38), SRSF2 (+0.88/+0.87), RELA (+0.36/+0.44).
- **Candidate drivers** (TF targets enriched in the set, q < 0.05, and the TF tracks its targets in both donor sets; overlap, fold): cortex seed:chromatin_transcription_regulators: SP1 (8, 4.5x), MYB (3, 13.1x) | cortex seed:epilepsy_dee: REST (5, 62.9x) | cortex seed:synaptic_and_channels: REST (6, 56.6x) | cortex module:CX01: E2F4 (30, 10.7x), E2F1 (24, 5.3x), E2F3 (10, 9.7x), FOXM1 (8, 9.5x), E2F2 (7, 9.3x) | cortex module:CX02: REST (8, 8.2x) | cortex module:CX03: SOX2 (9, 14.6x), SP1 (29, 2.9x), JUN (18, 3.0x), SP3 (12, 3.9x), NFKB2 (17, 2.9x) | cortex module:CX04: E2F4 (29, 14.9x), E2F1 (34, 10.8x), E2F2 (13, 25.0x), E2F3 (12, 16.8x), FOXM1 (8, 13.7x) | cortex module:CX09: NEUROG2 (4, 72.8x) | cortex module:CXw01: REST (12, 11.3x) | cortex module:CXw02: E2F4 (29, 11.7x), E2F1 (25, 6.2x), E2F3 (10, 11.0x), FOXM1 (8, 10.7x), E2F2 (7, 10.6x) | cortex module:CXw03: SP1 (36, 3.5x), SOX2 (9, 14.2x), NFKB2 (22, 3.7x), NFKB1 (22, 3.2x), JUN (20, 3.2x) | cortex module:CXw04: E2F4 (29, 22.3x), E2F1 (31, 14.7x), E2F2 (13, 37.5x), E2F3 (12, 25.2x), TFDP1 (6, 28.9x) | cortex module:CXw08: NEUROG2 (4, 84.0x) | human_dev seed:asd_high_confidence: REST (3, 20.2x) | human_dev seed:chromatin_transcription_regulators: MBD2 (3, 27.6x), SP1 (8, 4.6x) | human_dev seed:epilepsy_dee: REST (6, 73.1x), JUND (5, 9.6x).

Method, limitations and output files: [results/30_tf_regulators/SUMMARY.md](results/30_tf_regulators/SUMMARY.md)

## 31. Mutation intolerance (gnomAD LOEUF) and Mendelian disease genes in lists and modules

**Question.** Are NDD lists, sub-modules and modules enriched for loss-of-function-intolerant and Mendelian disease genes beyond expression-matched genes, and do constrained genes vary less between donors, connect more within lists, or change more with age?

- **Gene sets more constrained than matched genes (mean LOEUF; null)**: cortex module:CXw05 0.79 vs 0.97; cortex module:CX05 0.81 vs 0.99; human_dev seed:epilepsy_dee 0.47 vs 0.93; cortex seed:epilepsy_dee 0.47 vs 0.96; human_dev seed:synaptic_and_channels 0.46 vs 0.87; human_dev seed:chromatin_transcription_regulators 0.36 vs 0.81; human_dev seed:asd_high_confidence 0.44 vs 0.84; cortex seed:synaptic_and_channels 0.46 vs 0.93; cortex seed:chromatin_transcription_regulators 0.36 vs 0.88; human_dev module:HD03 0.73 vs 0.92; cortex seed:asd_high_confidence 0.44 vs 0.91; cortex module:CXw01 0.74 vs 0.92; human_dev seed:id_dd_dominant 0.21 vs 0.82; cortex module:CX02 0.71 vs 0.93; cortex seed:id_dd_dominant 0.21 vs 0.88.
- **Gene sets richer in Mendelian disease genes (share; null)**: cortex module:CX04 0.37 vs 0.27; human_dev module:HD05 0.35 vs 0.20; cortex module:CX02 0.35 vs 0.27; cortex module:CXw04 0.42 vs 0.28; human_dev seed:synaptic_and_channels 0.92 vs 0.28; cortex seed:synaptic_and_channels 0.92 vs 0.27; human_dev seed:epilepsy_dee 1.00 vs 0.26; cortex seed:epilepsy_dee 1.00 vs 0.25; human_dev seed:chromatin_transcription_regulators 1.00 vs 0.30; cortex seed:id_dd_dominant 1.00 vs 0.29; human_dev seed:id_dd_dominant 1.00 vs 0.31; cortex seed:chromatin_transcription_regulators 1.00 vs 0.29; human_dev seed:asd_high_confidence 0.92 vs 0.29; cortex seed:asd_high_confidence 0.92 vs 0.28.
- **LOEUF vs between-donor variability percentile (15)** (Spearman; median over classes / strata, range): median +0.01 (-0.09 to +0.12; n = 27).
- **LOEUF vs |age-trend Z| (03)** (Spearman; median over classes / strata, range): median -0.05 (-0.08 to +0.10; n = 12).

Method, limitations and output files: [results/31_constraint_and_disease/SUMMARY.md](results/31_constraint_and_disease/SUMMARY.md)

## 32. Ligand-receptor signalling between cell types over development

**Question.** Which cell types could signal to which through ligand-receptor pairs, how does that potential change with age in both donor sets, and where do NDD genes act as ligands or receptors?

- **cortex: busiest sender -> receiver routes** (pairs expressed per age, mean of chemistries): Neuroblast -> Radial glia 174, Neuron -> Neuronal IPC 170, Neuron -> Radial glia 169, Neuroblast -> Neuronal IPC 169, Neuroblast -> Neuroblast 162, Neuroblast -> Glioblast 161.
- **human_dev: busiest sender -> receiver routes** (pairs expressed per age, mean of chemistries): Glioblast -> Glioblast 327, Fibroblast -> Glioblast 320, Fibroblast -> Fibroblast 315, Radial glia -> Glioblast 314, Glioblast -> Radial glia 314, Glioblast -> Oligo 303.
- **cortex: replicated changes in signalling potential** (174 of 1605 pair x route tests): Neuroblast -> Radial glia: 16 rise / 7 fall; Neuroblast -> Neuronal IPC: 4 rise / 5 fall; Neuron -> Radial glia: 5 rise / 5 fall; Neuroblast -> Neuroblast: 10 rise / 4 fall; Radial glia -> Neuroblast: 7 rise / 4 fall; Neuroblast -> Neuron: 14 rise / 3 fall; Neuron -> Neuronal IPC: 5 rise / 3 fall; Neuronal IPC -> Neuronal IPC: 5 rise / 3 fall.
- **cortex: replicated changes involving NDD genes** (ligand -> receptor, sender -> receiver; rho v2/v3): CNTN2 -> CNTNAP2 (Neuroblast -> Neuronal IPC) falls with age (-0.96/-0.89); CNTN2 -> CNTNAP2 (Neuroblast -> Neuroblast) falls with age (-0.93/-0.94); CNTN2 -> CNTNAP2 (Neuroblast -> Neuron) falls with age (-0.86/-1.00).
- **human_dev: replicated changes in signalling potential** (1802 of 10094 pair x route tests): Radial glia -> Radial glia: 26 rise / 45 fall; Neuronal IPC -> Radial glia: 21 rise / 31 fall; Neuroblast -> Radial glia: 28 rise / 31 fall; Neuron -> Radial glia: 31 rise / 29 fall; Radial glia -> Neuronal IPC: 30 rise / 29 fall; Radial glia -> Vascular: 22 rise / 27 fall; Radial glia -> Neuron: 30 rise / 26 fall; Radial glia -> Neuroblast: 39 rise / 24 fall.
- **human_dev: replicated changes involving NDD genes** (ligand -> receptor, sender -> receiver; rho v2/v3): NLGN1 -> NRXN1 (Radial glia -> Radial glia) rises with age (+0.96/+1.00); NLGN1 -> NRXN1 (Neuroblast -> Neuroblast) rises with age (+0.94/+0.98); NLGN1 -> NRXN1 (Neuroblast -> Neuron) rises with age (+0.96/+0.93); NLGN1 -> NRXN1 (Neuroblast -> Neuronal IPC) rises with age (+0.93/+0.98); TAFA2 -> NRXN1 (Radial glia -> Neuronal IPC) rises with age (+0.93/+0.98); NLGN1 -> NRXN1 (Neuronal IPC -> Radial glia) rises with age (+0.85/+1.00); NLGN1 -> NRXN1 (Neuronal IPC -> Neuronal IPC) rises with age (+0.93/+0.95); TAFA2 -> NRXN1 (Neuroblast -> Neuronal IPC) rises with age (+0.90/+0.95); NLGN1 -> NRXN1 (Radial glia -> Glioblast) rises with age (+0.92/+0.96); NLGN1 -> NRXN1 (Neuroblast -> Radial glia) rises with age (+0.86/+0.98) (+97 more).

Method, limitations and output files: [results/32_cell_communication/SUMMARY.md](results/32_cell_communication/SUMMARY.md)

## 33. Gene-list symbols recovered through HGNC, and entries still missing

**Question.** Skipped: no gene lists.

- **Not run**: no gene lists. Fetch annotations with running_scripts/fetch_annotations.sh.

Method, limitations and output files: [results/33_symbol_rescue/SUMMARY.md](results/33_symbol_rescue/SUMMARY.md)

## 34. External validation: cortex age trends vs BrainSpan neocortex

**Question.** Do genes and gene lists that change with age in cortex change the same way in an independent atlas (BrainSpan neocortex) over the overlapping ages?

- **BrainSpan**: 11 donors, 8-16 pcw, neocortical samples averaged per donor; 8,527 genes shared with cortex.
- **Gene-level agreement**: Spearman of our combined Z with BrainSpan rho +0.26; of 1378 genes with a replicated whole-cortex age trend, 68% change the same way in BrainSpan (50% by chance), 72% of the 216 also nominal there.
- **Gene sets moving with age in BrainSpan** (mean rho BrainSpan vs ours; same direction?): seed:s_phase -0.39 vs -0.55 (yes); seed:synaptic_and_channels +0.38 vs +0.67 (yes); module:CX01 -0.22 vs -0.33 (yes); module:CX02 +0.31 vs +0.67 (yes); module:CX04 -0.21 vs -0.51 (yes); module:CX05 +0.17 vs +0.53 (yes); module:CX07 +0.14 vs +0.50 (yes); module:CXw01 +0.31 vs +0.70 (yes); module:CXw02 -0.21 vs -0.31 (yes); module:CXw04 -0.22 vs -0.52 (yes); module:CXw05 +0.14 vs +0.64 (yes); seed:g2m_phase -0.21 vs -0.35 (yes); seed:epilepsy_dee +0.23 vs +0.31 (yes).

Method, limitations and output files: [results/34_brainspan_validation/SUMMARY.md](results/34_brainspan_validation/SUMMARY.md)


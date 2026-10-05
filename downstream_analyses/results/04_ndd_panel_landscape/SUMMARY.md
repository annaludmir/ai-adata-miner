# NDD gene panels across cell classes and age

_Generated 2026-10-05 08:03 UTC by `downstream_analyses/04_ndd_panel_landscape.py` from `csv_exports/`._

## Question

Which cell classes preferentially express each neurodevelopmental-disorder panel, consistently across donors and chemistries -- and are NDD genes enriched among genes that change with age?

## Inputs

- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/10_markers/specificity_cell_class.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v2/17_gsea/gsea_results.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v3/10_markers/specificity_cell_class.csv`
- `cortex__v3/11_panels/panel_coverage.csv`
- `cortex__v3/17_gsea/gsea_results.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/10_markers/specificity_cell_class.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v2/17_gsea/gsea_results.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v3/10_markers/specificity_cell_class.csv`
- `human_dev__v3/11_panels/panel_coverage.csv`
- `human_dev__v3/17_gsea/gsea_results.csv`
- `results/03_age_trends_within_cell_class/age_trends_combined.csv`

## Method

- Panels: the NDD group of the seed panels (panels/README.md), restricted to genes in the pseudobulk export.
- A: log2 TMM-CPM per (class, age point), genes >= 5 CPM, Z-scored across columns; panel score = mean Z. T = mean over age points of (class - mean of other classes) -- age points are donors. Null: 5,000 random sets, one gene drawn from each panel gene's expression decile. Effect in null SDs; v2/v3 by signed Stouffer, BH per dataset.
- B: overlap with 03's replicated trends; p from expression-matched random sets drawn from the genes tested in 03 (plain hypergeometric reported alongside -- it over-calls highly expressed panels).
- C: 17_gsea results joined across chemistries. D: 10_markers top class per gene.

## Key findings

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

## Limitations

- Seed panels are short hand-picked lists, not SFARI/DDG2P releases; a missing enrichment can be a panel-coverage problem (11_panels/panel_coverage.csv).
- The matched null controls for expression level, not for gene length or connectivity; long neuronal genes can still look neuron-enriched for that reason.
- Gene-set nulls ignore gene-gene correlation, so p-values are optimistic for co-regulated panels; the replication requirement is the main safeguard.
- human_dev classes pool brain regions (see 03).

## What would strengthen this

- Swap in full SFARI / DDG2P / Epi25 releases via panels/*.csv and re-run 11, 17 and this.
- Add gene-length-matched nulls (gene_id_map has Start/End).
- Repeat part A at cluster resolution to localise panels to specific neuron types.

## Output files

- `panel_class_preference_per_stratum.csv` -- Panel x class: mean Z difference vs other classes over age points, matched-null p
- `panel_class_preference_combined.csv` -- v2 x v3 combined panel preference per class; tier replicated / supported
- `panel_age_trend_enrichment.csv` -- NDD panel genes among replicated rising/falling genes of 03; p vs expression-matched random sets (hypergeometric shown for comparison)
- `gsea_panel_enrichment_by_chemistry.csv` -- 17_gsea enrichments significant in either chemistry, with the other alongside
- `ndd_gene_top_class.csv` -- Per NDD gene: most specific cell class and tau in each chemistry
- `panel_preference_cortex.png` -- NDD panel x cell class preference, cortex
- `panel_preference_human_dev.png` -- NDD panel x cell class preference, human_dev

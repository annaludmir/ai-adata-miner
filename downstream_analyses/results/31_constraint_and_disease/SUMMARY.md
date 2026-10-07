# Mutation intolerance (gnomAD LOEUF) and Mendelian disease genes in lists and modules

_Generated 2026-10-07 11:40 UTC by `downstream_analyses/31_constraint_and_disease.py` from `csv_exports/`._

## Question

Are NDD lists, sub-modules and modules enriched for loss-of-function-intolerant and Mendelian disease genes beyond expression-matched genes, and do constrained genes vary less between donors, connect more within lists, or change more with age?

## Inputs

- `annotations/constraint/gnomad.v4.1.constraint_metrics.tsv`
- `annotations/disease/genes_to_disease.txt`
- `cortex__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `results/03_age_trends_within_cell_class/age_trends_combined.csv`
- `results/08_coexpression_modules/modules.csv`
- `results/15_donor_variability/gene_variability.csv`

## Method

- gnomAD v4.1 constraint (LOEUF; MANE / canonical transcript); HPO gene-to-disease (Mendelian = OMIM / Orphanet genes annotated MENDELIAN); symbols mapped to each dataset.
- A: per set, mean LOEUF, share in the genome's most constrained LOEUF decile, Mendelian share, each vs 2,000 random sets matched on expression decile; BH per dataset x statistic.
- B: Spearman of LOEUF with 15's variability percentile, 07's connectivity and 03's |combined Z|.

## Key findings

- **Gene sets more constrained than matched genes (mean LOEUF; null)**: cortex module:CXw05 0.79 vs 0.97; cortex module:CX05 0.81 vs 0.99; human_dev seed:epilepsy_dee 0.47 vs 0.93; cortex seed:epilepsy_dee 0.47 vs 0.96; human_dev seed:synaptic_and_channels 0.46 vs 0.87; human_dev seed:chromatin_transcription_regulators 0.36 vs 0.81; human_dev seed:asd_high_confidence 0.44 vs 0.84; cortex seed:synaptic_and_channels 0.46 vs 0.93; cortex seed:chromatin_transcription_regulators 0.36 vs 0.88; human_dev module:HD03 0.73 vs 0.92; cortex seed:asd_high_confidence 0.44 vs 0.91; cortex module:CXw01 0.74 vs 0.92; human_dev seed:id_dd_dominant 0.21 vs 0.82; cortex module:CX02 0.71 vs 0.93; cortex seed:id_dd_dominant 0.21 vs 0.88.
- **Gene sets richer in Mendelian disease genes (share; null)**: cortex module:CX04 0.37 vs 0.27; human_dev module:HD05 0.35 vs 0.20; cortex module:CX02 0.35 vs 0.27; cortex module:CXw04 0.42 vs 0.28; human_dev seed:synaptic_and_channels 0.92 vs 0.28; cortex seed:synaptic_and_channels 0.92 vs 0.27; human_dev seed:epilepsy_dee 1.00 vs 0.26; cortex seed:epilepsy_dee 1.00 vs 0.25; human_dev seed:chromatin_transcription_regulators 1.00 vs 0.30; cortex seed:id_dd_dominant 1.00 vs 0.29; human_dev seed:id_dd_dominant 1.00 vs 0.31; cortex seed:chromatin_transcription_regulators 1.00 vs 0.29; human_dev seed:asd_high_confidence 0.92 vs 0.29; cortex seed:asd_high_confidence 0.92 vs 0.28.
- **LOEUF vs between-donor variability percentile (15)** (Spearman; median over classes / strata, range): median +0.01 (-0.09 to +0.12; n = 27).
- **LOEUF vs |age-trend Z| (03)** (Spearman; median over classes / strata, range): median -0.05 (-0.08 to +0.10; n = 12).

## Limitations

- Seed panels and many user lists were built from disease genes, so their Mendelian share is circular; the informative cases are modules and sub-modules found without lists.
- LOEUF is unreliable for short genes (few expected LoF variants); its upper bound is conservative there.

## What would strengthen this

- Weight list genes by constraint in 06/07, or test constrained and unconstrained halves separately.

## Output files

- `set_constraint.csv` -- Per gene set x statistic: mean LOEUF / constrained-decile share / Mendelian share vs expression-matched random sets; BH per dataset x statistic
- `constraint_relations.csv` -- Spearman of LOEUF with step-3 gene properties (variability, connectivity, age-trend strength)
- `gene_constraint.csv` -- Per gene: expression level (log2 CPM + 1), LOEUF, pLI, missense Z, Mendelian flag

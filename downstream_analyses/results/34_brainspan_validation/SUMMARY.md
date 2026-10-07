# External validation: cortex age trends vs BrainSpan neocortex

_Generated 2026-10-07 11:40 UTC by `downstream_analyses/34_brainspan_validation.py` from `csv_exports/`._

## Question

Do genes and gene lists that change with age in cortex change the same way in an independent atlas (BrainSpan neocortex) over the overlapping ages?

## Inputs

- `annotations/brainspan/columns_metadata.csv`
- `annotations/brainspan/expression_matrix.csv`
- `annotations/brainspan/rows_metadata.csv`
- `cortex__v2/09_pseudobulk/age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/age__pseudobulk_counts.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- BrainSpan RNA-seq (Gencode v10), neocortical samples of donors 8-16 pcw, log2(RPKM + 1) averaged per donor, genes >= 1 RPKM; Spearman with age, permutation p; matched via Ensembl id.
- Ours: whole-cortex pseudobulk per age point, log2 TMM-CPM, genes >= 5 CPM; Spearman per chemistry, v2 x v3 signed Stouffer.
- Lists: mean BrainSpan rho vs 2,000 random genes matched on BrainSpan expression decile; BH.

## Key findings

- **BrainSpan**: 11 donors, 8-16 pcw, neocortical samples averaged per donor; 8,527 genes shared with cortex.
- **Gene-level agreement**: Spearman of our combined Z with BrainSpan rho +0.26; of 1378 genes with a replicated whole-cortex age trend, 68% change the same way in BrainSpan (50% by chance), 72% of the 216 also nominal there.
- **Gene sets moving with age in BrainSpan** (mean rho BrainSpan vs ours; same direction?): seed:s_phase -0.39 vs -0.55 (yes); seed:synaptic_and_channels +0.38 vs +0.67 (yes); module:CX01 -0.22 vs -0.33 (yes); module:CX02 +0.31 vs +0.67 (yes); module:CX04 -0.21 vs -0.51 (yes); module:CX05 +0.17 vs +0.53 (yes); module:CX07 +0.14 vs +0.50 (yes); module:CXw01 +0.31 vs +0.70 (yes); module:CXw02 -0.21 vs -0.31 (yes); module:CXw04 -0.22 vs -0.52 (yes); module:CXw05 +0.14 vs +0.64 (yes); seed:g2m_phase -0.21 vs -0.35 (yes); seed:epilepsy_dee +0.23 vs +0.31 (yes).

## Limitations

- Bulk tissue mixes cell types, and so does our whole-cortex pseudobulk: agreement here includes shifts in cell-type composition with age, which is the point of a like-for-like comparison but not a within-cell-type validation.
- BrainSpan has few prenatal donors in this window (about one per age point) and starts at 8 pcw; our cortex reaches back to 5.5 pcw.
- BrainSpan dissections are cortical areas, ours are EMX1-lineage cells; interneurons and non-neural cells are in BrainSpan only.

## What would strengthen this

- Deconvolve BrainSpan with our cell-class profiles to compare within-class trends.

## Output files

- `brainspan_donors.csv` -- BrainSpan donors used: age and neocortical samples averaged
- `gene_trends_ours_vs_brainspan.csv` -- Per gene: our whole-cortex combined trend and BrainSpan neocortex rho with age
- `agreement.csv` -- Gene-level agreement between our cortex trends and BrainSpan
- `set_trends_ours_vs_brainspan.csv` -- Per gene set: mean rho with age in BrainSpan (vs matched null) and in our whole-cortex pseudobulk

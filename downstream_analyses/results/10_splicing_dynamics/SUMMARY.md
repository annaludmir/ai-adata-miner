# Splicing dynamics: where and when genes are being switched on or off

_Generated 2026-10-05 10:57 UTC by `downstream_analyses/10_splicing_dynamics.py` from `csv_exports/`._

## Question

In which cell classes are genes and gene sets being switched on or off (nascent vs mature RNA), and inside a class, is their unspliced share rising or falling with age alongside their expression?

## Inputs

- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v2/12_splicing/cell_class__spliced_counts.csv`
- `cortex__v2/12_splicing/cell_class__unspliced_counts.csv`
- `cortex__v2/12_splicing/cell_class_x_age__spliced_counts.csv`
- `cortex__v2/12_splicing/cell_class_x_age__unspliced_counts.csv`
- `cortex__v3/12_splicing/cell_class__spliced_counts.csv`
- `cortex__v3/12_splicing/cell_class__unspliced_counts.csv`
- `cortex__v3/12_splicing/cell_class_x_age__spliced_counts.csv`
- `cortex__v3/12_splicing/cell_class_x_age__unspliced_counts.csv`
- `results/03_age_trends_within_cell_class/age_trends_combined.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- Cortex only (human_dev has no spliced/unspliced layers). Pseudobulk spliced and unspliced counts from script 12; groups with >= 50 cells; genes with U+S >= 20 in a group.
- Log unspliced/spliced ratio log((U+0.5)/(S+0.5)), each group's median over genes subtracted (its overall unspliced level shifts for technical reasons); comparing a gene with itself across groups cancels its intron content.
- A gene is compared only across groups where it is clearly on (>= 10 CPM and >= 25% of its peak; for age trends, within the class's own ages): where a gene is nearly off, its residual reads are disproportionately unspliced. The remaining coupling between log U/S and expression (a pooled within-gene slope) is subtracted.
- A: induction score = centred ratio in a class minus the gene's mean over classes.
- B: Spearman of the centred ratio with age per class (>= 5 ages, complete profiles), exact permutation p, v2/v3 signed Stouffer, BH per class, tiered; joined with 03's expression trend.
- Gene sets: mean score vs 2,000 random sets matched on expression level x baseline ratio (5 x 5 bins); v2/v3 combined, BH, tiered.

## Key findings

- **Induction scores reproduce across donor sets** (Spearman v2 vs v3 per class): Glioblast 0.59, Neuroblast 0.57, Neuron 0.71, Neuronal IPC 0.53, Radial glia 0.54. Low values would mean the unspliced signal is noise at this depth.
- **Where gene sets are being switched on or off along the lineage** (replicated; + = more nascent RNA than their mRNA explains) -- seed:s_phase being switched off in Radial glia (-3.2/-5.4 null SDs); seed:s_phase being switched on in Neuronal IPC (+4.5/+2.3 null SDs); seed:g2m_phase being switched off in Radial glia (-4.8/-2.0 null SDs); seed:g2m_phase being switched on in Glioblast (+4.0/+2.0 null SDs); seed:epilepsy_dee being switched on in Neuronal IPC (+1.9/+1.9 null SDs); seed:asd_high_confidence being switched on in Radial glia (+1.7/+1.9 null SDs). Modules: CXw04 being switched on in Neuronal IPC; CXw02 being switched off in Radial glia; CXw01 being switched on in Radial glia; CXw02 being switched on in Glioblast; CXw01 being switched on in Glioblast; CX02 being switched on in Neuronal IPC; CX02 being switched on in Radial glia; CX01 being switched on in Glioblast.
- **Genes by dynamics within classes** (unspliced-share trend and expression trend both replicated, or expression flat): Radial glia: 2 induction ahead of expression, 3 shut-down ahead of expression.
- **Gene sets whose unspliced share changes with age inside a class** (replicated) -- seed:s_phase unspliced share rises in Radial glia (mean rho +0.33/+0.28); seed:g2m_phase unspliced share rises in Neuronal IPC (mean rho +0.24/+0.38); seed:asd_high_confidence unspliced share falls in Neuroblast (mean rho -0.13/-0.16); seed:chromatin_transcription_regulators unspliced share falls in Neuroblast (mean rho -0.11/-0.16). Read with 06/09: a list rising with age whose unspliced share also rises is being actively induced; one whose unspliced share is flat is accumulating mRNA.
- **Level coupling removed**: across groups, a gene's log U/S falls as its expression rises (pooled within-gene slope per log2 CPM: v2 cell_class -0.125, v3 cell_class -0.196, v2 cell_class_x_age +0.010, v3 cell_class_x_age -0.035); scores above are net of it, and genes are compared only where clearly on.
- **Overall unspliced level vs age** (removed before every test above; Spearman): Glioblast +0.54/-0.80, Neuroblast +0.36/+0.71, Neuron +0.39/+1.00, Neuronal IPC +0.36/+0.64, Radial glia +0.43/+0.64.

## Limitations

- Pseudobulk ratios average over cells in different states; this is a population measure of induction, not per-cell velocity.
- Unspliced counts include intron-retaining and nuclear transcripts; the group-level centring removes global shifts but not gene-specific changes in retention.
- Age points are donors (5-7 per chemistry in cortex).
- Cortex only; no replication in human_dev is possible from these files.

## What would strengthen this

- Per-cell velocity (scVelo) on the cortex layers, to place induction along differentiation pseudotime.
- Intron-length-aware matching for the gene-set nulls.

## Output files

- `induction_by_cell_class.csv` -- Per gene x cell class: centred log U/S minus the gene's mean over classes (+ = being switched on there), per chemistry
- `set_induction_by_cell_class_per_stratum.csv` -- Gene set x class: mean induction score vs matched random sets, per chemistry
- `set_induction_by_cell_class.csv` -- Gene set x class: induction combined over chemistries; tier replicated / supported
- `unspliced_age_trends_per_stratum.csv` -- Per gene x class x chemistry: Spearman of centred log U/S with age (exact permutation p)
- `unspliced_age_trends.csv` -- Per gene x class: unspliced-share age trend combined over chemistries, joined with 03's expression trend into a dynamics category
- `set_unspliced_age_trends_per_stratum.csv` -- Gene set x class: mean unspliced-share age trend vs matched random sets, per chemistry
- `set_unspliced_age_trends.csv` -- Gene set x class: unspliced-share age trend combined over chemistries; tiered
- `group_unspliced_levels.csv` -- Per class x age: overall unspliced level (technical and biological), removed before testing
- `splicing_dynamics_cortex.png` -- Gene-set induction by class and unspliced-share age trends

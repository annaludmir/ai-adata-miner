# Cell-class composition across age

_Generated 2026-10-05 06:04 UTC by `downstream_analyses/02_composition_vs_age.py` from `csv_exports/`._

## Question

Which cell classes expand or shrink with developmental age, consistently in two independent donor sets (v2 and v3 chemistry)?

## Inputs

- `cortex__v2/01_overview/sample_summary.csv`
- `cortex__v2/02_composition/counts_cell_class_by_sample.csv`
- `cortex__v3/01_overview/sample_summary.csv`
- `cortex__v3/02_composition/counts_cell_class_by_sample.csv`
- `human_dev__v2/01_overview/sample_summary.csv`
- `human_dev__v2/02_composition/counts_cell_class_by_sample.csv`
- `human_dev__v3/01_overview/sample_summary.csv`
- `human_dev__v3/02_composition/counts_cell_class_by_sample.csv`
- `results/01_data_audit/donor_sex_inferred.csv`

## Method

- Replicate unit: donor. Cells summed per donor within a tissue scope; donors with < 200 cells in the scope dropped; scopes need >= 5 donors in each chemistry.
- Classes below 0.5% of the scope merged into 'Other'; centred log-ratio with pseudocount 0.5.
- Spearman rho of CLR vs age per chemistry; two-sided permutation p, exact over all donor orderings (ties in age handled by permuting the observed ages).
- v2 and v3 combined by signed Stouffer (weights sqrt(n donors)); BH within each dataset x scope. Replicated = same direction, combined q < 0.05 and each chemistry nominally significant (one-sided p < 0.05); supported = combined q < 0.05 but one chemistry weak.
- Sex check: Spearman of CLR vs male indicator, reported next to each trend.

## Key findings

- **cortex / all cells** (7 + 7 donors), replicated: Radial glia down (rho v2 -0.79, v3 -0.95, q = 0.0025); Neuron up (rho v2 +0.79, v3 +0.88, q = 0.0046). Supported by the combined test but weak in one chemistry: Neuroblast down (rho v2 -0.71, v3 -0.59).
- **human_dev / Cerebellum** (6 + 6 donors): no class trend replicates in both chemistries. Strongest: Glioblast (rho v2 +0.93, v3 +0.64, combined p = 0.0093).
- **human_dev / Midbrain** (8 + 7 donors), replicated: Glioblast up (rho v2 +0.95, v3 +0.99, q = 4.3e-05); Radial glia down (rho v2 -0.82, v3 -0.95, q = 9.9e-04); Neuroblast down (rho v2 -0.69, v3 -0.88, q = 0.0069). Supported by the combined test but weak in one chemistry: Other up (rho v2 +0.53, v3 +0.99).
- **human_dev / Telencephalon** (7 + 6 donors), replicated: Radial glia down (rho v2 -0.81, v3 -0.93, q = 0.012).
- **Sex does not explain the replicated trends.** Where a class tracks donor sex, it does so in one chemistry only, or in opposite directions (e.g. cortex/all cells Radial glia: rho with male +0.79 in v2, -0.29 in v3; cortex/all cells Neuron: rho with male -0.79 in v2, +0.29 in v3; human_dev/Midbrain Neuroblast: rho with male +0.62 in v2, -0.14 in v3), while the age trend keeps its direction in both.
- **Trends that flip between chemistries** (significant in one, opposite sign in the other) -- not replicated, possibly donor- or window-specific: human_dev/Cerebellum Other; human_dev/Cerebellum Neuron; human_dev/Midbrain Neuron.

## Limitations

- v2 spans ~6-10 pcw and v3 ~5-14 pcw with a gap at 7-11.5, so 'replicated' means monotonic across two different windows -- a rise-then-fall would fail.
- Composition depends on how each sample was dissected and dissociated; a trend can reflect changing dissection practice with age, which donor-level data cannot rule out.
- cortex and human_dev share donors, so their agreement is not independent evidence.
- Spearman tests monotonic trends only; n = 4-9 donors per chemistry limits power.

## What would strengthen this

- Add donors at the ages one chemistry lacks, so trends are tested on one window.
- Use per-sample dissection metadata to model composition with dissection as a covariate.

## Output files

- `donor_composition.csv` -- Per donor: cell-class count, fraction and CLR
- `composition_age_trends.csv` -- Per stratum: Spearman of CLR vs age across donors, exact permutation p
- `composition_trends_replicated.csv` -- v2 and v3 trends combined (signed Stouffer); tier: replicated = both chemistries nominal + combined q < 0.05, supported = combined only
- `composition_cortex_all_cells.png` -- Cell-class fraction per donor vs age, cortex / all cells
- `composition_human_dev_Cerebellum.png` -- Cell-class fraction per donor vs age, human_dev / Cerebellum
- `composition_human_dev_Midbrain.png` -- Cell-class fraction per donor vs age, human_dev / Midbrain
- `composition_human_dev_Telencephalon.png` -- Cell-class fraction per donor vs age, human_dev / Telencephalon

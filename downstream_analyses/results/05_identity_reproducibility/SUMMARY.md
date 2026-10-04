# Cell-class identity across chemistries and files

_Generated 2026-10-04 13:15 UTC by `downstream_analyses/05_identity_reproducibility.py` from `csv_exports/`._

## Question

Does each cell class carry the same expression identity in v2 and v3, and in cortex and human_dev -- i.e. can class-level results be compared across strata?

## Inputs

- `cortex__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `cortex__v2/10_markers/top_markers_cell_class.csv`
- `cortex__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `cortex__v3/10_markers/top_markers_cell_class.csv`
- `human_dev__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v2/10_markers/top_markers_cell_class.csv`
- `human_dev__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v3/10_markers/top_markers_cell_class.csv`

## Method

- Per stratum: cell-class pseudobulk -> log2 TMM-CPM. Identity profile = class minus mean of the other classes shared by both strata in the comparison.
- Pearson correlation of profiles over the 2000 genes with the largest identity contrast on both sides; margin = same-class r minus best other-class r.
- Top-100 marker Jaccard from 10_markers/top_markers_cell_class.csv.

## Key findings

- **cortex__v2 vs cortex__v3** (different chemistry, independent donors): same-class r 0.96-0.97 (median 0.96) vs best other class 0.80 at most; top-100 marker Jaccard median 0.39. Closest neighbours (margin < 0.2): Glioblast (r 0.96 vs Radial glia 0.77); Radial glia (r 0.96 vs Glioblast 0.80).
- **human_dev__v2 vs human_dev__v3** (different chemistry, independent donors): same-class r 0.82-0.98 (median 0.95) vs best other class 0.88 at most; top-100 marker Jaccard median 0.47. Closest neighbours (margin < 0.2): Neuroblast (r 0.98 vs Neuron 0.88); Neuron (r 0.98 vs Neuroblast 0.84).
- **cortex__v2 vs human_dev__v2** (same donors, different file): same-class r 0.90-0.93 (median 0.91) vs best other class 0.79 at most; top-100 marker Jaccard median 0.14. Closest neighbours (margin < 0.2): Glioblast (r 0.90 vs Radial glia 0.79); Neuroblast (r 0.93 vs Neuron 0.79); Neuron (r 0.91 vs Neuroblast 0.72); Radial glia (r 0.93 vs Glioblast 0.77).
- **cortex__v3 vs human_dev__v3** (same donors, different file): same-class r 0.87-0.91 (median 0.89) vs best other class 0.79 at most; top-100 marker Jaccard median 0.09. Closest neighbours (margin < 0.2): Neuroblast (r 0.91 vs Neuron 0.79); Neuron (r 0.89 vs Neuroblast 0.74).
- **Class identity reproduces across chemistries**: a class matches itself first in 16/16 chemistry comparisons, so v2 and v3 results about the same class describe the same cell type. The smallest margins fall between adjacent lineages (radial glia / glioblast, neuroblast / neuron), as expected along a differentiation continuum. Marker lists overlap less than profiles agree, because v2's lower sensitivity reorders the tail of each list.
- **Across files, profiles agree but marker lists do not** (median Jaccard 0.14): human_dev classes span the whole brain, so their top markers include regional genes that cortex cells never express. Compare the files by profile, not by marker list.

## Limitations

- Class-level pseudobulk pools all donors in a stratum, so this tests the label's programme, not donor-to-donor variability.
- cortex vs human_dev share donors; their agreement bounds processing and regional pooling differences, not biological replication.
- Few shared classes (5) make the off-diagonal a coarse reference.

## What would strengthen this

- Repeat at cluster level with a cluster-to-cluster matching to find sub-types that do not reproduce across chemistry.

## Output files

- `identity_correlation__cortex__v2__cortex__v3.csv` -- Correlation of class identity profiles, cortex__v2 vs cortex__v3
- `identity_correlation__human_dev__v2__human_dev__v3.csv` -- Correlation of class identity profiles, human_dev__v2 vs human_dev__v3
- `identity_correlation__cortex__v2__human_dev__v2.csv` -- Correlation of class identity profiles, cortex__v2 vs human_dev__v2
- `identity_correlation__cortex__v3__human_dev__v3.csv` -- Correlation of class identity profiles, cortex__v3 vs human_dev__v3
- `identity_reproducibility.csv` -- Per class: same-class r, best other-class r, marker Jaccard
- `identity_correlation.png` -- Class x class identity correlation per comparison

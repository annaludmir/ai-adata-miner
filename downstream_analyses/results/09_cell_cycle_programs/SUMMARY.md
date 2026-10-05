# Cell-cycle programmes: proliferation over development and a phase map of genes

_Generated 2026-10-05 06:37 UTC by `downstream_analyses/09_cell_cycle_programs.py` from `csv_exports/`._

## Question

How does proliferation change over development within progenitor types; which genes and gene groups follow proliferation, and which lean to S or G2/M; and do list-level age trends in progenitors survive once proliferation-linked genes are set aside?

## Inputs

- `cortex__v2/03_cellcycle/cycle_scores_by_cluster_ClustersSurprise.csv`
- `cortex__v2/03_cellcycle/phase_fractions_by_cluster_ClustersSurprise.csv`
- `cortex__v2/03_cellcycle/proliferation_trajectory.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/03_cellcycle/cycle_scores_by_cluster_ClustersSurprise.csv`
- `cortex__v3/03_cellcycle/phase_fractions_by_cluster_ClustersSurprise.csv`
- `cortex__v3/03_cellcycle/proliferation_trajectory.csv`
- `human_dev__v2/03_cellcycle/cycle_scores_by_cluster_cluster_id.csv`
- `human_dev__v2/03_cellcycle/proliferation_trajectory.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/03_cellcycle/cycle_scores_by_cluster_cluster_id.csv`
- `human_dev__v3/03_cellcycle/proliferation_trajectory.csv`
- `results/08_coexpression_modules/modules.csv`

## Method

- A: per stratum x cell class, cycling fraction (cortex; cells not in Non-cycling or Post-M), mean cell-cycle score, and G1/S/G2M shares among cycling cells (where 03 exports phase fractions per class x age), vs age over points with >= 50 cells; Spearman, exact permutation p, v2/v3 signed Stouffer, BH per dataset, tiered.
- B: cluster pseudobulks as in 07/08 (log2 TMM-CPM). Proliferation association = Spearman with the cluster's cycling fraction (cortex) or mean cycle score (human_dev). S-vs-G2/M lean (cortex) = partial Spearman with S/(S+G2M) controlling for cycling fraction, over clusters with >= 5% of cells in S or G2/M. Gene classes need |rho| >= 0.3 (lean: >= 0.2) in both chemistries; cycle-independent = |rho| < 0.2 in both.
- Set profiles: mean association of a set's genes vs 2,000 random sets matched member by member on mean level x spread; v2/v3 combined, BH per dataset, tiered.
- C: 06's list age-coordination test in Radial glia, Neuronal IPC and Glioblast, restricted to cycle-independent genes, null drawn from cycle-independent genes of matched expression level.

## Key findings

- **cortex: proliferation over development in progenitors** -- fraction cycling falls with age in Radial glia (0.88->0.763 in v2, 0.746->0.6 in v3; replicated).
- **human_dev: no replicated proliferation trend in progenitors.**
- **G1 lengthening not testable yet**: phase fractions per cell class x age are exported by script 03 from this version on; re-run stage 1 to add them.
- **cortex gene phase map** (10,439 genes): 2323 proliferative (251 lean S, 397 lean G2/M), 2342 anti-proliferative, 2292 cycle-independent.
- **human_dev gene phase map** (12,018 genes): 3129 proliferative, 2903 anti-proliferative, 2655 cycle-independent.
- **Sanity check against the seed phase panels** (partly circular: phase calls come from such genes): cortex g2m_phase: lean -16.3/-18.2 null SDs (+ = S), 0 S / 47 G2/M genes; cortex s_phase: lean +12.1/+12.0 null SDs (+ = S), 33 S / 0 G2/M genes.
- **cortex: gene sets tied to proliferation, replicated** -- seed:g2m_phase proliferative, leans G2/M (rho +0.58/+0.65 vs +0.07 matched; 0 S-lean, 47 G2/M-lean genes); seed:s_phase proliferative, leans S (rho +0.75/+0.74 vs +0.04 matched; 33 S-lean, 0 G2/M-lean genes); seed:asd_high_confidence anti-proliferative (rho -0.30/-0.28 vs +0.01 matched); seed:chromatin_transcription_regulators anti-proliferative (rho -0.25/-0.22 vs +0.00 matched; 0 S-lean, 1 G2/M-lean genes); seed:epilepsy_dee anti-proliferative (rho -0.30/-0.24 vs +0.03 matched); seed:id_dd_dominant anti-proliferative (rho -0.14/-0.15 vs +0.00 matched); seed:synaptic_and_channels anti-proliferative (rho -0.53/-0.44 vs +0.04 matched).
- **cortex: proliferative modules with a replicated phase lean** -- CX01 leans G2/M; CX04 leans S; CXw02 leans G2/M; CXw04 leans S.
- **human_dev: gene sets tied to proliferation, replicated** -- seed:g2m_phase proliferative (rho +0.81/+0.78 vs +0.04 matched); seed:s_phase proliferative (rho +0.75/+0.71 vs +0.06 matched); seed:asd_high_confidence anti-proliferative (rho -0.26/-0.26 vs -0.03 matched); seed:epilepsy_dee anti-proliferative (rho -0.24/-0.24 vs +0.03 matched); seed:synaptic_and_channels anti-proliferative (rho -0.48/-0.45 vs -0.04 matched).

## Limitations

- Cluster pseudobulks pool cells of all phases, so gene-phase associations are between clusters with different phase mixes -- an ecological measure, not per-cell phase expression.
- Phase calls are derived from cell-cycle marker genes; canonical markers' phase lean is partly circular and serves only as a check.
- human_dev has no per-cell phase calls, so it has no G1/S/G2M shares or phase lean.
- Age points are donors (5-9 per chemistry); trends are across that many people.

## What would strengthen this

- Export per-cell-phase pseudobulks (cell class x phase) in stage 2 for a direct, non-ecological phase profile of every gene.
- Add an explicit G1 duration estimate (e.g. from S-phase fraction under a steady-state model) per progenitor type and age.

## Output files

- `proliferation_trajectories_per_stratum.csv` -- Per stratum x class: proliferation metric vs age (exact permutation p)
- `proliferation_trajectories_combined.csv` -- v2 x v3 combined proliferation trends per class and metric; tier replicated / supported
- `gene_phase_map_per_stratum.csv` -- Per gene x stratum: Spearman with cluster proliferation; partial Spearman with the S share of S+G2M cells (cortex), controlling for proliferation
- `gene_phase_map.csv` -- Per gene: proliferation class (|rho| >= 0.3 in both chemistries) and, for proliferative genes, S or G2/M lean (|partial rho| >= 0.2 in both; cortex)
- `set_phase_profile_per_stratum.csv` -- Per gene set x stratum: mean proliferation association / S-vs-G2M lean vs matched null
- `set_phase_profile.csv` -- Per gene set: proliferation association and S-vs-G2/M lean (v2/v3 combined, tiered), with counts of proliferative, S-, G2/M-leaning, anti-proliferative genes
- `proliferation_cortex.png` -- Progenitor proliferation vs age, cortex
- `proliferation_human_dev.png` -- Progenitor proliferation vs age, human_dev
- `phase_profile_cortex.png` -- Gene-set proliferation association vs S/G2M lean, cortex

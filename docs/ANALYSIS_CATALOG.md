# Analysis catalogue

What can be asked of these two datasets, what steps 1–2 now extract, and what is
deliberately left for step 3. Written for whoever picks up the repo next.

## The two datasets are not interchangeable

| | `cortex` | `human_dev` |
|---|---|---|
| cells × genes | 297,927 × 33,538 | 1,665,937 × 59,459 |
| scope | cortex / forebrain only, EMX1 lineage, QC-passed | whole brain, 10 regions, week 5 removed |
| `var` index | gene **symbol** | **versioned** Ensembl accession |
| cell classes | 5 (progenitor→neuron axis only) | 12 (incl. immune, vascular, erythrocyte) |
| clusterings | `Clusters`, `ClustersModularity`, `ClustersSurprise`, `leiden_scVI`, `louvain` | `cluster_id` (617) |
| latent spaces | `PCA`, `UMAP`, `UMAP3D`, `TSNE`, `X_scVI`, `nhoods` (Milo) | `Factors` + `varm['Loadings']`, `X_Embedding` |
| layers | `spliced`, `unspliced`, `ambiguous`, `matrix` | none |
| cycle phase label | `PostM` | `Post-M` |

### Chemistry is a confound, not a covariate

Both datasets mix 10x v2 and v3, and in `human_dev` the two cover **almost
disjoint age ranges** — v2 roughly 6–10 pcw, v3 roughly 5–5.5 and 11.5–14, with
only a handful of ages in both. A gene that "rises with development" in pooled
data may simply be one that v3 captures better.

So **the pipeline stratifies by default**: `--chemistry each` runs v2 and v3
separately into `csv_exports/<dataset>__v2/` and `__v3/`. `--chemistry all`
pools them into `csv_exports/<dataset>/`.

Stratifying removes the confound but costs coverage — neither chemistry spans
development alone. Script 18 quantifies that trade-off on the real file before
anything is interpreted; read `18_chemistry/age_chemistry_overlap.csv` first.

Three further traps this pipeline handles, each of which silently corrupts results if missed:

1. **`human_dev.var['gene_symbol']` holds accessions, not symbols.** The real
   symbols are in `var['Gene']`. `lib.io_utils.gene_frame` detects and corrects this.
2. **`Age`, `total_genes` and `total_UMIs` are categorical *strings* in `human_dev`.**
   A naive `.mean()` fails or sorts lexically. `coerce_numeric` forces float.
3. **The two files share no gene ids.** Symbols vs versioned accessions. The join
   key is the version-stripped accession (script 13).

## Implemented — steps 1 & 2

Everything below is extracted to CSV. Scripts 00–08 and 13 never touch `X`;
only 09 and 12 stream the matrix.

| # | Script | Analysis |
|---|---|---|
| 00 | `00_inventory_and_schema` | Ground-truth column/element inventory; regenerates the schema from the file and flags disagreement with `schemas/*.json` |
| 01 | `01_obs_overview_and_qc` | Cell counts per level of every categorical; sample and donor tables; QC distributions per class/age/region; global QC quantiles; **MAD-based sample outlier detection** |
| 02 | `02_composition` | Cell-class and cycle-phase composition across age, region, subregion, dissection, donor, sample, chemistry, cluster — as counts, fractions, **log2(observed/expected)** and **standardised χ² residuals**; Shannon/Simpson/entropy diversity; Cramér's V ranking of which grouping restructures composition most |
| 03 | `03_cellcycle` | Phase composition and continuous cycle scores per group; **proliferation trajectory** (cell class × post-conception week) — when each progenitor class exits the cycle |
| 04 | `04_cluster_profiles` | Per-cluster purity and entropy against biology *and* batch; `batch_dominated_flag` where a cluster is >90% one donor; **pairwise ARI/NMI between the alternative clusterings** |
| 05 | `05_covariate_confounds` | Cramér's V between all covariates + **explicit full-nesting detection**; η² of each QC metric per covariate; crosstabs for the pairs that matter. Age↔chemistry is nested in these data, so "developmental" effects are partly chemistry effects |
| 06 | `06_embeddings` | Per-group centroids and dispersion in every stored embedding; pairwise centroid distances; **kNN batch-mixing test** (observed vs expected same-donor neighbours); subsampled per-cell coordinates for plotting |
| 07 | `07_factor_modules` | Gene modules from `varm['Loadings']`: top-loading genes, module concentration, pairwise Jaccard overlap, and **hypergeometric enrichment of each NDD panel per module with BH-FDR** |
| 08 | `08_factor_activity` | Mean/median activity of every latent dimension per group; τ specificity across cell classes; factor×factor correlation; **factor×QC correlation flagging technical axes** (depth/mito/doublet only; cell-cycle tracking flagged separately) |
| 09 | `09_pseudobulk` | One streaming pass → summed raw counts (for edgeR/DESeq2), mean log1p(CP10K), detection fraction and CPM, per cell class / class×age / class×region / sample / donor / region / age / each clustering |
| 10 | `10_marker_specificity` | From 09's CSVs: τ, entropy specificity, log2FC vs rest, detection difference; ranked markers per group; gene co-expression matrix (the stand-in for `Loadings` on `cortex`) |
| 11 | `11_gene_panels` | Panel coverage (present vs absent — not the same as "not expressed"); tidy panel expression per group; z-scored panel scores; **marker-panel label check** against the supplied `CellClass` |
| 12 | `12_splicing_layers` | `cortex` only: spliced/unspliced counts per gene per group and the **unspliced fraction**, a directional proxy for genes being switched on or off |
| 13 | `13_cross_dataset_keys` | Gene id map, shared-gene table on the version-stripped accession, label overlap per categorical role |
| 14 | `14_normalization_diagnostics` | **TMM** scaling factors (Robinson & Oshlack) and **MA diagnostics** per group; per-group log2-CPM quantiles for box/violin plots. Flags groups whose median M departs from 0 — composition bias that CPM cannot fix |
| 15 | `15_sample_relationships` | **Correlation matrices** (Pearson + Spearman), **hierarchical clustering** under all three linkages with merge heights, and **PCA** with variance explained, sample scores and gene loadings |
| 16 | `16_expression_patterns` | **Z-scored K-means** gene clustering with a K sweep reporting homogeneity, separation and silhouette (suggested K = max silhouette, flagged when flat), per-cluster mean ± SD profiles, and the Z-score matrix for heatmaps |
| 19 | `19_rg_subtypes` | **Streams X (stage 2).** Scores every radial-glia cell for outer (HOPX, PTPRZ1, FAM107A, TNC, MOXD1, LIFR) vs ventricular (FBXO32, CTGF/CCN2, CYR61/CCN1, PALLD, PDGFD) radial-glia programmes; exports sub-type calls per age and donor, phase composition per sub-type x age (also within UMI quintiles, and threshold-free by score quintile), and pseudobulk per sub-type x age |
| 18 | `18_chemistry_comparability` | **Pooled by design.** Measures the v2/v3 overlap: cells, donors and age span per chemistry, and per level of every covariate whether both chemistries are present with enough cells to compare. Flags donor nested in chemistry |
| 17 | `17_gsea_panels` | **GSEA** of each panel against each group's ranked gene list: weighted running-sum ES, gene-set permutation p, BH-FDR, and leading-edge genes |

## Queued for step 3 — CSV-only, no h5ad needed

**Done** (see `downstream_analyses/REPORT.md`): composition vs age with donor as
the replicate unit (02); monotonic within-class age trends replicated across
chemistries (03; transient patterns not yet); NDD panel cell-class preference
against an expression-matched null, and panel overlap with age trends (04);
cortex vs human_dev and v2 vs v3 concordance (03, 05); user gene lists, covering
coverage, overlap, cell-class preference and age coordination (06); gene-list
co-expression coherence, cross-validated sub-modules and candidate members (07);
WGCNA-style co-expression modules with preservation across donor sets, hub
genes, age dynamics, list enrichment and seed-panel labels (08); cell-cycle
programs, covering proliferation and G1/S/G2M shares over development, a gene
phase map, gene-set phase profiles, and list age trends with proliferation
genes set aside (09); splicing dynamics from script 12's spliced/unspliced
pseudobulks, covering where genes and lists are being switched on or off along
the lineage, and unspliced-share age trends joined with expression trends (10);
lists at sub-type resolution, by fine cluster and along neuron sub-type axes
(11); sub-type abundance over age within classes (12); regional differences
within cell types (13); sex differences, exploratory (14); between-donor
variability against matched genes (15); a cross-donor-set age clock (16);
human_dev's latent factors (17); co-expression rewiring between early and late
clusters (18); robustness to gene length, QC metrics and dissociation stress
(19); agreement of CellClass with the files' other annotations (20, needs stage
1 from this version for its cross-tabulations).

**Part B (new extraction passes) -- built**: class x region x age and class x
phase (x age) pseudobulks in 09 (B1, B2; analyses 21, 22); per-cell programme
scores and within-cell co-expression in stage-2 script 22 (B3, B4; analyses 23,
24); lineage pseudotime from the stored latent space in stage-1 script 20 and
pseudotime pseudobulks in 09 (B5; analysis 25); Milo neighbourhood x donor
counts in stage-1 script 21 (B6; analysis 26); a tRG flag in 19 and neuron /
glial state scores in 22 (B7; analysis 27); numeric QC thresholds in the
exclusion rules, `exclusions_strict.csv` and `AIM_STRICT=true` (B8; analysis 28).

**Part C (external annotations) -- built**: `running_scripts/fetch_annotations.sh`
downloads GO / Reactome / KEGG (Enrichr), Lambert TFs and CollecTRI targets,
gnomAD v4.1 constraint, HPO gene-to-disease, OmniPath ligand-receptor pairs,
HGNC and BrainSpan into `AIM_ANNOTATIONS`. Analyses: functional enrichment
(29), TF regulators (30), constraint and disease genes (31), cell-cell
signalling (32), HGNC symbol rescue (33; also used by every list analysis),
BrainSpan validation (34). OMIM itself is licensed and not fetched.

**Neuron types -- built**: stage-2 script 23 calls every neuron and neuroblast
excitatory or inhibitory from transmitter and lineage genes; analysis 35 tracks
the inhibitory share over age within regions, compares the types at matched
region x age, and splits list age trends by type. cortex (EMX1 lineage) holds
essentially no inhibitory neurons (GAD1 in 0.2% of neurons), so the
comparisons are human_dev's.

**Knowledge-guided gene programmes (Spectra) -- built**: script 24 draws a
stratified ~25k-cell subsample per dataset x chemistry into `AIM_WORK`;
script 25 fits Spectra (Kunes et al. 2023) with a prior of the user lists, seed
NDD and cell-cycle panels, 29 core GO processes and class markers, plus free
factors (`running_scripts/slurm_spectra.sh`, one stratum per job, cpu or gpu).
Analysis 36 does not trust a factor's label, since a prior-steered factor can
echo its prior regardless of the data: it pairs v2 and v3 factors by gene
weights, tests each factor's top genes for co-expression in the other
chemistry's cluster pseudobulks (independent donors) against matched genes,
tests whether genes Spectra added to a list co-vary with the list genes it kept,
tests the NDD lists for enrichment among the top genes of a second, core-prior
fit's programmes (GO, cell cycle and markers only, so no NDD set shaped them;
against the expression-matched expectation), tracks programme activity
over age per class -- within each region as well as pooled -- and checks new
factors against the 08 modules. List programmes reproduce their lists almost
whole (the prior holds them together), so Spectra does not prune lists; the
informative programmes for the lists are those of the core-prior fit (in the
list-guided fit the lists dominate the global gene-gene graph, which pulls even
the free factors into them). A third fit adds STRING network modules
(experimental / curated-database edges only, no text mining) to the core prior;
36 tests it the same way and asks which STRING modules are rich in NDD genes and
also form a programme that replicates across chemistries and holds in held-out
donors. Next
candidate: a prior-free comparison (scHPF / cNMF).

**Data facts found along the way**: cortex's fine clusterings (`Clusters`,
`ClustersModularity`, `ClustersSurprise`) were computed per chemistry and reuse
labels for different cells, so "cluster 12" in v2 is not cluster 12 in v3;
`leiden_scVI`, `louvain` and human_dev's `cluster_id` are shared
(`_common.cluster_label_identity`). Males are few (2–4 all-male age points per
chemistry), which limits any sex analysis. The two files' cell-cycle labels
differ: human_dev calls non-cycling cells "Post-M" (its neurons are ~98% Post-M
and it has almost no "Non-cycling"), while cortex uses "Non-cycling" and keeps
"Post-M" for a small separate group; "cycling" = G1 + S + G2M works for both.

These all run off `csv_exports/`. Discover inputs via `csv_exports/<dataset>/_manifest.csv`.

**Composition statistics**
- Differential abundance across age/region with **donor as the replicate unit** —
  cells are not independent; the n is donors (~15–26), not cells (~2M).
  Dirichlet-multinomial or a `propeller`-style logit-transformed t-test.
- Composition entropy trajectory over development.

**Differential expression**
- Pseudobulk DE on `sample__pseudobulk_counts.csv` with edgeR/DESeq2,
  `~ age + chemistry` — never a per-cell test, which inflates significance by
  treating cells as replicates.
- Age-spline models of expression per cell class on `cell_class_x_age__mean_lognorm.csv`;
  classify genes as monotonic up/down, transient, or flat.

**Gene programmes**
- WGCNA-style co-expression modules from `gene_coexpression_*.csv` (script 16 gives
  K-means patterns; WGCNA's soft-thresholded topological overlap is the step up).
- Cross-dataset module preservation: do `human_dev` factor modules reproduce in `cortex`?
- TF-target inference by intersecting module membership with TF annotation.
- **GO enrichment** of the script-16 clusters — needs a GO annotation file in `panels/`;
  the hypergeometric machinery (script 07) and GSEA (script 17) are already built.
- **Bi-clustering**, for patterns shared by only a subset of groups — script 16's
  K-means requires a gene to behave consistently across *all* groups.

**Disease genetics — the point of the repo**
- **EWCE-style specificity-weighted enrichment**: are NDD genes enriched in the
  specificity distribution of particular cell classes, rather than just "expressed"?
  Uses `specificity_*.csv` and needs a matched-expression null.
- Constraint-weighted analysis (pLI/LOEUF) once gnomAD is added to `panels/`.
- Temporal windows of NDD gene expression — which developmental week carries peak risk-gene expression, per cell class.
- Overlap between ASD, ID and epilepsy panels' cellular profiles.

**Cross-dataset**
- Concordance of cortex vs `human_dev`-telencephalon pseudobulk on shared genes;
  disagreement localises annotation differences rather than biology.

**Technical**
- Variance partitioning of pseudobulk profiles into age / region / donor / chemistry.

## Needs the h5ad again — cannot come from these CSVs

State this explicitly so step 3 does not attempt it from tables:

- **Per-cell** differential expression, and any per-cell gene signature score.
- **Trajectory / pseudotime** inference (DPT, PAGA, Palantir) — needs the cell × gene matrix.
- **Full RNA velocity** (scVelo) — script 12 gives only aggregated unspliced fractions.
- **Milo differential abundance** using `cortex.obsm['nhoods']` — a cells × neighbourhoods
  sparse matrix; scripts 06/08 skip it deliberately. A dedicated script could aggregate it.
- **Cell–cell communication** (CellPhoneDB/CellChat) — needs per-cell ligand/receptor.
- Re-clustering, integration, doublet detection, ambient-RNA correction.

## Statistical cautions that apply to everything downstream

1. **Cells are not replicates.** Donors are. With 15 (`cortex`) and 26 (`human_dev`)
   donors, any p-value computed over ~2M cells is meaningless.
2. **Age is nested in chemistry** (script 05 quantifies it). Developmental effects
   and v2→v3 chemistry effects are not fully separable in these data.
3. **Composition is compositional.** Fractions sum to 1, so one cell type expanding
   forces all others down. Use the log2(O/E) tables, not raw fractions, for enrichment.
4. **Absent ≠ not expressed.** `panel_coverage.csv` distinguishes a gene missing from
   the file from a gene measured at zero.
5. **Seed gene panels are not authoritative.** See `panels/README.md`.
6. **RPKM/FPKM is deliberately never computed.** It corrects the transcript-length
   bias of full-length bulk protocols. Both these datasets are 10x UMI data counted
   from the 3' end, where that bias does not exist — dividing by gene length would
   *introduce* an artefact. Use CPM/CP10K within a sample and TMM between groups.
7. **Never compare a v2 group against a v3 group** and call the difference
   biological. Within `<dataset>__v2/` everything is chemistry-matched by
   construction; across the two folders it is not. For a contrast that must
   span chemistries, restrict to the ages script 18 marks `comparable`.
8. **GSEA p-values here come from gene-set permutation**, which does not preserve
   gene–gene correlation and is therefore anti-conservative for co-regulated panels.
   The ranking of hypotheses is sound; the absolute p is optimistic. Sample-label
   permutation would need per-donor replicates and belongs in step 3.

# Gene programmes and NDD genes in the developing brain

_As of 11 October 2026; analysis 36 and scripts 24-26._

## Summary

Neurodevelopmental-disorder (NDD) genes gather in two gene programmes of the developing brain, and otherwise do not form programmes of their own. We asked whether NDD gene lists form gene programmes when programmes are learned from single-cell data, with and without prior knowledge of gene networks. We fitted four kinds of programme model to two developing-brain datasets, each split by 10x chemistry into two donor-disjoint replicates, and kept only results seen in both.

- **Microcephaly genes sit in the progenitor proliferation programme** (human\_dev): 9.3–13.2× more than expression level predicts, in two different fits and both chemistries (BLM, CENPJ, CIT, PLK4, POC1A, STIL).
- **Synaptic NDD genes sit in a glutamate-receptor network module** (cortex): a STRING module of DLG4, GRIN1, GRIN2A, GRIA1/2, CACNG2 and GRIK2 forms a programme that replicates across chemistries and holds in held-out donors.
- **AD GWAS genes sit in microglia** (human\_dev, 10.3×/14.0×; BIN1, SPI1, MAF): known biology, recovered with no prior at all, which shows the test can find real signal.

Beyond these, the NDD lists are not concentrated in any programme more than their expression level predicts. When the lists themselves were used as the prior, Spectra simply reproduced them, so that fit cannot show list support.

## Data and design

Every fit used the same cells, and every result had to hold in two independent sets of donors.

- **Two datasets.** *cortex* is a dorsal (EMX1-lineage) cortex atlas: radial glia, IPCs, neuroblasts, neurons, glioblasts. *human\_dev* is a whole-brain developmental atlas with more classes, including erythrocytes, immune cells, vascular cells, fibroblasts and neural crest.
- **Chemistry as replicate.** Each dataset was split into 10x v2 and v3 cells. The two strata have disjoint donors, so v2 and v3 serve as two replicates, and each fit was run separately on each.
- **One subsample per stratum.** Up to 2,500 cells per cell class (classes with at least 200 cells), at most 25,000 cells, fixed seed: 12,500 cells per cortex stratum and about 24,000 per human\_dev stratum. Models used the 3,000 most variable genes, plus each fit's prior-set genes.
- **Gene lists.** The user NDD lists (ASD SFARI tiers, ASD\_HC, autism\_strong, ASD syndromic, microcephaly, BP and AD GWAS) and curated NDD panels (epilepsy/DEE, ID/DD, chromatin regulators, synaptic genes). GWAS lists count one gene per locus.

## The four programme fits

Only the three fits whose prior holds no NDD list can test the lists; the list-guided fit can only show what the prior imposes. Spectra (Kunes et al. 2023) learns programmes steered by a gene-gene graph built from prior gene sets, plus free programmes with no prior; cNMF (Kotliar et al. 2019) uses no prior at all.

| Fit | Prior | Programmes per stratum (cortex / human\_dev) | Can show | Cannot show |
| --- | --- | --- | --- | --- |
| Spectra, list-guided | NDD lists and panels, 29 core GO processes, cell-cycle panels, class markers | 60 / 71–72 | Which prior programmes the data support | Whether the lists cohere: the lists dominate the graph and even free programmes are pulled into them |
| Spectra, core prior | 29 GO processes, cell-cycle panels, class markers; 15 free global programmes | 56 / 67–68 | NDD lists in programmes they did not shape | Programmes outside the core biology the prior names |
| Spectra, STRING prior | Core prior plus 60 STRING network modules (experimental or curated evidence ≥ 700, no text mining) | 116 / 127–128 | Network modules rich in NDD genes that the data confirm as programmes | NDD genes a module already contains, as data-driven additions |
| cNMF, no prior | None; K chosen by stability, the same K for both chemistries | 25 / 30 | Programmes found with no prior knowledge | Fine programmes in rare cell types below cNMF's resolution |

Spectra runs used a plateau learning-rate rule in place of the package's own. The package's rule stopped one fit at epoch 431 of 5,000 and left the others at a learning rate of 0.001 from about epoch 150. All final fits converged, after 2,900–3,900 epochs.

## How programmes were judged

A programme's label proves nothing: a programme steered by a gene set can echo that set whatever the data say. Each programme was therefore judged on data the fit never saw.

1. **Replication.** v2 and v3 programmes are paired when each is the other's best match by cosine similarity of gene weights (≥ 0.5).
2. **Held-out support.** The mean pairwise Spearman correlation of a programme's top 50 genes across the *other* chemistry's cluster pseudobulks (other donors), against 500 random gene sets matched on expression level and spread.
3. **NDD enrichment.** For each list, the number of its genes among a programme's top 50 against the expected number. The expectation gives each list gene the chance that genes of similar expression level have of being a top gene, since NDD genes tend to be long and highly expressed in neurons; the p-value is exact (Poisson-binomial). For a programme steered by a GO set, marker panel or STRING module, that set's genes are left out first. Results are BH-corrected per dataset and count only when q < 0.05 in both chemistries.
4. **Age trends.** Spearman correlation of mean programme activity with age within each cell class, pooled over regions and within each region (human\_dev sampled different regions at different ages), combined across chemistries.

## Findings

Three NDD-gene results replicate in both chemistries; four more are significant in one chemistry and point the same way in the other.

**Replicated (q < 0.05 in v2 and v3)**

| Programme | Dataset | Fit | NDD set | Fold over expectation (v2 / v3) | Genes shared by both chemistries |
| --- | --- | --- | --- | --- | --- |
| Progenitor proliferation (free; RRM2, CDK1, MELK, ASPM, TOP2A) | human\_dev | Core prior | Microcephaly | 9.7 / 9.3 | BLM, CENPJ, CIT, H4C3, PLK4, POC1A, STIL |
| Same programme | human\_dev | STRING prior | Microcephaly | 13.2 / 10.9 | adds BRIP1, GMNN, KIF11 |
| STRING module M074, glutamate receptors and PSD (DLG4, GRIN1, GRIN2A/B, GRIA1–4, CACNG2, GRIK2) | cortex | STRING prior | Synaptic panel / autism\_strong / ASD\_SFARI\_2 | 47 / 14 / 7 (the module's own content) | DLG4, GRIN1, GRIN2A, GRIA1, GRIA2, CACNG2, GRIK2 |
| Microglia (free; P2RY12, AIF1, GPR34) | human\_dev | cNMF | AD GWAS | 10.3 / 14.0 | BIN1, MAF, SPI1 |

The proliferation programme is the most robust: two priors, both chemistries, and its top genes co-express in held-out donors at +56 to +63 SD. The synapse module's NDD content is a property of the network; what the data add is that the module forms a programme that replicates and holds in held-out donors (+13.9 SD in v2, +3.5 SD in v3).

**Suggestive (one chemistry)**

- Chromatin regulators in cortex: programmes seeded by two STRING chromatin modules pull in further autism genes beyond the modules' own (AUTS2, ATRX, CHD2, CHD7, KDM5C, SETD2, SATB1); v3 only.
- Syndromic ASD genes in a human\_dev cortical-neuron programme (cNMF; BCL11A, NFIB, SOX5, ZBTB18): q 0.07 / 0.03.
- Microcephaly genes in proliferation in cNMF (human\_dev, 5.3× / 8.4×, v3 significant), and in cortex a DNA-replication programme with the Meier–Gorlin genes CDT1, GINS2, GMNN and MCM7 (3.7× in both, not significant).
- autism\_strong in a free human\_dev neuron programme marking midbrain and hindbrain inhibitory neurons (LHX1, GATA3, SOX14; MYT1L, DSCAM, SCN2A, SHANK2): 4.0× in v3, 1.8× in v2.

**Programmes in general.** Most programmes are real: 41–72 of each Spectra fit's and 23–28 of each cNMF fit's are supported in held-out donors. The Spectra core and STRING fits recover about half of cNMF's cortex programmes and two thirds of its human\_dev ones; what they miss is mostly rare cell types (multiciliated cells, Cajal–Retzius cells, cortical hem, astrocytes) and stress or dissociation programmes (FOS, EGR1, HSPA1A). Neither replicated NDD programme changes with age consistently across chemistries.

## What is circular, and caveats

The list-guided fit cannot test the lists, and several early results were artefacts of the method rather than biology.

- **Lists reproduce themselves.** In the list-guided fit, 35 of 50 list programmes had top 50 genes made only of list genes, and 92% of list genes outweighed 95% of non-list genes. Spectra does not prune a list it is given.
- **Free programmes are not free under a big prior.** Every global programme, free ones included, is fitted to one gene-gene graph; ASD\_SFARI\_2 alone puts about 150,000 gene pairs in it. Two "free" programmes ended up as SFARI lists (top-gene overlap 0.90 and 1.00), so their enrichment was built in. This is why the core- and STRING-prior fits exist.
- **Overlap with the steering set.** A GO-steered programme looks enriched for an NDD panel that shares its genes (synaptic panel in the GO synaptic-transmission programme, about 11×). The test now leaves out the steering set's genes.
- **GWAS loci.** Counting each gene of a GWAS locus made the chr6 histone cluster look like 8–18 BP GWAS genes in a cortex histone programme. GWAS lists now count one gene per locus.

Caveats:

- Only each programme's top 50 genes are tested, so a list spread thinly across a programme would be missed.
- Each fit uses a subsample of about 12,500–25,000 cells per stratum, so rare cell types contribute few cells, and per-region age points can rest on a few dozen cells.
- human\_dev sampled different brain regions at different ages, so age trends pooled over regions can follow sampling; within-region trends have 5–7 age points per chemistry.
- STRING modules were built from experimental and curated-database evidence only; using text mining would bias them towards genes already studied with NDDs.

## Where the results live

Everything is in the ai-adata-miner repository: analysis 36 (`downstream_analyses/36_spectra_programs.py`) reads the fits and writes `results/36_spectra_programs/`, with a SUMMARY.md and one CSV per table.

| Step | Script | Output |
| --- | --- | --- |
| Cell subsample | `scripts_generated/24_spectra_input.py` | `aim_work/spectra/<dataset>__<chemistry>/input.h5ad` |
| Spectra fits | `scripts_generated/25_spectra_fit.py` (`SPECTRA_PRIOR=full`, `core` or `string`) | `csv_exports/<stratum>/25_spectra/`, `25_spectra_core/`, `25_spectra_string/` |
| cNMF fit | `scripts_generated/26_cnmf.py`, then `--finalise` for one K per dataset | `csv_exports/<stratum>/26_cnmf/` |
| Tests and comparison | `downstream_analyses/36_spectra_programs.py` (step 3) | `core_*`, `string_*`, `cnmf_*` CSVs; `string_module_ndd_content.csv`; `cnmf_vs_spectra.csv` |

On powerslurm, `running_scripts/slurm_spectra.sh` and `slurm_cnmf.sh` run one stratum per job; their headers and `running_scripts/README.md` give the exact commands. STRING v12 files come from `running_scripts/fetch_annotations.sh`.

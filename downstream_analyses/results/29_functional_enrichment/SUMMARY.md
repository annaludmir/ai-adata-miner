# Functional enrichment of modules, sub-modules, age-trend genes and phase classes

_Generated 2026-10-07 11:39 UTC by `downstream_analyses/29_functional_enrichment.py` from `csv_exports/`._

## Question

Which biological processes, components, functions and pathways do the modules, list sub-modules, age-trend genes, proliferation classes and user lists carry?

## Inputs

- `annotations/gmt/GO_Biological_Process_2023.gmt`
- `annotations/gmt/GO_Cellular_Component_2023.gmt`
- `annotations/gmt/GO_Molecular_Function_2023.gmt`
- `annotations/gmt/KEGG_2021_Human.gmt`
- `annotations/gmt/Reactome_2022.gmt`
- `cortex__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `cortex__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `human_dev__v3/09_pseudobulk/cell_class__pseudobulk_counts.csv`
- `results/03_age_trends_within_cell_class/age_trends_combined.csv`
- `results/08_coexpression_modules/modules.csv`
- `results/09_cell_cycle_programs/gene_phase_map.csv`

## Method

- Libraries: GO BP / CC / MF 2023, Reactome 2022, KEGG 2021 (Enrichr .gmt, fetch_annotations.sh); term genes mapped to each dataset's symbols.
- Universe: expressed genes (mean CPM >= 1); terms of 10-500 universe genes; hypergeometric p, BH within each query group.

## Key findings

- **Libraries**: GO_Biological_Process_2023 (5,407 terms), GO_Cellular_Component_2023 (474 terms), GO_Molecular_Function_2023 (1,147 terms), KEGG_2021_Human (320 terms), Reactome_2022 (1,818 terms); query groups tested: 54.
- **Modules -- top terms** (overlap, fold, q): cortex CX01 (G2/M phase): Cell Cycle, Mitotic R-HSA-69278 (80, 9.1x, q 1.7e-52); M Phase R-HSA-68886 (63, 9.9x, q 2.2e-42); Mitotic Prometaphase R-HSA-68877 (46, 13.7x, q 8.4e-37) | cortex CX02: Neuron Projection (35, 5.1x, q 5.7e-12); Chemical Synaptic Transmission (24, 8.0x, q 5.7e-12); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.7x, q 8.7e-12) | cortex CX03 (radial glia markers): Regulation Of Epithelial Cell Proliferation (9, 13.6x, q 7.2e-05); Collagen-Containing Extracellular Matrix (13, 6.3x, q 3.6e-04); Extracellular Matrix Organization R-HSA-1474244 (12, 6.2x, q 8.4e-04) | cortex CX04 (S phase): DNA Metabolic Process (39, 11.8x, q 1.4e-27); Cell Cycle, Mitotic R-HSA-69278 (46, 7.6x, q 7.1e-25); DNA-templated DNA Replication (23, 26.8x, q 2.7e-24) | cortex CX05: Cardiac Conduction R-HSA-5576891 (7, 10.6x, q 0.016); Muscle Contraction R-HSA-397014 (8, 8.0x, q 0.016) | cortex CXw01: Neuronal System R-HSA-112316 (31, 6.1x, q 1.8e-12); Neuron Projection (36, 4.8x, q 9.5e-12); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.1x, q 6.2e-11) | cortex CXw02 (G2/M phase): Cell Cycle, Mitotic R-HSA-69278 (75, 9.7x, q 4.1e-51); M Phase R-HSA-68886 (59, 10.5x, q 4.9e-41); Mitotic Prometaphase R-HSA-68877 (44, 14.8x, q 1.3e-36) | cortex CXw03 (radial glia markers): Collagen-Containing Extracellular Matrix (16, 7.5x, q 1.7e-06); Regulation Of Epithelial Cell Proliferation (9, 13.3x, q 4.4e-05); Negative Regulation Of Cell Differentiation (13, 6.6x, q 1.3e-04) | cortex CXw04 (S phase): DNA Metabolic Process (35, 15.9x, q 1.8e-29); Cell Cycle, Mitotic R-HSA-69278 (41, 10.1x, q 1.2e-27); Cell Cycle Checkpoints R-HSA-69620 (32, 15.0x, q 3.4e-26) | cortex CXw05: Neuron Projection (14, 5.5x, q 9.0e-04); Central Nervous System Development (9, 6.4x, q 0.026); Inorganic Cation Import Across Plasma Membrane (5, 15.4x, q 0.027) | human_dev HD01 (G2/M phase; S phase): Cell Cycle, Mitotic R-HSA-69278 (102, 9.8x, q 3.4e-72); Cell Cycle Checkpoints R-HSA-69620 (67, 12.6x, q 7.3e-53); Mitotic Prometaphase R-HSA-68877 (48, 12.0x, q 1.1e-35) | human_dev HD02: Collagen-Containing Extracellular Matrix (17, 4.7x, q 5.0e-04); Hemostasis R-HSA-109582 (24, 3.3x, q 5.0e-04); Vesicle (16, 4.7x, q 5.0e-04) | human_dev HD03: Neuron Projection (41, 6.0x, q 3.5e-17); Transmission Across Chemical Synapses R-HSA-112315 (24, 7.6x, q 1.9e-11); Neuronal System R-HSA-112316 (27, 5.9x, q 1.3e-10) | human_dev HD05: Cilium Movement (9, 96.8x, q 5.6e-13); Cilium (12, 19.6x, q 1.2e-09); Axoneme Assembly (5, 66.6x, q 1.4e-05) | human_dev HDw01 (G2/M phase; S phase): Cell Cycle, Mitotic R-HSA-69278 (79, 10.5x, q 2.4e-57); Cell Cycle Checkpoints R-HSA-69620 (56, 14.5x, q 3.9e-47); M Phase R-HSA-68886 (54, 10.0x, q 3.6e-36) | human_dev HDw02: Cilium (14, 20.5x, q 1.2e-11); Cilium Movement (8, 77.1x, q 1.2e-10); Axoneme Assembly (6, 71.6x, q 2.4e-07) | (+1 more groups).
- **Age trends -- top terms** (overlap, fold, q): cortex Neuroblast: up: Ventricular Septum Development (4, 30.3x, q 0.035); Cohesin Loading Onto Chromatin R-HSA-2470946 (3, 59.0x, q 0.035); Outflow Tract Septum Morphogenesis (3, 49.2x, q 0.042) | cortex Neuronal IPC: down: Transcriptional Regulation By RUNX3 R-HSA-8878159 (5, 16.0x, q 0.036); Ubiquitin-dependent Degradation Of Cyclin D R-HSA-75815 (4, 23.2x, q 0.036); Hh Mutants Abrogate Ligand Secretion R-HSA-5387390 (4, 22.2x, q 0.036) | cortex Radial glia: down: Glycolysis / Gluconeogenesis (8, 12.7x, q 3.2e-04); S Phase R-HSA-69242 (14, 5.7x, q 3.2e-04); Synthesis Of DNA R-HSA-69239 (12, 6.7x, q 3.2e-04) | cortex Radial glia: up: Microtubule Binding (21, 7.3x, q 4.5e-09); Tubulin Binding (22, 5.8x, q 7.3e-08); Microtubule Cytoskeleton (22, 5.4x, q 1.7e-07) | human_dev Erythrocyte: up: RNA Polymerase III Transcription Termination R-HSA-73980 (4, 53.4x, q 0.0038); RNA Polymerase III Abortive And Retractive Initiation R-HSA-749476 (4, 28.1x, q 0.028) | human_dev Immune: up: Regulation Of Actin Filament Polymerization (5, 20.4x, q 0.02) | human_dev Neuroblast: up: Neuronal System R-HSA-112316 (44, 2.4x, q 1.5e-04); Axon Guidance (25, 2.8x, q 0.0042); Axonogenesis (29, 2.5x, q 0.0052) | human_dev Neuron: up: Selenocysteine Synthesis R-HSA-2408557 (28, 4.3x, q 5.2e-08); Peptide Chain Elongation R-HSA-156902 (26, 4.2x, q 2.2e-07); Selenoamino Acid Metabolism R-HSA-2408522 (29, 3.8x, q 2.2e-07) | human_dev Neuronal IPC: up: Neuron Projection Development (27, 3.0x, q 0.0012); Regulation Of Dendrite Development (11, 5.7x, q 0.0031); Cytoskeleton (51, 1.9x, q 0.0097) | human_dev Radial glia: up: Tubulin Binding (30, 2.4x, q 0.031) | human_dev Vascular: down: Purine-Containing Compound Biosynthetic Process (4, 35.8x, q 0.018); Cholesterol Biosynthesis R-HSA-191273 (4, 26.0x, q 0.034); Sterol Biosynthetic Process (4, 23.9x, q 0.034).
- **Proliferations -- top terms** (overlap, fold, q): cortex anti-proliferative: Neuron Projection (169, 2.3x, q 1.8e-26); Neuronal System R-HSA-112316 (116, 2.4x, q 1.2e-18); Transmission Across Chemical Synapses R-HSA-112315 (85, 2.6x, q 2.6e-16) | cortex proliferative: Cell Cycle, Mitotic R-HSA-69278 (226, 2.5x, q 1.2e-43); M Phase R-HSA-68886 (147, 2.2x, q 3.2e-21); Mitotic Prometaphase R-HSA-68877 (94, 2.7x, q 2.3e-20) | cortex proliferative, G2/M-leaning: Cell Cycle, Mitotic R-HSA-69278 (97, 6.3x, q 8.3e-48); M Phase R-HSA-68886 (79, 7.0x, q 4.1e-42); Mitotic Prometaphase R-HSA-68877 (54, 9.1x, q 3.3e-34) | cortex proliferative, S-leaning: DNA Metabolic Process (49, 9.2x, q 4.4e-30); DNA-templated DNA Replication (29, 21.0x, q 2.1e-28); DNA Strand Elongation R-HSA-69190 (19, 33.7x, q 1.9e-23) | human_dev anti-proliferative: Neuron Projection (188, 2.1x, q 9.1e-27); Transmission Across Chemical Synapses R-HSA-112315 (105, 2.6x, q 8.0e-23); Neuronal System R-HSA-112316 (135, 2.3x, q 8.0e-23) | human_dev proliferative: Cell Cycle, Mitotic R-HSA-69278 (247, 2.1x, q 1.7e-34); DNA Metabolic Process (143, 2.2x, q 4.6e-23); Processing Of Capped Intron-Containing Pre-mRNA R-HSA-72203 (136, 2.2x, q 1.0e-22).

## Limitations

- Terms overlap heavily (GO's hierarchy), so several top terms often describe one signal.
- Enrichment ignores expression level within the universe; highly expressed housekeeping processes (ribosome, translation) can surface for groups of highly expressed genes.

## What would strengthen this

- Collapse redundant GO terms (e.g. by semantic similarity) for a shorter summary.

## Output files

- `enrichment.csv` -- Per query group x term: overlap, fold, hypergeometric p, BH q within the query (rows with q < 0.25 and >= 3 genes)
- `top_terms.csv` -- Top 5 terms (q < 0.05) per query group

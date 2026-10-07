# Ligand-receptor signalling between cell types over development

_Generated 2026-10-07 11:40 UTC by `downstream_analyses/32_cell_communication.py` from `csv_exports/`._

## Question

Which cell types could signal to which through ligand-receptor pairs, how does that potential change with age in both donor sets, and where do NDD genes act as ligands or receptors?

## Inputs

- `annotations/lr/intercell_ligand_receptor.tsv`
- `annotations/lr/omnipath_interactions.tsv`
- `cortex__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `cortex__v2/11_panels/panel_coverage.csv`
- `cortex__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`
- `human_dev__v2/11_panels/panel_coverage.csv`
- `human_dev__v3/09_pseudobulk/cell_class_x_age__pseudobulk_counts.csv`

## Method

- Pairs: OmniPath interactions, ligand (secreted / membrane) -> receptor (plasma membrane), each role annotated by >= 4 intercell resources; symbols mapped to each dataset.
- Score = sqrt(ligand CPM in sender x receptor CPM in receiver), CPM from TMM over all class x age pseudobulks of a stratum; expressed = both >= 10 CPM.
- Age trend per pair x route expressed at >= half of >= 5 shared age points: Spearman, exact permutation; v2 x v3 signed Stouffer, BH per dataset, tiered.

## Key findings

- **cortex: busiest sender -> receiver routes** (pairs expressed per age, mean of chemistries): Neuroblast -> Radial glia 174, Neuron -> Neuronal IPC 170, Neuron -> Radial glia 169, Neuroblast -> Neuronal IPC 169, Neuroblast -> Neuroblast 162, Neuroblast -> Glioblast 161.
- **human_dev: busiest sender -> receiver routes** (pairs expressed per age, mean of chemistries): Glioblast -> Glioblast 327, Fibroblast -> Glioblast 320, Fibroblast -> Fibroblast 315, Radial glia -> Glioblast 314, Glioblast -> Radial glia 314, Glioblast -> Oligo 303.
- **cortex: replicated changes in signalling potential** (174 of 1605 pair x route tests): Neuroblast -> Radial glia: 16 rise / 7 fall; Neuroblast -> Neuronal IPC: 4 rise / 5 fall; Neuron -> Radial glia: 5 rise / 5 fall; Neuroblast -> Neuroblast: 10 rise / 4 fall; Radial glia -> Neuroblast: 7 rise / 4 fall; Neuroblast -> Neuron: 14 rise / 3 fall; Neuron -> Neuronal IPC: 5 rise / 3 fall; Neuronal IPC -> Neuronal IPC: 5 rise / 3 fall.
- **cortex: replicated changes involving NDD genes** (ligand -> receptor, sender -> receiver; rho v2/v3): CNTN2 -> CNTNAP2 (Neuroblast -> Neuronal IPC) falls with age (-0.96/-0.89); CNTN2 -> CNTNAP2 (Neuroblast -> Neuroblast) falls with age (-0.93/-0.94); CNTN2 -> CNTNAP2 (Neuroblast -> Neuron) falls with age (-0.86/-1.00).
- **human_dev: replicated changes in signalling potential** (1802 of 10094 pair x route tests): Radial glia -> Radial glia: 26 rise / 45 fall; Neuronal IPC -> Radial glia: 21 rise / 31 fall; Neuroblast -> Radial glia: 28 rise / 31 fall; Neuron -> Radial glia: 31 rise / 29 fall; Radial glia -> Neuronal IPC: 30 rise / 29 fall; Radial glia -> Vascular: 22 rise / 27 fall; Radial glia -> Neuron: 30 rise / 26 fall; Radial glia -> Neuroblast: 39 rise / 24 fall.
- **human_dev: replicated changes involving NDD genes** (ligand -> receptor, sender -> receiver; rho v2/v3): NLGN1 -> NRXN1 (Radial glia -> Radial glia) rises with age (+0.96/+1.00); NLGN1 -> NRXN1 (Neuroblast -> Neuroblast) rises with age (+0.94/+0.98); NLGN1 -> NRXN1 (Neuroblast -> Neuron) rises with age (+0.96/+0.93); NLGN1 -> NRXN1 (Neuroblast -> Neuronal IPC) rises with age (+0.93/+0.98); TAFA2 -> NRXN1 (Radial glia -> Neuronal IPC) rises with age (+0.93/+0.98); NLGN1 -> NRXN1 (Neuronal IPC -> Radial glia) rises with age (+0.85/+1.00); NLGN1 -> NRXN1 (Neuronal IPC -> Neuronal IPC) rises with age (+0.93/+0.95); TAFA2 -> NRXN1 (Neuroblast -> Neuronal IPC) rises with age (+0.90/+0.95); NLGN1 -> NRXN1 (Radial glia -> Glioblast) rises with age (+0.92/+0.96); NLGN1 -> NRXN1 (Neuroblast -> Radial glia) rises with age (+0.86/+0.98) (+97 more).

## Limitations

- Expression of a ligand and its receptor is potential, not signalling: protein levels, processing, localisation and spatial proximity are not measured.
- Complex receptors (multiple subunits) are reduced to pairwise interactions.
- Within a class, a trend can come from sub-type mix (as everywhere in 03).

## What would strengthen this

- Spatial data to check which sender-receiver pairs are neighbours.

## Output files

- `sender_receiver_overview.csv` -- Per stratum x sender x receiver: mean number of pairs with both partners >= 10 CPM per age point
- `interaction_age_trends_per_stratum.csv` -- Per stratum x pair x sender x receiver: Spearman of log2 score with age
- `interaction_age_trends_combined.csv` -- v2 x v3 combined; BH per dataset; tier; NDD gene flag
- `ndd_ligands_receptors.csv` -- NDD list / seed-panel genes that are ligands or receptors, with partners

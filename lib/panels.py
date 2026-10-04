"""Gene panels used to focus the expression extractions.

PROVENANCE WARNING
------------------
These are *seed* lists curated from domain knowledge so the pipeline produces
something useful on first run.  They are NOT authoritative releases.  Before
publishing anything, replace them with versioned downloads:
  * SFARI Gene       -> https://gene.sfari.org  (score 1/2, syndromic flag)
  * DDG2P / G2P      -> https://www.ebi.ac.uk/gene2phenotype
  * Epi25 / ClinGen epilepsy panels
  * gnomAD pLI / LOEUF constraint for the constraint-weighted analyses

`load_panels()` merges any JSON/CSV in panels/ over these defaults, so adding a
real SFARI export needs no code change.

Marker panels are deliberately included alongside the disease panels: they let
the extraction sanity-check the provided CellClass labels instead of trusting
them, which matters because the two datasets were annotated separately.
"""
from __future__ import annotations

import json
from pathlib import Path

import config

PANEL_DIR = config.REPO_ROOT / "panels"

# ---------------------------------------------------------------------------
# Neurodevelopmental disorder risk genes (seed lists -- see warning above)
# ---------------------------------------------------------------------------
NDD_PANELS: dict[str, list[str]] = {
    "asd_high_confidence": [
        "CHD8", "SCN2A", "SYNGAP1", "ADNP", "FOXP1", "POGZ", "ARID1B", "SHANK3",
        "DYRK1A", "GRIN2B", "PTEN", "TBR1", "ANK2", "CHD2", "KDM5B", "MED13L",
        "SETD5", "KMT2C", "ASH1L", "DSCAM", "TRIP12", "CTNNB1", "SHANK2",
        "NRXN1", "RELN", "MECP2", "TSC1", "TSC2", "NF1", "FMR1", "CNTNAP2",
        "SLC6A1", "STXBP1", "GIGYF1", "DEAF1", "CREBBP", "SRCAP", "WAC",
    ],
    "id_dd_dominant": [
        "ARID1B", "ANKRD11", "KMT2A", "KMT2D", "EP300", "CREBBP", "SETD5",
        "DYRK1A", "MED13L", "SATB2", "TCF4", "EHMT1", "SMARCA2", "SMARCB1",
        "ASXL1", "ASXL3", "KAT6A", "KAT6B", "NSD1", "CHD7", "DDX3X", "PURA",
        "GATAD2B", "PACS1", "PPP2R5D", "SON", "HNRNPU", "ZEB2", "FOXG1",
    ],
    "epilepsy_dee": [
        "SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3", "KCNT1", "STXBP1",
        "CDKL5", "PCDH19", "GABRA1", "GABRB3", "GABRG2", "GRIN1", "GRIN2A",
        "GRIN2B", "DEPDC5", "TSC1", "TSC2", "SLC2A1", "ARX", "SPTAN1", "WWOX",
    ],
    "chromatin_transcription_regulators": [
        "CHD8", "CHD2", "CHD7", "ARID1B", "ARID1A", "SMARCA4", "SMARCB1",
        "KMT2A", "KMT2C", "KMT2D", "KDM5B", "KDM6A", "KDM6B", "SETD5", "NSD1",
        "EP300", "CREBBP", "ASH1L", "EZH2", "DNMT3A", "DNMT3B", "TET2",
        "HDAC4", "SIN3A", "BCL11A", "BCL11B", "MECP2", "ADNP", "POGZ",
    ],
    "synaptic_and_channels": [
        "SYNGAP1", "SHANK2", "SHANK3", "NRXN1", "NRXN2", "NLGN3", "NLGN4X",
        "DLG4", "GRIN2B", "GRIA1", "GRIA2", "CACNA1A", "CACNA1C", "CACNA1E",
        "SCN1A", "SCN2A", "SCN8A", "KCNQ2", "STXBP1", "SYN1", "SNAP25",
        "SLC6A1", "GABRB2", "CNTNAP2", "ANK2", "ANK3",
    ],
}

# ---------------------------------------------------------------------------
# Cell-identity markers -- used to validate the supplied CellClass labels
# ---------------------------------------------------------------------------
MARKER_PANELS: dict[str, list[str]] = {
    "radial_glia":   ["VIM", "SOX2", "PAX6", "HES1", "HES5", "SLC1A3", "FABP7",
                      "NES", "GLI3", "TNC", "PTPRZ1", "HOPX", "CRYAB", "MOXD1"],
    "neuronal_ipc":  ["EOMES", "NEUROG1", "NEUROG2", "PPP1R17", "NEUROD4",
                      "ASCL1", "GADD45G", "SSTR2", "PENK"],
    "neuroblast":    ["NEUROD2", "NEUROD6", "SOX11", "SOX4", "DCX", "STMN2",
                      "NHLH1", "TBR1", "BCL11B"],
    "neuron":        ["SNAP25", "SYT1", "RBFOX3", "MAP2", "STMN2", "NEFL",
                      "GAP43", "SLC17A7", "GAD1", "GAD2", "DLX1", "DLX2"],
    "glioblast_opc": ["OLIG1", "OLIG2", "PDGFRA", "SOX10", "EGFR", "APOE",
                      "AQP4", "GFAP", "S100B", "CSPG4"],
    "oligo":         ["MBP", "PLP1", "MOG", "MAG", "CNP", "SOX10"],
    "immune_microglia": ["AIF1", "C1QA", "C1QB", "P2RY12", "CX3CR1", "PTPRC",
                         "SPI1", "CSF1R"],
    "vascular_endothelial": ["CLDN5", "PECAM1", "FLT1", "RGS5", "PDGFRB",
                             "COL4A1", "DCN", "COL1A1"],
    "erythrocyte":   ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "ALAS2"],
    "neural_crest":  ["SOX10", "FOXD3", "TFAP2A", "TFAP2B", "PHOX2B"],
}

# ---------------------------------------------------------------------------
# Cell cycle (Tirosh et al. 2016 core sets) and regional patterning TFs
# ---------------------------------------------------------------------------
CELL_CYCLE_PANELS: dict[str, list[str]] = {
    "s_phase": [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
        "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2",
        "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
        "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
        "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
    ],
    "g2m_phase": [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
        "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4",
        "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
        "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK", "CDC25C",
        "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
        "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE",
        "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ],
}

PATTERNING_PANELS: dict[str, list[str]] = {
    "dorsoventral_telencephalon": ["EMX1", "EMX2", "PAX6", "LHX2", "NR2F1",
                                   "GSX2", "NKX2-1", "DLX1", "DLX2", "OLIG2",
                                   "SHH", "GLI3", "WNT3A", "RSPO3", "TTR"],
    "anteroposterior": ["FOXG1", "SIX3", "OTX1", "OTX2", "EN1", "EN2", "GBX2",
                        "HOXA2", "HOXB2", "HOXB4", "IRX3", "PAX2", "PAX7"],
    "cortical_layers": ["RELN", "TBR1", "BCL11B", "FEZF2", "SATB2", "CUX1",
                        "CUX2", "POU3F2", "RORB", "FOXP2", "NR4A2", "LHX6"],
    "cortical_hem_signalling": ["WNT3A", "WNT2B", "WNT5A", "BMP4", "BMP7",
                                "RSPO1", "RSPO2", "RSPO3", "TTR", "LMX1A"],
}

ALL_PANEL_GROUPS: dict[str, dict[str, list[str]]] = {
    "ndd": NDD_PANELS,
    "marker": MARKER_PANELS,
    "cell_cycle": CELL_CYCLE_PANELS,
    "patterning": PATTERNING_PANELS,
}


USER_LIST_GROUP = "user_lists"
_GENE_COLUMNS = ("gene", "Gene", "genes", "symbol", "Symbol", "gene_symbol", "gene_name")


def read_gene_list(path) -> list[str]:
    """One gene list file -> unique genes in file order.

    CSV/TSV: the first column named like a gene column ('gene', 'symbol', ...),
    otherwise the first column. Anything else: one gene per line (first
    tab/comma field), '#' lines skipped. Same conventions as ndd_gene_modules.
    """
    import pandas as pd
    path = Path(path)
    if path.suffix.lower() in (".csv", ".tsv"):
        df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")
        col = next((c for c in _GENE_COLUMNS if c in df.columns), df.columns[0])
        raw = df[col].tolist()
    else:
        raw = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                raw.append(line.replace(",", "\t").split("\t")[0])
    seen, out = set(), []
    for g in raw:
        g = str(g).strip()
        if g and g.lower() not in ("nan", "none") and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def load_user_lists(folder=None) -> dict[str, list[str]]:
    """Every gene list in config.GENE_LISTS_DIR (or `folder`), keyed by file stem."""
    folder = Path(folder or config.GENE_LISTS_DIR)
    if not folder.is_dir():
        return {}
    lists = {}
    for path in sorted(folder.iterdir()):
        if path.suffix.lower() in (".csv", ".tsv", ".txt") and not path.name.startswith("."):
            genes = read_gene_list(path)
            if genes:
                lists[path.stem] = genes
    return lists


def load_panels() -> dict[str, dict[str, list[str]]]:
    """Built-in panels, user gene lists, and any file in panels/ merged on top.

    User gene lists (config.GENE_LISTS_DIR) form the group 'user_lists'.
    A JSON file in panels/ must map panel-name -> list of symbols; it is filed
    under the group named by its stem (e.g. panels/ndd.json overrides 'ndd').
    A CSV in panels/ must have columns `panel,gene`.
    """
    panels = {g: {k: list(v) for k, v in d.items()} for g, d in ALL_PANEL_GROUPS.items()}
    user = load_user_lists()
    if user:
        panels[USER_LIST_GROUP] = user
    if not PANEL_DIR.exists():
        return panels

    for path in sorted(PANEL_DIR.glob("*.json")):
        group = path.stem
        data = json.loads(path.read_text())
        panels.setdefault(group, {}).update({k: list(v) for k, v in data.items()})
    for path in sorted(PANEL_DIR.glob("*.csv")):
        import pandas as pd
        df = pd.read_csv(path)
        if {"panel", "gene"} <= set(df.columns):
            group = path.stem
            grouped = df.groupby("panel")["gene"].apply(lambda s: sorted(set(s.astype(str))))
            panels.setdefault(group, {}).update(grouped.to_dict())
    return panels


def panel_long_frame():
    """All panels as a tidy group/panel/gene table."""
    import pandas as pd
    rows = [
        {"panel_group": group, "panel": name, "gene": gene}
        for group, d in load_panels().items()
        for name, genes in d.items()
        for gene in genes
    ]
    return pd.DataFrame(rows).drop_duplicates()


def all_panel_genes() -> set[str]:
    return {g for d in load_panels().values() for genes in d.values() for g in genes}

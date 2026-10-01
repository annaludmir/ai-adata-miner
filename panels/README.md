# Replacing the seed gene panels

`lib/panels.py` ships *seed* lists curated from domain knowledge so the pipeline
produces something on first run. They are not authoritative. Drop real exports
here and they override the defaults with no code change:

| file | format | overrides |
|---|---|---|
| `ndd.json` | `{"panel_name": ["GENE1", ...]}` | the `ndd` panel group |
| `marker.json` | same | the `marker` group |
| `<anything>.csv` | columns `panel,gene` | the group named by the file stem |

Recommended real sources:

- **SFARI Gene** — <https://gene.sfari.org> (filter to score 1/2 + syndromic)
- **DDG2P / G2P** — <https://www.ebi.ac.uk/gene2phenotype>
- **ClinGen / Epi25** epilepsy panels
- **gnomAD constraint** (pLI, LOEUF) for constraint-weighted analyses

Use gene **symbols** matching `var['Gene']`; the pipeline maps to Ensembl
accessions itself via `lib/io_utils.gene_frame`.

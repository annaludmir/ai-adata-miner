# Prompts

System prompts that drive step 3 — AI-generated analyses over the CSV layer.

- `analysis_generators/00_system_context.md` — shared context: what CSVs exist,
  how to discover them via the manifest, and the statistical rules that must hold.
  Prepend to everything.
- `analysis_generators/01_generate_analysis.md` — template for generating one
  analysis. Substitute `{{QUESTION}}` and `{{DATASETS}}`.

Candidate questions are listed under "Queued for step 3" in
[`../docs/ANALYSIS_CATALOG.md`](../docs/ANALYSIS_CATALOG.md).

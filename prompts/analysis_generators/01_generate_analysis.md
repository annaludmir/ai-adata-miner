# Prompt: generate one downstream analysis (step 3)

{{INCLUDE 00_system_context.md}}

## Task

Write and run a single, self-contained analysis answering:

> {{QUESTION}}

Target dataset(s): {{DATASETS}}

## Procedure

1. Read `csv_exports/<dataset>/_manifest.csv` and list the exact CSVs you will use.
   If the question cannot be answered from what exists, say so and stop — propose
   which extraction script would need extending instead of fabricating a table.
2. State your approach in three sentences, including the replicate unit and the
   confounders you are controlling for.
3. Write `downstream_analyses/<nn>_<slug>.py`. It must be re-runnable, deterministic
   (seed everything), and write results to `downstream_analyses/results/<nn>_<slug>/`.
4. Run it. If it fails, fix it and run again — do not report untested code.
5. Write `downstream_analyses/results/<nn>_<slug>/SUMMARY.md`:
   - the question
   - the tables used (exact paths)
   - method, with the replicate unit stated explicitly
   - findings, with effect sizes and n(donors)
   - **limitations**, naming any confound from `05_confounds/`
   - what would be needed to strengthen the result

## Quality bar

A finding is only reported if it survives: the donor-level replicate rule, a
check against `confound_warnings.csv`, and a sanity check that it is not driven
by one donor or one sample (`04_clusters/cluster_profile_*.csv` has the flags).
Prefer one defensible result over five suggestive ones.

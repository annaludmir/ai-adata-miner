#!/usr/bin/env python3
"""Collect every results/<nn>_<slug>/SUMMARY.md into one REPORT.md.

The report carries each analysis's key findings and links to its full summary,
so one file answers "what did step 3 find?" and every claim stays traceable to
the script and tables behind it.
"""
from __future__ import annotations

import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _common as C

REPORT = C.RESULTS.parent / "REPORT.md"

PREAMBLE = """\
Read this first. Three facts about the data constrain every result below
(details in 01):

1. **Donors are the replicate unit, and each donor has one age.** Within a
   chemistry there are 7-15 donors, mostly one per age. A trend over age is a
   trend across those few people.
2. **v2 and v3 are disjoint donor sets covering different ages.** They are
   analysed separately and used to replicate each other; nothing here pools
   them. "Replicated" means the same direction and nominal significance in
   both, plus a combined FDR < 0.05.
3. **cortex is a subset of the human_dev donors.** Agreement between the two
   files is not independent evidence. cortex is the primary dataset for
   within-cell-type questions; human_dev classes pool brain regions whose
   sampling changes with age.
"""


def github_anchor(heading: str) -> str:
    """The id GitHub gives a markdown heading, so the contents links work."""
    return re.sub(r"[^a-z0-9 -]", "", heading.lower()).replace(" ", "-")


def section(text: str, name: str) -> str:
    m = re.search(rf"^## {re.escape(name)}\n(.*?)(?=^## |\Z)", text, flags=re.S | re.M)
    return m.group(1).strip() if m else ""


def main() -> None:
    summaries = sorted(p for p in C.RESULTS.glob("[0-9][0-9]_*/SUMMARY.md"))
    if not summaries:
        sys.exit(f"ERROR: no results/*/SUMMARY.md under {C.RESULTS} -- run the analyses first")
    lines = ["# Step 3 report: downstream analyses", "",
             f"_Built {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} from "
             f"{len(summaries)} analyses over `{C.EXPORTS.name}/`._", "",
             PREAMBLE, "## Contents", ""]
    parsed = []
    for p in summaries:
        text = p.read_text()
        title = text.splitlines()[0].lstrip("# ").strip()
        parsed.append((p, title, text))
        heading = f"{p.parent.name[:2]}. {title}"
        lines.append(f"- [{heading}](#{github_anchor(heading)})")
    lines.append("")
    for p, title, text in parsed:
        rel = p.relative_to(C.RESULTS.parent)
        lines += [f"## {p.parent.name[:2]}. {title}", "",
                  f"**Question.** {section(text, 'Question')}", "",
                  section(text, "Key findings"), "",
                  f"Method, limitations and output files: [{rel}]({rel})", ""]
    REPORT.write_text("\n".join(lines) + "\n")
    C.log(f"wrote {REPORT} ({len(summaries)} analyses)")


if __name__ == "__main__":
    main()

"""Shared helpers for gene-programme fits (Spectra, script 25; cNMF, script 26)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from lib.io_utils import Manifest


def write_group_scores(man: Manifest, obs: pd.DataFrame, scores: np.ndarray, factors: list[str],
                       subdir: str) -> None:
    """Mean cell score per factor in each class, class x age, class x region x age and donor."""
    S = pd.DataFrame(scores, index=obs.index, columns=factors)
    cls = obs.cell_class.astype(str)
    groups = {"cell_class": cls}
    if "age_pcw" in obs:
        groups["cell_class_x_age"] = cls + " | " + obs.age_pcw.astype(str)
        if "region" in obs:
            groups["cell_class_x_region_x_age"] = cls + " | " + obs.region.astype(str) + " | " + obs.age_pcw.astype(str)
    if "donor" in obs:
        groups["donor"] = obs.donor.astype(str)
    long = []
    for gname, lab in groups.items():
        m = S.groupby(lab.to_numpy()).mean()
        n = lab.value_counts()
        t = m.reset_index(names="group").melt(id_vars="group", var_name="factor", value_name="mean_score")
        t.insert(0, "grouping", gname)
        t["n_cells"] = t.group.map(n).astype(int)
        long.append(t)
    man.write(pd.concat(long, ignore_index=True), "factor_scores_by_group",
              "Mean cell score per factor per class, class x age, class x region x age and donor", subdir=subdir)

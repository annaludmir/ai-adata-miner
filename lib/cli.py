"""Shared argument parsing and bootstrap for every extraction script."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def bootstrap() -> None:
    """Make `import config` / `import lib...` work when run as a script."""
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))


bootstrap()

import config  # noqa: E402
from lib.io_utils import log, require_file  # noqa: E402


def build_parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--dataset", default="all",
                   choices=[*config.DATASETS, "all"],
                   help="which dataset to process")
    p.add_argument("--chunk-size", type=int, default=config.CHUNK_SIZE,
                   help="cells per chunk when streaming X")
    p.add_argument("--limit-cells", type=int, default=None,
                   help="process only the first N cells -- use for a fast smoke "
                        "test before committing to a full run")
    p.add_argument("--overwrite", action="store_true",
                   help="rewrite CSVs that already exist (default: rewrite anyway; "
                        "flag kept for symmetry with --skip-existing)")
    p.add_argument("--skip-existing", action="store_true",
                   help="skip a dataset whose outputs for this script already exist")
    p.add_argument("--chemistry", default="each",
                   choices=["each", "all", "v2", "v3"],
                   help="10x chemistry stratification. 'each' (the default) runs v2 "
                        "and v3 separately into <dataset>__v2/ and <dataset>__v3/; "
                        "'all' pools them into <dataset>/. Age and chemistry are "
                        "heavily confounded in these data, so pooled results mix a "
                        "developmental effect with the v2->v3 switch")
    return p


def chemistry_values(args) -> list[str | None]:
    """Chemistries to iterate. None means 'pool everything'."""
    choice = getattr(args, "chemistry", "each")
    if choice == "all":
        return [None]
    if choice == "each":
        return list(config.CHEMISTRIES)
    return [choice]


def dataset_variants(args):
    """Yield (dataset_key, chemistry, namespace) for every unit of work.

    The namespace is what the CSVs are filed under: 'cortex' when pooled,
    'cortex__v2' when stratified. Keeping the chemistry in the folder name --
    rather than only in a column -- means a downstream analysis cannot
    accidentally pool the two by globbing.
    """
    for key in selected_datasets(args):
        for chem in chemistry_values(args):
            yield key, chem, config.namespace(key, chem)


def selected_datasets(args) -> list[str]:
    return list(config.DATASETS) if args.dataset == "all" else [args.dataset]


def resolve_h5ad(key: str) -> Path:
    return require_file(config.dataset(key)["h5ad"])


def banner(script: str, key: str, chem: str | None = None) -> None:
    meta = config.dataset(key)
    suffix = f" | chemistry={chem}" if chem else " | chemistry=pooled"
    log(f"=== {script} | dataset={key}{suffix} ({meta['label']})")

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
    return p


def selected_datasets(args) -> list[str]:
    return list(config.DATASETS) if args.dataset == "all" else [args.dataset]


def resolve_h5ad(key: str) -> Path:
    return require_file(config.dataset(key)["h5ad"])


def banner(script: str, key: str) -> None:
    meta = config.dataset(key)
    log(f"=== {script} | dataset={key} ({meta['label']})")

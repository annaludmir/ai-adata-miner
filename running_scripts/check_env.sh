#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=FAIL
#SBATCH --job-name=aim_checkenv
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:20:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# Does AIM_ENV actually satisfy this pipeline? Run this before anything else.
#
# Light enough for the login node, and submittable if you would rather not run
# anything there:
#
#   bash check_env.sh            # login node, ~1 minute
#   sbatch check_env.sh          # as a job
#   AIM_ENV=/path/to/other/env bash check_env.sh
#
# It checks the five packages the pipeline needs, resolves anndata's read_elem
# (which has moved twice between versions), imports the repo's own modules, and
# then does the real test: opens both h5ad files and reads their shape and a
# slice of .obs. An env that passes this can run the pipeline.

set -uo pipefail

AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
AIM_ENV="${AIM_ENV:-/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new}"
export AI_ADATA_DATA_ROOT="${AI_ADATA_DATA_ROOT:-/miridan-data/annaludmir/ndd_gene_modules/data}"

module load mamba/mamba-1.5.8 2>/dev/null || true
mamba activate "$AIM_ENV" 2>/dev/null || true
cd "$AIM_ROOT" || { echo "FAIL: AIM_ROOT=$AIM_ROOT does not exist"; exit 1; }

echo "env : $AIM_ENV"
echo "repo: $AIM_ROOT"
echo "data: $AI_ADATA_DATA_ROOT"
echo

mamba run -p "$AIM_ENV" python -u - <<'PY'
import importlib, sys, traceback
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
fail = []

print(f"python {sys.version.split()[0]}  ({sys.executable})")
print()
print("--- required packages ---")
# The pipeline deliberately does NOT need scanpy; these five are the whole list.
REQUIRED = {"numpy": "1.24", "pandas": "2.0", "scipy": "1.10",
            "h5py": "3.8", "anndata": "0.9"}
for mod, minimum in REQUIRED.items():
    try:
        m = importlib.import_module(mod)
        try:
            from importlib.metadata import version as _v
            ver = _v(mod)
        except Exception:
            ver = getattr(m, "__version__", "?")
        def tup(v):
            out = []
            for p in str(v).split("."):
                digits = "".join(c for c in p if c.isdigit())
                if not digits:
                    break
                out.append(int(digits))
            return tuple(out)
        ok = tup(ver) >= tup(minimum) if ver != "?" else True
        print(f"  {'ok  ' if ok else 'OLD '} {mod:<10} {ver:<12} (need >= {minimum})")
        if not ok:
            fail.append(f"{mod} {ver} < {minimum}")
    except Exception as e:
        print(f"  MISSING {mod:<10} ({e})")
        fail.append(f"{mod} missing")

print()
print("--- anndata read_elem (moved twice between versions) ---")
try:
    from anndata.io import read_elem                      # >= 0.11
    print("  ok   anndata.io.read_elem")
except ImportError:
    try:
        from anndata.experimental import read_elem        # 0.8 - 0.10
        print("  ok   anndata.experimental.read_elem")
    except ImportError:
        try:
            from anndata._io.specs import read_elem       # fallback
            print("  ok   anndata._io.specs.read_elem")
        except ImportError as e:
            print(f"  FAIL no read_elem anywhere ({e})")
            fail.append("anndata read_elem unavailable")

print()
print("--- repo modules ---")
for mod in ("config", "lib.io_utils", "lib.aggregate", "lib.stats_utils",
            "lib.bulk_stats", "lib.panels", "lib.cli"):
    try:
        importlib.import_module(mod)
        print(f"  ok   {mod}")
    except Exception as e:
        print(f"  FAIL {mod}: {e}")
        fail.append(f"import {mod}")

print()
print("--- the real test: can we open the data? ---")
try:
    import config
    from lib.io_utils import XReader, read_obs, read_var
    for key, meta in config.DATASETS.items():
        path = Path(meta["h5ad"])
        if not path.exists():
            print(f"  MISSING {key}: {path}")
            fail.append(f"{key} h5ad not found")
            continue
        with XReader(path) as xr:
            shape, enc = xr.shape, xr.encoding
        var = read_var(path)
        obs = read_obs(path)
        agree = (shape[0] == meta["expected_n_obs"] and shape[1] == meta["expected_n_vars"])
        print(f"  ok   {key:<10} X {shape[0]:>9,} x {shape[1]:>6,} ({enc})  "
              f"obs {obs.shape[1]} cols, var {var.shape[1]} cols"
              f"{'' if agree else '   <-- SHAPE DIFFERS FROM schemas/'}")
        if not agree:
            print(f"       schema expects {meta['expected_n_obs']:,} x "
                  f"{meta['expected_n_vars']:,}; update config.py or re-extract")
except Exception:
    traceback.print_exc()
    fail.append("reading h5ad")

print()
print("=" * 62)
if fail:
    print("NOT READY -- fix these first:")
    for f in fail:
        print(f"  - {f}")
    print("\nIf packages are missing, see running_scripts/README.md -> 'Creating the environment'")
    sys.exit(1)
print("READY -- this environment can run the pipeline.")
PY
rc=$?
echo "Python exit code: $rc"
exit $rc

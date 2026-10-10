#!/bin/bash
# Download the external annotation files step-3 part C reads (analyses 29-34).
#
# Light downloads (~0.3 GB), so it runs on the LOGIN node -- no sbatch needed:
#
#   ./fetch_annotations.sh                       # into /miridan-data/annaludmir/aim_annotations
#   AIM_ANNOTATIONS=/elsewhere ./fetch_annotations.sh
#   FORCE=true ./fetch_annotations.sh            # re-download files that already exist
#
# Every file records where and when it came from in MANIFEST.tsv. Nothing here
# is committed to git (sizes and licences: Enrichr libraries are free for
# academic use; gnomAD, HGNC, HPO and OmniPath are open; BrainSpan is free with
# attribution). OMIM is licensed and not fetched -- HPO's gene-to-disease table
# carries OMIM and Orphanet disease ids openly.

set -uo pipefail

AIM_ROOT="${AIM_ROOT:-/miridan-data/annaludmir/ai-adata-miner}"
if [[ -n "${AIM_ANNOTATIONS:-}" ]]; then
  OUT="$AIM_ANNOTATIONS"
elif [[ -d /miridan-data/annaludmir ]]; then
  OUT=/miridan-data/annaludmir/aim_annotations
else
  OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/annotations"
fi
FORCE="${FORCE:-false}"
mkdir -p "$OUT"/{gmt,tf,constraint,disease,lr,hgnc,brainspan}
MANIFEST="$OUT/MANIFEST.tsv"
[[ -f "$MANIFEST" ]] || printf 'file\turl\tdownloaded_utc\tbytes\n' > "$MANIFEST"
FAILED=()

ENRICHR="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
OMNIPATH="https://omnipathdb.org"

fetch() {   # fetch <relative path> <url> [min bytes]
  local rel="$1" url="$2" min="${3:-1000}" dest="$OUT/$1"
  if [[ -s "$dest" && "$FORCE" != "true" ]]; then
    echo "  have  $rel"
    return
  fi
  if curl -sSL --fail --retry 3 --retry-delay 5 --max-time 1800 -o "$dest.part" "$url"; then
    local bytes; bytes=$(wc -c < "$dest.part" | tr -d ' ')
    if (( bytes >= min )); then
      mv "$dest.part" "$dest"
      printf '%s\t%s\t%s\t%s\n' "$rel" "$url" "$(date -u '+%F %T')" "$bytes" >> "$MANIFEST"
      echo "  ok    $rel ($bytes bytes)"
      return
    fi
    echo "  SHORT $rel ($bytes bytes < $min)"
  else
    echo "  FAIL  $rel"
  fi
  rm -f "$dest.part"
  FAILED+=("$rel")
}

echo "annotations -> $OUT"
echo "C1 functional gene sets (Enrichr .gmt)"
for lib in GO_Biological_Process_2023 GO_Cellular_Component_2023 GO_Molecular_Function_2023 \
           Reactome_2022 KEGG_2021_Human; do
  fetch "gmt/${lib}.gmt" "${ENRICHR}${lib}" 10000
done

echo "C2 transcription factors (Lambert et al. 2018) and their targets (CollecTRI)"
fetch tf/TF_names_v_1.01.txt "http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt" 5000
fetch tf/collectri.tsv "${OMNIPATH}/interactions?datasets=collectri&genesymbols=yes&format=tsv" 100000

echo "C3 constraint (gnomAD v4.1) and disease genes (HPO)"
fetch constraint/gnomad.v4.1.constraint_metrics.tsv \
  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" 1000000
fetch disease/genes_to_disease.txt "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt" 100000

echo "C4 ligand-receptor (OmniPath interactions + intercell roles)"
fetch lr/omnipath_interactions.tsv "${OMNIPATH}/interactions?datasets=omnipath,ligrecextra&genesymbols=yes&format=tsv" 100000
fetch lr/intercell_ligand_receptor.tsv "${OMNIPATH}/intercell?categories=ligand,receptor&format=tsv" 100000

echo "C5 HGNC symbols (previous / alias)"
fetch hgnc/hgnc_complete_set.txt "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" 1000000

echo "C6 BrainSpan developmental transcriptome (RNA-seq, Gencode v10 genes)"
fetch brainspan/genes_matrix_csv.zip "http://www.brainspan.org/api/v2/well_known_file_download/267666525" 1000000
if [[ -s "$OUT/brainspan/genes_matrix_csv.zip" && ! -s "$OUT/brainspan/expression_matrix.csv" ]]; then
  (cd "$OUT/brainspan" && unzip -o -q genes_matrix_csv.zip) && echo "  unzipped BrainSpan" \
    || { echo "  FAIL  unzip BrainSpan"; FAILED+=("brainspan unzip"); }
fi

echo "C7 STRING v12.0 human protein links (per-channel scores) and protein names (Spectra STRING prior)"
STRING="https://stringdb-downloads.org/download"
fetch string/9606.protein.links.detailed.v12.0.txt.gz \
  "${STRING}/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz" 100000000
fetch string/9606.protein.info.v12.0.txt.gz "${STRING}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" 1000000

echo
if (( ${#FAILED[@]} )); then
  echo "FAILED: ${FAILED[*]} -- rerun later; analyses that need them say so and skip"
  exit 1
fi
echo "all annotation files present in $OUT"

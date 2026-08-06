#!/usr/bin/env bash

set -euo pipefail

OUT_DIR="mges_scanner/db"
TMP_DIR=$(mktemp -d -t mge_hmms-XXXXXX)
FINAL_DB="${OUT_DIR}/mge_signatures.hmm"
PFAM_FTP_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

# Array of Pfam Accessions
PFAM_ACCS=(
  "PF00872"
  "PF01609"
  "PF01526"
  "PF13358"
  "PF02371"
  "PF01797"
  "PF00589"
  "PF01051"
  "PF01446"
  "PF08704"
  "PF01719"
)

cleanup() {
  echo "[INFO] Cleaning up temporary files..."
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

echo "[1/4] Checking required dependencies..."
for cmd in curl gunzip hmmpress hmmfetch grep; do
  if ! command -v "$cmd" &> /dev/null; then
    echo "[ERROR] Required tool '$cmd' is not installed." >&2
    exit 1
  fi
done

mkdir -p "${OUT_DIR}"

echo "[2/4] Downloading Pfam-A database archive..."
curl -sSfL "${PFAM_FTP_URL}" -o "${TMP_DIR}/Pfam-A.hmm.gz"

echo "[3/4] Extracting target MGE profiles..."
gunzip "${TMP_DIR}/Pfam-A.hmm.gz"

# Build regex string to match accession headers with any version decimal (e.g., ACC   PF00872.24)
ACC_REGEX=$(printf "|ACC   %s\\." "${PFAM_ACCS[@]}")
ACC_REGEX="${ACC_REGEX#|}"

# Parse full versioned accessions from the uncompressed Pfam-A.hmm
VERSIONED_ACCS_FILE="${TMP_DIR}/versioned_accs.txt"
grep -E "${ACC_REGEX}" "${TMP_DIR}/Pfam-A.hmm" | awk '{print $2}' > "${VERSIONED_ACCS_FILE}"

# Index and extract exact versioned entries
hmmpress "${TMP_DIR}/Pfam-A.hmm"
hmmfetch -f "${TMP_DIR}/Pfam-A.hmm" "${VERSIONED_ACCS_FILE}" > "${FINAL_DB}"

echo "[4/4] Pressing final MGE database..."
rm -f "${FINAL_DB}".h3*
hmmpress "${FINAL_DB}"

echo "-------------------------------------------------------------------"
echo "✅ Success! MGE HMM database built successfully at:"
echo "   ${FINAL_DB}"
echo "-------------------------------------------------------------------"

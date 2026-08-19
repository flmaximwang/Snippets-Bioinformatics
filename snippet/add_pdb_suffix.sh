#!/bin/bash
# Append a .pdb suffix to every file in the given directory (or directories).
# Idempotent: files whose names already end in .pdb are left untouched.
#
# Use case: MASTER (create_master_db.sh) decompress/filter steps name the output
# files pdbXXXX (no extension), e.g. /mnt/data/public/MASTER/{pdb,pdb_filtered}.
# This adds the .pdb extension in place so downstream tools that require it work.
#
# Usage:
#   bash add_pdb_suffix.sh <dir> [<dir> ...]
#
# Example:
#   bash add_pdb_suffix.sh /mnt/data/public/MASTER/pdb /mnt/data/public/MASTER/pdb_filtered

set -euo pipefail

[ "$#" -ge 1 ] || { echo "usage: $0 <dir> [<dir> ...]" >&2; exit 1; }

total_renamed=0
total_skipped=0

for dir in "$@"; do
    [ -d "$dir" ] || { echo "SKIP (not a dir): $dir" >&2; continue; }
    renamed=0; skipped=0
    while IFS= read -r -d '' f; do
        if [[ "$f" == *.pdb ]]; then
            skipped=$((skipped + 1))
            continue
        fi
        mv -- "$f" "$f.pdb"
        renamed=$((renamed + 1))
    done < <(find "$dir" -maxdepth 1 -type f -print0)

    echo "$dir: renamed $renamed, already-.pdb skipped $skipped"
    total_renamed=$((total_renamed + renamed))
    total_skipped=$((total_skipped + skipped))
done

echo "DONE: total renamed $total_renamed, skipped $total_skipped"

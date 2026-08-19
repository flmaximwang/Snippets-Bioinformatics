#!/bin/bash
# Build a MASTER structural-search database from the local RCSB PDB archive.
# Input:  <rcsb>/{aa}/pdbXXXX.ent.gz   (gzipped PDB, divided layout)
# Output: <out>/pdb/     decompressed plain-text PDB
#         <out>/pds/     createPDS target PDS files (the actual database)
#         <out>/pdb_list.txt     input list for createPDS (absolute paths)
#         <out>/target_list.txt  --targetList for master (absolute paths)
#
# createPDS only reads *uncompressed* PDB text, so gz files must be decompressed
# first. Steps are idempotent: rerunning skips already-done work.
#
# Usage:
#   bash create_master_db.sh --rcsb <SRC> --out <OUT> \
#                            [--no-cleanPDB] [--dCut X --dStep X --phiStep X --psiStep X] \
#                            {decompress|build|verify|all}
#
# Required (no default): --rcsb, --out, --createpds
# Expert (default):       --dCut 25.0 --dStep 5.0 --phiStep 10.0 --psiStep 10.0
# cleanPDB is ON by default; pass --no-cleanPDB to disable.
#
# Deps: createPDS + master (MASTER v1.6, grigoryanlab.org).

set -euo pipefail

# ---------------- expert params (defaults) ----------------
DCUT=25.0; DSTEP=5.0; PHISTEP=10.0; PSISTEP=10.0   # must match between DB and query
CLEAN_PDB="--cleanPDB"                              # ON by default; --no-cleanPDB clears it

# createPDS is resolved via PATH (a MASTER install provides it, e.g. /usr/local/bin).

# ---------------- required (no default) ----------------
RCSB_DIR=""
MASTERDIR=""

ACTION=""

usage() {
    cat >&2 <<EOF
Usage: $0 --rcsb <SRC> --out <OUT> [options] {decompress|build|verify|all}

Required:
  --rcsb <SRC>       source gz PDB archive (divided layout, e.g. .../RCSB)
  --out <OUT>        output root for this database

Expert (defaults):
  --dCut <X>         CA-CA distance table upper limit   (default 25.0)
  --dStep <X>        CA-CA distance bin size            (default 5.0)
  --phiStep <X>      phi angle bin size                 (default 10.0)
  --psiStep <X>      psi angle bin size                 (default 10.0)

Other:
  --no-cleanPDB      disable --cleanPDB (default: ON)

Action:
  decompress         gunzip .ent.gz -> out/pdb, write pdb_list.txt
  build              createPDS pdb_list.txt -> out/pds, write target_list.txt
  verify             compare PDS count vs target_list count
  all                decompress + build + verify

Notes:
  createPDS must be on PATH (MASTER install, e.g. /usr/local/bin).
EOF
}

# ---------------- parse named options ----------------
while [ $# -gt 0 ]; do
    case "$1" in
        --rcsb)        RCSB_DIR="$2"; shift 2 ;;
        --out)         MASTERDIR="$2"; shift 2 ;;
        --no-cleanPDB) CLEAN_PDB=""; shift ;;
        --dCut)        DCUT="$2"; shift 2 ;;
        --dStep)       DSTEP="$2"; shift 2 ;;
        --phiStep)     PHISTEP="$2"; shift 2 ;;
        --psiStep)     PSISTEP="$2"; shift 2 ;;
        decompress|build|verify|all) ACTION="$1"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------------- validate required ----------------
err=0
[ -n "$RCSB_DIR" ]   || { echo "missing required: --rcsb" >&2; err=1; }
[ -n "$MASTERDIR" ]  || { echo "missing required: --out" >&2; err=1; }
[ -n "$ACTION" ]     || { echo "missing action: decompress|build|verify|all" >&2; err=1; }
if [ "$err" -ne 0 ]; then usage; exit 1; fi

# verify createPDS is reachable via PATH
if ! command -v createPDS >/dev/null 2>&1; then
    echo "createPDS not found on PATH (MASTER v1.6 provides it, e.g. /usr/local/bin)" >&2
    exit 1
fi
CREATEPDS="$(command -v createPDS)"

PDB_DIR="$MASTERDIR/pdb"
PDS_DIR="$MASTERDIR/pds"
PDB_LIST="$MASTERDIR/pdb_list.txt"
TARGET_LIST="$MASTERDIR/target_list.txt"
LOG="$MASTERDIR/create_master_db.log"

mkdir -p "$MASTERDIR" "$PDB_DIR" "$PDS_DIR"

# Resolve MASTERDIR to absolute so relative --out canonicalizes here. All of
# PDB_DIR/PDS_DIR/PDB_LIST/TARGET_LIST and the paths inside pdb_list.txt are
# derived from MASTERDIR and thus become absolute too. This keeps them valid
# after build() cds into PDS_DIR (makePDS runs there).
MASTERDIR="$(cd "$MASTERDIR" && pwd)"
PDB_DIR="$MASTERDIR/pdb"
PDS_DIR="$MASTERDIR/pds"
PDB_LIST="$MASTERDIR/pdb_list.txt"
TARGET_LIST="$MASTERDIR/target_list.txt"
LOG="$MASTERDIR/create_master_db.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# ---------------- action: decompress ----------------
decompress() {
    local NJOBS; NJOBS=${NJOBS:-$(nproc)}     # 这台 nproc=32,SSD 可直接跑满
    # bash -c 会新开子 shell,不继承普通 shell 变量,须 export 才能读到路径
    export PDB_DIR RCSB_DIR
    log "Decompressing $RCSB_DIR/*/*.ent.gz -> $PDB_DIR/  ($NJOBS 并行,流式)"
    mkdir -p "$PDB_DIR"

    # find 流式输出:边遍历边喂给 xargs,不等 glob 全部展开
    # 已存在则跳过:worker 内用 [ -f ] 判断,避免重复解压
    find "$RCSB_DIR" -name '*.ent.gz' -print0 | \
    xargs -0 -P "$NJOBS" -n1 bash -c '
        f="$1"
        id=$(basename "$f" .ent.gz)
        [ -f "$PDB_DIR/$id" ] && exit 0          # 目标已存在 → 跳过
        gunzip -c "$f" > "$PDB_DIR/$id"
    ' _

    local n; n=$(find "$RCSB_DIR" -name '*.ent.gz' | wc -l)
    log "Processed $n gz entries (含已存在跳过)."

    log "Writing $PDB_LIST"
    find "$PDB_DIR" -maxdepth 1 -type f | sort > "$PDB_LIST"
    log "  $(wc -l < "$PDB_LIST") uncompressed PDB files listed."
}

# ---------------- action: build (convert to PDS) ----------------
build() {
    if [ ! -s "$PDB_LIST" ]; then
        log "ERROR: $PDB_LIST empty/missing. Run: $0 --rcsb ... --out ... --createpds ... decompress first." >&2
        exit 1
    fi
    log "Converting PDB -> PDS (target) into $PDS_DIR/ (cleanPDB='${CLEAN_PDB:-off}')"
    # createPDS writes each .pds in cwd, naming from the PDB base name + .pds.
    # Run from $PDS_DIR so outputs land here.
    ( cd "$PDS_DIR" && \
      "$CREATEPDS" --type target --pdbList "$PDB_LIST" \
                   --dCut "$DCUT" --dStep "$DSTEP" \
                   --phiStep "$PHISTEP" --psiStep "$PSISTEP" \
                   $CLEAN_PDB )
    log "Writing $TARGET_LIST"
    find "$PDS_DIR" -maxdepth 1 -type f -name '*.pds' | sort > "$TARGET_LIST"
    log "  $(wc -l < "$TARGET_LIST") PDS target files."
}

# ---------------- action: verify counts ----------------
verify() {
    local npds nlist
    npds=$(find "$PDS_DIR" -maxdepth 1 -name '*.pds' | wc -l)
    nlist=$(wc -l < "$TARGET_LIST" 2>/dev/null || echo 0)
    echo "PDS files on disk   : $npds"
    echo "lines in target_list: $nlist"
    [ "$npds" -eq "$nlist" ] && [ "$npds" -gt 0 ] && echo "OK: consistent" \
        || echo "MISMATCH: re-run build"
}

# ---------------- dispatch ----------------
case "$ACTION" in
    decompress) decompress ;;
    build)      build ;;
    verify)     verify ;;
    all)        decompress; build; verify ;;
    *)          usage; exit 1 ;;   # unreachable (validated above)
esac

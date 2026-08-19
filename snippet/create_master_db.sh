#!/bin/bash
# Build a MASTER structural-search database from the local RCSB PDB archive.
# Input:  <rcsb>/{aa}/pdbXXXX.ent.gz   (gzipped PDB, divided layout)
# Output: <out>/pdb/            decompressed plain-text PDB (*.pdb)
#         <out>/pdb_filtered/   PDBs containing >=1 MASTER-supported residue
#         <out>/pds/            createPDS target PDS files (the actual database)
#         <out>/pdb_filtered_list.txt  createPDS --pdbList (absolute paths)
#         <out>/target_list.txt createPDS --pdsList & master --targetList (absolute paths)
#
# createPDS only reads *uncompressed* PDB text, so gz files must be decompressed
# first. Steps are idempotent: rerunning skips already-done work.
#
# A note on paths: --rcsb and --out are canonicalized to ABSOLUTE paths up front,
# so every path derived from them (PDB_DIR / PDB_FILTERED / PDS_DIR /
# PDB_FILTERED_LIST / TARGET_LIST / the *.pds outputs) and every path written
# INSIDE pdb_filtered_list.txt / target_list.txt
# is absolute. This keeps everything valid regardless of cwd changes — build()
# cds into PDS_DIR so createPDS drops .pds there, and because all paths are
# absolute the cd is harmless.
#
# Usage:
#   bash create_master_db.sh --rcsb <SRC> --out <OUT> \
#                            [--no-cleanPDB] [--dCut X --dStep X --phiStep X --psiStep X] \
#                            {decompress|filter|gen_list|build|verify|all}
#
# Required (no default): --rcsb, --out
# Expert (default):       --dCut 25.0 --dStep 5.0 --phiStep 10.0 --psiStep 10.0
# cleanPDB is ON by default; pass --no-cleanPDB to disable.
#
# Deps: createPDS + master (MASTER v1.6, grigoryanlab.org).

set -euo pipefail

# ---------------- expert params (defaults) ----------------
DCUT=25.0; DSTEP=5.0; PHISTEP=10.0; PSISTEP=10.0   # must match between DB and query
CLEAN_PDB="--cleanPDB"                              # ON by default; --no-cleanPDB clears it

# Residue names createPDS accepts with --cleanPDB (CreateOptions::setLegalAA):
# 20 natural + HIS protonation variants + MSE + modified CSO/HIP/PTR/SEC/SEP/TPO.
MASTER_SUPPORTED="ALA ARG ASN ASP CYS GLN GLU GLY HIS HSC HSD HSE HSP ILE LEU LYS MET MSE PHE PRO SER THR TRP TYR VAL CSO HIP PTR SEC SEP TPO"

# createPDS is resolved via PATH (a MASTER install provides it, e.g. /usr/local/bin).

# ---------------- required (no default) ----------------
RCSB_DIR=""
MASTERDIR=""

ACTION=""

usage() {
    cat >&2 <<EOF
Usage: $0 --rcsb <SRC> --out <OUT> [options] {decompress|filter|gen_list|build|verify|all}

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
  decompress         gunzip .ent.gz -> out/pdb/*.pdb
  filter             keep PDBs w/ any MASTER-supported residue -> out/pdb_filtered
  gen_list           write pdb_filtered_list.txt + target_list.txt (before build)
  build              createPDS (via both lists) -> out/pds/*.pds
  verify             compare PDS count vs target_list count
  all                decompress + filter + gen_list + build + verify

Notes:
  createPDS must be on PATH (MASTER install, e.g. /usr/local/bin).
  All paths are canonicalized to absolute before use.
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
        decompress|filter|gen_list|build|verify|all) ACTION="$1"; shift ;;
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

# ---------------- canonicalize to absolute paths ----------------
# Both inputs become absolute, so every derived path — PDB_DIR, PDB_FILTERED,
# PDB_FILTERED_LIST, TARGET_LIST, the *.pds outputs, and the lines written
# inside pdb_filtered_list.txt / target_list.txt — is absolute. Relative --rcsb/--out work regardless of cwd.
RCSB_DIR="$(cd "$RCSB_DIR" && pwd)"
MASTERDIR="$(cd "$MASTERDIR" && pwd)"

PDB_DIR="$MASTERDIR/pdb"
PDB_FILTERED="$MASTERDIR/pdb_filtered"
PDB_FILTERED_LIST="$MASTERDIR/pdb_filtered_list.txt"
PDS_DIR="$MASTERDIR/pds"
TARGET_LIST="$MASTERDIR/target_list.txt"
LOG="$MASTERDIR/create_master_db.log"

mkdir -p "$MASTERDIR" "$PDB_DIR" "$PDB_FILTERED" "$PDS_DIR"

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
        out="$PDB_DIR/$(basename "$f" .ent.gz).pdb"
        [ -f "$out" ] && exit 0                  # 目标已存在 → 跳过
        gunzip -c "$f" > "$out"
    ' _

    local n; n=$(find "$RCSB_DIR" -name '*.ent.gz' | wc -l)
    log "Processed $n gz entries (含已存在跳过)."
}

# ---------------- action: filter (keep PDBs with >=1 MASTER-supported residue) ----------------
filter() {
    if [ ! -d "$PDB_DIR" ]; then
        log "ERROR: $PDB_DIR missing. Run: $0 --rcsb ... --out ... decompress first." >&2
        exit 1
    fi
    local NJOBS; NJOBS=${NJOBS:-$(nproc)}
    # bash -c 会新开子 shell,不继承普通 shell 变量,须 export 才能读到路径与残基表
    export PDB_FILTERED MASTER_SUPPORTED
    log "Filtering $PDB_DIR -> $PDB_FILTERED (keep PDBs containing any MASTER-supported residue)"
    mkdir -p "$PDB_FILTERED"

    find "$PDB_DIR" -maxdepth 1 -type f -print0 | \
    xargs -0 -P "$NJOBS" -n1 bash -c '
        f="$1"
        id=$(basename "$f")
        [ -f "$PDB_FILTERED/$id" ] && exit 0          # 已拷过 → 跳过(幂等)
        # 保留 iff 任一 ATOM/HETATM 残基名(第18-20列)落在 MASTER 支持表内。
        # awk 的退出码只决定是否 cp,worker 自身必须始终以 0 退出,
        # 否则被过滤文件会让 xargs 返回 123 并在 set -e 下中断脚本。
        if awk '\''
            BEGIN { n=split(ENVIRON["MASTER_SUPPORTED"],a," "); for(i=1;i<=n;i++) allowed[a[i]]=1 }
            /^(ATOM  |HETATM)/ { if (substr($0,18,3) in allowed) { f=1; exit } }
            END { exit (f ? 0 : 1) }
        '\'' "$f"; then
            cp "$f" "$PDB_FILTERED/$id"
        fi
        exit 0
    ' _

    local n; n=$(find "$PDB_FILTERED" -maxdepth 1 -type f | wc -l)
    log "  $n PDBs kept in $PDB_FILTERED."
}

# ---------------- action: gen_list (write the 2 lists createPDS needs) ----------------
# Runs BEFORE build. Writes:
#   pdb_filtered_list.txt -> createPDS --pdbList (input PDBs, absolute paths)
#   target_list.txt       -> createPDS --pdsList (output .pds paths into out/pds,
#                           one per input), and doubles as master's --targetList.
gen_list() {
    if [ ! -d "$PDB_FILTERED" ]; then
        log "ERROR: $PDB_FILTERED missing. Run: $0 ... filter first." >&2
        exit 1
    fi
    mkdir -p "$PDS_DIR"
    # PDB_FILTERED is absolute -> lines written here are absolute.
    log "Writing $PDB_FILTERED_LIST (createPDS --pdbList inputs)"
    find "$PDB_FILTERED" -maxdepth 1 -type f -name '*.pdb' | sort > "$PDB_FILTERED_LIST"
    # Derive each output .pds path: same basename, in PDS_DIR, .pdb -> .pds.
    log "Writing $TARGET_LIST (createPDS --pdsList outputs / master --targetList)"
    awk -v fd="$PDB_FILTERED/" -v pds="$PDS_DIR" '
        { base=$0; sub("^" fd, "", base); sub(/\.pdb$/, "", base); print pds "/" base ".pds" }
    ' "$PDB_FILTERED_LIST" > "$TARGET_LIST"
    log "  $(wc -l < "$PDB_FILTERED_LIST") PDB inputs, $(wc -l < "$TARGET_LIST") PDS outputs."
}

# ---------------- action: build (convert to PDS) ----------------
build() {
    if [ ! -s "$PDB_FILTERED_LIST" ] || [ ! -s "$TARGET_LIST" ]; then
        log "ERROR: $PDB_FILTERED_LIST or $TARGET_LIST missing/empty. Run: $0 ... gen_list first." >&2
        exit 1
    fi
    log "Converting PDB -> PDS (target) into $PDS_DIR/ (cleanPDB='${CLEAN_PDB:-off}')"
    # Output .pds paths come from --pdsList (target_list.txt), so createPDS writes
    # them into $PDS_DIR regardless of cwd.
    ( cd "$PDS_DIR" && \
      "$CREATEPDS" --type target \
                   --pdbList "$PDB_FILTERED_LIST" --pdsList "$TARGET_LIST" \
                   --dCut "$DCUT" --dStep "$DSTEP" \
                   --phiStep "$PHISTEP" --psiStep "$PSISTEP" \
                   $CLEAN_PDB )
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
    filter)     filter ;;
    gen_list)   gen_list ;;
    build)      build ;;
    verify)     verify ;;
    all)        decompress; filter; gen_list; build; verify ;;
    *)          usage; exit 1 ;;   # unreachable (validated above)
esac

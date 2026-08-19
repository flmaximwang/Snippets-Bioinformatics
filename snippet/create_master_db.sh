#!/bin/bash
# Build a MASTER structural-search database from the local RCSB PDB archive.
# Input:  <rcsb>/{aa}/pdbXXXX.ent.gz   (gzipped PDB, divided layout)
# Output: <out>/pdb/            decompressed plain-text PDB (*.pdb)
#         <out>/pdb_filtered/   PDBs containing >=1 MASTER-supported residue
#         <out>/pds/            createPDS target PDS files (the actual database)
#         <out>/pdb_filtered_list.txt  createPDS --pdbList (incremental: only PDBs missing a .pds)
#         <out>/build_pds_list.txt     createPDS --pdsList (incremental outputs, for build)
#         <out>/target_list.txt        full .pds list for master --targetList (via gen_target_list)
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
#                            {decompress|filter|gen_build_list|build|gen_target_list|verify|all}
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
Usage: $0 --rcsb <SRC> --out <OUT> [options] {decompress|filter|gen_build_list|build|gen_target_list|verify|all}

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
  gen_build_list     write incremental lists (PDBs missing a .pds) -> build uses them
  build              createPDS (via incremental lists) -> out/pds/*.pds
  gen_target_list    write full .pds list out/target_list.txt (for master --targetList)
  verify             compare PDS count vs target_list count
  all                decompress + filter + gen_build_list + build + gen_target_list + verify

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
        decompress|filter|gen_build_list|build|gen_target_list|verify|all) ACTION="$1"; shift ;;
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
BUILD_PDS_LIST="$MASTERDIR/build_pds_list.txt"
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

# ---------------- action: gen_build_list (incremental; run BEFORE build) ----------------
# Scans the .pds already in $PDS_DIR and writes lists containing ONLY the PDBs that
# do not yet have a corresponding .pds (i.e. the missing work). build() consumes these,
# so re-running build only processes what's still missing and never touches finished .pds.
#   pdb_filtered_list.txt -> createPDS --pdbList (PDBs missing a .pds, absolute paths)
#   build_pds_list.txt    -> createPDS --pdsList (their output .pds paths, absolute)
gen_build_list() {
    if [ ! -d "$PDB_FILTERED" ]; then
        log "ERROR: $PDB_FILTERED missing. Run: $0 ... filter first." >&2
        exit 1
    fi
    mkdir -p "$PDS_DIR"
    # PDB_FILTERED is absolute -> lines written here are absolute.
    log "Writing incremental build lists (PDBs missing a .pds) to $PDS_DIR/scan"
    # Collect basenames of .pds already on disk into a set.
    # Then, for every filtered PDB, keep it only if its .pds is NOT in the set.
    # awk substring removal keeps the awk one-pass; PDB_FILTERED basenames map
    # 1:1 to .pds basenames (both from the same <base>, dir differs).
    find "$PDB_FILTERED" -maxdepth 1 -type f -name '*.pdb' -printf '%f\n' | sort | \
    awk -v pf="$PDB_FILTERED" -v pds="$PDS_DIR" -v ofile="$PDB_FILTERED_LIST" -v ofile2="$BUILD_PDS_LIST" '
        function pdsName(b) { sub(/\.pdb$/, ".pds", b); return pds "/" b }
        BEGIN {
            # load existing .pds basenames
            while (("find " pds " -maxdepth 1 -type f -name \"*.pds\" -printf \"%f\\n\"" | getline b) > 0)
                have[b] = 1
            close("")
        }
        { if (!(($0 in have) || (pdsName($0) in have))) {
              print pf "/" $0 > ofile
              print pdsName($0) > ofile2
          } }
    '
    log "  $(wc -l < "$PDB_FILTERED_LIST" 2>/dev/null || echo 0) PDB inputs pending, $(wc -l < "$BUILD_PDS_LIST" 2>/dev/null || echo 0) PDS outputs pending."
}

# ---------------- action: gen_target_list (complete .pds list; run AFTER build) ----------------
# Writes the FULL list of every .pds on disk into $TARGET_LIST, for master --targetList.
# This is the database index for searching, independent of build bookkeeping.
gen_target_list() {
    mkdir -p "$PDS_DIR"
    log "Writing $TARGET_LIST (full .pds list for master --targetList)"
    find "$PDS_DIR" -maxdepth 1 -type f -name '*.pds' | sort > "$TARGET_LIST"
    log "  $(wc -l < "$TARGET_LIST") PDS entries."
}

# ---------------- action: build (convert to PDS) ----------------
build() {
    if [ ! -s "$PDB_FILTERED_LIST" ] || [ ! -s "$BUILD_PDS_LIST" ]; then
        log "ERROR: $PDB_FILTERED_LIST or $BUILD_PDS_LIST missing/empty. Run: $0 ... gen_build_list first." >&2
        exit 1
    fi
    log "Converting PDB -> PDS (target) into $PDS_DIR/ (cleanPDB='${CLEAN_PDB:-off}')"
    # Output .pds paths come from --pdsList (build_pds_list.txt), so createPDS writes
    # them into $PDS_DIR regardless of cwd.
    ( cd "$PDS_DIR" && \
      "$CREATEPDS" --type target \
                   --pdbList "$PDB_FILTERED_LIST" --pdsList "$BUILD_PDS_LIST" \
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
    gen_build_list) gen_build_list ;;
    build)      build ;;
    gen_target_list) gen_target_list ;;
    verify)     verify ;;
    all)        decompress; filter; gen_build_list; build; gen_target_list; verify ;;
    *)          usage; exit 1 ;;   # unreachable (validated above)
esac

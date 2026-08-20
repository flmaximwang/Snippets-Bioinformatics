#!/bin/bash
# Run a MASTER structural-motif SEARCH (the `master` binary only) against a built
# PDS database, with piece-wise parallelization and resume.
#
# This is the query side. It does NOT build a DB (that's create_master_db.sh) and
# does NOT convert PDB -> PDS (that's createPDS). It only takes a query .pds and a
# database --targetList, and runs the `master` search binary against it.
#
# Split-by-count = the foundation of resume:
#   --targetList is divided into consecutive pieces of --target-num-per-list
#   structures each (last piece may be smaller). Each piece is one `master
#   --targetList <piece>` job. A `.done.<i>` sentinel is written only after a piece
#   finishes successfully, so a re-run of the same command skips every already-done
#   piece and only searches the unfinished ones, then re-merges the piece outputs.
#
# --njobs is ONLY a concurrency cap: it controls how many pieces are searched at the
# same time (to utilize all cores), NOT the split size.
#
# A fingerprint of (targetList content + --target-num-per-list) is stored in
# out/.chunks/. If either changes, all `.done.*` sentinels are invalidated so a
# stale "done" can never wrongly skip work that now needs doing.
#
# Usage:
#   bash master_search.sh --query <q.pds> --targetList <list> --rmsdCut <X> \
#                         --out <dir> [options]
#
# Required:
#   --query <q.pds>          query structure already converted to PDS
#   --targetList <list>      file, one target .pds path per line (the database)
#   --rmsdCut <x>            RMSD cutoff (Angstrom) for a match
#   --out <dir>              output root (per-piece files + final merged results)
#
# Options / passthrough to `master`:
#   --target-num-per-list <N> structures per piece (default 10000); sets the resume granularity
#   --njobs <N>               how many pieces to search concurrently (default 8)
#   --bbRMSD                  full-backbone RMSD instead of CA-only
#   --topN <N>                keep best N matches (per piece, not global — see Notes)
#   --minN <N>                return at least N best matches (per piece)
#   --gapLen <s>              gap restraints, e.g. '1-10;0-3'
#   --outType <t>             match|full|wgap (region written by structOut)
#   --seqOut <file>           also write match sequences, merged to this file
#   --no-structs              do not write match structures (default: writes $out/structs)
#   --structOut <dir>         match structure dir  (default: $out/structs)
#   --matchOut <file>         final merged match-address file (default: $out/match.txt)
#   --master <path>           master binary (default: found via PATH)
#   --force                   reprocess all pieces even if done
#
# Notes:
#   --topN/--minN are applied PER PIECE, not globally across the whole database.
#   For a strict global top-N, set --topN per piece of equal size, or merge and
#   re-rank by RMSD yourself. This is a documented limitation of the split design.
#
# Deps: master (MASTER, grigoryanlab.org, e.g. /usr/local/bin).

set -euo pipefail

# ---------------- config ----------------
MASTER=""
QUERY=""
TARGETLIST=""
RMSDCUT=""
OUT=""
PER_LIST=10000          # --target-num-per-list
NJOBS=8                 # concurrency cap only
BBRMSD=""
TOPN=""
MINN=""
GAPLEN=""
OUTTYPE=""
SEQOUT=""
STRUCT_OUT=""
NO_STRUCTS=""
MATCH_OUT=""
FORCE=""

usage() {
    cat >&2 <<EOF
Usage: $0 --query <q.pds> --targetList <list> --rmsdCut <X> --out <dir> [options]

Required:
  --query <q.pds>          query structure already converted to PDS
  --targetList <list>      file, one target .pds path per line (the database)
  --rmsdCut <x>            RMSD cutoff (Angstrom) for a match
  --out <dir>              output root

Options:
  --target-num-per-list <N> structures per piece (default 10000); the resume granularity
  --njobs <N>               pieces searched concurrently (default 8; cap only, not split size)
  --bbRMSD                  full-backbone RMSD instead of CA-only
  --topN <N>                keep best N matches (per piece)
  --minN <N>                at least N best matches (per piece)
  --gapLen <s>              gap restraints, e.g. '1-10;0-3'
  --outType <t>             match|full|wgap
  --seqOut <file>           also write match sequences, merged to this file
  --no-structs              do not write match structures
  --structOut <dir>         match structure dir  (default: \$out/structs)
  --matchOut <file>         final merged match-address file (default: \$out/match.txt)
  --master <path>           master binary (default: via PATH)
  --force                   reprocess all pieces even if done
EOF
}

# ---------------- parse args ----------------
while [ $# -gt 0 ]; do
    case "$1" in
        --query)       QUERY="$2"; shift 2 ;;
        --targetList)  TARGETLIST="$2"; shift 2 ;;
        --rmsdCut)     RMSDCUT="$2"; shift 2 ;;
        --out)         OUT="$2"; shift 2 ;;
        --target-num-per-list) PER_LIST="$2"; shift 2 ;;
        --njobs)       NJOBS="$2"; shift 2 ;;
        --master)      MASTER="$2"; shift 2 ;;
        --topN)        TOPN="$2"; shift 2 ;;
        --minN)        MINN="$2"; shift 2 ;;
        --gapLen)      GAPLEN="$2"; shift 2 ;;
        --outType)     OUTTYPE="$2"; shift 2 ;;
        --seqOut)      SEQOUT="$2"; shift 2 ;;
        --structOut)   STRUCT_OUT="$2"; shift 2 ;;
        --matchOut)    MATCH_OUT="$2"; shift 2 ;;
        --bbRMSD)      BBRMSD=1; shift ;;
        --no-structs)  NO_STRUCTS=1; shift ;;
        --force)       FORCE=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------------- validate ----------------
err=0
[ -n "$QUERY" ]      || { echo "missing required: --query" >&2; err=1; }
[ -n "$TARGETLIST" ] || { echo "missing required: --targetList" >&2; err=1; }
[ -n "$RMSDCUT" ]    || { echo "missing required: --rmsdCut" >&2; err=1; }
[ -n "$OUT" ]        || { echo "missing required: --out" >&2; err=1; }
if [ "$err" -ne 0 ]; then usage; exit 1; fi

if ! [[ "$NJOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "invalid --njobs: '$NJOBS' (must be a positive integer)" >&2; exit 1
fi
if ! [[ "$PER_LIST" =~ ^[1-9][0-9]*$ ]]; then
    echo "invalid --target-num-per-list: '$PER_LIST' (must be a positive integer)" >&2; exit 1
fi

# master binary
if [ -z "$MASTER" ]; then
    command -v master >/dev/null 2>&1 || { echo "master not found on PATH (MASTER install, e.g. /usr/local/bin)" >&2; exit 1; }
    MASTER="$(command -v master)"
fi
[ -x "$MASTER" ] || { echo "master not executable: $MASTER" >&2; exit 1; }

# query must be a .pds produced by createPDS --type query
[ -f "$QUERY" ] || { echo "query PDS not found: $QUERY" >&2; exit 1; }
case "$QUERY" in
    *.pds) : ;;
    *) echo "query must be a .pds file (run createPDS --type query first)" >&2; exit 1 ;;
esac

[ -f "$TARGETLIST" ] || { echo "targetList not found: $TARGETLIST" >&2; exit 1; }

# ---------------- canonicalize ----------------
OUT="$(cd "$(dirname "$OUT")" && pwd)/$(basename "$OUT")"  # absolute
CHUNKDIR="$OUT/.chunks"
LOG="$OUT/master_search.log"
mkdir -p "$OUT"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# defaults for optional outputs
[ -z "$MATCH_OUT" ] && MATCH_OUT="$OUT/match.txt"
if [ -z "$NO_STRUCTS" ]; then
    [ -z "$STRUCT_OUT" ] && STRUCT_OUT="$OUT/structs"
    mkdir -p "$STRUCT_OUT"
fi

total=$(wc -l < "$TARGETLIST")
if [ "$total" -eq 0 ]; then
    log "targetList is empty: nothing to search."
    : > "$MATCH_OUT"
    exit 0
fi

# ---------------- fingerprint: invalidate sentinels when split input changes ----------------
# Split depends only on (targetList content) + PER_LIST, so NJOBS must NOT be in the
# fingerprint (changing --njobs must never force a redo; it only changes concurrency).
newfp=$( { cat "$TARGETLIST"; echo "per_list=$PER_LIST"; } | sha256sum | cut -d' ' -f1 )
if [ -f "$CHUNKDIR/fingerprint" ] && [ "$(cat "$CHUNKDIR/fingerprint")" = "$newfp" ] \
   && [ -z "$FORCE" ]; then
    log "Resuming: --targetList / --target-num-per-list unchanged since previous run."
else
    rm -rf "$CHUNKDIR"
    mkdir -p "$CHUNKDIR"
    printf '%s\n' "$newfp" > "$CHUNKDIR/fingerprint"
    [ -n "$FORCE" ] && log "Forcing reprocess of all pieces."
fi

# number of pieces = ceil(total / PER_LIST)
npieces=$(( (total + PER_LIST - 1) / PER_LIST ))

log "Searching $total targets against query $QUERY (pieces=$npieces, per-piece<=$PER_LIST, njobs=$NJOBS cap, rmsdCut=$RMSDCUT)"

# build piece lists (deterministic consecutive split of the given order)
for i in $(seq 1 "$npieces"); do
    start=$(( (i-1)*PER_LIST + 1 ))
    end=$(( start + PER_LIST - 1 )); [ "$end" -gt "$total" ] && end=$total
    sed -n "${start},${end}p" "$TARGETLIST" > "$CHUNKDIR/targets.$i"
done

# ---------------- run pieces (njobs at once) with resume ----------------
run_piece() {
    # $1 = piece index. A done piece is skipped (this is the resume base).
    local i="$1" cl
    cl="$CHUNKDIR/targets.$i"
    if [ -f "$CHUNKDIR/.done.$i" ]; then
        log "  piece $i: already done (resume skip)"
        return 0
    fi
    # per-piece match/seq go under .chunks/, structs share the final --structOut dir
    # (pieces are disjoint target sets => no filename collisions).
    local cmd=( "$MASTER" --query "$QUERY" --targetList "$cl" --rmsdCut "$RMSDCUT" )
    [ -n "$BBRMSD" ]  && cmd+=( --bbRMSD )
    [ -n "$TOPN" ]    && cmd+=( --topN "$TOPN" )
    [ -n "$MINN" ]    && cmd+=( --minN "$MINN" )
    [ -n "$GAPLEN" ]  && cmd+=( --gapLen "$GAPLEN" )
    [ -n "$OUTTYPE" ] && cmd+=( --outType "$OUTTYPE" )
    cmd+=( --matchOut "$CHUNKDIR/match.$i" )
    [ -n "$SEQOUT" ]  && cmd+=( --seqOut "$CHUNKDIR/seq.$i" )
    if [ -n "$STRUCT_OUT" ]; then
        mkdir -p "$STRUCT_OUT"
        cmd+=( --structOut "$STRUCT_OUT" )
    fi
    if "${cmd[@]}" > "$CHUNKDIR/run.$i.log" 2>&1; then
        local nm=0; [ -f "$CHUNKDIR/match.$i" ] && nm=$(wc -l < "$CHUNKDIR/match.$i")
        : > "$CHUNKDIR/.done.$i"
        log "  piece $i: OK ($nm matches)"
    else
        local rc=$?
        log "  piece $i FAILED (exit $rc); tail:"
        tail -n 15 "$CHUNKDIR/run.$i.log" | sed 's/^/    /'
        return 1
    fi
}

# Rolling pool: at most NJOBS pieces in flight, then wait the oldest before
# launching the next. `active` holds the pids of currently-running pieces; each
# pid is waited exactly once (drain after the loop).
failed=0
active=()
for i in $(seq 1 "$npieces"); do
    ( run_piece "$i" ) &
    active+=( $! )
    if [ "${#active[@]}" -ge "$NJOBS" ]; then
        if ! wait "${active[0]}"; then failed=$((failed+1)); fi
        active=( "${active[@]:1}" )
    fi
done
for p in "${active[@]}"; do
    if ! wait "$p"; then failed=$((failed+1)); fi
done

# ---------------- merge finished pieces ----------------
: > "$MATCH_OUT"
for i in $(seq 1 "$npieces"); do
    if [ -f "$CHUNKDIR/.done.$i" ] && [ -f "$CHUNKDIR/match.$i" ]; then
        cat "$CHUNKDIR/match.$i" >> "$MATCH_OUT"
    fi
done
log "  merged $(wc -l < "$MATCH_OUT") matches -> $MATCH_OUT"

if [ -n "$SEQOUT" ]; then
    : > "$SEQOUT"
    for i in $(seq 1 "$npieces"); do
        if [ -f "$CHUNKDIR/.done.$i" ] && [ -f "$CHUNKDIR/seq.$i" ]; then
            cat "$CHUNKDIR/seq.$i" >> "$SEQOUT"
        fi
    done
    log "  merged seqs -> $SEQOUT"
fi

if [ "$failed" -gt 0 ]; then
    log "WARNING: $failed pieces failed. Re-run the same command to resume (only unfinished/failed pieces reprocess)."
    exit 2
fi
log "All $npieces pieces done."

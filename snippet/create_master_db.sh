#!/bin/bash
# Build a MASTER structural-search database from the local RCSB PDB archive.
#
# Layout mirrors the RCSB *divided* archive: every result is indexed by the
# middle 2 chars of the 4-char PDB id, i.e. pdb100d -> idx "00".
#
#   Input:   <rcsb>/{aa}/pdbXXXX.ent.gz          (gzipped PDB, divided)
#   Output:  <out>/pdb/{aa}/pdbXXXX.pdb          (decompressed, divided)
#            <out>/pdb_filtered/{aa}/pdbXXXX.pdb (kept single-model PDBs, divided;
#                                                 empty {aa} dirs are created so the
#                                                 index set is complete)
#            <out>/pdb_filtered/{aa}/pdbXXXX/pdbXXXX_m{n}.pdb
#                                                 (kept multi-ensemble PDBs: split per
#                                                 MODEL block into a per-id subdir;
#                                                 incomplete blocks are dropped)
#            <out>/pds/{aa}/pdbXXXX.pds          (createPDS target, divided, mirroring
#                                                 the pdb_filtered layout; the actual DB)
#            <out>/pdb_filtered_list.txt  createPDS --pdbList (absolute, incremental: only
#                                                 PDBs missing a .pds)
#            <out>/build_pds_list.txt     createPDS --pdsList (their absolute .pds paths)
#            <out>/target_list.txt        full .pds list for master --targetList
#
# createPDS only reads *uncompressed* PDB text, so gz must be decompressed first.
# Steps are idempotent: rerunning skips already-done work.
#
# Incremental skip is decided per index {aa} (see gen_build_list): an index whose
# pdb_filtered/{aa} is empty, or whose .pdb count equals its .pds count with no
# stray files, is treated as complete and skipped; otherwise it is descended into
# and each .pdb is compared individually.
#
# --rcsb / --out are canonicalized to ABSOLUTE paths up front, so every derived
# path and every line written inside the lists is absolute.
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
  decompress         gunzip .ent.gz -> out/pdb/{aa}/*.pdb
  filter             keep PDBs w/ any MASTER-supported residue -> out/pdb_filtered/{aa}
  gen_build_list     write incremental lists (PDBs missing a .pds, per-index) -> build uses them
  build              createPDS (via incremental lists) -> out/pds/{aa}/*.pds
  gen_target_list    write full .pds list out/target_list.txt (for master --targetList)
  verify             compare PDS count vs target_list count
  all                decompress + filter + gen_build_list + build + gen_target_list + verify

Notes:
  createPDS must be on PATH (MASTER install, e.g. /usr/local/bin).
  All paths are canonicalized to absolute before use. All output is divided by
  the middle 2 chars of the PDB id (RCSB-style).
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

# Index of a pdbXXXX(.pdb/.pds/.ent.gz) basename: middle 2 chars of the 4-char id.
#   "pdb100d" -> id "100d" -> idx "00"
idx_of() { local b="$1"; b="${b%.*}"; b="${b#pdb}"; echo "${b:1:2}"; }

# ---------------- action: decompress ----------------
decompress() {
    local NJOBS; NJOBS=${NJOBS:-$(nproc)}
    # bash -c 会新开子 shell,不继承普通 shell 变量,须 export 才能读到路径
    export PDB_DIR RCSB_DIR
    log "Decompressing $RCSB_DIR/*/*.ent.gz -> $PDB_DIR/{aa}/  ($NJOBS 并行,流式)"
    mkdir -p "$PDB_DIR"

    # find 流式输出:边遍历边喂给 xargs,不等 glob 全部展开
    # 已存在则跳过:worker 内用 [ -f ] 判断,避免重复解压
    find "$RCSB_DIR" -name '*.ent.gz' -print0 | \
    xargs -0 -P "$NJOBS" -n1 bash -c '
        f="$1"
        base="$(basename "$f" .ent.gz)"            # pdb100d
        idx="$(idx_of "$base")"                     # 00   (function exported below)
        [ -z "$idx" ] && { echo "bad basename: $f" >&2; exit 0; }
        out="$PDB_DIR/$idx/$base.pdb"
        [ -f "$out" ] && exit 0                      # 目标已存在 → 跳过
        mkdir -p "$PDB_DIR/$idx"
        gunzip -c "$f" > "$out"
    ' _
    # NOTE: idx_of is defined above; make it visible to the bash -c subshells.
    # (bash -c starts a fresh shell, so export the function before the pipeline.)
    local n; n=$(find "$RCSB_DIR" -name '*.ent.gz' | wc -l)
    log "Processed $n gz entries (含已存在跳过)."
}
# export idx_of so the decompress/filter bash -c workers can call it
export -f idx_of

# filter_one <pdb_file> <out_dir>: keep-or-split a single decompressed PDB.
#   * no MODEL record or exactly 1 MODEL block -> single model: keep as-is under
#     out_dir/{base}.pdb iff it has >=1 complete residue (CA,C,N,O on same chain:resSeq:iCode).
#   * >=2 MODEL blocks -> multi-ensemble: split each block into its own file,
#     written flat inside out_dir/{base}/pdbXXXX_m{n}.pdb (n = 1-based MODEL index),
#     keeping only blocks that contain >=1 complete residue. This keeps the {aa}
#     index dirs unchanged (each multi-model PDB gets one per-id subdir, no new
#     top-level split), so gen_build_list's depth-1 index scan and the mirror-image
#     fast-path count stay valid.
# Always exits 0 (worker must not abort xargs under set -e).
filter_one() {
    local f="$1" OUT="$2" base nmodels
    base="$(basename "$f" .pdb)"
    nmodels=$(grep -c '^MODEL' "$f" 2>/dev/null || true)
    [ -z "$nmodels" ] && nmodels=0

    if [ "$nmodels" -le 1 ]; then
        # single model: keep flat, unchanged name
        if awk -v supported="$MASTER_SUPPORTED" '
            BEGIN { n=split(supported,a," "); for(i=1;i<=n;i++) allowed[a[i]]=1 }
            /^(ATOM  |HETATM)/ {
                rn = substr($0,18,3)
                if (!(rn in allowed)) next
                key = substr($0,22,1) ":" substr($0,23,4) ":" substr($0,27,1)
                atom = substr($0,13,4); gsub(/ /,"",atom)
                if      (atom=="CA") hasCA[key]=1
                else if (atom=="C")  hasC[key]=1
                else if (atom=="N")  hasN[key]=1
                else if (atom=="O")  hasO[key]=1
            }
            END { for (k in hasCA) { if (hasC[k] && hasN[k] && hasO[k]) exit 0 } exit 1 }
        ' "$f"; then
            cp "$f" "$OUT/$base.pdb"
        fi
        return 0
    fi

    # multi-ensemble: split per MODEL block into $OUT/$base/, keep complete blocks.
    # Streaming version: each line is written straight to the per-model file
    # (print > fname) instead of accumulating the whole block in a buf string.
    # mawk's `buf = buf $0` is O(n^2) (immutable strings, full re-copy per append),
    # which turns large ensembles (e.g. 76 MB / ~10 models) into hours of work.
    # Write-first then rm incomplete blocks at finalize -> O(n), same output.
    mkdir -p "$OUT/$base"
    awk -v out="$OUT/$base" -v base="$base" -v supported="$MASTER_SUPPORTED" '
        BEGIN {
            ns=split(supported,a," "); for(i=1;i<=ns;i++) allowed[a[i]]=1
            thisn=0
        }
        /^MODEL/ {
            if (thisn>0) finalize()
            if (match($0,/[0-9]+/)) thisn=substr($0,RSTART,RLENGTH)+0
            else thisn++
            delete ca; delete c; delete n; delete o
            fname = sprintf("%s/%s_m%d.pdb", out, base, thisn)
        }
        /^(ATOM  |HETATM)/ {
            if (thisn==0) next
            rn=substr($0,18,3)
            if (rn in allowed) {
                key=substr($0,22,1) ":" substr($0,23,4) ":" substr($0,27,1)
                atom=substr($0,13,4); gsub(/ /,"",atom)
                if      (atom=="CA") ca[key]=1
                else if (atom=="C")  c[key]=1
                else if (atom=="N")  n[key]=1
                else if (atom=="O")  o[key]=1
            }
            print $0 > fname
        }
        /^ENDMDL/ { if (thisn>0) { finalize(); thisn=0 } }
        END { if (thisn>0) finalize() }
        function finalize(   k,ok){
            ok=0
            for (k in ca) if (c[k]&&n[k]&&o[k]) { ok=1; break }
            if (ok) {
                print "END" >> fname
            } else {
                system(sprintf("rm -f %s", fname))
            }
            close(fname)
        }
    ' "$f"
    # completion sentinel: written ONLY after the whole file is fully processed
    # (including the all-models-incomplete case, which drops everything).
    # filter() skips a multi-model file iff this sentinel exists, so an
    # interrupted run (partial _m* written, no sentinel) is re-processed next
    # time -> no silent data loss. It is neither .pdb nor .pds, so it does not
    # perturb the gen_build_list count fast-path.
    : > "$OUT/$base/.done"
    return 0
}
export -f filter_one

# ---------------- action: filter (keep PDBs with >=1 MASTER-supported residue) ----------------
filter() {
    if [ ! -d "$PDB_DIR" ]; then
        log "ERROR: $PDB_DIR missing. Run: $0 --rcsb ... --out ... decompress first." >&2
        exit 1
    fi
    local NJOBS; NJOBS=${NJOBS:-$(nproc)}
    export PDB_FILTERED MASTER_SUPPORTED
    log "Filtering $PDB_DIR/{aa} -> $PDB_FILTERED/{aa} (single keep / multi split; keep >=1 full MASTER residue: CA, C, N, O)"
    mkdir -p "$PDB_FILTERED"

    # ensure an index dir exists for every index present in pdb/ (even if empty),
    # so gen_build_list has a complete index universe to scan.
    find "$PDB_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r idir; do
        mkdir -p "$PDB_FILTERED/$(basename "$idir")"
    done

    # worker: single model -> $idx/base.pdb (flat, unchanged name);
    # multi-model -> $idx/base/*.pdb (per-id subdir, _m{n}) + a .done sentinel.
    # Idempotent: skip if the single flat file exists, or the multi-model .done
    # sentinel exists (written only after that file was fully processed, so an
    # interrupted run is re-processed -> no silent data loss). awk's exit code
    # only decides keep/copy; the worker always exits 0, otherwise filtered files
    # would make xargs return 123 and abort under set -e.
    find "$PDB_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.pdb' -print0 | \
    xargs -0 -P "$NJOBS" -n1 bash -c '
        f="$1"
        idx="$(basename "$(dirname "$f")")"
        id="$(basename "$f")"
        base="${id%.pdb}"
        [ -f "$PDB_FILTERED/$idx/$id" ] && exit 0                    # single already done
        [ -f "$PDB_FILTERED/$idx/$base/.done" ] && exit 0            # multi already done (sentinel)
        filter_one "$f" "$PDB_FILTERED/$idx"
        exit 0
    ' _

    local n; n=$(find "$PDB_FILTERED" -type f -name '*.pdb' | wc -l)
    log "  $n PDB model files kept in $PDB_FILTERED."
}

# ---------------- action: gen_build_list (incremental; run BEFORE build) ----------------
# Incremental skip is decided PER INDEX {aa}:
#   - index dir empty                        -> nothing to build (complete)
#   - .pdb count == .pds count             -> complete (fast path, skip)
#   - otherwise descend: compare each .pdb individually against its .pds
# Writes pdb_filtered_list.txt (missing .pdb inputs) and build_pds_list.txt
# (their .pds outputs). Lines are absolute paths.
gen_build_list() {
    if [ ! -d "$PDB_FILTERED" ]; then
        log "ERROR: $PDB_FILTERED missing. Run: $0 ... filter first." >&2
        exit 1
    fi
    mkdir -p "$PDS_DIR"
    : > "$PDB_FILTERED_LIST"
    : > "$BUILD_PDS_LIST"
    log "Writing incremental build lists (per-index, PDBs missing a .pds) to $PDS_DIR"

    find "$PDB_FILTERED" -mindepth 1 -maxdepth 1 -type d | sort | while read -r idir; do
        idx="$(basename "$idir")"
        # counts are RECURSIVE now: single models sit flat in {aa}, multi-model
        # splits live in {aa}/<id>/ subdirs. Counting all depths keeps the fast
        # path consistent with the mirror-image pds/ layout.
        npdb=$(find "$idir" -type f -name '*.pdb' | wc -l)
        # pds/{idx} 可能还不存在: find 会因 No such file 返回非零,
        # 在 set -euo pipefail 下会让 $( ) 赋值失败并中断整个 while 管道,
        # 故追加 || true,让"目录不存在"计为 0,而非致命错误。
        npds=$(find "$PDS_DIR/$idx" -type f -name '*.pds' 2>/dev/null | wc -l || true)

        # 空索引 → 无需建库
        [ "$npdb" -eq 0 ] && continue
        # 快速路径: .pdb 与 .pds 个数一致 → 判定该索引完整并跳过。
        # 哨兵 .done 既非 .pdb 也非 .pds,不参与任一计数,故不影响快路径;
        # 安全前提是管线纯增量(不删文件),否则由下方逐文件下钻兜底。
        if [ "$npdb" -eq "$npds" ]; then
            continue
        fi

        # 否则下钻逐个比较(递归):为每个尚缺 .pds 的 .pdb 写一对绝对路径。
        # .pds 位置镜像 .pdb 的相对路径:pdb_filtered/{aa}/... -> pds/{aa}/...
        find "$idir" -type f -name '*.pdb' | while read -r f; do
            rel="${f#$PDB_FILTERED/}"                     # {aa}/.../pdbXXXX.pdb (结构可含子目录)
            pout="$PDS_DIR/${rel%.pdb}.pds"
            if [ ! -e "$pout" ]; then
                mkdir -p "$(dirname "$pout")"
                printf '%s\n' "$f"     >> "$PDB_FILTERED_LIST"
                printf '%s\n' "$pout"  >> "$BUILD_PDS_LIST"
            fi
        done
    done

    log "  $(wc -l < "$PDB_FILTERED_LIST") PDB inputs pending, $(wc -l < "$BUILD_PDS_LIST") PDS outputs pending."
}

# ---------------- action: gen_target_list (complete .pds list; run AFTER build) ----------------
gen_target_list() {
    mkdir -p "$PDS_DIR"
    log "Writing $TARGET_LIST (full .pds list for master --targetList)"
    find "$PDS_DIR" -type f -name '*.pds' | sort > "$TARGET_LIST"
    log "  $(wc -l < "$TARGET_LIST") PDS entries."
}

# ---------------- action: build (convert to PDS) ----------------
build() {
    if [ ! -s "$PDB_FILTERED_LIST" ] || [ ! -s "$BUILD_PDS_LIST" ]; then
        log "ERROR: $PDB_FILTERED_LIST or $BUILD_PDS_LIST missing/empty. Run: $0 ... gen_build_list first." >&2
        exit 1
    fi
    log "Converting PDB -> PDS (target) into $PDS_DIR/{aa}/ (cleanPDB='${CLEAN_PDB:-off}')"
    # Output .pds paths come from --pdsList (build_pds_list.txt), so createPDS writes
    # them into $PDS_DIR/{aa} regardless of cwd. Parent dirs were created by gen_build_list.
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
    npds=$(find "$PDS_DIR" -type f -name '*.pds' | wc -l)
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

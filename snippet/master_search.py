#!/usr/bin/env python3
"""Run a MASTER structural-motif SEARCH (the `master` binary only) against a built
PDS database, with piece-wise parallelization, resume, and a tqdm progress bar.

This is the query side. It does NOT build a DB (that's create_master_db.sh) and
does NOT convert PDB -> PDS (that's createPDS). It only takes a query .pds and a
database --targetList, and runs the `master` search binary against it.

Design (kept from the original master_search.sh):

Split-by-count = the foundation of resume:
  --targetList is divided into consecutive pieces of --chunk-size
  structures each (last piece may be smaller). Each piece is one
  `master --targetList <piece>` job. A `.done.<i>` sentinel is written only after
  a piece finishes successfully, so a re-run of the same command skips every
  already-done piece and only searches the unfinished ones, then re-merges the
  piece outputs.

--njobs is ONLY a concurrency cap: it controls how many pieces are searched at the
same time (to utilize all cores), NOT the split size.

A fingerprint of (targetList content + --chunk-size) is stored in
out/.chunks/. If either changes, all `.done.*` sentinels are invalidated so a
stale "done" can never wrongly skip work that now needs doing.

Output layout (differs from the bash version -- fixes two bugs):

  OUT/
    match.txt                  merged match-address file (default --matchOut)
    seqs.txt                   merged match sequences      (if --seqOut)
    structs/                   final --structOut dir (default OUT/structs)
      piece.<i>/match1.pdb ... per-piece subdir  <-- no cross-piece overwrite
    pieces/
      piece.<i>/match.txt      per-piece match (result, NOT hidden in .chunks)
      piece.<i>/seq.txt        per-piece seq  (if --seqOut)
    .chunks/                   WORK/metadata only: fingerprint, targets.<i>,
                               run.<i>.log, done.<i> sentinels. No result files.

Bugs fixed vs. the bash original:
  1. Structures no longer overwrite each other. `master` names every piece's
     output `match1.pdb, match2.pdb, ...` (numbered from 1, see Search.cpp
     renameStruct). Running pieces concurrently into ONE shared --structOut dir
     made piece N clobber piece M's `matchN.pdb`. Fix: each piece writes into its
     own --structOut/piece.<i>/ subdir.
  2. Per-piece match files no longer land in the hidden .chunks work dir. They
     now live under OUT/pieces/piece.<i>/ (real results), while .chunks holds
     only resume metadata (fingerprint, sentinels, logs).

Usage:
  master_search.py --query <q.pds> --targetList <list> --rmsdCut <X> \
                   --out <dir> [options]

Required:
  --query <q.pds>          query structure already converted to PDS
  --targetList <list>      file, one target .pds path per line (the database)
  --rmsdCut <x>            RMSD cutoff (Angstrom) for a match
  --out <dir>              output root (per-piece files + final merged results)

Options / passthrough to `master`:
  --chunk-size <N>          structures per piece (default 10000); sets the resume granularity
  --njobs <N>               how many pieces to search concurrently (default 8)
  --bbRMSD                  full-backbone RMSD instead of CA-only
  --topN <N>                keep best N matches (global top-N after merge, see Notes)
  --minN <N>                return at least N best matches (per piece)
  --gapLen <s>              gap restraints, e.g. '1-10;0-3'
  --outType <t>             match|full|wgap (region written by structOut)
  --seqOut <file>           also write match sequences, merged to this file
  --no-structs              do not write match structures (default: writes $out/structs)
  --structOut <dir>         match structure dir  (default: $out/structs)
  --matchOut <file>         final merged match-address file (default: $out/match.txt)
  --master <path>           master binary (default: found via PATH)
  --force                   reprocess all pieces even if done

Notes:
  --topN is a GLOBAL top-N: it is forwarded to each piece (so a piece never
  carries more than N candidates, keeping per-piece files small) and then
  re-applied after the merge, where all pieces' matches are sorted by RMSD
  ascending and truncated to the overall best N. This is exact — a global top-N
  match is always inside some piece's per-piece top-N (a piece that discards a
  match m has N matches with RMSD <= m, so m can never rank in the global top-N).
  --minN stays per piece (it widens the cutoff so a piece returns at least N
  matches). The merged match.txt (and seqs.txt, aligned to it) is RMSD-sorted.
  Adding an output flag (e.g. --seqOut) after pieces already ran does NOT
  invalidate sentinels -> done pieces emit nothing for it; use --force to rebuild.

Deps: master (MASTER, grigoryanlab.org, e.g. /usr/local/bin); tqdm (optional).
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

_LOG = None  # set in main; used by log()


def log(msg):
    line = "[%s] %s" % (datetime.now().strftime("%F %T"), msg)
    if _LOG is not None:
        with open(_LOG, "a") as f:
            f.write(line + "\n")
    if tqdm is not None:
        tqdm.write(line)
    else:
        print(line, file=sys.stderr)


def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="master_search.py",
        description="Parallel, resumable MASTER structural-motif search (query side only).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--query", required=True, help="query structure already converted to PDS")
    p.add_argument("--targetList", required=True, help="file, one target .pds path per line")
    p.add_argument("--rmsdCut", required=True, help="RMSD cutoff (Angstrom) for a match")
    p.add_argument("--out", required=True, help="output root")
    p.add_argument("--chunk-size", dest="chunk_size", type=int, default=10000,
                   help="structures per piece; the resume granularity")
    p.add_argument("--njobs", type=int, default=8,
                   help="pieces searched concurrently (cap only, not split size)")
    p.add_argument("--master", default=None, help="master binary (default: via PATH)")
    p.add_argument("--topN", type=int, default=None, help="keep best N matches (per piece)")
    p.add_argument("--minN", type=int, default=None, help="return at least N best matches (per piece)")
    p.add_argument("--gapLen", default=None, help="gap restraints, e.g. '1-10;0-3'")
    p.add_argument("--outType", default=None, help="match|full|wgap (region written by structOut)")
    p.add_argument("--seqOut", default=None, help="also write match sequences, merged to this file")
    p.add_argument("--structOut", default=None, help="match structure dir (default: $out/structs)")
    p.add_argument("--matchOut", default=None, help="final merged match-address file (default: $out/match.txt)")
    p.add_argument("--bbRMSD", action="store_true", help="full-backbone RMSD instead of CA-only")
    p.add_argument("--no-structs", action="store_true", help="do not write match structures")
    p.add_argument("--force", action="store_true", help="reprocess all pieces even if done")
    return p.parse_args(argv)


def validate(a):
    if a.njobs < 1:
        sys.exit("invalid --njobs: %r (must be a positive integer)" % a.njobs)
    if a.chunk_size < 1:
        sys.exit("invalid --chunk-size: %r (must be a positive integer)" % a.chunk_size)
    if a.master is None:
        a.master = shutil.which("master")
    if a.master is None:
        sys.exit("master not found on PATH (MASTER install, e.g. /usr/local/bin)")
    if not os.access(a.master, os.X_OK):
        sys.exit("master not executable: %s" % a.master)
    if not os.path.isfile(a.query):
        sys.exit("query PDS not found: %s" % a.query)
    if not a.query.endswith(".pds"):
        sys.exit("query must be a .pds file (run createPDS --type query first)")
    if not os.path.isfile(a.targetList):
        sys.exit("targetList not found: %s" % a.targetList)


def count_lines(path):
    n = 0
    with open(path, "rb") as f:
        for _ in f:
            n += 1
    return n


def read_last_lines(path, n):
    try:
        with open(path, "rb") as f:
            data = f.read()
        lines = data.decode("utf-8", "replace").splitlines()
        return lines[-n:]
    except OSError:
        return []


def fingerprint(targetlist, per_list):
    h = hashlib.sha256()
    with open(targetlist, "rb") as f:
        h.update(f.read())
    h.update(("per_list=%d" % per_list).encode())
    return h.hexdigest()


def write_piece_targets(targetlist, total, per_list, work_dir):
    """Write .chunks/targets.<i> for each piece; return npieces."""
    npieces = (total + per_list - 1) // per_list
    for i in range(1, npieces + 1):
        start = (i - 1) * per_list
        end = min(start + per_list, total)
        piece_list = os.path.join(work_dir, "targets.%d" % i)
        with open(targetlist) as src, open(piece_list, "w") as dst:
            for _ in range(start):  # skip to this piece's first line
                src.readline()
            for _ in range(end - start):
                dst.write(src.readline())
    return npieces


def build_cmd(a, piece_list, match_out, seq_out, struct_dir):
    cmd = [a.master, "--query", a.query, "--targetList", piece_list,
           "--rmsdCut", a.rmsdCut]
    if a.bbRMSD:
        cmd.append("--bbRMSD")
    # --topN is forwarded per piece AND re-applied globally after the merge.
    # Both are safe: a global top-N match is always inside some piece's top-N
    # (proof: if a piece discards m there are N matches with RMSD <= m in it, so
    # m cannot rank in the global top-N). Per-piece trimming only shrinks the
    # candidate set / file sizes. --minN widens the cutoff per piece so a piece
    # returns at least N matches.
    if a.topN is not None:
        cmd += ["--topN", str(a.topN)]
    if a.minN is not None:
        cmd += ["--minN", str(a.minN)]
    if a.gapLen:
        cmd += ["--gapLen", a.gapLen]
    if a.outType:
        cmd += ["--outType", a.outType]
    cmd += ["--matchOut", match_out]
    if seq_out:
        cmd += ["--seqOut", seq_out]
    if struct_dir:
        cmd += ["--structOut", struct_dir]
    return cmd


def run_piece(a, i, work_dir, pieces_root, struct_root):
    """Search piece i. Returns (ok, msg). Writes .done sentinel on success."""
    piece_dir = os.path.join(pieces_root, "piece.%d" % i)
    os.makedirs(piece_dir, exist_ok=True)
    piece_list = os.path.join(work_dir, "targets.%d" % i)
    match_out = os.path.join(piece_dir, "match.txt")
    seq_out = os.path.join(piece_dir, "seq.txt") if a.seqOut else None
    struct_dir = None
    if not a.no_structs:
        struct_dir = os.path.join(struct_root, "piece.%d" % i)
        os.makedirs(struct_dir, exist_ok=True)
    run_log = os.path.join(work_dir, "run.%d.log" % i)
    cmd = build_cmd(a, piece_list, match_out, seq_out, struct_dir)
    try:
        with open(run_log, "wb") as rf:
            rc = subprocess.run(cmd, stdout=rf, stderr=subprocess.STDOUT).returncode
    except OSError as e:
        return False, "spawn error: %s" % e
    if rc != 0:
        tail = read_last_lines(run_log, 15)
        return False, "exit %d; log tail:\n%s" % (rc, "\n".join("    " + t for t in tail))
    nm = count_lines(match_out) if os.path.isfile(match_out) else 0
    with open(os.path.join(work_dir, "done.%d" % i), "w"):
        pass
    return True, "%d matches" % nm


def _rmsd(line):
    """Leading RMSD of a match line (first whitespace-separated token)."""
    try:
        return float(line.split()[0])
    except (ValueError, IndexError):
        return float("inf")


def merge_from_done(a, npieces, work_dir, pieces_root, match_out, seq_out):
    """Merge done pieces, sort globally by RMSD ascending, keep overall top-N.

    Collects (match, seq) pairs in piece order, sorts by the match line's RMSD
    (stable on ties), truncates to --topN if set, and writes BOTH files in the
    same sorted order so a match line and its sequence stay aligned.
    """
    rows = []  # (rmsd, order, match_line, seq_line_or_None)
    order = 0
    for i in range(1, npieces + 1):
        if not os.path.isfile(os.path.join(work_dir, "done.%d" % i)):
            continue
        pm = os.path.join(pieces_root, "piece.%d" % i, "match.txt")
        ps = os.path.join(pieces_root, "piece.%d" % i, "seq.txt") if seq_out else None
        if not os.path.isfile(pm):
            continue
        if ps:
            with open(pm) as mf, open(ps) as sf:
                for ml in mf:
                    rows.append((_rmsd(ml), order, ml, sf.readline()))
                    order += 1
        else:
            with open(pm) as mf:
                for ml in mf:
                    rows.append((_rmsd(ml), order, ml, None))
                    order += 1
    rows.sort(key=lambda t: (t[0], t[1]))
    if a.topN and len(rows) > a.topN:
        rows = rows[:a.topN]
    with open(match_out, "w") as mf:
        for _r, _o, ml, _sl in rows:
            mf.write(ml)
    log("  merged %d matches (RMSD-sorted%s) -> %s" % (
        len(rows), ", top %d" % a.topN if a.topN else "", match_out))
    if seq_out:
        with open(seq_out, "w") as sf:
            for _r, _o, _ml, sl in rows:
                if sl is not None:
                    sf.write(sl)
        log("  merged %d seqs (aligned to matches) -> %s" % (len(rows), seq_out))


def main(argv):
    global _LOG
    a = parse_args(argv)
    validate(a)

    a.out = os.path.abspath(a.out)
    work_dir = os.path.join(a.out, ".chunks")
    pieces_root = os.path.join(a.out, "pieces")
    match_out = a.matchOut or os.path.join(a.out, "match.txt")
    if a.no_structs:
        struct_root = None
    else:
        struct_root = a.structOut or os.path.join(a.out, "structs")
        os.makedirs(struct_root, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(pieces_root, exist_ok=True)
    _LOG = os.path.join(a.out, "master_search.log")

    total = count_lines(a.targetList)
    if total == 0:
        log("targetList is empty: nothing to search.")
        open(match_out, "w").close()
        return 0

    # ---- fingerprint: invalidate sentinels when the split input changes ----
    # Split depends only on (targetList content) + PER_LIST, so NJOBS must NOT be
    # in the fingerprint (changing --njobs must never force a redo).
    newfp = fingerprint(a.targetList, a.chunk_size)
    fp_file = os.path.join(work_dir, "fingerprint")
    resume = False
    if os.path.isfile(fp_file) and open(fp_file).read().strip() == newfp and not a.force:
        resume = True
        log("Resuming: --targetList / --chunk-size unchanged since previous run.")
    if not resume:
        # Clear the whole work area + per-piece results on a changed split (or --force).
        shutil.rmtree(work_dir, ignore_errors=True)
        shutil.rmtree(pieces_root, ignore_errors=True)
        os.makedirs(work_dir, exist_ok=True)
        os.makedirs(pieces_root, exist_ok=True)
        with open(fp_file, "w") as f:
            f.write(newfp + "\n")
        # remove stale per-piece struct subdirs from a previous split
        if struct_root:
            for name in os.listdir(struct_root):
                if name.startswith("piece."):
                    shutil.rmtree(os.path.join(struct_root, name), ignore_errors=True)
        if a.force:
            log("Forcing reprocess of all pieces.")

    npieces = write_piece_targets(a.targetList, total, a.chunk_size, work_dir)
    log("Searching %d targets against query %s (pieces=%d, per-piece<=%d, "
        "njobs=%d cap, rmsdCut=%s)" % (total, a.query, npieces,
        a.chunk_size, a.njobs, a.rmsdCut))

    # ---- run pieces (njobs at once) with resume + tqdm progress ----
    failed = 0
    pbar = tqdm(total=npieces, desc="searching pieces", unit="piece", leave=False) if tqdm else None
    todo = []
    for i in range(1, npieces + 1):
        if os.path.isfile(os.path.join(work_dir, "done.%d" % i)):
            log("  piece %d: already done (resume skip)" % i)
            if pbar:
                pbar.update(1)
        else:
            todo.append(i)
    if todo:
        with ProcessPoolExecutor(max_workers=a.njobs) as ex:
            futs = {ex.submit(run_piece, a, i, work_dir, pieces_root, struct_root): i
                    for i in todo}
            for fut in as_completed(futs):
                i = futs[fut]
                try:
                    ok, msg = fut.result()
                except Exception as e:  # noqa: BLE001
                    ok, msg = False, "unexpected error: %r" % e
                if ok:
                    log("  piece %d: OK (%s)" % (i, msg))
                else:
                    log("  piece %d FAILED: %s" % (i, msg))
                    failed += 1
                if pbar:
                    pbar.update(1)
    if pbar:
        pbar.close()

    # ---- merge finished pieces ----
    merge_from_done(a, npieces, work_dir, pieces_root, match_out, a.seqOut)

    if failed > 0:
        log("WARNING: %d pieces failed. Re-run the same command to resume "
            "(only unfinished/failed pieces reprocess)." % failed)
        return 2
    log("All %d pieces done." % npieces)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

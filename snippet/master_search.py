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
    match.txt                  merged match-address file, global top-N, RMSD-sorted
    seqs.txt                   merged match sequences, aligned to match.txt (if --seqOut)
    structs/                   final --structOut dir (default OUT/structs)
      match<n>.pdb             globally ranked n-th match (match1.pdb = overall best)
    pieces/
      piece.<i>/match.txt      per-piece match (result, NOT hidden in .chunks)
      piece.<i>/seq.txt        per-piece seq  (if --seqOut)
      piece.<i>/structs/       per-piece match structures match1.pdb, ... (no
                               cross-piece overwrite, numbered within the piece)
    .chunks/                   WORK/metadata only: fingerprint, targets.<i>,
                               run.<i>.log, done.<i> sentinels. No result files.

Bugs fixed vs. the bash original:
  1. Structures no longer overwrite each other. `master` names every piece's
     output `match1.pdb, match2.pdb, ...` (numbered from 1, see Search.cpp
     renameStruct). Running pieces concurrently into ONE shared --structOut dir
     made piece N clobber piece M's `matchN.pdb`. Fix: each piece writes into its
     own pieces/piece.<i>/structs/ subdir, and the final structs/match<n>.pdb are
     renumbered by GLOBAL rank at merge time.
  2. Per-piece match files no longer land in the hidden .chunks work dir. They
     now live under OUT/pieces/piece.<i>/ (real results), while .chunks holds
     only resume metadata (fingerprint, sentinels, logs).

Usage (two steps: search then merge):

  master_search.py search --query <q.pds> --targetList <list> --rmsdCut <X> \
                          --out <dir> [search options]
  master_search.py merge --out <dir> [merge options]

`search` builds the piece lists, runs one `master` job per piece (parallel,
resumable via .done sentinels) and writes the per-piece results under
OUT/pieces/. `merge` combines the finished pieces into the final OUT/match.txt,
OUT/seqs.txt and OUT/structs/ (RMSD-sorted, global top-N, optional
--no-breakage). `merge` can be re-run any number of times on the same OUT (it
only reads the finished pieces), so you can change --merged-topN / --no-breakage
without re-searching.

Search options (required: --query --targetList --rmsdCut --out):
  --chunk-size <N>          structures per piece (default 10000); sets the resume granularity
  --njobs <N>               how many pieces to search concurrently (default 8)
  --bbRMSD                  full-backbone RMSD instead of CA-only
  --topN <N>                keep best N matches per piece (forwarded to master; 0 = no limit)
  --minN <N>                return at least N best matches (per piece)
  --gapLen <s>              gap restraints, e.g. '1-10;0-3'
  --outType <t>             match|full|wgap (region written by structOut)
  --seqOut <file>           also write per-piece match sequences (for the later merge)
  --no-structs              do not write match structures
  --master <path>           master binary (default: found via PATH)
  --force                   reprocess all pieces even if done

Merge options (required: --out):
  --merged-topN <N>         final top-N kept after merge + --no-breakage (default: same as --topN)
  --topN <N>                default for --merged-topN when not given
  --no-breakage             drop matches whose adjacent segments are not really
                            peptide-bonded (C-to-N distance of boundary residues >= 1.4 A)
  --seqOut <file>           also write match sequences, merged to this file
  --no-structs              do not write match structures
  --structOut <dir>         match structure dir  (default: $out/structs)
  --matchOut <file>         final merged match-address file (default: $out/match.txt)

Notes:
  --topN is a PER-PIECE cap: it is forwarded to each master piece so a piece
  never carries more than N candidates (keeping per-piece files small). The final
  output keeps --merged-topN matches (default: --topN) after the merge, where all
  pieces' matches are sorted by RMSD ascending, optionally filtered by
  --no-breakage, then truncated. Per-piece topN never loses a global-top candidate:
  a piece that discards a match m has N matches with RMSD <= m, so m cannot rank
  in the global top-N (so keep --topN >= --merged-topN, or set --topN 0).
  --minN stays per piece (it widens the cutoff so a piece returns at least N).
  The merged match.txt (and seqs.txt, aligned to it) is RMSD-sorted.

  --no-breakage requires the per-piece match structures (i.e. not --no-structs).
  It measures, for each adjacent pair of query segments, the distance between the
  C of seg_i's last matched residue and the N of seg_{i+1}'s first matched residue
  in the target structure, and drops any match where that distance is >= 1.4 A
  (i.e. the two segments are NOT directly peptide-bonded — a real gap). NOTE this
  keeps only DIRECTLY-bonded (contiguous) matches: if a loop residue sits between
  the two segments, that boundary C-N is not a direct bond and the match is
  dropped.  Adding an output flag (e.g. --seqOut) after pieces already ran does
  NOT invalidate sentinels -> done pieces emit nothing for it; use --force to rebuild.

Deps: master (MASTER, grigoryanlab.org, e.g. /usr/local/bin); tqdm (optional).
"""
import argparse
import hashlib
import math
import os
import re
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


def build_parser():
    p = argparse.ArgumentParser(
        prog="master_search.py",
        description="Parallel, resumable MASTER structural-motif search (query side only). "
                    "Two steps: `search` runs the piece-wise search, `merge` merges results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", metavar="{search,merge}")
    sub.required = True

    # ---- search: build pieces + run master (no merge) ----
    ps = sub.add_parser("search",
                        help="run the piece-wise MASTER search (no merge)",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ps.add_argument("--query", required=True, help="query structure already converted to PDS")
    ps.add_argument("--targetList", required=True, help="file, one target .pds path per line")
    ps.add_argument("--rmsdCut", required=True, help="RMSD cutoff (Angstrom) for a match")
    ps.add_argument("--out", required=True, help="output root")
    ps.add_argument("--chunk-size", dest="chunk_size", type=int, default=10000,
                    help="structures per piece; the resume granularity")
    ps.add_argument("--njobs", type=int, default=8,
                    help="pieces searched concurrently (cap only, not split size)")
    ps.add_argument("--master", default=None, help="master binary (default: via PATH)")
    ps.add_argument("--topN", type=int, default=100,
                    help="keep best N matches per piece (forwarded to master; 0 = no limit)")
    ps.add_argument("--minN", type=int, default=None, help="return at least N best matches (per piece)")
    ps.add_argument("--gapLen", default=None, help="gap restraints, e.g. '1-10;0-3'")
    ps.add_argument("--outType", default=None, help="match|full|wgap (region written by structOut)")
    ps.add_argument("--seqOut", default=None,
                    help="also write per-piece match sequences (for the later merge)")
    ps.add_argument("--bbRMSD", action="store_true", help="full-backbone RMSD instead of CA-only")
    ps.add_argument("--no-structs", action="store_true", help="do not write match structures")
    ps.add_argument("--force", action="store_true", help="reprocess all pieces even if done")
    ps.set_defaults(func=cmd_search)

    # ---- merge: combine finished pieces into final outputs ----
    pm = sub.add_parser("merge",
                        help="merge finished piece results into final outputs",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pm.add_argument("--out", required=True, help="output root written by `search`")
    pm.add_argument("--topN", type=int, default=100,
                    help="default for --merged-topN when it is not given")
    pm.add_argument("--merged-topN", type=int, default=None,
                    help="final top-N kept after merge (and --no-breakage). Default: same as --topN")
    pm.add_argument("--matchOut", default=None,
                    help="final merged match-address file (default: $out/match.txt)")
    pm.add_argument("--seqOut", default=None,
                    help="also write match sequences, merged to this file")
    pm.add_argument("--structOut", default=None,
                    help="match structure dir (default: $out/structs)")
    pm.add_argument("--no-structs", action="store_true", help="do not write match structures")
    pm.add_argument("--no-breakage", action="store_true",
                    help="drop matches whose adjacent segments are not really peptide-bonded "
                         "(C-to-N distance >= 1.4 A)")
    pm.set_defaults(func=cmd_merge)

    return p


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


def run_pieces(a, npieces, work_dir, pieces_root):
    """Run unfinished pieces (njobs at once) with resume + tqdm. Returns failed count."""
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
            futs = {ex.submit(run_piece, a, i, work_dir, pieces_root): i for i in todo}
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
    return failed


def cmd_search(a):
    """`search`: build pieces, run master on each, write per-piece results + .done."""
    global _LOG
    validate(a)
    a.out = os.path.abspath(a.out)
    work_dir = os.path.join(a.out, ".chunks")
    pieces_root = os.path.join(a.out, "pieces")
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(pieces_root, exist_ok=True)
    _LOG = os.path.join(a.out, "master_search.log")

    total = count_lines(a.targetList)
    if total == 0:
        log("targetList is empty: nothing to search.")
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
        if a.force:
            log("Forcing reprocess of all pieces.")

    npieces = write_piece_targets(a.targetList, total, a.chunk_size, work_dir)
    log("Searching %d targets against query %s (pieces=%d, per-piece<=%d, "
        "njobs=%d cap, rmsdCut=%s)" % (total, a.query, npieces,
        a.chunk_size, a.njobs, a.rmsdCut))

    failed = run_pieces(a, npieces, work_dir, pieces_root)
    if failed > 0:
        log("WARNING: %d pieces failed. Re-run `search` to resume "
            "(only unfinished/failed pieces reprocess)." % failed)
        return 2
    log("All %d pieces done. Run `master_search.py merge --out %s` to merge results."
        % (npieces, a.out))
    return 0


def cmd_merge(a):
    """`merge`: combine finished pieces into final match.txt / seqs.txt / structs/."""
    global _LOG
    a.out = os.path.abspath(a.out)
    work_dir = os.path.join(a.out, ".chunks")
    pieces_root = os.path.join(a.out, "pieces")
    match_out = a.matchOut or os.path.join(a.out, "match.txt")
    struct_root = None if a.no_structs else (a.structOut or os.path.join(a.out, "structs"))
    os.makedirs(a.out, exist_ok=True)
    if struct_root:
        os.makedirs(struct_root, exist_ok=True)
    _LOG = os.path.join(a.out, "master_search.log")

    npieces = 0
    if os.path.isdir(work_dir):
        for name in os.listdir(work_dir):
            if name.startswith("done."):
                npieces = max(npieces, int(name.split(".")[1]))
    if npieces == 0:
        log("No finished pieces under %s. Run `master_search.py search --out %s` first."
            % (a.out, a.out))
        return 1

    if a.no_breakage:
        any_struct = False
        if os.path.isdir(pieces_root):
            for i in range(1, npieces + 1):
                sd = os.path.join(pieces_root, "piece.%d" % i, "structs")
                if os.path.isdir(sd) and os.listdir(sd):
                    any_struct = True
                    break
        if not any_struct:
            log("--no-breakage needs per-piece match structures, but none were found. "
                "Re-run `search` without --no-structs.")
            return 1

    merge_from_done(a, npieces, work_dir, pieces_root, match_out, a.seqOut, struct_root)
    return 0


def main(argv):
    p = build_parser()
    args = p.parse_args(argv)
    return args.func(args)


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


def run_piece(a, i, work_dir, pieces_root):
    """Search piece i. Returns (ok, msg). Writes .done sentinel on success."""
    piece_dir = os.path.join(pieces_root, "piece.%d" % i)
    os.makedirs(piece_dir, exist_ok=True)
    piece_list = os.path.join(work_dir, "targets.%d" % i)
    match_out = os.path.join(piece_dir, "match.txt")
    seq_out = os.path.join(piece_dir, "seq.txt") if a.seqOut else None
    struct_dir = None
    if not a.no_structs:
        # per-piece structs live INSIDE the piece dir: pieces/piece.<i>/structs/
        struct_dir = os.path.join(piece_dir, "structs")
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


PEPTIDE_BOND = 1.4  # Angstrom: C-to-N peptide bond is ~1.33 A


def _parse_segs(match_line):
    """Return [(start,end), ...] residue ranges of each segment from a match line."""
    return [(int(a), int(b)) for a, b in re.findall(r"\((\d+),(\d+)\)", match_line)]


def _residue_list(pdb_path):
    """Parse a struct PDB into an ordered list of [chain, resSeq, {atomname:(x,y,z)}]."""
    residues = []
    cur = None
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        name = line[12:16].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        if cur is None or cur[0] != chain or cur[1] != resseq:
            cur = [chain, resseq, {}]
            residues.append(cur)
        cur[2][name] = (x, y, z)
    return residues


def _dist(p, q):
    return math.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2)


def _has_breakage(struct_path, segs):
    """True if any two adjacent segments are not really peptide-bonded.

    For each adjacent pair (seg_i, seg_{i+1}) measure the distance between the C
    of seg_i's last matched residue and the N of seg_{i+1}'s first matched residue
    (coordinates come from the per-piece struct file). If that C-N distance is
    >= PEPTIDE_BOND the two segments are not directly bonded (a real gap: loop
    residue in between, chain break, or cross-chain) -> breakage.
    """
    if not os.path.isfile(struct_path):
        return True
    residues = _residue_list(struct_path)
    if not residues:
        return True
    n = len(residues)
    lengths = [b - a + 1 for a, b in segs]
    gaps = [segs[i + 1][0] - segs[i][1] - 1 for i in range(len(segs) - 1)]
    # If the struct file holds the segments plus their gap residues (outType wgap)
    # its residue count is sum(lengths)+sum(gaps); otherwise (outType match) it is
    # just sum(lengths) and the gap residues are absent. Coordinates are real in
    # both cases, so only the position arithmetic changes.
    if sum(lengths) + sum(gaps) == n:
        eff_gaps = gaps
    else:
        eff_gaps = [0] * len(gaps)
    cum = 0
    for i in range(len(segs) - 1):
        last_idx = cum + lengths[i]
        first_idx = cum + lengths[i] + eff_gaps[i] + 1
        if last_idx > n or first_idx > n:
            return True
        c = residues[last_idx - 1][2].get("C")
        nn = residues[first_idx - 1][2].get("N")
        if c is None or nn is None:
            return True
        if _dist(c, nn) >= PEPTIDE_BOND:
            return True
        cum += lengths[i] + eff_gaps[i]
    return False


def merge_from_done(a, npieces, work_dir, pieces_root, match_out, seq_out, struct_root):
    """Merge done pieces into the global top-N, RMSD-sorted.

    For each done piece, collect every (match line, seq line) and remember which
    piece + local index it came from (local index maps 1:1 to that piece's
    per-piece struct file match<local_index+1>.pdb). Sort all rows by the match
    line's RMSD ascending (stable on ties) and truncate to --topN if set, then:

      - match.txt      <- merged match lines (top-N, RMSD-sorted)
      - seqs.txt       <- merged seq lines, aligned to the match lines
      - structs/<n>.pdb <- per-piece struct file renamed by GLOBAL rank n, so
                           structs/match1.pdb is the overall best match.
    """
    rows = []  # (rmsd, order, match_line, seq_line_or_None, piece_i, local_idx)
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
                for lidx, ml in enumerate(mf):
                    rows.append((_rmsd(ml), order, ml, sf.readline(), i, lidx))
                    order += 1
        else:
            with open(pm) as mf:
                for lidx, ml in enumerate(mf):
                    rows.append((_rmsd(ml), order, ml, None, i, lidx))
                    order += 1
    rows.sort(key=lambda t: (t[0], t[1]))
    if a.no_breakage:
        before = len(rows)
        kept = []
        for row in rows:
            segs = _parse_segs(row[2])
            spath = os.path.join(pieces_root, "piece.%d" % row[4], "structs",
                                 "match%d.pdb" % (row[5] + 1))
            if not _has_breakage(spath, segs):
                kept.append(row)
        rows = kept
        log("  no-breakage: kept %d / %d (dropped %d with real gaps)" %
            (len(rows), before, before - len(rows)))
    n_final = a.merged_topN if a.merged_topN is not None else a.topN
    if n_final and len(rows) > n_final:
        rows = rows[:n_final]

    with open(match_out, "w") as mf:
        for _r, _o, ml, _sl, _i, _l in rows:
            mf.write(ml)
    log("  merged %d matches (RMSD-sorted%s) -> %s" % (
        len(rows), ", top %d" % n_final if n_final else "", match_out))
    if seq_out:
        with open(seq_out, "w") as sf:
            for _r, _o, _ml, sl, _i, _l in rows:
                if sl is not None:
                    sf.write(sl)
        log("  merged %d seqs (aligned to matches) -> %s" % (len(rows), seq_out))

    if struct_root:
        # renumber per-piece structs by GLOBAL rank
        for name in os.listdir(struct_root):
            if name.startswith("match") and name.endswith(".pdb"):
                os.remove(os.path.join(struct_root, name))
        n_struct = 0
        for rank, (_r, _o, _ml, _sl, pi, li) in enumerate(rows, start=1):
            src = os.path.join(pieces_root, "piece.%d" % pi, "structs",
                               "match%d.pdb" % (li + 1))
            if os.path.isfile(src):
                shutil.copyfile(src, os.path.join(struct_root, "match%d.pdb" % rank))
                n_struct += 1
        log("  merged %d structs (global rank) -> %s" % (n_struct, struct_root))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

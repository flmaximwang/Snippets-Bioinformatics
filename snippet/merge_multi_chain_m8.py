#!/usr/bin/env python3
"""merge_multi_chain_m8.py — merge a ColabFold unpaired-mode multi-chain
pdb70.m8 into a single-chain coordinate system.

ColabFold unpaired mode searches each chain of a multi-chain query as a
separate MMseqs2 query (ids 101, 102, ...). The resulting pdb70.m8 therefore
carries qstart/qend coordinates relative to each fragment chain. When the
fragments are really segments of one chain (split only to improve per-segment
template coverage), the hits must be renumbered onto the merged chain and
re-sorted so that a downstream top-N template picker finds hits for every
fragment.

The per-chain query id and sequence length are read from the a3m files: the
first header line of each a3m gives the query id, the first block gives the
query sequence. The a3m arguments are order-sensitive: they define the chain
order in the merged chain.

Renumbering:  new qstart/qend = old qstart/qend + cumulative length of the
preceding chains.  tstart/tend, evalue, bits and cigar are left untouched.
The query id column of every row is rewritten to --query-name.

Sorting: within each chain, hits are sorted by bits descending (e-value
ascending, then qstart ascending as tie-breaks); the chains are then
interleaved round-robin by rank (chain0-best, chain1-best, ..., chain0-2nd,
chain1-2nd, ...) so that any prefix of the output contains hits from every
chain.  Chains without hits (e.g. a fragment with no template hit) drop out
of the interleaving automatically.

Usage:
    merge_multi_chain_m8.py \\
        --m8 unpaired/pdb70.m8 \\
        --a3m unpaired/unpaired_0.a3m \\
        --a3m unpaired/unpaired_1.a3m \\
        --a3m unpaired/unpaired_2.a3m \\
        --query-name v8_0_0_split \\
        -o unpaired/pdb70.merged.m8

Pure stdlib — runs with any python3.
"""

import argparse
import sys
from collections import OrderedDict

# MMseqs2 convertalis m8 column indexes (0-based), format:
# query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits cigar
QID, TARGET, FIDENT, ALNLEN, MISMATCH, GAPOPEN = 0, 1, 2, 3, 4, 5
QSTART, QEND, TSTART, TEND = 6, 7, 8, 9
EVALUE, BITS, CIGAR = 10, 11, 12
N_COLS = 13


def read_query(path):
    """Return (query_id, length) from the first block of a ColabFold a3m.

    Only the first '>' block is used: it holds the query sequence.  Everything
    after it (database hits, or duplicated query blocks) is ignored, so the
    length is that of the query chain alone.
    """
    qid = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if qid is not None:
                    break  # first block complete
                qid = line[1:].split()[0]
            elif qid is not None:
                seq.append(line.strip())
    if qid is None:
        raise ValueError(f"{path}: no query header found")
    return qid, len("".join(seq).replace("-", ""))


def parse_m8(path):
    """Return (rows, skipped): rows are lists of the 13 m8 fields."""
    rows, skipped = [], 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < N_COLS:
                skipped += 1
                sys.stderr.write(f"WARN: skipping malformed m8 line: {line}\n")
                continue
            try:
                int(cols[QSTART])
                int(cols[QEND])
                float(cols[EVALUE])
                float(cols[BITS])
            except ValueError:
                skipped += 1
                sys.stderr.write(f"WARN: skipping non-numeric m8 line: {line}\n")
                continue
            rows.append(cols)
    return rows, skipped


def main():
    ap = argparse.ArgumentParser(
        description="Merge ColabFold unpaired multi-chain pdb70.m8 into "
        "single-chain coordinates, renumbered and rank-interleaved."
    )
    ap.add_argument("--m8", required=True, help="unpaired-mode pdb70.m8 (hits of all chains)")
    ap.add_argument(
        "--a3m",
        required=True,
        action="append",
        metavar="FILE",
        help="per-chain a3m, in merged-chain order (repeatable)",
    )
    ap.add_argument("-o", "--out", required=True, help="output merged m8 path")
    ap.add_argument(
        "--query-name",
        default="merged",
        help="query id written for every merged row (default: merged)",
    )
    args = ap.parse_args()

    # chain order = a3m argument order; offset = cumulative preceding length
    chains = []  # (qid, length, offset)
    offset = 0
    for path in args.a3m:
        qid, length = read_query(path)
        chains.append((qid, length, offset))
        offset += length
    id_to_chain = {qid: (i, off) for i, (qid, _, off) in enumerate(chains)}
    total = offset

    rows, skipped = parse_m8(args.m8)
    groups = OrderedDict((i, []) for i in range(len(chains)))
    unmatched = 0
    for cols in rows:
        chain = id_to_chain.get(cols[QID])
        if chain is None:
            unmatched += 1
            continue
        i, off = chain
        cols[QID] = args.query_name
        cols[QSTART] = str(int(cols[QSTART]) + off)
        cols[QEND] = str(int(cols[QEND]) + off)
        groups[i].append(cols)
    if unmatched:
        sys.stderr.write(
            f"WARN: dropped {unmatched} rows whose query id matches no a3m\n"
        )

    # within each chain: bits desc, evalue asc, qstart asc (stable)
    for i in groups:
        groups[i].sort(
            key=lambda r: (-float(r[BITS]), float(r[EVALUE]), int(r[QSTART]))
        )

    # round-robin interleave by rank across chains (in a3m argument order)
    out_rows = []
    rank = 0
    while True:
        progressed = False
        for i in groups:
            if rank < len(groups[i]):
                out_rows.append(groups[i][rank])
                progressed = True
        if not progressed:
            break
        rank += 1

    with open(args.out, "w") as f:
        for cols in out_rows:
            f.write("\t".join(cols) + "\n")

    print(
        f"chains: {[(qid, length) for qid, length, _ in chains]}  "
        f"merged length: {total}\n"
        f"rows: input {len(rows)}  written {len(out_rows)}  "
        f"skipped {skipped}  unmatched {unmatched}\n"
        f"per chain: "
        + ", ".join(f"{qid}:{len(groups[i])}" for i, (qid, _, _) in enumerate(chains))
        + f"\noutput: {args.out}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

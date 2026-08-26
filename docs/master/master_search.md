# master_search.py — MASTER structural-motif query runner

`master_search.py` runs the MASTER search binary (`master`) against a built PDS
database, **query-side only**. It does not build the database (that's
`snippet/create_master_db.sh`) and does not convert PDB→PDS (that's `createPDS`).
It takes a query `.pds` and a `--targetList`, divides the list into pieces, and
runs one `master --targetList <piece>` job per piece — in parallel, with
per-piece resume. Rewritten from `master_search.sh` in Python (tqdm progress bar);
the CLI is identical.

Two design facts that drive every parameter below:

1. **The split is the resume base.** `--chunk-size` chops the
   `--targetList` into consecutive pieces. Each piece gets a `.done.<i>` sentinel
   after it finishes; re-running the same command reuses every finished piece and
   only searches the unfinished ones. Nothing about the split itself is
   re-computed on a resume once the sentinel exists.
2. **`--njobs` is only a concurrency cap.** It controls *how many* pieces run at
   the same time (to use all cores), never *how big* a piece is.

Files produced under `--out`:

```
out/
  match.txt          merged match-address file — global top-N (default 100), RMSD-sorted
  seqs.txt           merged match sequences, aligned to match.txt (only if --seqOut)
  structs/           globally ranked match structures (unless --no-structs)
    match<n>.pdb     n-th best match overall (match1.pdb = overall best)
  pieces/
    piece.<i>/match.txt    piece i's matches       (before merge)
    piece.<i>/seq.txt      piece i's sequences     (before merge, if --seqOut)
    piece.<i>/structs/     piece i's own structures: match1.pdb, match2.pdb, ...
  master_search.log  run log
  .chunks/           work/metadata only — NO result files
    targets.<i>      piece i's slice of the targetList
    run.<i>.log      piece i's master stdout/stderr
    done.<i>         sentinel: piece i finished
    fingerprint      hash of (targetList + --chunk-size)
```

`master` names every piece's output structures `match1.pdb, match2.pdb, ...`
(numbered from 1, regardless of the target — see `Search.cpp` `renameStruct`).
Writing all pieces into ONE flat `--structOut` dir would make concurrent pieces
clobber each other's `matchN.pdb`, so each piece writes into its own
`pieces/piece.<i>/structs/` subdir. At merge time those per-piece files are
renumbered by GLOBAL rank into the flat `structs/match<n>.pdb`, so
`structs/match1.pdb` is the overall best match (rank n matches `match.txt` line n).

---

## Full usage

```
python3 master_search.py \
  --query <q.pds> --targetList <list> --rmsdCut <X> --out <dir> \
  [--chunk-size N] [--njobs N] \
  [--bbRMSD] [--topN N] [--minN N] [--gapLen S] [--outType T] \
  [--seqOut FILE] [--no-structs] [--structOut DIR] [--matchOut FILE] \
  [--master PATH] [--force]
```

**Prerequisite — the query must already be PDS.** `master` cannot read a raw PDB.
Convert before calling:

```
createPDS --type query --pdb in.pdb --pds in.pds \
  --dCut 25.0 --dStep 5.0 --phiStep 10.0 --psiStep 10.0
```

`dCut`/`dStep`/`phiStep`/`psiStep` must match the values used when the database
was built (see `create_master_db.sh`). If you pass a non-`.pds` `--query`, the
script errors out and prints this hint.

---

## Required arguments

### `--query <file>`

The query structure, already converted to PDS (`createPDS --type query`).

- **Input:** `/path/to/query.pds` (from `createPDS --type query --pdb ...`).
- **Effect:** This exact file is passed unchanged as `master --query` to every
  piece. The script only checks that it exists and ends in `.pds`.
- **Error if:** the file is missing → `query PDS not found`; the path does not
  end in `.pds` → `query must be a .pds file (run createPDS --type query first)`.

```
--query query/1akha.pds
```

### `--targetList <file>`

The database: a text file, **one absolute/relative `.pds` path per line**. This
is the file that gets split into pieces.

- **Input:** e.g. `target_list.txt` produced by `create_master_db.sh`'s
  `gen_target_list` step, or any hand-written list of target PDS files.
- **Effect:** split into `ceil(lines / --chunk-size)` consecutive slices
  (preserving the file's order — it is never sorted). Each slice becomes one
  piece's `master --targetList`.
- **Error if:** missing → `targetList not found`. Empty → logs
  `nothing to search` and exits 0 with an empty `match.txt`.

```
--targetList /mnt/data/public/PDB_db/target_list.txt
```

### `--rmsdCut <x>`

RMSD cutoff in Ångstrom; a target "matches" only if its backbone RMSD to the
query is at or below this value.

- **Input:** a non-negative float, e.g. `3.0`.
- **Effect:** passed as `master --rmsdCut`. Smaller → stricter (fewer matches,
  slower cutoff but same completeness); larger → looser.
- **Example effect:** searching the homeobox `1akha` query against a database
  containing `1b72a` at `--rmsdCut 3.0` yields one match at RMSD 0.619:

```
--rmsdCut 3.0
# out/match.txt:
 0.61906 targets/1b72a.pds [(13,61)]
```

### `--out <dir>`

Root directory for all outputs (merged results and the `.chunks/` work area).

- **Input:** a (possibly non-existent) directory path, e.g. `--out searches/a`.
- **Effect:** created if needed; every derived path is canonicalized to absolute.
  Merge targets default to `$out/match.txt`, `$out/seqs.txt`, `$out/structs/`.
- **Resume meaning:** re-running with the **same** `--out` (and unchanged list +
  per-list) is what makes the `.done` sentinels reusable.

```
--out searches/a
```

---

## Splitting & concurrency

### `--chunk-size <N>`

**Structures per piece.** The full `--targetList` is divided into consecutive
pieces of up to N entries each (the last piece may be smaller). This is the
**resume granularity**: a whole piece is either done or not.

- **Input:** a positive integer. Default `10000`.
- **Effect on piece count:** `number_of_pieces = ceil(lines / N)`.

```
# 5 targets (1azw,1b72a,1ju3,1l7a,1pq5), --chunk-size 2
-> pieces=3 : [1azw,1b72a] [1ju3,1l7a] [1pq5]
```

- **Effect on resume:** smaller N → more, smaller pieces → finer resume (a crash
  loses less work) but more pieces to manage. Larger N → fewer, bigger pieces →
  coarser resume, less merge overhead.
- **On a full archive** the default 10000 keeps per-piece memory/log size tame
  while preserving meaningful resume granularity.

### `--njobs <N>`

**Concurrency cap** — how many pieces are searched at the same time. It does NOT
change the split.

- **Input:** a positive integer. Default `8`.
- **Effect:** at most N `master` processes run concurrently (rolling pool: when N
  are in flight it waits for the oldest before launching the next). `master` is
  single-core per process, so set N ≈ your core count to use the machine.
- **Not part of the fingerprint:** changing `--njobs` alone never forces a redo —
  the script only re-looks at pieces that aren't `.done` yet.

```
--njobs 8      # 8 pieces in-flight on an 8+ core box
--njobs 1      # fully serial (equivalent to one piece at a time)
```

---

## `master` search-criteria options (passed through per piece)

These are forwarded verbatim to each `master --targetList <piece>` call. Most
apply **per piece**; the exception is `--topN`, which is also re-applied globally
at merge time (see below).

### `--bbRMSD`

Flag. Switch the metric from the default CA-only RMSD to full-backbone RMSD
(uses N, CA, C, O).

- **Input:** just the presence of the flag.
- **Effect:** `master` gets `--bbRMSD`; matches are scored/ranked on all four
  backbone atoms instead of CA alone. Use when backbone geometry beyond the Cα
  trace matters.

```
--bbRMSD
```

### `--topN <N>`

Keep only the best N matches (by the search metric). It is **global**: applied
per piece (so a piece never carries more than N candidates, keeping per-piece
files small) **and** re-applied after the merge, where all pieces' matches are
sorted by RMSD ascending and truncated to the overall best N.

- **Input:** a non-negative integer (default `100`; `0` = no limit).
- **Effect:** forwarded as `master --topN` (also speeds search: it lets `master`
  lower the effective threshold once N matches are found), then re-applied
  globally at merge time. The merged `match.txt` is RMSD-sorted.
- **Why both are safe (exactness):** a global top-N match is always inside some
  piece's per-piece top-N — a piece that discards a match `m` has N matches with
  RMSD ≤ `m`, so `m` can never rank in the global top-N. Per-piece trimming only
  shrinks the candidate set, never drops a global-top candidate.

```sh
--rmsdCut 3.0 --topN 5000     # overall best 5000 across the whole database
```

### `--minN <N>`

Return at least the N best matches (by the metric) regardless of `--rmsdCut`.
`--minN` stays **per piece** (it widens the cutoff so a piece returns at least N
matches); it does not change the final global top-N.

- **Input:** a non-negative integer; if both `--minN` and `--topN` are given,
  `--minN` must be ≤ `--topN`. Default `0`.
- **Effect:** passed as `master --minN`. Guarantees a minimum number of hits even
  when few structures fall within the cutoff.

```sh
--rmsdCut 1.5 --minN 10       # at least 10 best matches, even above 1.5 A
```

### `--gapLen <s>`

Constrain the length of sequence that maps between adjacent query segments (only
relevant for multi-segment queries).

- **Input:** one `min-max` per gap, semicolon-separated. A query with `k` disjoint
  segments has `k-1` gaps. A single `min` (e.g. `3`) means exactly that length.
- **Effect:** passed as `master --gapLen`. Example: a three-segment query where
  gap 1 (between segments 1–2) must be 1–10 residues and gap 2 (2–3) must be 0–3:

```
--gapLen '1-10;0-3'
```

Important: with `--gapLen`, MASTER assumes the segment order in the matches is the
same as the order in the query PDB. Without it, MASTER considers all orderings.

---

## Output options

### `--matchOut <file>`

Final **merged** match-address file (one line per match).

- **Input:** an output path. Default `$out/match.txt`.
- **Effect:** after all pieces finish (or resume), each `.done` piece's
  `pieces/piece.<i>/match.txt` is concatenated in piece order into this file.
  Every line records a match: RMSD, target entry, and aligned residues.
- **Example line:**
  ```
   0.61906 targets/1b72a.pds [(13,61)]
  ```
  RMSD in Å · the target PDS path · the aligned residue interval. This file is
  also `master --matchIn`-compatible for reprocessing a previous search.

### `--seqOut <file>`

Also write match sequences, **merged** to this file. Off by default.

- **Input:** an output path, e.g. `--seqOut out/seqs.txt`.
- **Effect:** each `.done` piece's `pieces/piece.<i>/seq.txt` is concatenated
  into the file (same order as matches). One line per match: RMSD then the
  sequence.
- **Example line:**
  ```
   0.61906 PHE THR THR ARG GLN LEU THR GLU LEU ...
  ```
- **Resume caveat:** per-piece `seq.txt` is only written on a run where
  `--seqOut` is active. If you first search without `--seqOut` (pieces now
  `.done`) and then re-run *with* `--seqOut`, those pieces are skipped and
  produce no `seq.txt`, so the merged seqs file comes back empty — add `--force`
  to regenerate.

### `--structOut <dir>`

Directory for match **structures** in PDB format (one PDB per match, written
optimally superimposed onto the query over the matching region).

- **Input:** an output dir. Default `$out/structs`. This is the **final, flat**
  merged-structures dir — per-piece structures are written to
  `pieces/piece.<i>/structs/` (see layout above).
- **Effect:** `master` numbers every piece's output `match1.pdb, match2.pdb, ...`
  starting from 1 (it does **not** name files after the target — see `Search.cpp`
  `renameStruct`). Writing all pieces into ONE flat dir would make concurrent
  pieces overwrite each other's `matchN.pdb`, so each piece writes into its own
  `pieces/piece.<i>/structs/` subdir. At merge time the per-piece files are
  renumbered by GLOBAL rank into `--structOut/match<n>.pdb` (rank n = the n-th
  lowest-RMSD match, matching `match.txt` line n). The region written is governed
  by `--outType`.
- **Example produced file:** `structs/match1.pdb` (the overall best match).

### `--outType <t>`

Which region of each match is written by `--structOut` (and affects `--seqOut`).

- **Input:** `match` (default), `full`, or `wgap`.
  - `match` — just the residues that align onto the query.
  - `full` — the entire database entry.
  - `wgap` — the matching region plus any residues mapping between query segments
    (only valid when `--gapLen` is given).
- **Effect:** passed as `master --outType`; the structs/seqs reflect that region.

```
--outType full        # save whole target structures, not just the aligned patch
```

### `--no-structs`

Flag. Disable writing match structures entirely (e.g. when you only want
`match.txt`).

- **Input:** just the presence of the flag; clearing it restores the default.
- **Effect:** `--structOut` is not passed to `master`; no `structs/` output.

```
--no-structs --rmsdCut 3.0    # fast: only the merged address/rank file
```

---

## Environment / misc

### `--master <path>`

Explicit path to the `master` binary (default: resolved via `PATH`).

- **Input:** a path to an executable, e.g. `--master /opt/master/1.5/bin/master`.
- **Effect:** used for every piece; errors out if missing or not executable.
  Without it the script uses `command -v master` (linstall adds it to
  `/usr/local/bin`, which is already on `PATH`).

### `--force`

Flag. Reprocess every piece even if `.done` sentinels exist.

- **Input:** just the presence of the flag.
- **Effect:** clears the whole `.chunks/` work area and the fingerprint, then
  re-runs every piece from scratch. Needed when piece output flags change after
  pieces are already done (e.g. adding `--seqOut`), or when you suspect a `.done`
  piece's outputs are corrupt.

```
--force --seqOut out/seqs.txt    # regenerate match sequences for all pieces
```

---

## Worked example

Build a 3-piece search and show the piece split:

```
# 11 targets -> --chunk-size 4 => 3 pieces: [4][4][3]
python3 snippet/master_search.py \
  --query query/1akha.pds --targetList db/target_list.txt --rmsdCut 3.0 \
  --out searches/a --njobs 3 --chunk-size 4 \
  --bbRMSD --topN 5000 --seqOut searches/a/seqs.txt
```

Log:

```
Searching 11 targets ... (pieces=3, per-piece<=4, njobs=3 cap, rmsdCut=3.0)
  piece 1: OK (2 matches)
  piece 2: OK (5 matches)
  piece 3: OK (1 match)
  merged 8 matches -> searches/a/match.txt
  merged seqs -> searches/a/seqs.txt
```

You are interrupted after pieces 1–2 finish. Rerun the **same command**:

```
Resuming: --targetList / --chunk-size unchanged since previous run.
  piece 1: already done (resume skip)
  piece 2: already done (resume skip)
  piece 3: OK (1 match)
  merged 8 matches -> searches/a/match.txt
```

Only piece 3 was searched again; 1–2 were reused, and the merge is recomputed from
the union of all `.done` pieces.

---

## Parameter reference table

| Parameter | Default | Type | Purpose |
|---|---|---|---|
| `--query` | — | required | query PDS (`createPDS --type query`) |
| `--targetList` | — | required | DB list of target PDS files, split into pieces |
| `--rmsdCut` | — | required | RMSD cutoff (Å) for a match |
| `--out` | — | required | output root |
| `--chunk-size` | `10000` | int | structures per piece = resume granularity |
| `--njobs` | `8` | int | max concurrent pieces (cap only) |
| `--bbRMSD` | off | flag | full-backbone RMSD metric |
| `--topN` | `100` | int | keep best N matches globally (0 = no limit) |
| `--minN` | `0` | int | return ≥N best matches (per piece) |
| `--gapLen` | — | str | gap-length restraints, `'min-max;...'` |
| `--outType` | `match` | str | `match`\|`full`\|`wgap` region for structs/seqs |
| `--seqOut` | off | file | merged match-sequence output |
| `--structOut` | `$out/structs` | dir | match-structure output dir |
| `--no-structs` | off | flag | skip writing match structures |
| `--matchOut` | `$out/match.txt` | file | merged match-address output |
| `--master` | `PATH` | path | explicit `master` binary |
| `--force` | off | flag | reprocess all pieces |

---

## Notes & limitations

- **`--topN`/`--minN` are per piece**, not global (see `--topN`).
- **`--seqOut` merge is empty on a pure resume** unless the original run already
  used `--seqOut` (see `--seqOut`).
- The script is deliberately **`master`-only**; it never calls `createPDS` or
  `parsePDS`. Convert your query with `createPDS` first.
- `--rmsdCut` completeness guarantee: MASTER's default `--rmsdMode 0` returns
  *all* matches within the cutoff; that guarantee is preserved because each piece
  is an independent provable search of a disjoint target subset.
- Verified against MASTER v1.5 binaries (`master`, `createPDS`) at
  `/usr/local/bin`.

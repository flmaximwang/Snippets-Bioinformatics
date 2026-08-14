#!/usr/bin/env python3
"""aa_to_dna.py — codon-optimize protein sequences with dnachisel.

Takes either a raw amino-acid string or a protein FASTA file as the positional
argument, translates each sequence to a coding DNA, then runs a dnachisel
DnaOptimizationProblem with:

  - EnforceTranslation(translation=<aa>)      keep the exact amino-acid sequence
  - AvoidPattern("<base>" * n) for ACGT       forbid homopolymer runs of n or
                                              more identical nucleotides
  - CodonOptimize(species=<species>)          maximize codon usage for the
                                              target organism (default: e_coli,
                                              i.e. codon-usage-table-driven CAI
                                              optimization, method
                                              "use_best_codon")

Homopolymer and k-mer constraints:
    AvoidPattern("<base>" * max_repeat) forbids homopolymer runs of
    max_repeat or more nucleotides; --max-repeat sets it directly (default
    7 — the value is passed as-is, no scan).
    UniquifyAllKmers(k) forbids repeated k-mers of length k (including
    reverse-complement homologies, dnachisel's default). Smaller k is
    harder to satisfy, so k is scanned upward from 7 and the FIRST k for
    which the design converges is kept (tightest achievable);
    --max-kmer-len is the scan upper bound (default 25, hard error if
    nothing works by then). With a fixed --seed the whole scan is
    deterministic.

GC-content constraint (optional):
    --gc-content LOWER_UPPER (fractions, e.g. 0.3_0.7) adds dnachisel's
    EnforceGCContent with the GC fraction in [LOWER, UPPER]; --gc-window N
    restricts the check to sliding N-bp windows (whole sequence when
    omitted). The initial sequence is pre-optimized with a GC-only dnachisel
    problem so the constraint only has to be kept during the final design.
    Without --gc-content the GC limit is off.

Usage:
    aa_to_dna.py MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFP -o opt.fa
    aa_to_dna.py proteins.fasta --species e_coli --seed 1 -o opt.fa
    aa_to_dna.py proteins.fasta --gc-content 0.35_0.65 --gc-window 50 -o opt.fa
    aa_to_dna.py proteins.fasta --max-repeat 30 --seed 7 --fmt csv -o opt.csv
    aa_to_dna.py MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFP --fmt txt > opt.txt

Output formats (--fmt / -o suffix): fasta (default; 60-char wrapped
">name desc [meta]" records), csv (header row "name,desc,dna", one record
per row), txt (a "name desc [meta]" line followed by the full unwrapped DNA
line, blank line between records). With -o the format is inferred from the
output suffix (.fa/.fasta/.fna/.fas, .csv, .txt — an unrecognized suffix
falls back to an explicit --fmt, otherwise errors); without -o, --fmt
decides (default fasta).

The output is DNA-only (no stop codon appended: the coding sequence
codes exactly for the input amino acids). Codon tables come from
biorazer.database.codon_usage (python-codon-tables). Requires: dnachisel,
numpy, biorazer (testing env: `mamba activate BioRazer`).
"""

import argparse
import csv
import os
import random
import re
import sys

import numpy as np
from dnachisel import (AvoidPattern, DnaOptimizationProblem, EnforceGCContent,
                       EnforceTranslation, NoSolutionError, UniquifyAllKmers)
from dnachisel.builtin_specifications import CodonOptimize
from biorazer.database.codon_usage import get_codon_usage_table_by_aa

START_KMER = 7          # k-mer uniqueness scan begins here (smaller = harder)
START_REPEAT = 7        # default homopolymer ceiling (fixed, not scanned)
LINE_WIDTH = 60         # FASTA sequence wrap
FORMAT_BY_SUFFIX = {    # -o suffix -> output format
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta", ".fas": "fasta",
    ".csv": "csv",
    ".txt": "txt",
}


def build_codon_tables(species):
    """Fetch the species' {aa: {codon: freq}} table from biorazer's
    database and derive: the initial codon per amino acid (highest
    frequency), the synonymous codon groups (used for the GC reachability
    diagnostic) and the valid amino-acid alphabet. CodonOptimize consumes
    the same underlying python-codon-tables data, so the initial sequence
    and the final objective agree on the codon repertoire."""
    try:
        by_aa = get_codon_usage_table_by_aa(species)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    init_codons = {aa: max(codons, key=codons.get)
                   for aa, codons in by_aa.items()}
    synonymous = {aa: sorted(codons) for aa, codons in by_aa.items()}
    return init_codons, synonymous, set(by_aa)


def read_input(arg):
    """Return [(name, desc, aa), ...] from a raw AA string or a FASTA file.

    A path to an existing file is parsed as FASTA; anything else is treated
    as a single amino-acid sequence (whitespace stripped, upper-cased).
    """
    if os.path.isfile(arg):
        records = []
        name = desc = None
        seq_lines = []
        with open(arg) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if name is not None:
                        records.append((name, desc, "".join(seq_lines)))
                    fields = line[1:].split(None, 1)
                    name = fields[0]
                    desc = fields[1] if len(fields) > 1 else ""
                    seq_lines = []
                elif line and not line.startswith(";"):
                    seq_lines.append(line.strip())
        if name is not None:
            records.append((name, desc, "".join(seq_lines)))
        if not records:
            sys.exit(f"ERROR: no FASTA record found in {arg!r}")
        return records
    aa = "".join(arg.split()).upper()
    return [("seq1", "", aa)]


def check_aa(aa, valid_aa):
    bad = sorted(set(aa) - valid_aa)
    if bad:
        sys.exit(
            f"ERROR: invalid amino-acid character(s): {bad!r} "
            f"(allowed: {''.join(sorted(valid_aa - {'*'}))} plus '*')"
        )
    if not aa:
        sys.exit("ERROR: empty amino-acid sequence")


def parse_gc_range(s):
    """'LOWER_UPPER' -> (lo, hi): GC fraction bounds in [0, 1], lo < hi."""
    parts = s.split("_")
    if len(parts) != 2:
        sys.exit(f"ERROR: --gc-content expects LOWER_UPPER (e.g. 0.3_0.7), "
                 f"got {s!r}")
    try:
        lo, hi = (float(p) for p in parts)
    except ValueError:
        sys.exit(f"ERROR: --gc-content values must be numbers, got {s!r}")
    if not (0.0 <= lo < hi <= 1.0):
        sys.exit(f"ERROR: --gc-content needs 0 <= LOWER < UPPER <= 1 "
                 f"(fractions, e.g. 0.3_0.7), got {s!r}")
    return lo, hi


def _gc_count(seq):
    """Number of G/C nucleotides."""
    return seq.count("G") + seq.count("C")


def _window_gc_bounds(aa, synonymous, s, e):
    """Optimistic local reachable GC fraction bounds for DNA window [s, e):
    min/max G/C count over the synonymous codons of the window's amino
    acids, weighted by how much of each codon overlaps the window."""
    lo = hi = 0
    for p in range(s // 3, (e + 2) // 3):
        if p >= len(aa):
            break
        ov = min(e, 3 * p + 3) - max(s, 3 * p)
        if ov <= 0:
            continue
        gcs = [_gc_count(c) for c in synonymous[aa[p]]]
        lo += min(gcs) * ov / 3
        hi += max(gcs) * ov / 3
    return lo / (e - s), hi / (e - s)


def gc_preoptimize(dna, aa, lo, hi, window, synonymous):
    """Push the GC content into [lo, hi] with dnachisel's own local search:
    a small problem whose only objective is EnforceGCContent(mini, maxi,
    window), constrained by EnforceTranslation. This is fast (location-wise
    exhaustive enumeration accepts improving variants directly, which the
    constraint solver does not) and uses the exact same sliding-window
    semantics as the later design problem. The objective is best-effort, so
    the result is verified window-by-window here; an unreachable target is
    diagnosed with the sequence's synonymous-codon GC range.
    """
    problem = DnaOptimizationProblem(
        sequence=dna,
        constraints=[EnforceTranslation(translation=aa)],
        objectives=[EnforceGCContent(mini=lo, maxi=hi, window=window)],
        logger=None,
    )
    problem.optimize()
    fixed = str(problem.sequence)
    n = len(fixed)
    if window is not None and window < n:
        spans = [(i, i + window) for i in range(0, n - window + 1)]
    else:
        spans = [(0, n)]
    offenders = [
        i for i, (s, e) in enumerate(spans)
        if not (lo <= _gc_count(fixed[s:e]) / (e - s) <= hi)
    ]
    if offenders:
        gmin = sum(min(_gc_count(c) for c in synonymous[a]) for a in aa)
        gmax = sum(max(_gc_count(c) for c in synonymous[a]) for a in aa)
        if window is None:
            sys.exit(
                f"ERROR: GC range [{lo}, {hi}] not reachable for this "
                f"amino-acid sequence: synonymous-codon GC range is "
                f"[{gmin / n:.3f}, {gmax / n:.3f}]"
            )
        # With windows, the global range is not the binding constraint: a
        # window's OWN amino-acid composition caps its GC. Report the worst
        # offending window's local reachable range so the user sees why.
        worst = None
        worst_gap = -1.0
        for idx in offenders:
            s, e = spans[idx]
            wlo, whi = _window_gc_bounds(aa, synonymous, s, e)
            gap = max(lo - whi, wlo - hi)  # > 0: target unreachable here
            if gap > worst_gap:
                worst_gap, worst = gap, (s, e, wlo, whi)
        assert worst is not None  # offenders is non-empty here
        ws, we, wlo, whi = worst
        sys.exit(
            f"ERROR: GC range [{lo}, {hi}] not reachable with window="
            f"{window}: {len(offenders)} offending window(s) remain; "
            f"global synonymous-codon GC range is [{gmin / n:.3f}, "
            f"{gmax / n:.3f}], but the worst window [{ws},{we}) has a local "
            f"reachable GC of [{wlo:.3f}, {whi:.3f}] — its amino-acid "
            f"composition cannot reach [{lo}, {hi}] (try a larger "
            f"--gc-window or loosen --gc-content)"
        )
    return fixed


def max_homopolymer_run(dna):
    """Longest run of one identical nucleotide in dna."""
    return max(len(m.group()) for m in re.finditer(r"(.)\1+", dna))


def repeated_kmers(dna, min_len, max_len):
    """All "maximal" repeated k-mers with length in [min_len, max_len]:
    a repeated k-mer that also occurs inside a longer reported repeat is
    skipped (nested repeats are implied by their parent, so scanning from
    long to short keeps only the informative ones). Returns
    {(length, kmer): count} — a k-mer repeated r times counts once.
    """
    result = {}
    covered = []  # sequences of already-reported (longer) repeats
    for L in range(max_len, min_len - 1, -1):
        seen = {}
        for i in range(len(dna) - L + 1):
            kmer = dna[i:i + L]
            seen[kmer] = seen.get(kmer, 0) + 1
        for kmer, cnt in seen.items():
            if cnt > 1 and not any(kmer in longer for longer in covered):
                result[(L, kmer)] = cnt
                covered.append(kmer)
    return result


def optimize_aa(aa, species, max_repeat, max_kmer_len, seed, gc_range,
                gc_window, init_codons, synonymous):
    """Scan k = START_KMER..max_kmer_len and return (k, dna, max_run) of
    the first k for which the design converges with all constraints
    passing. Smaller k is harder (UniquifyAllKmers uniqueness), so the
    first success is the tightest achievable; the homopolymer ceiling
    (AvoidPattern at max_repeat) is fixed and never scanned.

    dnachisel raises NoSolutionError when the local search cannot satisfy
    the constraints for the current k — that is the "try a larger k" signal
    of the scan, so it is caught and the scan continues.
    """
    init_dna = "".join(init_codons[x] for x in aa)
    # A window >= the sequence length would be an empty constraint for
    # dnachisel's sliding windows — normalize it to the global check so the
    # pre-optimization and the final constraint agree (per-record: records
    # in one FASTA may differ in length).
    if gc_window is not None and gc_window >= len(init_dna):
        gc_window = None
    if gc_range is not None:
        lo, hi = gc_range
        init_dna = gc_preoptimize(init_dna, aa, lo, hi, gc_window, synonymous)
    for k in range(START_KMER, max_kmer_len + 1):
        print(f"  k={k}: trying ...", file=sys.stderr)
        constraints = [EnforceTranslation(translation=aa)]
        if gc_range is not None:
            lo, hi = gc_range
            constraints.append(EnforceGCContent(mini=lo, maxi=hi,
                                                window=gc_window))
        constraints += [AvoidPattern(base * max_repeat) for base in "ACGT"]
        constraints.append(UniquifyAllKmers(k=k))
        problem = DnaOptimizationProblem(
            sequence=init_dna,
            constraints=constraints,
            objectives=[CodonOptimize(species=species)],
            logger=None,
        )
        try:
            problem.optimize()
        except NoSolutionError:
            pass
        if problem.all_constraints_pass():
            dna = str(problem.sequence)
            return k, dna, max_homopolymer_run(dna)
        print(f"  k={k}: not satisfiable", file=sys.stderr)
    sys.exit(
        f"ERROR: {species} design with unique k-mers not found in kmer "
        f"scan {START_KMER}..{max_kmer_len}; raise --max-kmer-len or "
        f"change --seed"
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage:")[1],
    )
    ap.add_argument("sequence_or_fasta", metavar="SEQ_OR_FASTA",
                    help="amino-acid string (e.g. MKTAY...) or path to a "
                         "protein FASTA file")
    ap.add_argument("-o", "--output", metavar="OUT",
                    help="output file (default: stdout); the format is "
                         "inferred from the suffix (.fa/.fasta/.fna/.fas, "
                         ".csv, .txt); an unrecognized suffix falls back to "
                         "--fmt")
    ap.add_argument("--fmt", choices=("fasta", "csv", "txt"), default=None,
                    help="output format when writing to stdout, or fallback "
                         "for an unrecognized -o suffix (default: fasta)")
    ap.add_argument("--species", default="e_coli",
                    help="target species codon-usage table for dnachisel "
                         "(default: e_coli)")
    ap.add_argument("--max-repeat", type=int, default=START_REPEAT, metavar="N",
                    help="forbid homopolymer runs of N or more nucleotides "
                         "via AvoidPattern (fixed value, no scan; "
                         "default: 7)")
    ap.add_argument("--max-kmer-len", type=int, default=25, metavar="N",
                    help="upper bound of the k-mer uniqueness scan: try "
                         "UniquifyAllKmers(k) for k = 7..N and keep the "
                         "first k that works (smaller k is harder; "
                         "default: 25)")
    ap.add_argument("--seed", type=int, default=42,
                    help="random seed for reproducibility (default: 42)")
    ap.add_argument("--gc-content", metavar="LOWER_UPPER", default=None,
                    help="constrain GC content to the [LOWER, UPPER] fraction "
                         "range (e.g. 0.3_0.7) via dnachisel EnforceGCContent "
                         "(default: off)")
    ap.add_argument("--gc-window", type=int, default=None, metavar="N",
                    help="check the GC content per sliding N-bp window "
                         "instead of the whole sequence; requires "
                         "--gc-content")
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    gc_range = parse_gc_range(args.gc_content) if args.gc_content else None
    if args.gc_window is not None:
        if args.gc_window < 1:
            sys.exit(f"ERROR: --gc-window must be >= 1 bp, got {args.gc_window}")
        if gc_range is None:
            sys.exit("ERROR: --gc-window requires --gc-content")
    if args.max_kmer_len < START_KMER:
        sys.exit(f"ERROR: --max-kmer-len must be >= {START_KMER}, "
                 f"got {args.max_kmer_len}")
    if args.max_repeat < 1:
        sys.exit(f"ERROR: --max-repeat must be >= 1, got {args.max_repeat}")

    # Output format: an -o suffix wins when recognized; an explicit --fmt is
    # the fallback for stdout or for an unrecognized suffix; default fasta.
    fmt = args.fmt
    if args.output:
        ext = os.path.splitext(args.output)[1].lower()
        fmt = FORMAT_BY_SUFFIX.get(ext, fmt)
        if fmt is None:
            sys.exit(
                f"ERROR: cannot infer output format from -o suffix {ext!r}; "
                f"use .fa/.fasta/.fna/.fas, .csv or .txt, or pass --fmt"
            )
    fmt = fmt or "fasta"

    records = read_input(args.sequence_or_fasta)
    init_codons, synonymous, valid_aa = build_codon_tables(args.species)
    for _, _, aa in records:
        check_aa(aa, valid_aa)

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        writer = None
        if fmt == "csv":
            writer = csv.writer(out)
            writer.writerow(["name", "desc", "dna"])
        for name, desc, aa in records:
            print(f"[{name}] optimizing {len(aa)} aa for {args.species} "
                  f"(homopolymer<{args.max_repeat}, kmer scan from "
                  f"{START_KMER}) ...", file=sys.stderr)
            k, dna, max_run = optimize_aa(aa, args.species, args.max_repeat,
                                          args.max_kmer_len, args.seed,
                                          gc_range, args.gc_window,
                                          init_codons, synonymous)
            gc_frac = (dna.count("G") + dna.count("C")) / len(dna)
            reps = repeated_kmers(dna, START_KMER, args.max_kmer_len)
            meta = (f"[dnachisel:{args.species} homopolymer<{args.max_repeat} "
                    f"kmer_unique={k} repeats={len(reps)} "
                    f"max_run={max_run}nt]")
            if fmt == "fasta":
                out.write(">" + " ".join(x for x in (name, desc, meta) if x)
                          + "\n")
                for i in range(0, len(dna), LINE_WIDTH):
                    out.write(dna[i:i + LINE_WIDTH] + "\n")
            elif fmt == "csv":
                assert writer is not None  # created above when fmt == "csv"
                writer.writerow([name, desc, dna])
            else:  # txt
                out.write(" ".join(x for x in (name, desc, meta) if x)
                          + "\n" + dna + "\n\n")
            print(f"[{name}] done: {len(aa)} aa -> {len(dna)} bp | "
                  f"GC {gc_frac:.1%} | avoid homopolymer runs >= "
                  f"{args.max_repeat} nt | unique {k}-mers | actual max run "
                  f"{max_run} nt", file=sys.stderr)
            # report every repeated k-mer with length in the scan range
            # [START_KMER, max_kmer_len] (lengths >= k are unique by
            # construction, so repeats, if any, are shorter than k)
            by_len = {}
            for (L, kmer), cnt in reps.items():
                by_len.setdefault(L, []).append((kmer, cnt))
            if reps:
                summary = " | ".join(
                    f"{L}-mer x{len(v)}" for L, v in sorted(by_len.items()))
                print(f"[{name}] repeated k-mers "
                      f"({START_KMER}..{args.max_kmer_len} nt): "
                      f"{len(reps)} distinct [{summary}]", file=sys.stderr)
                for L in sorted(by_len):
                    for kmer, cnt in sorted(by_len[L]):
                        print(f"  {L}-mer {kmer} x{cnt}", file=sys.stderr)
            else:
                print(f"[{name}] no repeated k-mers in "
                      f"{START_KMER}..{args.max_kmer_len} nt range",
                      file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
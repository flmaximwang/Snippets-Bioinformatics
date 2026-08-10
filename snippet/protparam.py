#!/Applications/BioRazer/bin/python3
"""
protparam — local CLI for protein physico-chemical parameters.
Calculates the same metrics as Expasy ProtParam, backed by Biopython.

Usage:
  protparam SEQUENCE
  protparam -f FILE.fasta
  cat FILE.fasta | protparam

Output: Tabular summary of molecular weight, pI, extinction coefficient,
instability index, aliphatic index, GRAVY, amino acid composition, etc.
"""

import sys
import argparse
import io
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq


# ── helpers ──────────────────────────────────────────────────────────────

def parse_sequences(source: str, fmt: str = "fasta"):
    """Yield (header, seq_str) from a string or file-like."""
    from Bio import SeqIO

    if fmt == "fasta":
        handle = io.StringIO(source) if isinstance(source, str) else source
        for record in SeqIO.parse(handle, "fasta"):
            yield record.id or record.name, str(record.seq).upper()
    else:
        # plain sequence: one or more lines, skip blank/comments
        lines = source.splitlines() if isinstance(source, str) else source.read().splitlines()
        seq_parts = []
        header = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if seq_parts:
                    yield header or "seq", "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.upper().replace(" ", "").replace("\t", ""))
        if seq_parts:
            yield header or "seq", "".join(seq_parts)


def format_results(header: str, seq: str) -> str:
    """Compute and format ProtParam results for one sequence."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    analysed = ProteinAnalysis(seq)

    mw = analysed.molecular_weight()
    pi = analysed.isoelectric_point()
    gravy = analysed.gravy()
    aromaticity = analysed.aromaticity()

    # molar extinction coefficient
    # Biopython returns (reduced, oxidized) tuple
    ec = analysed.molar_extinction_coefficient()
    ec_reduced, ec_oxidized = ec

    # instability index
    ii = analysed.instability_index()
    stable = "stable (II < 40)" if ii < 40 else "unstable (II >= 40)"

    # aliphatic index
    ai = _aliphatic_index(seq)

    # secondary structure
    helix, turn, sheet = analysed.secondary_structure_fraction()

    # amino acid composition
    aa_counts = analysed.count_amino_acids()
    length = len(seq)

    # charge at pH 7
    try:
        charge_pH7 = analysed.charge_at_pH(7.0)
    except Exception:
        charge_pH7 = None

    lines = []
    lines.append("=" * 56)
    lines.append(f"  ProtParam — {header}")
    lines.append("=" * 56)
    lines.append(f"  Sequence length:      {length} aa")
    lines.append(f"  Molecular weight:     {mw:.2f} Da")
    lines.append(f"  Theoretical pI:       {pi:.2f}")
    lines.append(f"  Charge at pH 7:       {charge_pH7:+.2f}" if charge_pH7 is not None else "")
    lines.append(f"  Extinction coeff:")
    lines.append(f"    Reduced Cys:        {ec_reduced} M⁻¹·cm⁻¹")
    lines.append(f"    Oxidized Cys:       {ec_oxidized} M⁻¹·cm⁻¹")
    lines.append(f"  Instability index:    {ii:.2f} — {stable}")
    lines.append(f"  Aliphatic index:      {ai:.2f}")
    lines.append(f"  GRAVY:                {gravy:.3f}")
    lines.append(f"  Aromaticity:          {aromaticity:.3f}")
    lines.append(f"  Secondary structure:")
    lines.append(f"    Helix:              {helix:.2%}")
    lines.append(f"    Turn:               {turn:.2%}")
    lines.append(f"    Sheet:              {sheet:.2%}")
    lines.append(f"  Amino acid composition:")
    lines.append(f"    {'AA':>3s}  {'Count':>5s}  {'%':>7s}")
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        cnt = aa_counts.get(aa, 0)
        pct = 100.0 * cnt / length if length else 0
        lines.append(f"    {aa:>3s}  {cnt:5d}  {pct:6.2f}%")
    lines.append("-" * 56)
    return "\n".join(l for l in lines if l)


def _aliphatic_index(seq: str) -> float:
    """Aliphatic index: relative volume of Ala, Val, Ile, Leu side-chains.
    Same formula as ExPASy ProtParam.
    """
    residues = len(seq)
    if residues == 0:
        return 0.0
    counts = {aa: seq.count(aa) for aa in "AVIL"}
    # Ala: 1.0, Val: 2.9, Ile: 3.9, Leu: 3.9
    return (
        counts.get("A", 0) * 1.0
        + counts.get("V", 0) * 2.9
        + (counts.get("I", 0) + counts.get("L", 0)) * 3.9
    ) * 100.0 / residues


# ── main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="protparam — local protein physico-chemical parameter calculator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  protparam MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRS\n"
            "  protparam -f my_seqs.fasta\n"
            "  cat my_seqs.fasta | protparam\n"
        ),
    )
    parser.add_argument(
        "sequence", nargs="?", default=None,
        help="Protein sequence in single-letter code (or piped via stdin)"
    )
    parser.add_argument(
        "-f", "--file", metavar="FASTA", default=None,
        help="FASTA file with one or more sequences"
    )
    parser.add_argument(
        "--plain", action="store_true", default=False,
        help="Input is a plain sequence line, not FASTA"
    )
    args = parser.parse_args()

    sequences = []

    # 1) FASTA file
    if args.file:
        with open(args.file) as fh:
            for hdr, seq in parse_sequences(fh, fmt="fasta"):
                sequences.append((hdr, seq))

    # 2) Inline argument
    elif args.sequence:
        seq = args.sequence.upper()
        # strip non-standard chars
        clean = "".join(c for c in seq if c in "ACDEFGHIKLMNPQRSTVWY")
        sequences.append(("seq", clean))

    # 3) Stdin
    elif not sys.stdin.isatty():
        data = sys.stdin.read()
        # detect FASTA
        if data.lstrip().startswith(">"):
            for hdr, seq in parse_sequences(data, fmt="fasta"):
                sequences.append((hdr, seq))
        else:
            # plain text, treat as one sequence (skip whitespace/numbers)
            clean = "".join(c.upper() for c in data if c.isalpha() and c.upper() in "ACDEFGHIKLMNPQRSTVWY")
            sequences.append(("stdin", clean))

    else:
        parser.print_help()
        sys.exit(1)

    # Compute and print
    results = []
    for hdr, seq in sequences:
        if not seq:
            results.append(f"  [WARNING] {hdr}: empty sequence – skipped")
            continue
        results.append(format_results(hdr, seq))
    print("\n\n".join(results))

    # Warnings
    for hdr, seq in sequences:
        if seq and ("X" in seq or "B" in seq or "Z" in seq or "J" in seq or "O" in seq or "U" in seq):
            print(f"  [WARNING] {hdr}: contains ambiguous residues (B/Z/J/O/U/X) — "
                  "results may be inaccurate", file=sys.stderr)


if __name__ == "__main__":
    main()

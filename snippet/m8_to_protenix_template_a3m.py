#!/usr/bin/env python3
"""m8_to_protenix_template_a3m.py — convert a ColabFold pdb70.m8 hit table to a
Protenix-compatible hmmsearch-style template a3m.

Protenix reads template hits from `proteinChain.templatesPath`, which must be an
a3m (hmmsearch style) or hhr file. Each entry needs:

    >{pdb}_{chain}/{h_start}-{h_end} [subseq from] mol:protein length:{len}
    <query-width masked alignment, 60 chars/line>

The hit *name* (pdb_id + author chain id) is what Protenix uses to select the
chain from its mmCIF structures; the sequence in the a3m only feeds the
pre-filter (align ratio / duplicate / length checks), since Protenix re-aligns
the query to the mmCIF chain sequence with kalign afterwards.

The sequence line spans the whole query with the alignment embedded:
columns before q_start and after q_end are masked with '-' and the aligned
window q_start..q_end carries the letters, hmmsearch-style:

    M -> uppercase hit residue (aligned column)
    I -> '-' (query residue with no hit residue; gap in the target)
    D -> lowercase hit residue (hit insert between query columns)

So the parser's aligned_cols (uppercase non-'-' characters) equals the true
number of aligned columns — emitting the full template chain instead would
report the whole chain as aligned (align_ratio inflated) and the description
span would not match the actual template region.

Usage:
    m8_to_protenix_template_a3m.py --m8 pdb70.m8 --query query.fa \
        --cif-dir templates/ -o templates.hmmsearch.a3m

Pure stdlib — runs with any python3.
"""

import argparse
import os
import re
import sys
from collections import OrderedDict

# 3-letter -> 1-letter, standard amino acids + common modified residues.
# Anything unknown maps to 'X' (same convention as Protenix's parser).
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K",
    "ASX": "B", "GLX": "Z", "UNK": "X",
}

LINE_WIDTH = 60


def parse_m8(path):
    """Return [(pdb_id, chain, qstart, qend, hstart, hend, cigar), ...].

    Coordinates are m8 columns 7-10, 1-based inclusive; hits may repeat.
    """
    hits = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 13:
                sys.stderr.write(f"WARN: skipping malformed m8 line: {line}\n")
                continue
            target = cols[1]
            m = re.match(r"^([a-zA-Z0-9]{4})_(.+)$", target)
            if not m:
                sys.stderr.write(f"WARN: skipping unparseable target {target!r}\n")
                continue
            try:
                qstart, qend, hstart, hend = (int(cols[i]) for i in (6, 7, 8, 9))
            except ValueError:
                sys.stderr.write(f"WARN: skipping bad coordinates in m8 line: {line}\n")
                continue
            hits.append(
                (m.group(1).lower(), m.group(2), qstart, qend, hstart, hend, cols[12])
            )
    return hits


def read_query_seq(path_or_seq):
    """Return the query sequence from a FASTA/a3m file (first block) or raw string."""
    if os.path.exists(path_or_seq):
        qid = None
        seq = []
        with open(path_or_seq) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if qid is not None:
                        break  # first block complete
                    qid = line[1:].split()[0]
                elif qid is not None:
                    seq.append(line.strip())
        if qid is None:
            raise ValueError(f"{path_or_seq}: no '>' header found")
        return "".join(seq).replace("-", "")
    return path_or_seq  # raw sequence


def generate_masked_seq(query_seq, hit_seq, q_start, q_end, h_start, h_end, cigar):
    """Reconstruct the full-query-width masked alignment for one m8 hit.

    The output spans the WHOLE query: columns before q_start and after q_end
    are masked with '-'; the aligned window q_start..q_end is rendered
    hmmsearch-style:

        M -> uppercase hit residue (aligned column)
        I -> '-' (query residue with no hit residue; gap in the target)
        D -> lowercase hit residue (hit insert between query columns)

    The window rendering is M+I+D characters long (D inserts are extra
    characters), so the output length is len(query) + D. Consistencies
    enforced (violation -> ValueError): M+I == q span and M+D == h span,
    where the spans are q_end-q_start+1 and h_end-h_start+1.
    """
    if q_start < 1 or q_end > len(query_seq) or q_start > q_end:
        raise ValueError(
            f"query window {q_start}-{q_end} outside query length {len(query_seq)}"
        )
    if h_start < 1 or h_end > len(hit_seq) or h_start > h_end:
        raise ValueError(
            f"hit window {h_start}-{h_end} outside hit length {len(hit_seq)}"
        )

    ops = re.findall(r"(\d+)([MID])", cigar)
    if not ops:
        raise ValueError(f"unparseable CIGAR {cigar!r}")

    m = sum(int(n) for n, op in ops if op == "M")
    i = sum(int(n) for n, op in ops if op == "I")
    d = sum(int(n) for n, op in ops if op == "D")
    if m + i != q_end - q_start + 1 or m + d != h_end - h_start + 1:
        raise ValueError(
            f"CIGAR {cigar} (M={m} I={i} D={d}) inconsistent with "
            f"query span {q_end - q_start + 1} / hit span {h_end - h_start + 1}"
        )

    out = []
    hi = h_start - 1  # 0-based cursor into hit_seq
    for n, op in ops:
        n = int(n)
        if op == "M":
            out.append(hit_seq[hi : hi + n].upper())
            hi += n
        elif op == "I":
            out.append("-" * n)
        else:  # D
            out.append(hit_seq[hi : hi + n].lower())
            hi += n
    # frame the window in the full query: non-template regions are '-'
    return "-" * (q_start - 1) + "".join(out) + "-" * (len(query_seq) - q_end)


def extract_chain_sequences(cif_path):
    """Parse the _atom_site loop of an mmCIF and return {auth_chain_id: sequence}.

    Prefers `_atom_site.auth_asym_id` for chain identity (matches the m8 target
    names); falls back to `_atom_site.label_asym_id` if auth is missing.
    Only group_PDB == 'ATOM' rows of model 1 are used, ordered by residue
    position (label_seq_id, else auth_seq_id + insertion code).
    """
    seqs = {}  # chain -> OrderedDict[(pos_key)] -> resname
    with open(cif_path) as f:
        lines = f.readlines()

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line != "loop_":
            i += 1
            continue
        # collect column names of this loop
        i += 1
        cols = []
        while i < n and lines[i].strip().startswith("_"):
            cols.append(lines[i].strip())
            i += 1
        if not any(c.startswith("_atom_site.") for c in cols):
            continue  # not the atom_site loop; skip its data rows
        # data rows end at '#', 'loop_', or another '_'-prefixed line
        while i < n:
            row = lines[i].strip()
            if row == "" or row == "#" or row.startswith("_"):
                break
            if row == "loop_":
                break
            tokens = row.split()
            if len(tokens) < len(cols):
                tokens += ["."] * (len(cols) - len(tokens))
            rec = dict(zip(cols, tokens))
            i += 1
            if rec.get("_atom_site.group_PDB") != "ATOM":
                continue
            model = rec.get("_atom_site.pdbx_PDB_model_num", ".")
            if model not in ("1", ".", "?"):
                continue
            chain = rec.get("_atom_site.auth_asym_id", ".")
            if chain in (".", "?"):
                chain = rec.get("_atom_site.label_asym_id", ".")
            if chain in (".", "?"):
                continue
            res = rec.get("_atom_site.label_comp_id", ".")
            if res in (".", "?"):
                continue
            # residue position key
            lseq = rec.get("_atom_site.label_seq_id", ".")
            if lseq not in (".", "?"):
                pos = ("L", int(lseq))
            else:
                aseq = rec.get("_atom_site.auth_seq_id", ".")
                ins = rec.get("_atom_site.pdbx_PDB_ins_code", ".")
                if aseq in (".", "?"):
                    continue
                pos = ("A", int(aseq), ins if ins not in (".", "?") else "")
            seqs.setdefault(chain, OrderedDict())[pos] = res
        break  # _atom_site loop consumed; nothing else needed

    out = {}
    for chain, residues in seqs.items():
        letters = []
        for res in residues.values():
            letters.append(THREE_TO_ONE.get(res, "X"))
        out[chain] = "".join(letters)
    return out


def wrap(seq, width=LINE_WIDTH):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(
        description="Convert ColabFold pdb70.m8 template hits to a "
        "Protenix hmmsearch-style template a3m."
    )
    ap.add_argument("--m8", required=True, help="pdb70.m8 hit table (MMseqs2 format)")
    ap.add_argument(
        "--query",
        required=True,
        help="query sequence: FASTA/a3m file (first block) or raw AA string",
    )
    ap.add_argument(
        "--cif-dir",
        required=True,
        help="directory with template structures: tries {pdb}_{chain}.cif, then {pdb}.cif",
    )
    ap.add_argument("-o", "--out", required=True, help="output a3m path")
    args = ap.parse_args()

    query_seq = read_query_seq(args.query)
    hits = parse_m8(args.m8)
    if not hits:
        sys.stderr.write("ERROR: no hits parsed from m8 file\n")
        return 1

    # cache parsed sequences per file (multiple chains can live in one file)
    cif_cache = {}
    written, skipped = 0, 0
    with open(args.out, "w") as out:
        for pdb, chain, qstart, qend, hstart, hend, cigar in hits:
            seq = None
            for fname in (f"{pdb}_{chain}.cif", f"{pdb}.cif"):
                path = os.path.join(args.cif_dir, fname)
                if os.path.exists(path):
                    if path not in cif_cache:
                        cif_cache[path] = extract_chain_sequences(path)
                    seq = cif_cache[path].get(chain)
                    if seq is None:
                        sys.stderr.write(
                            f"WARN: chain {chain!r} not found in {fname} "
                            f"(available: {sorted(cif_cache[path])})\n"
                        )
                    break
            if seq is None:
                sys.stderr.write(f"WARN: skipping {pdb}_{chain}: no structure/chain\n")
                skipped += 1
                continue
            try:
                masked = generate_masked_seq(
                    query_seq, seq, qstart, qend, hstart, hend, cigar
                )
            except ValueError as e:
                sys.stderr.write(f"WARN: skipping {pdb}_{chain}: {e}\n")
                skipped += 1
                continue
            header = (
                f">{pdb}_{chain}/{hstart}-{hend} [subseq from] "
                f"mol:protein length:{len(seq)}"
            )
            out.write(header + "\n" + wrap(masked) + "\n")
            written += 1

    print(
        f"m8 hits: {len(hits)}  written: {written}  skipped: {skipped}\n"
        f"output: {args.out}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

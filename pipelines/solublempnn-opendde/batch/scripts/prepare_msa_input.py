#!/usr/bin/env python3
"""prepare_msa_input.py - 取 SolubleMPNN seqs/*.fa 的第 2 条记录(第一条设计序列),
按模板切分位置切成 3 片段 fasta 作为 colabfold-msa 输入。

切分位置与模板一致: 前 frag1_len 个残基 : 中间 mid_len 个残基 : 其余。
中间片段直接采用设计序列自身残基, 不做任何替换。

stdout 输出三行供上层脚本解析:
    DESIGNED <完整设计序列>      (记录到 summary.csv 用)
    NORMALIZED <完整设计序列>    (opendde.json sequence 用, 与 MSA query 一致)
    SPLIT <frag1>:<mid>:<frag2>
"""
import argparse
import sys
from pathlib import Path


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", required=True, help="SolubleMPNN 输出 seqs/*.fa")
    p.add_argument("-o", "--output", required=True, help="输出 3 片段 fasta (colabfold-msa 输入)")
    p.add_argument("--record-name", default="candidate-8.n.0")
    p.add_argument("--frag1-len", type=int, default=79)
    p.add_argument("--mid-len", type=int, default=6)
    args = p.parse_args()

    records = parse_fasta(Path(args.input).read_text())
    if len(records) < 2:
        sys.exit(f"error: {args.input} 只有 {len(records)} 条记录, 需要 >= 2 条(第 1 条为参考, 第 2 条为第一条设计)")
    ref_seq = records[0][1].upper()
    header, designed = records[1]
    seq = designed.upper()

    if len(seq) != len(ref_seq):
        sys.exit(f"error: 设计序列长度 {len(seq)} != 参考序列长度 {len(ref_seq)}")

    cut1 = args.frag1_len
    cut2 = args.frag1_len + args.mid_len
    frag1 = seq[:cut1]
    mid = seq[cut1:cut2]
    frag2 = seq[cut2:]
    if len(mid) != args.mid_len or len(frag1) + len(mid) + len(frag2) != len(seq):
        sys.exit(f"error: 片段长度不匹配: {len(frag1)}+{len(mid)}+{len(frag2)} != {len(seq)}")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(f">{args.record_name}\n{frag1}:{mid}:{frag2}\n")

    print(f"DESIGNED {seq}")
    print(f"NORMALIZED {seq}")
    print(f"SPLIT {frag1}:{mid}:{frag2}")
    print(f"# record2 header={header}, len={len(seq)}", file=sys.stderr)


if __name__ == "__main__":
    main()
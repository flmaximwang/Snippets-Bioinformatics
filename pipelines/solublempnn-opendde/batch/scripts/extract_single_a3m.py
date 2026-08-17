#!/usr/bin/env python3
"""extract_single_a3m.py - 提取 paired.a3m 的第一条记录(query)为 single.a3m。

OpenDDE 的 pairedMsaPath 需要 query-only 的 a3m (README.md Step 2:
"Extract the first sequence (query) to single.a3m")。
"""
import argparse
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paired-msa", required=True, help="paired/paired.a3m 路径")
    p.add_argument("-o", "--output", required=True, help="输出 single.a3m 路径")
    args = p.parse_args()

    text = Path(args.paired_msa).read_text()
    for rec in text.split(">"):
        if not rec.strip():
            continue
        header, _, rest = rec.partition("\n")
        seq = "".join(rest.split())
        Path(args.output).write_text(">" + header.strip() + "\n" + seq + "\n")
        print(f"# query extracted: {len(seq)} aa -> {args.output}")
        return
    sys.exit(f"error: {args.paired_msa} 中没有记录")


if __name__ == "__main__":
    main()
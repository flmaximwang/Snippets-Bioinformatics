#!/usr/bin/env python3
"""make_opendde_json.py - 基于 opendde_template.json 生成单个 candidate 的 OpenDDE 输入。

按 AGENTS.md: paired.a3m 放 unpairedMsaPath, single.a3m 放 pairedMsaPath。
不填 templatesPath (本阶段不用 template)。sequence 必须与 MSA query 一致。
"""
import argparse
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", required=True, help="opendde_template.json 路径")
    p.add_argument("--sequence", required=True, help="全长序列(NORMALIZED, 与 MSA query 一致)")
    p.add_argument("--paired-msa", required=True, metavar="ABS_PATH", help="single.a3m 绝对路径")
    p.add_argument("--unpaired-msa", required=True, metavar="ABS_PATH", help="paired/paired.a3m 绝对路径")
    p.add_argument("--name", required=True, help="job name, 例如 candidate-8.0.0")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args()

    data = json.loads(Path(args.template).read_text())
    job = data[0]
    chain = job["sequences"][0]["proteinChain"]
    chain["sequence"] = args.sequence
    chain["pairedMsaPath"] = args.paired_msa
    chain["unpairedMsaPath"] = args.unpaired_msa
    job["name"] = args.name

    Path(args.output).write_text(json.dumps(data, indent=4, ensure_ascii=False) + "\n")
    print(f"# {args.name}: sequence={len(args.sequence)} aa, pairedMsaPath={args.paired_msa}, unpairedMsaPath={args.unpaired_msa}")


if __name__ == "__main__":
    main()
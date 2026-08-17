#!/usr/bin/env python3
"""check_candidate.py - OpenDDE 预测结果一致性检查。

每个 sample: 读 summary_confidence_*.json 的 ptm; 用现成 CLI
(renumber_protein_structure_file.py) 把 .cif 重编号成 .pdb
(--renumber-res 1_1,86_107, 与 v8-0-0 相同); biotite 叠合对齐参考结构,
计算 backbone (N/CA/C/O) RMSD。默认排除 linker 残基(80-85, 任意序列)。

按 ptm 降序检查, 第一个同时满足 ptm>=cutoff 且 rmsd<=cutoff 的 sample 即 ACCEPT。
stdout: 每 sample 一行 + 最终 VERDICT 行。
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as bs
from biorazer.structure.io.protein import Pdb_AtomArray

BB_ATOMS = ("N", "CA", "C", "O")


def parse_res_range(s: str) -> set[int]:
    m = re.match(r"^(\d+)-(\d+)$", s)
    if not m:
        sys.exit(f"error: 无法解析残基范围: {s} (期望如 80-85)")
    return set(range(int(m.group(1)), int(m.group(2)) + 1))


def backbone_rmsd(ref_path: str, mob_path: str, exclude_res: set[int]) -> tuple[float, int]:
    ref = Pdb_AtomArray(input_io=ref_path).read()
    mob = Pdb_AtomArray(input_io=mob_path).read()
    ref = ref[ref.chain_id == ref.chain_id[0]]
    mob = mob[mob.chain_id == mob.chain_id[0]]

    rid_r = bs.get_residues(ref)[0]
    rid_m = bs.get_residues(mob)[0]
    common = np.intersect1d(rid_r, rid_m)
    common = common[~np.isin(common, np.array(sorted(exclude_res)))]

    mask_ref = np.isin(ref.res_id, common) & np.isin(ref.atom_name, BB_ATOMS)
    mask_mob = np.isin(mob.res_id, common) & np.isin(mob.atom_name, BB_ATOMS)
    if mask_ref.sum() != mask_mob.sum():
        sys.exit(f"error: backbone 原子数不匹配 ref={mask_ref.sum()} vs mob={mask_mob.sum()}")

    aligned, _ = bs.superimpose(ref[mask_ref], mob[mask_mob])
    rmsd = float(bs.rmsd(ref[mask_ref], aligned))
    n_res = int(mask_ref.sum() / len(BB_ATOMS))
    return rmsd, n_res


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cand-dir", required=True, help="candidate 数据目录 (data/candidate-8.n.0)")
    p.add_argument("--name", required=True, help="opendde job name (candidate-8.n.0)")
    p.add_argument("--ref-pdb", required=True, help="参考结构 v8-0-0_mutate_conserved.pdb")
    p.add_argument("--renumber-script", required=True, help="renumber_protein_structure_file.py 路径")
    p.add_argument("--renumber-spec", default="1_1,86_107")
    p.add_argument("--ptm-cutoff", type=float, default=0.8)
    p.add_argument("--rmsd-cutoff", type=float, default=1.5)
    p.add_argument("--exclude-linker", default="80-85")
    args = p.parse_args()
    exclude = parse_res_range(args.exclude_linker)

    cand = Path(args.cand_dir)
    # predictions 位于 <cand-dir>/seed_*/predictions (cand-dir = CAND_DIR, 平铺)
    pred_dirs = sorted(cand.glob("seed_*/predictions"))
    if not pred_dirs:
        sys.exit(f"error: {cand} 下找不到 seed_*/predictions, OpenDDE 尚未运行?")
    pred_dir = pred_dirs[0]

    samples: dict[int, float] = {}
    for j in pred_dir.glob(f"{args.name}_summary_confidence_sample_*.json"):
        m = re.search(r"sample_(\d+)\.json$", j.name)
        if m is None:
            continue
        ptm = json.loads(j.read_text()).get("ptm")
        if ptm is not None:
            samples[int(m.group(1))] = float(ptm)
    if not samples:
        sys.exit(f"error: {pred_dir} 下没有 summary_confidence json")

    rows = []
    chosen = None
    for i in sorted(samples, key=lambda x: samples[x], reverse=True):
        cif = pred_dir / f"{args.name}_sample_{i}.cif"
        rpdb = pred_dir / f"{args.name}_sample_{i}_renumbered.pdb"
        if not cif.exists():
            rows.append((i, samples[i], None, 0, "no cif"))
            continue
        if not rpdb.exists():
            subprocess.run(
                [sys.executable, args.renumber_script, "-i", str(cif),
                 "-o", str(rpdb), "--renumber-res", args.renumber_spec],
                check=True, capture_output=True, text=True,
            )
        rmsd, n_res = backbone_rmsd(args.ref_pdb, str(rpdb), exclude)
        ok = samples[i] >= args.ptm_cutoff and rmsd <= args.rmsd_cutoff
        rows.append((i, samples[i], rmsd, n_res, "PASS" if ok else "fail"))
        if ok and chosen is None:
            chosen = (i, samples[i], rmsd)

    for i, ptm, rmsd, n_res, flag in rows:
        if rmsd is None:
            print(f"sample_{i}: ptm={ptm:.4f} rmsd=NA ({flag})")
        else:
            print(f"sample_{i}: ptm={ptm:.4f} rmsd={rmsd:.3f} ({n_res} res backbone, excl {args.exclude_linker}) ({flag})")
    if chosen is not None:
        print(f"VERDICT ACCEPT sample_{chosen[0]} ptm={chosen[1]:.4f} rmsd={chosen[2]:.3f}")
    else:
        print("VERDICT REJECT")


if __name__ == "__main__":
    main()
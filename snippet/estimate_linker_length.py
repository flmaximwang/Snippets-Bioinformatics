#!/usr/bin/env python3
"""estimate_linker_length.py — 根据两个链的末端几何估算 linker 长度范围 (三种物理状态)

输入: 一个含 2 条链 (默认 A, B) 的 PDB / CIF 文件。
输出: 基于链 A 的 C 端与链 B 的 N 端 (距离 + 朝向) 估算的 linker 残基数范围。

几何量定义:
  d_CN   = |C(A末) - N(B首)|               RPXDock --termini_max_dist 对应的缺口
  d_CC   = |CA(A末) - CA(B首)|             残基数估算的标准量 (Cα-Cα)
  出口方向 v_exit  = CA(A末-1) -> CA(A末)    链 A 末端 Cα 虚键方向 (走向)
  入口方向 v_entry = CA(B首) -> CA(B首+1)   链 B 起始 Cα 虚键方向 (走向)
  缺口方向 gap     = CA(A末) -> CA(B首)      即 d_CC 方向
  alpha = 夹角(v_exit, gap)                  A 端是否需要折返 (>90° 需要)
  beta  = 夹角(v_entry, gap)                 B 端是否需要折返 (>90° 需要)
  只依赖 Cα, 无需重建缺失原子; 单残基链端无虚键 -> 该端角度记 N/A (不计折返)。

linker 残基数为 L 时, 在 CA(A末) 与 CA(B首) 之间构成 L+1 个 Cα-Cα 键 (虚键)。

三种物理状态的估算 (同一缺口 d_CC):
  1) 完全伸展 (β 伸展态): 每键 ~3.8 Å
       L_ext = max(2, ceil(d_CC/3.8) - 1 + 2*n_turn)      # 硬下限 (链近直线)
       3.8 Å = trans 肽键 Cα-Cα 虚键长 / WLC contour 每残基长 [van Rosmalen 2017]
  2) 平均 (柔性 linker, 蠕虫链 WLC): 端到端均方距离
       <R^2> = 2*lp*Lc*(1 - (lp/Lc)*(1 - e^(-Lc/lp))),  Lc = (L+1)*3.8
       解 <R^2> = d_CC^2 求 Lc -> L_avg = Lc/3.8 - 1 + 2*n_turn  # 柔性 linker 典型值
       lp (persistence length) 默认 4.5 Å, 实验 FRET 校准 [van Rosmalen 2017]
  3) 蜷缩 (α-螺旋, 最紧凑规则结构): 轴向每残基 1.5 Å
       L_col = max(2, ceil(d_CC/1.5) - 1 + 2*n_turn)      # 紧凑上界
       1.5 Å = α-螺旋轴向每残基位移 (3.6 残基/圈 x 1.5 Å = 5.4 Å pitch)
       [Pauling 1951; Fitzkee & Rose 2004]

折返: 若 α 或 β > 90°, 该端需要折返, 三种状态统一各 +2 残基 (作转向)。
       末端若背离缺口, 无论柔性还是刚性 linker 都得先折返回来, 故统一加。

推荐: 三种物理状态给的是"同一缺口 d_CC 下不同 linker 形态所需的残基数",
      并非总有 L_ext < L_avg < L_col — 缺口很长时柔性线圈因端距~√n 增长,
      L_avg 会超过 L_col。各代表一种设计取向:
        L_ext  最短可行 (链拉直, contour 极限, 硬下限)
        L_avg  若用柔性 (Gly/Ser) 线圈 linker, 自然舒展时的典型值
        L_col  若刻意设计成刚性 α-螺旋/紧凑形式所需
      常规柔性 linker 实际可行区间取 [L_ext, L_avg]; 更短会绷紧 (熵代价高),
      L_avg 即热力学平均端距等于缺口的长度。

用法:
  python estimate_linker_length.py -i <file.pdb|file.cif> [--chain-a A] [--chain-b B]
      [--lp 4.5]   # WLC persistence length (Å); 柔性 Gly/Ser ~4.5, polyserine ~6.2

参考文献 (数值均已核对原文):
  [1] van Rosmalen, M.; Krom, M.; Merkx, M. Tuning the Flexibility of
      Glycine-Serine Linkers To Allow Rational Design of Multidomain Proteins.
      Biochemistry 2017, 56, 6565-6574. DOI: 10.1021/acs.biochem.7b00902
      (WLC 拟合柔性 Gly/Ser linker: lp = 4.5/4.8/6.2 Å 对应 GSSGSS/GSSSSSS/SSSSSS;
       b0 = 3.8 Å/残基 contour 长)
  [2] Fitzkee, N. C.; Rose, G. D. Reassessing random-coil statistics in unfolded
      proteins. Proc. Natl. Acad. Sci. U.S.A. 2004, 101, 12497-12502.
      DOI: 10.1073/pnas.0404236101
      (随机线圈标度律 Rg = R0*N^v, v 0.33 蜷缩 / 0.5 理想 / 0.6 良溶剂;
       引 "α-helical rod has a length of 1.50 Å per residue")
  [3] Pauling, L.; Corey, R. B.; Branson, H. R. The structure of proteins: two
      hydrogen-bonded helical configurations of the polypeptide chain.
      Proc. Natl. Acad. Sci. U.S.A. 1951, 37, 205-211. DOI: 10.1073/pnas.37.4.205
      (α-螺旋 3.6 残基/圈, 5.4 Å pitch -> 轴向每残基 1.5 Å)
"""

import argparse
import io
import math
import sys

import numpy as np

from biorazer.structure.io.protein import Pdb_AtomArray, Cif_AtomArray
import biotite.structure as bio_struc

RESPACE_CA = 3.8    # 全伸展 Cα-Cα 虚键长 / WLC contour 每残基长 (Å/残基) [1]
RESPACE_HELIX = 1.5  # α-螺旋轴向每残基位移 (Å/残基) [2,3]
DEFAULT_LP = 4.5     # WLC persistence length, 柔性 Gly/Ser linker 实验值 (Å) [1]


def load_array(path: str):
    """按扩展名读入结构, 返回第一个 model 的 AtomArray。"""
    if path.lower().endswith((".cif", ".mmcif")):
        return Cif_AtomArray(path, io.StringIO()).read()
    return Pdb_AtomArray(path, io.StringIO()).read()


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else None


def angle_deg(u, v):
    """两方向向量夹角 (度)。任一为 None 返回 None。"""
    if u is None or v is None:
        return None
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-12)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def residue_atoms(chain, res_id):
    """取某残基的原子坐标 dict {原子名: coord}, 同名字取首个 (altloc A)。"""
    sel = chain[chain.res_id == res_id]
    atoms = {}
    for atom in sel:
        name = atom.atom_name.strip()
        if name not in atoms:
            atoms[name] = atom.coord
    return atoms


def need(name, atoms, chain_id, res_id):
    if name not in atoms:
        sys.exit(f"[错误] 链 {chain_id} 残基 {res_id} 缺少原子 {name}, 无法估算")
    return atoms[name]


def wlc_msd(contour, lp):
    """WLC 端到端均方距离 <R^2> = 2*lp*Lc*(1 - (lp/Lc)*(1 - e^(-Lc/lp)))。"""
    x = contour / lp
    return 2.0 * lp * contour * (1.0 - (1.0 - math.exp(-x)) / x)


def wlc_linker_length(d_cc, lp):
    """求柔性 linker 残基数 L, 使 WLC 的 RMS 端到端距离 = d_CC。

    返回 L (>=2)。contour Lc = (L+1)*RESPACE_CA; 解 <R^2>=d_CC^2。
    """
    target = d_cc * d_cc
    lo = RESPACE_CA  # 至少 1 个虚键 (L=0), 但其 RMS 可能已 > d_CC
    # 上界: 从 1 键起倍增直至 msd >= target
    hi = lo
    while wlc_msd(hi, lp) < target and hi < 1e6:
        hi *= 2
    if wlc_msd(lo, lp) >= target:
        # 缺口小于单个键的 RMS -> 无法"平均"对齐, 取最小 2 残基
        return 2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if wlc_msd(mid, lp) < target:
            lo = mid
        else:
            hi = mid
    lc = 0.5 * (lo + hi)
    return max(2, int(round(lc / RESPACE_CA - 1)))


def main():
    ap = argparse.ArgumentParser(
        description="根据链 A C 端与链 B N 端的距离和朝向, 按三种物理状态估算 linker 长度范围"
    )
    ap.add_argument("-i", "--input", required=True, help="输入 PDB / CIF (含 A, B 两条链)")
    ap.add_argument("--chain-a", default="A", help="链 A (C 端侧), 默认 A")
    ap.add_argument("--chain-b", default="B", help="链 B (N 端侧), 默认 B")
    ap.add_argument("--lp", type=float, default=DEFAULT_LP,
                    help=f"WLC persistence length (Å), 默认 {DEFAULT_LP} (柔性 Gly/Ser); "
                         f"更少 Gly 的 linker 用更大值 (polyserine ~6.2)")
    args = ap.parse_args()
    if args.lp <= 0:
        sys.exit("[错误] --lp 必须为正")

    arr = load_array(args.input)
    arr = arr[bio_struc.filter_amino_acids(arr)]

    chains = sorted(set(str(c) for c in arr.chain_id))
    if args.chain_a not in chains:
        sys.exit(f"[错误] 链 {args.chain_a} 不存在, 可用链: {chains}")
    if args.chain_b not in chains:
        sys.exit(f"[错误] 链 {args.chain_b} 不存在, 可用链: {chains}")
    if args.chain_a == args.chain_b:
        sys.exit("[错误] 链 A 与链 B 不能相同")

    ca = arr[arr.chain_id == args.chain_a]
    cb = arr[arr.chain_id == args.chain_b]

    res_a_last = max(ca.res_id)   # 链 A 的 C 端残基
    res_b_first = min(cb.res_id)  # 链 B 的 N 端残基
    aa = residue_atoms(ca, res_a_last)
    bb = residue_atoms(cb, res_b_first)

    c_a = need("C", aa, args.chain_a, res_a_last)
    ca_a = need("CA", aa, args.chain_a, res_a_last)
    n_b = need("N", bb, args.chain_b, res_b_first)
    ca_b = need("CA", bb, args.chain_b, res_b_first)

    # ---- 距离 ----
    d_cn = float(np.linalg.norm(c_a - n_b))
    d_cc = float(np.linalg.norm(ca_a - ca_b))

    # ---- 朝向 (用 Cα 虚键方向, 只依赖 Cα, 免去重建缺失原子) ----
    #  A 端出口方向 = 链 A 最末 Cα 虚键: CA(n-1) -> CA(n)
    #  B 端入口方向 = 链 B 最首 Cα 虚键: CA(1) -> CA(2)
    #  缺口方向     = CA(A末) -> CA(B首)  (即 d_CC 方向)
    res_a = sorted(set(int(r) for r in ca.res_id))
    res_b = sorted(set(int(r) for r in cb.res_id))
    ca_prev = None
    if len(res_a) >= 2:
        rp = res_a[-2]
        hit = ca[ca.res_id == rp]
        if len(hit):
            ca_prev = hit.coord[0]
    ca_next = None
    if len(res_b) >= 2:
        rn = res_b[1]
        hit = cb[cb.res_id == rn]
        if len(hit):
            ca_next = hit.coord[0]

    v_exit = unit(ca_a - ca_prev) if ca_prev is not None else None
    v_entry = unit(ca_next - ca_b) if ca_next is not None else None
    gap = unit(ca_b - ca_a)          # 缺口方向: CA(A末) -> CA(B首)

    alpha = angle_deg(gap, v_exit) if v_exit is not None else None   # A 端转向角
    beta = angle_deg(gap, v_entry) if v_entry is not None else None  # B 端转向角

    # ---- 三种状态的残基数估算 ----
    n_turn = int((alpha is not None and alpha > 90)) + int((beta is not None and beta > 90))
    turn_extra = 2 * n_turn

    l_ext = max(2, math.ceil(d_cc / RESPACE_CA) - 1 + turn_extra)   # 完全伸展 (硬下限)
    l_avg = max(2, wlc_linker_length(d_cc, args.lp) + turn_extra)   # 平均 (柔性 WLC, 含折返)
    l_col = max(2, math.ceil(d_cc / RESPACE_HELIX) - 1 + turn_extra)  # 蜷缩 (α-螺旋上界)

    # ---- 输出 ----
    print(f"输入      : {args.input}")
    print(f"链 A (C端): 残基 {res_a_last}  {''.join(ca.res_name[ca.res_id == res_a_last][:1])}")
    print(f"链 B (N端): 残基 {res_b_first}  {''.join(cb.res_name[cb.res_id == res_b_first][:1])}")
    print()
    print("几何量:")
    print(f"  C(A末)-N(B首) 距离 d_CN  = {d_cn:6.2f} Å   (RPXDock --termini_max_dist 参考)")
    print(f"  CA(A末)-CA(B首) 距离 d_CC = {d_cc:6.2f} Å   (残基数估算标准量)")
    print(f"  A 端出口 vs 缺口夹角 α    = {('%.1f°' % alpha) if alpha is not None else 'N/A'}   {'朝向缺口, 无需折返' if (alpha is not None and alpha <= 90) else '背离缺口, 需要折返'}")
    print(f"  B 端入口 vs 缺口夹角 β    = {('%.1f°' % beta) if beta is not None else 'N/A'}   {'朝向缺口, 无需折返' if (beta is not None and beta <= 90) else '背离缺口, 需要折返'}")
    print()
    print("三种状态的 linker 残基数 (折返端数 = %d, 各 +%d):" % (n_turn, 2 * n_turn))
    print(f"  完全伸展   L_ext = {l_ext:4d}   每键 ~3.8 Å/残基 [1]          硬下限 (链近直线)")
    print(f"  平均(WLC)  L_avg = {l_avg:4d}   lp = {args.lp:.1f} Å 柔性 Gly/Ser [1]  柔性 linker 典型值")
    print(f"  蜷缩(α螺旋) L_col = {l_col:4d}   每残基 ~1.5 Å [2,3]        紧凑/刚性螺旋选项")
    print()
    print("推荐 (三种状态各代表一种设计取向, 长缺口下 L_avg 会超 L_col):")
    print(f"  完全伸展 L_ext = {l_ext} 残基   硬下限 (链拉直, contour 极限)")
    print(f"  柔性线圈 区间   = {l_ext} - {l_avg} 残基   (更短绷紧, L_avg = 热力学平均端距=缺口)")
    print(f"  刚性 α-螺旋     = {l_col} 残基   (刻意设计成紧凑/螺旋时所需)")
    print()
    print("提示: RPXDock --termini_max_dist 可参考 d_CN + 2~3 Å 设置;")
    print("      折返端额外残基 (+%d/端) 已统一计入三种状态。" % (2))


if __name__ == "__main__":
    main()

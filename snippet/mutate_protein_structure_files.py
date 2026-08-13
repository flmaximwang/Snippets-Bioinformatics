import argparse
import re
from pathlib import Path
import numpy as np
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.io.protein import Pdb_AtomArray, Cif_AtomArray, AtomArray_Pdb, AtomArray_Cif
from biorazer.structure.manipulation.modification import mutate_without_side_chains

def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", metavar="INPUT", help="输入 PDB/CIF 结构文件")
    p.add_argument("--mutation-spec", metavar="SPEC", help="突变 specification, 逗号分隔多个突变; 格式 <源1字母><res_id><目标1字母>, 如 A1M,N56E")
    p.add_argument("-o", "--output", metavar="OUTPUT", help="输出结构文件, 格式由扩展名决定 (默认 <input>_mutated.<ext>)")

    args = p.parse_args()
    return args

def read_structure(path: str):
    """按扩展名选择 PDB/CIF 读取器, 返回 biotite AtomArray。"""
    suffix = Path(path).suffix.lower()
    if suffix == ".pdb":
        return Pdb_AtomArray(input_io=path).read()
    if suffix in (".cif", ".mmcif"):
        return Cif_AtomArray(input_io=path).read()
    raise ValueError(f"unsupported input format: {path} (expect .pdb or .cif)")

def write_structure(array, path: str):
    """按扩展名选择 PDB/CIF 写入器。"""
    suffix = Path(path).suffix.lower()
    if suffix == ".pdb":
        AtomArray_Pdb(output_io=path).write(array)
    elif suffix in (".cif", ".mmcif"):
        AtomArray_Cif(output_io=path).write(array)
    else:
        raise ValueError(f"unsupported output format: {path} (expect .pdb or .cif)")

def main():
    args = parse_args()
    if not args.mutation_spec:
        raise ValueError("--mutation-spec is required (e.g. A1M,N56E)")

    array = read_structure(args.input)
    # mutate_without_side_chains 以 PDB 标准大写 res_name 校验源字母, 但其自身输出
    # 的是首字母大写 (如 "Ala"), 导致输出文件再次突变时校验失败; 读入后统一大写。
    array.res_name = np.char.upper(array.res_name)
    spec = [entry.strip() for entry in args.mutation_spec.split(",")]
    mutated = mutate_without_side_chains(array, spec)

    output = args.output or str(Path(args.input).with_name(Path(args.input).stem + "_mutated" + Path(args.input).suffix))
    write_structure(mutated, output)

    print_with_decoration("Mutated Structure")
    print(f"{args.input} -> {output}")
    print_decoration_line()
    for entry in spec:
        match = re.match(r"^[A-Za-z](\d+)[A-Za-z]$", entry)
        res_id = int(match.group(1)) if match else None
        if res_id is None:
            continue
        res_mask = mutated.res_id == res_id
        n_atoms = int(res_mask.sum())
        n_res = len(set(zip(mutated.chain_id[res_mask], mutated.res_id[res_mask], mutated.ins_code[res_mask])))
        if n_res == 1:
            chain = mutated.chain_id[res_mask][0]
            name = mutated.res_name[res_mask][0]
            print(f"{entry}: chain {chain} residue {res_id} -> {name}, {n_atoms} backbone atoms kept")
        else:
            print(f"{entry}: {n_res} residues matched (res_id {res_id} shared), {n_atoms} backbone atoms kept")
    print_decoration_line()

if __name__ == "__main__":
    main()
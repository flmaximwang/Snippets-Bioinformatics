import argparse, re
from pathlib import Path
import numpy as np
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.io.protein import Pdb_AtomArray, Cif_AtomArray, AtomArray_Pdb, AtomArray_Cif
import biotite.structure as bio_struc

def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", metavar="INPUT", help="输入 PDB/CIF 结构文件")
    p.add_argument("--renumber-res", metavar="SPEC", help="重新编号 specification, 逗号分隔多个 anchor; 多链格式 C_P_ID, 单链格式 P_ID (P 为 1-based 位置)")
    p.add_argument("-o", "--output", metavar="OUTPUT", help="输出结构文件, 格式由扩展名决定 (默认 <input>_renumbered.<ext>)")

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

def parse_renumber_spec(renumber_res: str | None, chain_ids: list[str]) -> dict:
    """解析 --renumber-res 为 {chain_id: [(position, res_id), ...]}。

    单链格式 P_ID, 多链格式 C_P_ID; P 为 1-based 位置。
    """
    if renumber_res is None:
        return {chain_id: [] for chain_id in chain_ids}

    anchors_tmp = renumber_res.split(",")
    if len(chain_ids) == 1:
        anchors = []
        for anchor_tmp in anchors_tmp:
            if not re.match(r"^\d+_\d+$", anchor_tmp):
                raise ValueError(f"Unsupported renumber_res spec: {anchor_tmp}")
            position, res_id = anchor_tmp.split("_")
            anchors.append((int(position), int(res_id)))
        return {chain_ids[0]: anchors}

    anchors_map = {chain_id: [] for chain_id in chain_ids}
    for anchor_tmp in anchors_tmp:
        if not re.match(r"^[^_]+_\d+_\d+$", anchor_tmp):
            raise ValueError(f"Unsupported renumber_res spec: {anchor_tmp}")
        chain_id, position, res_id = anchor_tmp.split("_")
        if chain_id not in anchors_map:
            raise ValueError(f"renumber_res references unknown chain: {chain_id}")
        anchors_map[chain_id].append((int(position), int(res_id)))
    return anchors_map

def resolve_renumber_res(sequence_length: int, anchors: list) -> list:
    """返回 1..sequence_length 每个位置的新 res_id。

    与 design_annotation.resolve_renumber_res 同逻辑, 仅 anchor 位置改为 1-based;
    anchors 为空时返回 1..N。
    """
    if not anchors:
        return [i for i in range(1, sequence_length + 1)]

    anchors_positions = list(map(lambda x: x[0], anchors))
    current_shift = 0
    result = []
    for i in range(1, sequence_length + 1):
        if i in anchors_positions:
            res_id = anchors[anchors_positions.index(i)][1]
            current_shift = res_id - i
        result.append(i + current_shift)
    if set(result).__len__() < result.__len__():
        raise ValueError(f"renumber spec causes duplicated res_ids")
    return result

def main():
    args = parse_args()
    array = read_structure(args.input)
    chain_ids = bio_struc.get_chains(array)
    anchors_map = parse_renumber_spec(args.renumber_res, chain_ids)

    new_res_ids = array.res_id.copy()
    new_ids_per_chain = {}
    old_ids_per_chain = {}
    for chain_id in chain_ids:
        chain_mask = array.chain_id == chain_id
        chain_array = array[chain_mask]
        old_ids, _ = bio_struc.get_residues(chain_array)
        old_ids_per_chain[chain_id] = list(old_ids)
        n_res = len(old_ids_per_chain[chain_id])
        anchors = anchors_map[chain_id]
        for position, _ in anchors:
            if not (1 <= position <= n_res):
                raise ValueError(f"chain {chain_id}: renumber_res position {position} out of range (1..{n_res})")
        new_ids = resolve_renumber_res(n_res, anchors)
        new_ids_per_chain[chain_id] = new_ids
        chain_atom_indices = np.where(chain_mask)[0]
        res_starts = bio_struc.get_residue_starts(chain_array)
        for k, new_id in enumerate(new_ids):
            start = res_starts[k]
            end = res_starts[k + 1] if k + 1 < len(res_starts) else len(chain_array)
            new_res_ids[chain_atom_indices[start:end]] = new_id

    array.res_id = new_res_ids
    array.ins_code[:] = ""

    output = args.output or str(Path(args.input).with_name(Path(args.input).stem + "_renumbered" + Path(args.input).suffix))
    write_structure(array, output)

    print_with_decoration("Renumbered Structure")
    print(f"{args.input} -> {output}")
    print_decoration_line()
    for chain_id in chain_ids:
        old_ids = old_ids_per_chain[chain_id]
        new_ids = new_ids_per_chain[chain_id]
        print(f"Chain {chain_id}: {len(new_ids)} residues, {old_ids[0]}..{old_ids[-1]} -> {new_ids[0]}..{new_ids[-1]}")
    print_decoration_line()

if __name__ == "__main__":
    main()

import argparse
from pathlib import Path

BRIEF_HEADER = "tid	node_count	idf	rmsd	matching_residues	query_residues"
PER_STRUCTURE_HEADER = "tid	idf	total_match_count	node_count	edge_count	max_node_cov	min_rmsd	nres	plddt	matching_residues	db_key	query_residues"


def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument(
        "-i",
        "--input",
        required=True,
        nargs="+",
        metavar="INPUT",
        help="Path to the result text file",
    )
    p.add_argument(
        "--per-structure",
        action="store_true",
        help="Process output with --per-structure option"
    )

    return p.parse_args()

def main():

    args = parse_args()

    paths: list[Path] = []
    for path in args.input:
        paths.append(Path(path))
        if not paths[-1].exists():
            raise FileNotFoundError(f"{paths[-1]} not found")
        if paths[-1].suffix == ".tsv":
            raise FileExistsError("Do not name result files with .tsv suffix")

    
    for path in paths:
        with open(path) as f:
            content = f.read()
        with open(path.with_suffix(".tsv"), "w") as f:
            if args.per_structure:
                f.write(PER_STRUCTURE_HEADER + "\n" + content)
            else:
                f.write(BRIEF_HEADER + "\n" + content)

if __name__ == "__main__":
    main()

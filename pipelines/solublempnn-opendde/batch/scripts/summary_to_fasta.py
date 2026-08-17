import argparse
from pathlib import Path

def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument("input", help="Path to the input csv")
    p.add_argument("output", nargs="?", help="Path to the output fasta")

    return p.parse_args()

def resolve_output(input, output):

    if not output is None:
        return output
    else:
        return Path(input).with_suffix(".fa")

def main():

    args = parse_args()
    result = {}
    with open(args.input) as f:
        for line in f:
            record_name, sequence = line.strip().split(",")
            result[record_name] = sequence
    output = resolve_output(args.input, args.output)
    with open(output, "w") as f:
        for key, value in result.items():
            f.write(f">{key}\n")
            f.write(f"{value}\n")

if __name__ == "__main__":
    main()
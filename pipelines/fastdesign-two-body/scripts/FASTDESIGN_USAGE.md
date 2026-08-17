# FastDesign Two-Body Workflow

## Assumption

`comparison_107_83.pdb` contains two fragments separated by a numbering gap:

- fragment 1: PDB residues `1-79`
- fragment 2: PDB residues `107-179`

This matches the existing demo files in this directory.

## Design Rules

- designable positions: `6,7,10,11,13,14,17,18,21,24,25,28,32,35,37,39,40,43,44,47,50,51,54,57,58,61,65`
- only these positions can change amino acid identity and repack
- all other residues are fully frozen
- the inter-fragment rigid-body jump remains movable during `FastDesign`

## Files

- `fastdesign_two_body.xml`: RosettaScripts protocol
- `run_fastdesign_batch.py`: preprocessing + multiprocessing launcher

The batch script no longer carries task defaults. Pass all task parameters
explicitly from `main.ipynb` or the command line.

## Example

```bash
cd /Volumes/ZHITAI-TiPlus7100-002/Repositories/PDE002_HCFHB/pipelines/pipeline-3.1.0/0_BackbonePlacement

mamba activate biorazer

python run_fastdesign_batch.py \
  --inputs comparison_107_83.pdb \
  --xml fastdesign_two_body.xml \
  --rosetta-dir /Users/maxim/Applications/Rosetta/rosetta.source.release-408 \
  --outdir fastdesign_runs \
  --nstruct 50 \
  --jobs 8 \
  --design-resids 6,7,10,11,13,14,17,18,21,24,25,28,32,35,37,39,40,43,44,47,50,51,54,57,58,61,65 \
  --fragment1 1-79 \
  --fragment2 107-179 \
  --suffix _fd \
  --extra-flags \
    -in:file:fullatom \
    -ex1 \
    -ex2aro
```

## Outputs

- prepared two-chain input: `fastdesign_runs/prepared_inputs/*.twochain.pdb`
- per-input decoys and scorefile: `fastdesign_runs/<input_name>/`
- Rosetta stdout/stderr logs: `rosetta.log`, `rosetta.err`

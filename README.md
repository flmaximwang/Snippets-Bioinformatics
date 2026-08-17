# Snippets-Bioinformatics
My snippets for bioinformatics

## Protenix2 docs (docs/protenix2)

Vendored from bytedance/Protenix: `docs/infer_json_format.md` + `examples/`
(no upstream git history; shallow sparse clone on sync).

Sync — overwrites same-named files only, never deletes anything:

    ./scripts/update_protenix2_docs.sh            # follow upstream main (GitHub current content)
    ./scripts/update_protenix2_docs.sh v2.0.0     # pin to an upstream tag
    ./scripts/update_protenix2_docs.sh <ref>      # any upstream branch/tag/commit

Current source is recorded in scripts/protenix2_upstream.ref
(after sync, review `git status` and commit manually).

## SolubleMPNN → colabfold-msa → OpenDDE pipeline (pipelines/solublempnn-opendde)

Batch design loop: SolubleMPNN generates one sequence → 3-fragment MSA via
colabfold-msa (--pair-mode paired, no templates) → OpenDDE prediction →
ptm/RMSD consistency check → append accepted candidates to `report/summary.csv`.
Iterates until TARGET_CANDIDATES accepted; idempotent per candidate, resumable.

Run from `pipelines/solublempnn-opendde/batch/`:

    bash scripts/main.sh

Deps: micromamba envs /opt/envs/{ProteinMPNN,BioRazer,OpenDDE}.
Workflow doc: `batch/scripts/README.md` (5 steps, full pipeline details).
`batch/scripts/renumber_protein_structure_file.py` is a symlink to
`snippet/renumber_protein_structure_files.py` (same file, no duplication).

## Dependencies

Get `biorazer` from Github.
# Snippets-Bioinformatics
My snippets for bioinformatics

## PyMOL alignment pitfalls (docs/pymol)

PyMOL 3.x `align`/`super`/`extra_fit` cannot align two selections that share an
object — the C layer builds the target residue list with the mobile object as
`exclude`, so any same-object alignment (even `align(obj and c. B, obj and c. A)`)
is `invalid selections for alignment`.

`extra_fit`'s `reference` argument must be a **bare object name** (it is compared
by string equality to remove the reference from its per-object loop). To align
onto a specific chain, exclude the reference from the mobile selection:

    extra_fit * and not ref_obj, ref_obj and c. A, method=align, cycles=5, cutoff=2.0

Same-object chain superposition: use `pair_fit` (`fit`/`align`/`super` all fail).

Full root-cause analysis + verified working command on RPXDock poses:
`docs/pymol/extra-fit-alignment.md`.

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

## Rosetta FastDesign two-body perturbation (pipelines/fastdesign-two-body)

Batch Rosetta FastDesign around an initial pose for a two-fragment assembly.
Prepares each input PDB as a trimmed single-chain model, writes a two-fragment
fold tree (inter-fragment jump movable), applies RigidBodyPerturbNoCenter, then
FastDesign with only the design positions' chi and the inter-fragment jump
movable — all other residues frozen.

Files in `pipelines/fastdesign-two-body/scripts/`:

- `run_fastdesign_batch.py` — trim + fold-tree + design-position mapping, then
  multiprocessing launcher for `rosetta_scripts`. Task parameters are all CLI
  args (no hardcoded defaults; `--dry-run` prints jobs without launching).
- `fastdesign_two_body.xml` — RosettaScripts protocol. Design positions are
  passed via `design_positions`/`aa_allowed`; each fragment-1 design resid maps
  to a fragment-2 partner at `resid + 100` (missing partners skipped with a
  warning).
- `FASTDESIGN_USAGE.md` — usage doc with a full example command.

Run from `pipelines/fastdesign-two-body/scripts/`:

    python run_fastdesign_batch.py \
      --inputs comparison_107_83.pdb \
      --xml fastdesign_two_body.xml \
      --rosetta-dir /path/to/rosetta.source.release-408 \
      --outdir fastdesign_runs \
      --nstruct 50 --jobs 8 \
      --design-resids 6,7,10,11,13,14,17,18,21,24,25,28,32,35,37,39,40,43,44,47,50,51,54,57,58,61,65 \
      --fragment1 1-79 --fragment2 107-179 \
      --suffix _fd --extra-flags -in:file:fullatom -ex1 -ex2aro

Deps: `biorazer` env (`biorazer_ex.apps.rosetta.execution.RosettaApp`) and a
Rosetta source release (runs its `rosetta_scripts` binary).

Outputs: `prepared_inputs/<stem>.prepared.pdb` + `<stem>.fold_tree.txt`,
per-decoy dirs with PDB + scorefile + `rosetta.log`/`rosetta.err`,
`run_config.json`, `run_summary.json`.

## MASTER structural-motif search — query side (snippet/master_search.sh)

Runs the `master` search binary against a built PDS database (query side only;
does NOT build the DB — that's `snippet/create_master_db.sh` — and does NOT do
PDB→PDS, that's `createPDS`). Given a query `.pds` and a `--targetList`, it
searches the database with piece-wise splitting + resume.

Parallelism & resume both come from **splitting `--targetList`**:

- `--target-num-per-list` (default 10000) divides the list into consecutive
  pieces of ≤ that many structures each; each piece is one `master --targetList
  <piece>` job.
- A `.done.<i>` sentinel in `out/.chunks/` is written only after a piece finishes,
  so re-running the same command skips every already-done piece (resume) and only
  searches the rest, then re-merges.
- `--njobs` is ONLY a concurrency cap (how many pieces run at once to use the
  cores), NOT the split size.
- A fingerprint of (targetList content + `--target-num-per-list`) invalidates all
  sentinels if either changes, so stale "done" can never skip needed work.

Per-piece match/seq go under `out/.chunks/`; match structures go into the shared
`--structOut` dir (pieces are disjoint target sets, so filenames never collide).
Final merged outputs: `out/match.txt`, optional `out/seqs.txt`, `out/structs/`.

    bash snippet/master_search.sh \
      --query query/1akha.pds --targetList target_list.txt --rmsdCut 3.0 \
      --out out --njobs 8 --target-num-per-list 10000 [--bbRMSD] [--topN N]

Query PDB → PDS must be done beforehand:
`createPDS --type query --pdb in.pdb --pds in.pds --dCut 25.0 --dStep 5.0
--phiStep 10.0 --psiStep 10.0` (dCut/dStep/phiStep/psiStep must match the DB
build). Note `--topN`/`--minN` apply per piece, not globally.

Full per-parameter reference (every flag with input + effect examples):
`docs/master/master_search.md`.

## Dependencies

Get `biorazer` from Github.
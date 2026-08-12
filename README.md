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
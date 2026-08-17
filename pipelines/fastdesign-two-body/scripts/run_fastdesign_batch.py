#!/usr/bin/env python3
"""Batch Rosetta FastDesign runner for assembling two backbone fragments.

The script prepares each input PDB as a trimmed single-chain model, then launches
`rosetta_scripts` in parallel via `biorazer_ex.apps.rosetta.execution.RosettaApp`.

Task-specific parameters are intentionally not hardcoded here. Pass them
explicitly from the notebook or command line.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

DEFAULT_AA_ALLOWED = "ADEFGHIKLMNQRSTVWY"
CANONICAL_AA_CODES = set("ACDEFGHIKLMNPQRSTVWY")
FRAGMENT2_DESIGN_OFFSET = 100
DEFAULT_FASTDESIGN_REPEATS = 4
DEFAULT_RB_TRANSLATION = 3.0
DEFAULT_RB_ROTATION = 8.0
ROSETTA_APP_RUN_MODE = "subprocess.run"


@dataclass(frozen=True)
class ResidueRecord:
    pdb_resseq: int
    pdb_icode: str
    chain_id: str
    lines: tuple[str, ...]


@dataclass(frozen=True)
class JobSpec:
    input_pdb: Path
    prepared_pdb: Path
    fold_tree_path: Path
    output_dir: Path
    rosetta_dir: Path
    decoy_index: int
    seed: int
    log_path: Path
    app_args: tuple[str, ...]


@dataclass(frozen=True)
class JobResult:
    name: str
    return_code: int
    elapsed_seconds: float
    error: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare trimmed inputs and batch-run Rosetta FastDesign.",
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="PDB files or directories containing PDB files.",
    )
    parser.add_argument(
        "--xml",
        required=True,
        help="RosettaScripts XML file.",
    )
    parser.add_argument(
        "--rosetta-dir",
        required=True,
        help="Rosetta source release root directory used by biorazer_ex RosettaApp.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for prepared inputs, logs, and decoys.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        required=True,
        help="Number of concurrent Rosetta jobs.",
    )
    parser.add_argument(
        "--nstruct",
        type=int,
        required=True,
        help="Decoys to generate per prepared input.",
    )
    parser.add_argument(
        "--design-resids",
        required=True,
        help="Comma-separated PDB residue numbers to design in fragment 1.",
    )
    parser.add_argument(
        "--aa-allowed",
        default=DEFAULT_AA_ALLOWED,
        metavar="AA_ALLOWED",
        help=(
            "Allowed canonical amino acids at design positions, as a one-letter string. "
            f"Default: {DEFAULT_AA_ALLOWED}"
        ),
    )
    parser.add_argument(
        "--fastdesign-repeats",
        type=int,
        default=DEFAULT_FASTDESIGN_REPEATS,
        help=f"FastDesign repeats value passed into the XML. Default: {DEFAULT_FASTDESIGN_REPEATS}",
    )
    parser.add_argument(
        "--rb-translation",
        type=float,
        default=DEFAULT_RB_TRANSLATION,
        help=(
            "Initial rigid-body translation perturbation magnitude in Angstrom "
            f"before FastDesign. Default: {DEFAULT_RB_TRANSLATION}"
        ),
    )
    parser.add_argument(
        "--rb-rotation",
        type=float,
        default=DEFAULT_RB_ROTATION,
        help=(
            "Initial rigid-body rotation perturbation magnitude in degrees "
            f"before FastDesign. Default: {DEFAULT_RB_ROTATION}"
        ),
    )
    parser.add_argument(
        "--fragment1",
        required=True,
        help="Inclusive PDB residue range for fragment 1.",
    )
    parser.add_argument(
        "--fragment2",
        required=True,
        help="Inclusive PDB residue range for fragment 2.",
    )
    parser.add_argument(
        "--suffix",
        required=True,
        help="Suffix used for Rosetta output filenames.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing prepared inputs and decoys.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without launching Rosetta.",
    )
    parser.add_argument(
        "--extra-flags",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra Rosetta flags appended verbatim after '--'.",
    )
    return parser.parse_args()


def parse_range(text: str) -> tuple[int, int]:
    try:
        start_text, end_text = text.split("-", 1)
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise ValueError(f"Invalid residue range: {text!r}") from exc
    if start > end:
        raise ValueError(f"Residue range start must be <= end: {text!r}")
    return start, end


def parse_design_resids(text: str) -> list[int]:
    values: list[int] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(int(token))
    if not values:
        raise ValueError("At least one design residue must be provided.")
    return values


def parse_aa_allowed(text: str) -> str:
    cleaned = "".join(text.split()).upper()
    if not cleaned:
        raise ValueError("AA_ALLOWED must contain at least one canonical amino acid.")

    invalid = sorted({aa for aa in cleaned if aa not in CANONICAL_AA_CODES})
    if invalid:
        raise ValueError(
            "AA_ALLOWED contains non-canonical amino acid codes: "
            + ",".join(invalid)
        )

    deduplicated = "".join(dict.fromkeys(cleaned))
    return deduplicated


def expand_inputs(paths: Iterable[str]) -> list[Path]:
    pdbs: list[Path] = []
    for raw_path in paths:
        path = Path(raw_path).expanduser().resolve()
        if path.is_dir():
            pdbs.extend(sorted(path.glob("*.pdb")))
        elif path.is_file():
            pdbs.append(path)
        else:
            raise FileNotFoundError(f"Input path does not exist: {path}")
    if not pdbs:
        raise FileNotFoundError("No PDB inputs found.")
    return pdbs


def collect_residues(pdb_path: Path) -> list[ResidueRecord]:
    residues: list[ResidueRecord] = []
    current_key: tuple[str, int, str] | None = None
    current_lines: list[str] = []

    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain_id = line[21]
            resseq = int(line[22:26])
            icode = line[26]
            key = (chain_id, resseq, icode)
            if current_key is None:
                current_key = key
            if key != current_key:
                prev_chain, prev_resseq, prev_icode = current_key
                residues.append(
                    ResidueRecord(
                        pdb_resseq=prev_resseq,
                        pdb_icode=prev_icode,
                        chain_id=prev_chain,
                        lines=tuple(current_lines),
                    )
                )
                current_lines = []
                current_key = key
            current_lines.append(line.rstrip("\n"))

    if current_key is not None:
        prev_chain, prev_resseq, prev_icode = current_key
        residues.append(
            ResidueRecord(
                pdb_resseq=prev_resseq,
                pdb_icode=prev_icode,
                chain_id=prev_chain,
                lines=tuple(current_lines),
            )
        )

    if not residues:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return residues


def write_trimmed_pdb(
    source_pdb: Path,
    prepared_pdb: Path,
    fragment1: tuple[int, int],
    fragment2: tuple[int, int],
) -> None:
    residues = collect_residues(source_pdb)
    frag1_residues = [
        residue
        for residue in residues
        if fragment1[0] <= residue.pdb_resseq <= fragment1[1]
    ]
    frag2_residues = [
        residue
        for residue in residues
        if fragment2[0] <= residue.pdb_resseq <= fragment2[1]
    ]

    if not frag1_residues:
        raise ValueError(
            f"{source_pdb.name}: no residues found in fragment1 range {fragment1}"
        )
    if not frag2_residues:
        raise ValueError(
            f"{source_pdb.name}: no residues found in fragment2 range {fragment2}"
        )

    prepared_pdb.parent.mkdir(parents=True, exist_ok=True)

    with source_pdb.open() as handle:
        header_lines = [
            line.rstrip("\n")
            for line in handle
            if line.startswith(
                ("HEADER", "TITLE", "REMARK", "CRYST1", "SCALE", "ORIGX", "MTRIX")
            )
        ]

    with prepared_pdb.open("w") as out:
        for line in header_lines:
            out.write(f"{line}\n")
        for residue in frag1_residues:
            for atom_line in residue.lines:
                out.write(f"{atom_line}\n")
        for residue in frag2_residues:
            for atom_line in residue.lines:
                out.write(f"{atom_line}\n")
        out.write("END\n")


def write_fold_tree_file(
    fold_tree_path: Path,
    fragment1_length: int,
    fragment2_length: int,
) -> None:
    if fragment1_length < 1 or fragment2_length < 1:
        raise ValueError("Both fragments must contain at least one residue.")

    fragment1_start = 1
    fragment1_end = fragment1_length
    fragment2_start = fragment1_end + 1
    fragment2_end = fragment1_end + fragment2_length
    fragment1_anchor = (fragment1_start + fragment1_end) // 2
    fragment2_anchor = (fragment2_start + fragment2_end) // 2

    lines = [
        "FOLD_TREE",
        f"EDGE {fragment1_anchor} {fragment1_start} -1",
        f"EDGE {fragment1_anchor} {fragment1_end} -1",
        f"EDGE {fragment1_anchor} {fragment2_anchor} 1",
        f"EDGE {fragment2_anchor} {fragment2_start} -1",
        f"EDGE {fragment2_anchor} {fragment2_end} -1",
    ]
    fold_tree_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def format_resnum_for_selector(residue: ResidueRecord) -> str:
    chain_id = residue.chain_id.strip()
    if chain_id:
        return f"{residue.pdb_resseq}{chain_id}"
    return str(residue.pdb_resseq)


def build_command(
    xml_path: Path,
    input_stem: str,
    prepared_pdb: Path,
    fold_tree_path: Path,
    output_dir: Path,
    design_positions: str,
    aa_allowed: str,
    fastdesign_repeats: int,
    rb_translation: float,
    rb_rotation: float,
    decoy_index: int,
    seed: int,
    suffix: str,
    overwrite: bool,
    extra_flags: list[str],
) -> tuple[str, ...]:
    parser_vars = [
        f"design_positions={design_positions}",
        f"aa_allowed={aa_allowed}",
        f"fastdesign_repeats={fastdesign_repeats}",
        f"fold_tree_file={fold_tree_path}",
        f"rb_translation={rb_translation}",
        f"rb_rotation={rb_rotation}",
    ]
    run_tag = f"{input_stem}_{decoy_index:04d}"
    cmd = [
        "-s",
        str(prepared_pdb),
        "-parser:protocol",
        str(xml_path),
        "-parser:script_vars",
        *parser_vars,
        "-nstruct",
        "1",
        "-out:path:all",
        str(output_dir),
        "-out:prefix",
        f"{run_tag}_",
        "-out:suffix",
        suffix,
        "-out:file:scorefile",
        str(output_dir / f"{run_tag}{suffix}.sc"),
        "-packing:use_input_sc",
        "-run:constant_seed",
        "1",
        "-run:jran",
        str(seed),
    ]
    if overwrite:
        cmd.append("-overwrite")
    if extra_flags:
        cmd.extend(extra_flags)
    return tuple(cmd)


def finalize_output_structure(
    output_dir: Path,
    input_stem: str,
    decoy_index: int,
) -> None:
    pdb_paths = sorted(output_dir.glob("*.pdb"))
    if not pdb_paths:
        return

    expected_path = output_dir / f"{input_stem}_{decoy_index:04d}.pdb"
    if len(pdb_paths) == 1:
        source_path = pdb_paths[0]
    else:
        matching_paths = [path for path in pdb_paths if path.name != expected_path.name]
        if len(matching_paths) != 1:
            raise RuntimeError(
                f"Expected exactly one Rosetta output PDB in {output_dir}, found: "
                f"{[path.name for path in pdb_paths]}"
            )
        source_path = matching_paths[0]

    if source_path == expected_path:
        return
    if expected_path.exists():
        expected_path.unlink()
    source_path.rename(expected_path)


def build_design_positions(
    residues: list[ResidueRecord],
    fragment1: tuple[int, int],
    fragment2: tuple[int, int],
    design_resids: list[int],
) -> tuple[str, list[int]]:
    fragment1_lookup = {
        residue.pdb_resseq: residue
        for residue in residues
        if fragment1[0] <= residue.pdb_resseq <= fragment1[1]
    }
    fragment2_lookup = {
        residue.pdb_resseq: residue
        for residue in residues
        if fragment2[0] <= residue.pdb_resseq <= fragment2[1]
    }

    selector_tokens: list[str] = []
    missing_fragment2: list[int] = []

    for resid in design_resids:
        residue1 = fragment1_lookup.get(resid)
        if residue1 is None:
            raise ValueError(f"Design residue {resid} not found in fragment1.")
        selector_tokens.append(format_resnum_for_selector(residue1))

        partner_resid = resid + FRAGMENT2_DESIGN_OFFSET
        residue2 = fragment2_lookup.get(partner_resid)
        if residue2 is None:
            missing_fragment2.append(partner_resid)
            continue
        selector_tokens.append(format_resnum_for_selector(residue2))

    deduplicated_tokens = list(dict.fromkeys(selector_tokens))
    return ",".join(deduplicated_tokens), missing_fragment2


def import_rosetta_app():
    try:
        from biorazer_ex.apps.rosetta.execution import RosettaApp
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Failed to import biorazer_ex. Activate the 'biorazer' environment first."
        ) from exc
    return RosettaApp


def format_preview_command(rosetta_dir: Path, app_args: Iterable[str]) -> str:
    cmd = [
        "RosettaApp",
        f"[app_dir={rosetta_dir}]",
        "rosetta_scripts",
        *app_args,
    ]
    return " ".join(str(part) for part in cmd)


def write_run_config(
    outdir: Path,
    args: argparse.Namespace,
    xml_path: Path,
    rosetta_dir: Path,
    fragment1: tuple[int, int],
    fragment2: tuple[int, int],
    design_resids: list[int],
    aa_allowed: str,
    jobs: list[JobSpec],
) -> None:
    config = {
        "inputs": [str(Path(path).expanduser().resolve()) for path in args.inputs],
        "xml": str(xml_path),
        "rosetta_dir": str(rosetta_dir),
        "rosetta_app_run_mode": ROSETTA_APP_RUN_MODE,
        "outdir": str(outdir),
        "jobs": args.jobs,
        "nstruct": args.nstruct,
        "design_resids": design_resids,
        "aa_allowed": aa_allowed,
        "fastdesign_repeats": args.fastdesign_repeats,
        "rb_translation": args.rb_translation,
        "rb_rotation": args.rb_rotation,
        "fragment1": {"start": fragment1[0], "end": fragment1[1]},
        "fragment2": {"start": fragment2[0], "end": fragment2[1]},
        "suffix": args.suffix,
        "overwrite": args.overwrite,
        "dry_run": args.dry_run,
        "extra_flags": args.extra_flags,
        "job_outputs": [
            {
                "input_pdb": str(job.input_pdb),
                "prepared_pdb": str(job.prepared_pdb),
                "fold_tree_path": str(job.fold_tree_path),
                "decoy_index": job.decoy_index,
                "seed": job.seed,
                "output_dir": str(job.output_dir),
                "log_path": str(job.log_path),
                "app_args": list(job.app_args),
            }
            for job in jobs
        ],
    }
    (outdir / "run_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def write_run_summary(
    outdir: Path,
    results: list[JobResult],
    batch_elapsed_seconds: float,
) -> None:
    summary = {
        "job_count": len(results),
        "success_count": sum(1 for result in results if result.return_code == 0),
        "failure_count": sum(1 for result in results if result.return_code != 0),
        "batch_elapsed_seconds": batch_elapsed_seconds,
        "jobs": [
            {
                "name": result.name,
                "return_code": result.return_code,
                "elapsed_seconds": result.elapsed_seconds,
                "error": result.error,
            }
            for result in results
        ],
    }
    (outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def make_jobs(args: argparse.Namespace) -> list[JobSpec]:
    xml_path = Path(args.xml).expanduser().resolve()
    if not xml_path.is_file():
        raise FileNotFoundError(f"RosettaScripts XML not found: {xml_path}")
    if not args.rosetta_dir:
        raise FileNotFoundError(
            "Rosetta directory is not set. Pass --rosetta-dir /path/to/rosetta.source.release-XXX"
        )
    rosetta_dir = Path(args.rosetta_dir).expanduser().resolve()
    if not rosetta_dir.is_dir():
        raise FileNotFoundError(f"Rosetta directory not found: {rosetta_dir}")
    input_pdbs = expand_inputs(args.inputs)
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prepared_dir = outdir / "prepared_inputs"
    fragment1 = parse_range(args.fragment1)
    fragment2 = parse_range(args.fragment2)
    design_resids = parse_design_resids(args.design_resids)
    aa_allowed = parse_aa_allowed(args.aa_allowed)

    jobs: list[JobSpec] = []
    for input_pdb in input_pdbs:
        residues = collect_residues(input_pdb)
        design_positions, missing_fragment2 = build_design_positions(
            residues=residues,
            fragment1=fragment1,
            fragment2=fragment2,
            design_resids=design_resids,
        )
        if missing_fragment2:
            print(
                f"[WARNING] {input_pdb.name}: fragment2 design positions not present and will be skipped: {missing_fragment2}",
                file=sys.stderr,
                flush=True,
            )

        prepared_pdb = prepared_dir / f"{input_pdb.stem}.prepared.pdb"
        fold_tree_path = prepared_dir / f"{input_pdb.stem}.fold_tree.txt"
        if prepared_pdb.exists() and not args.overwrite:
            pass
        else:
            write_trimmed_pdb(input_pdb, prepared_pdb, fragment1, fragment2)
        fragment1_length = sum(
            1 for residue in residues if fragment1[0] <= residue.pdb_resseq <= fragment1[1]
        )
        fragment2_length = sum(
            1 for residue in residues if fragment2[0] <= residue.pdb_resseq <= fragment2[1]
        )
        write_fold_tree_file(
            fold_tree_path=fold_tree_path,
            fragment1_length=fragment1_length,
            fragment2_length=fragment2_length,
        )

        struct_output_dir = outdir / input_pdb.stem
        struct_output_dir.mkdir(parents=True, exist_ok=True)
        for decoy_index in range(1, args.nstruct + 1):
            seed = len(jobs) + 1
            output_dir = struct_output_dir / f"{decoy_index:04d}"
            output_dir.mkdir(parents=True, exist_ok=True)
            log_path = output_dir / f"rosetta_{decoy_index:04d}.log"
            app_args = build_command(
                xml_path=xml_path,
                input_stem=input_pdb.stem,
                prepared_pdb=prepared_pdb,
                fold_tree_path=fold_tree_path,
                output_dir=output_dir,
                design_positions=design_positions,
                aa_allowed=aa_allowed,
                fastdesign_repeats=args.fastdesign_repeats,
                rb_translation=args.rb_translation,
                rb_rotation=args.rb_rotation,
                decoy_index=decoy_index,
                seed=seed,
                suffix=args.suffix,
                overwrite=args.overwrite,
                extra_flags=args.extra_flags,
            )
            jobs.append(
                JobSpec(
                    input_pdb=input_pdb,
                    prepared_pdb=prepared_pdb,
                    fold_tree_path=fold_tree_path,
                    output_dir=output_dir,
                    rosetta_dir=rosetta_dir,
                    decoy_index=decoy_index,
                    seed=seed,
                    log_path=log_path,
                    app_args=app_args,
                )
            )
    write_run_config(
        outdir=outdir,
        args=args,
        xml_path=xml_path,
        rosetta_dir=rosetta_dir,
        fragment1=fragment1,
        fragment2=fragment2,
        design_resids=design_resids,
        aa_allowed=aa_allowed,
        jobs=jobs,
    )
    return jobs


def _pick_rosetta_scripts_bin(rosetta_dir: Path) -> Path:
    bin_dir = rosetta_dir / "main" / "source" / "bin"
    candidates = sorted(bin_dir.glob("rosetta_scripts.default*"))
    if not candidates:
        raise FileNotFoundError(f"Could not find rosetta_scripts under {bin_dir}")

    non_default = [path for path in candidates if ".default." not in path.name]
    if non_default:
        return non_default[0]
    return candidates[0]


def run_job(job: JobSpec) -> JobResult:
    job_label = (
        f"{job.input_pdb.name} [{job.decoy_index}]"
        f" seed={job.seed}"
    )
    started_at = time.perf_counter()
    try:
        RosettaApp = import_rosetta_app()
        app = RosettaApp(
            app_dir=job.rosetta_dir,
            logger=f"rosetta.{job.input_pdb.stem}.{job.decoy_index:04d}",
        )
        if job.log_path.exists():
            job.log_path.unlink()
        app.set_default_logger_style(
            handler_types=("FileHandler",), file_path=job.log_path
        )
        app.bin = _pick_rosetta_scripts_bin(job.rosetta_dir)
        print(f"[START] {job_label}", flush=True)
        app.run(
            *job.app_args,
            cwd=job.output_dir,
            get_output=False,
            verbose=True,
            mode=ROSETTA_APP_RUN_MODE,
        )
        finalize_output_structure(
            output_dir=job.output_dir,
            input_stem=job.input_pdb.stem,
            decoy_index=job.decoy_index,
        )
        elapsed_seconds = time.perf_counter() - started_at
        return JobResult(
            name=job_label,
            return_code=0,
            elapsed_seconds=elapsed_seconds,
        )
    except Exception as exc:
        elapsed_seconds = time.perf_counter() - started_at
        return JobResult(
            name=job_label,
            return_code=1,
            elapsed_seconds=elapsed_seconds,
            error=str(exc),
        )


def main() -> int:
    args = parse_args()
    jobs = make_jobs(args)
    batch_started_at = time.perf_counter()
    outdir = Path(args.outdir).expanduser().resolve()

    print(f"Prepared {len(jobs)} Rosetta job(s).")

    if args.dry_run:
        for job in jobs:
            print(
                f"[DRY-RUN] {job.input_pdb.name} [{job.decoy_index}/{args.nstruct}] seed={job.seed}",
                flush=True,
            )
        return 0
    try:
        import_rosetta_app()
    except ModuleNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    worker_count = min(max(1, args.jobs), len(jobs))
    with mp.Pool(worker_count) as pool:
        results: list[JobResult] = []
        total = len(jobs)
        for completed, result in enumerate(pool.imap_unordered(run_job, jobs), start=1):
            results.append(result)
            status = "OK" if result.return_code == 0 else "ERROR"
            print(
                f"[{completed}/{total}] {status} {result.name} "
                f"elapsed={result.elapsed_seconds:.1f}s",
                flush=True,
            )
            if result.error:
                print(result.error, file=sys.stderr, flush=True)

    failures = [result.name for result in results if result.return_code != 0]
    batch_elapsed_seconds = time.perf_counter() - batch_started_at
    write_run_summary(
        outdir=outdir,
        results=results,
        batch_elapsed_seconds=batch_elapsed_seconds,
    )
    if failures:
        print("Failed jobs:")
        for name in failures:
            print(f"  - {name}")
        return 1

    print(f"All Rosetta jobs completed successfully in {batch_elapsed_seconds:.1f}s.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

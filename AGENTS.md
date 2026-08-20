# reverlor — AGENTS.md

**RE**peat **VER**ification using **LO**ng **R**eads.

Early-stage project. Single-author (Maksim Sikolenko). Work happens on the `dev` branch.

## Architecture

- `reverlor/` is a Python 3 package built with flit (`pyproject.toml`)
- Three entry points (all in `reverlor/`):
  - `reverlor.py` — full pipeline: calls `find_repeats_minimap2()` then `find_unresolved_repeats()`
  - `reverlor_find.py` — repeat search only
  - `reverlor_verify.py` — repeat verification only
- `src/` modules:
  - `find_repeats_minimap2.py` — self-contained repeat finder using `mappy` (Python minimap2 bindings). No subprocess calls.
  - `verify_repeats.py` — verification logic: checks whether long reads span predicted repeat regions using `pysam`
  - `bed_lib.py` — BED I/O, merging with `pybedtools`, data classes (`RepeatRegion`, `VerifyResult`)
  - `CoordIntersecter.py` — BAM coordinate intersection via `pysam`
  - `FindArgs.py`, `VerifyArgs.py`, `ReverlorArgs.py` — argparse + typed arg containers
  - `defaults.py` — default parameter values
  - `_version.py` — version string
  - `reverlor_logging.py` — logging setup with 4 verbosity levels
  - `util.py` — file-removal helpers
- `tests/` — pytest tests:
  - `test_CoordIntersecter.py` — coordinate intersection via `pysam`
  - `test_bed_lib.py` — BED I/O, `merge_features`
  - `test_find_repeats_minimap2.py` — repeat finder end-to-end (inserts known repeats, checks detection)
  - `test_util.py` — file/directory removal helpers

## Dependencies

All pip-installable:
- `mappy` — Python bindings for minimap2
- `pysam` — Python bindings for samtools/htslib
- `pybedtools` — Python bindings for bedtools
- `biopython` (`Bio.SeqIO`) — FASTA parsing in verification

## Running

Full pipeline:
```bash
reverlor <input.fasta> <input.bam> <output_dir>
  [--min-repeat-len 200] [--min-repeat-interval 100]
  [--span 5] [--shoulder-len 200]
  [--minimap-k 19] [--minimap-w 19] [--minimap-m 65]
  [--minimap-x {map-ont,lr:hq,map-hifi,map-pb,...}]
  [--samtools-f FLAG] [--samtools-F FLAG]
  [--keep-tmp] [--tmpdir PATH] [-t 1] [-v]
```

Find only:
```bash
reverlor_find <input.fasta> <output_dir> [options]
```

Verify only:
```bash
reverlor_verify <input.fasta> <input.bed> <input.bam> <output_dir> [options]
```

Alternative module invocation:
```bash
python3 -m reverlor <input.fasta> <input.bam> <output_dir> [options]
```

## State of the code

- Version 1.0.0 (2026-08-20).
- `find_repeats_minimap2.py` uses `mappy` directly (no external `minimap2` binary required).
- Verification uses `pysam` (no external `samtools` binary required).
- Merging uses `pybedtools` (requires `bedtools` on `$PATH` at runtime).
- Tests exist in `tests/` (pytest). Test data for `CoordIntersecter` and `find_repeats` is in `tests/data/`.
- No CI, no linter/formatter config, no pre-commit hooks.

## Code style peculiarities

Python code blocks must end with `# end if`, `# end for`, `# end while`, `# end with`, `# end def`, `# end class` etc. comments.

## Import convention

- Entry-point modules in `reverlor/` import from `src/` via relative imports: `from .src.Module import ...`
- Modules inside `src/` import from each other via relative imports: `from .Module import ...`
- Tests import via absolute path: `from reverlor.src.Module import ...`

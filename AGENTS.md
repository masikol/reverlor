# reverlor — AGENTS.md

**RE**peat **VER**ification using **LO**ng **R**eads.

Early-stage project (v0.0.a). Single-author (Maksim Sikolenko). Work happens on the `dev` branch.

## Architecture

- `reverlor/` is a flat Python 3 package (no `pyproject.toml`/`setup.py`/`setup.cfg`)
- Entry point: `reverlor/reverlor.py` — calls `reverlor_find()` then `reverlor_verify()`
- `reverlor_find` — parses args via `src/FindArgs.py`, finds repeats via `src/find_repeats` (module expected at `reverlor/src/find_repeats.py` but only `find_repeats_OLD.py` exists; the `_OLD` version works but is not importable under the canonical name)
- `reverlor_verify` — **stub** (`from src.reverlor_verify import ???` is a placeholder). No verification logic exists yet.
- `param_selection/` — shell scripts and test data (FASTA files for *Mycoplasma* strains)

## External dependencies (not pip-installable)

These must be on `$PATH` or passed via `--seqkit`/`--minimap2`/`--samtools`/`--bedtools`/`--con-hi`:
- `seqkit`, `minimap2`, `samtools`, `bedtools` — common bioinformatics tools
- `con-hi.py` — external script (separate project)

## Running

```bash
python3 reverlor/reverlor.py <input.fasta> <output_dir> [--window-size 1000] [--min-repeat-len 200]
```

## State of the code

- `reverlor_verify.py:33` has a broken placeholder import (`???`). Do **not** run the pipeline end-to-end; `reverlor_find` may work independently.
- No tests, no CI, no linter/formatter config, no pre-commit hooks.
- `reverlor/src/find_repeats_OLD.py` is the working implementation; it runs seqkit→minimap2→samtools→con-hi.py in a subprocess pipeline with both no-rotate and rotate modes, then merges results with bedtools.

## Code style peculiarities

Python code blocks must end with `# end if`, `# end for`, `# end while`, `# end with`, `# end def`, `# end class` etc. comments.

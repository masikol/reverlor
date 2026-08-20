# REVERLOR

**RE**peat **VER**ification using **LO**ng **R**eads

A bioinformatic tool that finds exact and inexact interspersed repeats in genomic sequences via [minimap2](https://github.com/lh3/minimap2) self-alignment, then verifies whether those repeats are spanned by long reads from a BAM alignment file. If less than *N* reads span a repeat, it is considered *unresolved*.

---

## Table of Contents

- [Motivation](#motivation)
- [Pipeline Overview](#pipeline-overview)
- [Limitations](#limitations)
- [Installation](#installation)
- [Usage](#usage)
  - [Quick Start](#quick-start)
  - [BAM File Preparation](#bam-file-preparation)
  - [Full Pipeline (reverlor)](#full-pipeline-reverlorpy)
  - [Find Only (reverlor_find)](#find-only-reverlor_findpy)
  - [Verify Only (reverlor_verify)](#verify-only-reverlor_verifypy)
- [Options Reference](#options-reference)
- [Output Files](#output-files)
- [Examples](#examples)

---

## Motivation

1. Accurate reconstruction of repeat regions (e.g., rRNA clusters, IS elements, prophages) is necessary for complete and contiguous genome assemblies.

2. Even long reads (Oxford Nanopore, PacBio) sometimes fail to span a repeat end-to-end.

3. A program that would find such repeats and verify them (check how many reads span them end-to-end) would be useful.

---

## Pipeline Overview

### Step 1: Repeat Finding (`reverlor_find`)

1. Self-align the reference FASTA using minimap2 (via `mappy` Python package).
2. Extract alignment coordinates as raw BED entries.
3. Merge overlapping/adjacent repeats using `pybedtools`.
4. Filter by minimum repeat length.
5. Output: `repeats.bed`.

### Step 2: Repeat Verification (`reverlor_verify`)

1. Read repeat regions from the BED file created at Step 1.
2. For each repeat, compute **“shoulder” coordinates** (flanking regions upstream and downstream).
3. Use `pysam` to find long reads spanning from the upstream shoulder position to the downstream one.
4. If fewer than `--span` reads span the repeat, consider it **unresolved**.
5. Output: `unresolved_repeats.bed`.

---

## Limitations

Reverlor reports tandem repeats as a single repeat.

Reverlor does not cluster repeat families.

Reverlor offers no explicit and straightforward control over minimum repeat sequence identity to be reported. It can though be controled with `--minimap-m` option, which is the minimum minimap2 chaining score. The lower `--minimap-m` is, the more dissimilar repeats shall be reported. Chaining score equals the approximate number of matching bases minus a concave gap penalty (see [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html)).

---

## Installation

Reverlor was tested using Python 3.12.3.

### 1. From PyPI

We recommend to install reverlor to a separate Python virtual environment.

#### 1.1 Using [uv](https://docs.astral.sh/uv/)

```bash
# Create a virtual environment in, for example, ./reverlor_venv directory
uv venv --python 3.12.3 ./reverlor_venv
# Activate the environment
source ./reverlor_venv/bin/activate
# Install
uv pip install reverlor
```

#### 1.2 Using pip

```bash
# Create a virtual environment in, for example, ./reverlor_venv directory
python3 -m venv ./reverlor_venv
# Activate the environment
source ./reverlor_venv/bin/activate
# Install
pip install reverlor
```

### 2. From source

```bash
# Get source code
git clone git@github.com:masikol/reverlor.git # Or download a release archive from https://github.com/masikol/reverlor/releases
cd reverlor/
# Create a virtual environment
uv venv --python 3.12.3 ./reverlor_venv
# Activate the environment
source ./reverlor_venv/bin/activate
# Install the “build” package
uv pip install --upgrade build
# Build reverlor package. This will produce a ./dist directory
python3 -m build --installer uv
# Install reverlor (assuming its version is 1.0.0)
uv pip install dist/reverlor-1.0.0-py3-none-any.whl

# If you want to test your installation, install pytest
uv pip install pytest==8.4.2
# Test your installation
python3 -m pytest tests
```
---

## Usage

### Quick Start

```bash
# Full pipeline: find repeats + verify with long reads
reverlor genome.fasta reads.bam output_dir/

# With verbose output
reverlor -vv genome.fasta reads.bam output_dir/

# With custom parameters
reverlor \
    --min-repeat-len 127 \
    --span 10 \
    --threads 4 \
    genome.fasta reads.bam output_dir/
```

### BAM File Preparation

For best performance, we strongly recommend to pass **sorted and indexed** BAM files to `reverlor`.

Here is a quick example on how to create such file:

```bash
minimap2 -a genome.fasta reads.fastq.gz \
    | samtools view -F 4 \
    | samtools sort -O BAM -o reads.sorted.bam

samtools index reads.sorted.bam
```

### Full Pipeline (`reverlor`)

Runs both repeat finding and verification in sequence.

```bash
reverlor [options] <fasta> <input_bam> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `fasta` | Path to input FASTA file (reference genome) |
| `input_bam` | Path to input BAM file (sorted, indexed long reads) |
| `output_dir` | Path to output directory |

---

### Find Only (`reverlor_find`)

Finds repeats without verification.

```bash
reverlor_find [options] <input_fasta> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `input_fasta` | Path to input FASTA file |
| `output_dir` | Path to output directory |

---

### Verify Only (`reverlor_verify`)

Verifies pre-existing repeat predictions against a **sorted and indexed** BAM file.

```bash
reverlor_verify [options] <input_fasta> <input_bed> <input_bam> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `input_fasta` | Path to input FASTA file |
| `input_bed` | Path to input BED file (repeat predictions) |
| `input_bam` | Path to **sorted and indexed** input BAM file |
| `output_dir` | Path to output directory |

---

## Options Reference

### Common Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-V, --version` | flag | -- | Show version and exit |
| `-v, --verbose` | count | 0 | Verbosity level: `-v` (warnings), `-vv` (info), `-vvv` (debug) |
| `-t, --threads` | int | 1 | Number of CPU threads (passed to minimap2) |
| `--keep-tmp` | flag | False | Keep temporary files |
| `--tmpdir` | str | system | Temporary directory path |

### Repeat Finding Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-l/--min-repeat-len` | int | 200 | Minimum repeat length to report |
| `-i/--min-repeat-interval` | int | 100 | Minimum interval between repeats (shorter ones get merged) |
| `-k/--minimap-k` | int | 19 | minimap2 k-mer length |
| `-w/--minimap-w` | int | 19 | minimap2 minimizer window size |
| `-m/--minimap-m` | int | 65 | minimap2 min chain score |
| `--minimap-x` | choice | None | minimap2 preset: `map-ont`, `lr:hq`, `map-hifi`, `map-pb`, `map-iclr`, `asm5`, `asm10`, `asm20` |

### Repeat Verification Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-s/--span` | int | 5 | Minimum number os spanning reads to consider a repeat resolved |
| `-u/--shoulder-len` | int | 200 | Shoulder length (bp) for coordinate checking |
| `--samtools-f` | int | [] | samtools `-f` flag (includes segments having these flags, repeatable) |
| `--samtools-F` | int | [256] | samtools `-F` flag (excludes segments having these flags, repeatable) |

#### Filter reads using SAM flags (`--samtools-f, --samtools-F`)

- These flags are applied at Step 2.3 (see the [Pipeline Overview](#pipeline-overview) above).
- Default `-F 256` excludes secondary alignments from verification.

- Add `--samtools-F 2048` to also exclude supplementary alignments:
```bash
reverlor \
    --samtools-F 256 \
    --samtools-F 2048 \
    genome.fasta \
    reads.sorted.bam \
    output/dir/
```

- Add `--samtools-f 16` to only include reverse-complementedly mapped reads.
```bash
reverlor \
    --samtools-F 256 \
    --samtools-f 16 \
    genome.fasta \
    reads.sorted.bam \
    output/dir/
```

---

## Output Files

Reverlor outputs BED files.

In `unresolved_repeats.bed`, the 5th column holds the number of spanning reads.

### Full Pipeline (`reverlor`)

| File | Description |
|------|-------------|
| `repeats.bed` | Detected repeat regions |
| `unresolved_repeats.bed` | Repeats with insufficient spanning coverage |

### Find Only (`reverlor_find`)

| File | Description |
|------|-------------|
| `repeats.bed` | Detected repeat regions |

### Verify Only (`reverlor_verify`)

| File | Description |
|------|-------------|
| `unresolved_reapeats.bed` | Repeats with insufficient spanning coverage |


---

## Examples

### Example 1: Basic usage on a bacterial genome

```bash
reverlor \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 2: Exclude secondary and supplementary alignments

```bash
reverlor \
    --samtools-F 256 \
    --samtools-F 2048 \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 3: Find repeats only, no verification (no BAM required)

```bash
reverlor_find \
    genome.fasta \
    output_dir/
```

### Example 4: Verify existing repeat predictions (`repeats.bed`)

```bash
reverlor_verify \
    --span 10 \
    --shoulder-len 300 \
    genome.fasta \
    repeats.bed \
    reads.sorted.bam \
    output_dir/
```

### Example 5: Debug mode with temporary files

```bash
reverlor \
    -vvv \
    --keep-tmp \
    --tmpdir /tmp/reverlor_debug \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 6: Full pipeline with verbose messages

```bash
reverlor \
    --min-repeat-len 200 \
    --min-repeat-interval 50 \
    --span 5 \
    --shoulder-len 200 \
    --samtools-F 256 \
    -t 8 \
    -vv \
    genome.fasta \
    reads.sorted.bam \
    output_full/
```
